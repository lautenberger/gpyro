PROGRAM GPYRO_GA

USE GA_IO
USE GA_SUBS
USE GPYRO_INIT
USE GPYRO_VARS

USE SCE_VARS
USE SCE_IO
USE SCE_INIT
USE SCE_UA
USE SCE_FUNCS
USE SCE_OBJFUNC
USE MPI

IMPLICIT NONE

INTEGER       :: J,IGEN,IERR,IGS,NGS2,ICASE
LOGICAL       :: LEX
CHARACTER(3) :: THREE
CHARACTER(15) :: SHELLSTR
CHARACTER(300) :: MESSAGE

!SCE MPI VARIABLES:
INTEGER  :: IRANK,IDEST,IINDIV,I,IGENE,K1,K2
INTEGER  :: ISTATUS(MPI_STATUS_SIZE)
REAL(EB), ALLOCATABLE, DIMENSION(:) :: SEND, RECV, XFTEMP
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: XTEMP
LOGICAL :: ALLRECEIVED,DONE

!SHC VARIABLES:
REAL(EB) :: PMUT,RVEC(2),SPAN,OLDINDIV(1:100),OLDFIT

!For random number generator:
INTEGER :: M
INTEGER, ALLOCATABLE, DIMENSION(:) :: K

! Initialize system clock
CALL SYSTEM_CLOCK(COUNT_RATE=CLOCK_COUNT_RATE)
CALL SYSTEM_CLOCK(COUNT_MAX=CLOCK_COUNT_MAX)

IGPYRO_TYPE = 3 !Material pyrolysis property estimation

CALL GET_OPERATING_SYSTEM

! MPI initialization
GA%MPI=.FALSE.
CALL MPI_INIT(IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,GA%IRANK,IERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,GA%NPROC,IERR)
IF (GA%NPROC .GT. 1) GA%MPI = .TRUE.

IF ((GA%MPI .AND. GA%IRANK .EQ. 0) .OR. (.NOT. GA%MPI)) THEN
   WRITE(*,*) 'Gpyro 0.8200'
   WRITE(*,*) 'Build date:  August 29, 2020'
   WRITE(*,*) 'Material property estimation program (gpyro_propest).'
   WRITE(*,*) 'https://github.com/reaxfire/gpyro'
   WRITE(*,*)
   WRITE(*,'(A,I4)') 'Number of cores (GA%NPROC): ', GA%NPROC
   WRITE(*,*)
ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
WRITE(*,'(A,I4,A)') 'Process ', GA%IRANK, ' initializing.'

! Read input files, allocate arrays, initialize variables
CALL READ_GPYRO
CALL READ_GA_INPUT_FILES 

CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

DO ICASE = 1, GPG%NCASES
   CALL ALLOCATE_GPYRO(GPG%IMESH(ICASE))
   CALL INIT_GPYRO_VARS(ICASE,GPG%IMESH(ICASE))
ENDDO

CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
CALL GA_INIT
CALL CHECK_GPYRO(GA%IRANK)
CALL CHECK_GPYRO_GA
CALL GET_EXPERIMENTAL_DATA 

! Initialize random number generator
CALL RANDOM_SEED(SIZE=M)
ALLOCATE(K(M))
K(:) = SCE%ISEED !+ IRANK !Make sure each process gets a different seed
CALL RANDOM_SEED(PUT=K(1:M))
DEALLOCATE(K)

CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

! Delete stop.stop file if it exists
IF ((GA%MPI .AND. GA%IRANK .EQ. 0) .OR. (.NOT. GA%MPI)) THEN
   INQUIRE(FILE='stop.stop',EXIST=LEX)
   IF (LEX) THEN
      SHELLSTR = TRIM(DELETECOMMAND) // ' ' // 'stop.stop'
      CALL SYSTEM(TRIM(SHELLSTR)) 
   ENDIF
ENDIF

IF (TRIM(GA%OPTIMIZATION_TYPE) .EQ. 'GA') THEN

   IF (.NOT. GA%MPI .OR. GA%MPI .AND. GA%IRANK .EQ. 0) THEN 
      IF (GA%RESTART) THEN
         CALL RESTART(2) !Read population from restart file
      ELSE
         CALL POPULATE(1,GA%NINDIV,1,GA%NGENE) !Set initial population
         CALL NORMALIZE_MASS_FRACTIONS
      ENDIF
      IF (GA%MPI) CALL SEND_POPULATION 
   ELSE
      IF (GA%MPI) CALL RECV_POPULATION
   ENDIF

! Begin main stepping loop 
   DO IGEN=1,GA%NGEN
      CALL QUALITY(IGEN)

      IF ((GA%MPI .AND. GA%IRANK .EQ. 0) .OR. (.NOT. GA%MPI)) THEN
         CALL SIMULATED_ANNEALING(IGEN)
         CALL RANK_INDIVIDUALS
         CALL DUMP_GA(IGEN)
         IF(IGEN .LT. GA%NGEN) CALL REPRODUCTION
         CALL DUMP_COPIES(IGEN)
         WRITE(*,'(A,I4,A,E11.3)') 'Generation ', IGEN,' fitness: ',GA%FIT(1,0)

         CALL RESTART(1) !Dump restart file
         INQUIRE(FILE='stop.stop',EXIST=LEX)
         IF (LEX) GA%INDIV(:,:) = -9D30 !This sends signal to shutdown
         IF (GA%MPI) CALL SEND_POPULATION
         IF (LEX) THEN
            MESSAGE='stop.stop file found, shutting down.'
            CALL SHUTDOWN_GRACEFULLY(MESSAGE)
         ENDIF
      ELSE
         IF (GA%MPI) THEN
            CALL RECV_POPULATION

            IF (GA%INDIV(1,1) .LT. -8D30) THEN
               MESSAGE='Slave process received shutdown signal, shutting down.'
               CALL SHUTDOWN_GRACEFULLY(MESSAGE)
            ENDIF 
         ENDIF
      ENDIF

   ENDDO

   WRITE(THREE,'(I3.3)') GA%IRANK
   MESSAGE='IRANK ' // THREE // ' shutting down'         
   CALL SHUTDOWN_GRACEFULLY(MESSAGE)

ENDIF

IF (TRIM(GA%OPTIMIZATION_TYPE) .EQ. 'SHC') THEN

   IGEN = 1

   IF (.NOT. GA%MPI .OR. GA%MPI .AND. GA%IRANK .EQ. 0) THEN 
      CALL POPULATE(1,GA%NINDIV,1,GA%NGENE) !Set initial population
      CALL NORMALIZE_MASS_FRACTIONS
      IF (GA%MPI) CALL SEND_POPULATION 
   ELSE
      IF (GA%MPI) CALL RECV_POPULATION
   ENDIF

!   WRITE(*,'(A,I4,A)') 'Process ', GA%IRANK, ' done initializing.'; WRITE(*,*)
      
   CALL QUALITY(IGEN)

   IF ((GA%MPI .AND. GA%IRANK .EQ. 0) .OR. (.NOT. GA%MPI)) THEN
      CALL RANK_INDIVIDUALS
      CALL DUMP_GA(IGEN)
      WRITE(*,'(A,I6,A,E11.3)') 'Generation ', IGEN,' fitness: ',GA%FIT(1,0)
   ELSE
      MESSAGE='Shutting down.'
      CALL SHUTDOWN_GRACEFULLY(MESSAGE)
   ENDIF
   
   GA%NINDIV = 1
   GA%MPI = .FALSE.

   DO IGEN=2, GA%NGEN

      OLDFIT = GA%FIT(1,0)
      OLDINDIV(1:GA%NGENE) = GA%INDIV(1,1:GA%NGENE)

!Mutate
      DO IGENE=1,GA%NGENE
         PMUT = GA%VAR(IGENE)%PMUT
         CALL RANDOM_ARRAY(RVEC,2)
         IF (RVEC(1) .LT. PMUT) THEN
            SPAN = GA%VAR(IGENE)%MAXVAL-GA%VAR(IGENE)%MINVAL
            GA%INDIV(1,IGENE) = GA%INDIV(1,IGENE) + (RVEC(2)-0.5D0)*GA%VAR(IGENE)%VMUTMAX*SPAN
         ENDIF
      ENDDO

      CALL QUALITY(IGEN)
            
      IF (GA%FIT(1,0) .LE. OLDFIT) THEN
         GA%FIT(1,0) = OLDFIT
         GA%INDIV(1,1:GA%NGENE) = OLDINDIV(1:GA%NGENE)
         IF (MOD(IGEN, 100) .EQ. 0) THEN 
            CALL DUMP_GA(IGEN)
            WRITE(*,'(A,I6,A,E11.3)') 'Generation ', IGEN,' fitness: ',GA%FIT(1,0)
         ENDIF
      ELSE
         CALL DUMP_GA(IGEN)
         WRITE(*,'(A,I6,A,E11.3)') 'Generation ', IGEN,' fitness: ',GA%FIT(1,0)
      ENDIF
      
      INQUIRE(FILE='stop.stop',EXIST=LEX)
      
      IF (LEX) THEN
         MESSAGE='stop.stop file found, shutting down.'   
         CALL SHUTDOWN_GRACEFULLY(MESSAGE)
      ENDIF

   ENDDO

   MESSAGE='Shutting down SHC'
   CALL SHUTDOWN_GRACEFULLY(MESSAGE)

ENDIF !SHC

IF (TRIM(GA%OPTIMIZATION_TYPE) .EQ. 'SCE') THEN

   DO IINDIV = 1, GA%NGENE
      SCE%BL(IINDIV) = GA%VAR(IINDIV)%MINVAL
      SCE%BU(IINDIV) = GA%VAR(IINDIV)%MAXVAL
   ENDDO

   IF (SCE%NOPT .NE. GA%NGENE) THEN
      WRITE(*,*) 'Error:  NOPT must be equal to GA%NGENE'
      STOP
   ENDIF
      
   CALL SCE_INITIALIZE(GA%MPI, GA%IRANK)

   ALLOCATE( RECV  (1:SCE%NOPT + SCE%NPT*(SCE%NOPT + 1)) ) 
   ALLOCATE( SEND  (1:SCE%NOPT + SCE%NPT*(SCE%NOPT + 1)) )
 
   ALLOCATE( XFTEMP(1:SCE%NPT) ) 
   ALLOCATE( XTEMP (1:SCE%NPT, 1:GA%NGENE) )
   
   !Determine objective function for initial randomly generated population
   IF (GA%MPI) THEN

      SCE%IRANK_STATUS(:) = 1

      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

      IF (GA%IRANK .EQ. 0) THEN
 
         IINDIV = 1
 
         DO WHILE (IINDIV .LE. SCE%NPT1)
                     
            IRANK = 0

            DO WHILE (IRANK .LT. GA%NPROC-1)
               IRANK = IRANK + 1
               IDEST = 0
               IF (SCE%IRANK_STATUS(IRANK) .EQ. 1) THEN
                  IDEST = IRANK
                  SCE%IRANK_STATUS(IRANK) = 0
                  IRANK = GA%NPROC + 1
               ENDIF
            ENDDO

            IF (IDEST .EQ. 0) THEN !No nodes available, so wait for one to finish

               CALL MPI_RECV(RECV, 3, MPI_DOUBLE_PRECISION, &
               MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, ISTATUS, IERR)

               SCE%IRANK_STATUS(INT(RECV(1))) = 1
               SCE%XF(INT(RECV(2)))=RECV(3)

            ELSE !Send individual to a free node
               SEND(1           ) = REAL(IDEST,EB)
               SEND(2           ) = REAL(IINDIV,EB)
               SEND(3:SCE%NOPT+2) = SCE%X(IINDIV,1:SCE%NOPT)

               SCE%ICALL = SCE%ICALL + 1
               
               CALL MPI_SEND(SEND, SCE%NOPT+2, MPI_DOUBLE_PRECISION, & 
                    IDEST, IINDIV, MPI_COMM_WORLD, IERR)
               
               IINDIV = IINDIV + 1
            ENDIF 

         ENDDO ! (IINDIV .LE. SCE%NPT1)

         ALLRECEIVED = .FALSE. 
         DO WHILE (.NOT. ALLRECEIVED) 
            ALLRECEIVED = .TRUE.
            DO IDEST = 1, GA%NPROC-1
               IF (SCE%IRANK_STATUS(IDEST) .EQ. 0) ALLRECEIVED = .FALSE.
            ENDDO

            IF (.NOT. ALLRECEIVED) THEN
               CALL MPI_RECV(RECV, 3, MPI_DOUBLE_PRECISION, &
               MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, ISTATUS, IERR)

               SCE%IRANK_STATUS(INT(RECV(1))) = 1
               SCE%XF(INT(RECV(2)))=RECV(3)
            
            ENDIF
         ENDDO
 
 !Send message that we're done
         SEND(:) = -1.        
         DO IDEST = 1, GA%NPROC-1 
            CALL MPI_SEND(SEND, SCE%NOPT+2, MPI_DOUBLE_PRECISION, & 
                 IDEST, IINDIV, MPI_COMM_WORLD, IERR)
         ENDDO

      ELSE !(GA%IRANK .EQ. 0) - THIS IS FOR SLAVE NODE

         IINDIV = 1
         DO WHILE (IINDIV .GT. 0) 
            CALL MPI_RECV(RECV, SCE%NOPT+2, MPI_DOUBLE_PRECISION, &
              0    , MPI_ANY_TAG, MPI_COMM_WORLD, ISTATUS, IERR)
            IDEST  = INT(RECV(1))
            IINDIV = INT(RECV(2))
            IF (IINDIV .GT. 0) THEN
               SCE%X(IINDIV,1:SCE%NOPT) = RECV(3:SCE%NOPT+2)
               SCE%xf(IINDIV) = SCE_OBJECTIVE_FUNCTION(SCE%nopt,SCE%x(IINDIV,:))

               SEND(1) = REAL(GA%IRANK,EB)
               SEND(2) = REAL(IINDIV,EB)
               SEND(3) = SCE%XF(IINDIV)
               
               CALL MPI_SEND(SEND, 3, MPI_DOUBLE_PRECISION, &
               0, IINDIV, MPI_COMM_WORLD, IERR)
            ENDIF
         ENDDO !IINDIV .GT. 0
      ENDIF
   ELSE !Serial 
      do IINDIV = 1, SCE%NPT1
         SCE%icall = SCE%icall + 1
         SCE%xf(IINDIV) = SCE_OBJECTIVE_FUNCTION(SCE%nopt,SCE%x(IINDIV,:))
      end do   
   ENDIF
   
!  ARRANGE THE POINTS IN ORDER OF INCREASING FUNCTION VALUE 
   IF (.NOT. GA%MPI .OR. (GA%MPI .AND. GA%IRANK .EQ. 0)) THEN
   
      call SCE_sort(SCE%npt1,SCE%NOPT,SCE%x,SCE%xf,SCE%NPT,SCE%NOPT)
               
      CALL SCE_BEST_WORST !Record best and worst points

!  COMPUTE THE PARAMETER RANGE FOR THE INITIAL POPULATION
      call SCE_parstt(SCE%npt1,SCE%NOPT,SCE%x,SCE%xnstd,SCE%bound,SCE%gnrng,SCE%ipcnvg,SCE%NPT,SCE%NOPT)

      CALL SCE_DUMP
      
      GA%FIT(:,:) = 0.
      GA%FIT(1,0) = -SCE%XF(1)

      GA%INDIV (:,:) = 0.
      GA%INDIV (1,1:GA%NGENE) = SCE%X(1,1:GA%NGENE)
      
      CALL DUMP_GA(SCE%nloop)

      DO WHILE (SCE%ICALL .LT. SCE%MAXN) !Start main loop
   
         SCE%nloop = SCE%nloop + 1

         write (*,*) 'Starting loop #: ', SCE%nloop, 'with fitness: ', -SCE%XF(1)
          
         do igs = 1, SCE%ngs1 !  BEGIN LOOP ON COMPLEXES
               
            IF (GA%MPI) THEN !MASTER NODE
            
            !Pack up XNSTD, X, and XF to send
               I = 1
               DO IGENE = 1, GA%NGENE
                  SEND(I) = SCE%XNSTD(IGENE)
                  I = I + 1
               ENDDO
               
               DO IINDIV = 1, SCE%NPT 
                  DO IGENE = 1, GA%NGENE
                     SEND(I) = SCE%X(IINDIV,IGENE)
                     I = I + 1
                  ENDDO
                  SEND(I) = SCE%XF(IINDIV)
                  I = I + 1
               ENDDO
                              
               CALL MPI_SEND(SEND, SCE%NOPT+SCE%NPT*(SCE%NOPT+1), MPI_DOUBLE_PRECISION, &
               IGS, 1, MPI_COMM_WORLD, IERR)               
   
            ELSE
               CALL SCE_OPTIMIZATION(IGS)
            ENDIF

         ENDDO

         IF (GA%MPI) THEN

            do igs = 1, SCE%ngs1 

               RECV(:) = 0.
               CALL MPI_RECV(RECV, 1+SCE%NPT*(SCE%NOPT+1), MPI_DOUBLE_PRECISION, &
               IGS, 2, MPI_COMM_WORLD, ISTATUS, IERR)               

!Unpack X and XF
               XTEMP(:,:) = 0D0
               XFTEMP(:)  = 0D0
               
               SCE%ICALL = SCE%ICALL + INT(RECV(1))
               I = 2
               DO IINDIV = 1, SCE%NPT 
                  DO IGENE = 1, GA%NGENE
                     XTEMP(IINDIV,IGENE) = RECV(I)
                     I = I + 1
                  ENDDO
                  XFTEMP(IINDIV) = RECV(I)
                  I = I + 1
               ENDDO

               do k1 = 1, SCE%npg
                  k2 = (k1-1) * SCE%ngs1 + igs
                  do j = 1, SCE%NOPT
                     SCE%x(k2,j) = XTEMP(K2,J)
                  end do
                  SCE%xf(k2) = XFTEMP(K2)
               end do
                              
            ENDDO !igs = 1, SCE%ngs1 

         ENDIF !GA%MPI

         call SCE_sort(SCE%npt1,SCE%NOPT,SCE%x,SCE%xf,SCE%NPT,SCE%NOPT) !  RE-SORT THE POINTS

         CALL SCE_BEST_WORST !Record best and worst points

!  TEST THE POPULATION FOR PARAMETER CONVERGENCE
         call SCE_parstt(SCE%npt1,SCE%NOPT,SCE%x,SCE%xnstd,SCE%bound,SCE%gnrng,SCE%ipcnvg,SCE%NPT,SCE%NOPT)

!  CHECK FOR COMPLEX NUMBER REDUCTION
         if (SCE%ngs1 .gt. SCE%mings) then
            ngs2 = SCE%ngs1
            SCE%ngs1 = SCE%ngs1 - 1
            SCE%npt1 = SCE%ngs1 * SCE%npg
            call SCE_comp(SCE%NOPT,SCE%npt1,SCE%ngs1,ngs2,SCE%npg,SCE%x,SCE%xf,SCE%cx,SCE%cf,SCE%NPT,SCE%NOPT)
         end if
      
         CALL SCE_DUMP
         
         GA%FIT(:,:) = 0.
         GA%FIT(1,0) = -SCE%XF(1)

         GA%INDIV (:,:) = 0.
         GA%INDIV (1,1:GA%NGENE) = SCE%X(1,1:GA%NGENE)
      
         CALL DUMP_GA(SCE%nloop)

         INQUIRE(FILE='stop.stop',EXIST=LEX)
         IF (LEX) SCE%IPCNVG = 1 !This sends signal to shutdown

!  IF POPULATION IS CONVERGED INTO A SUFFICIENTLY SMALL SPACE
         if (SCE%ipcnvg .eq. 1) THEN
            CALL SCE_DUMP_CONVERGED
            
            IF (GA%MPI) THEN
               SEND(:) = -9D30
               WRITE(*,*) 'MASTER NODE SENDING SHUTDOWN SIGNAL'
               DO IRANK = 1, GA%NPROC-1
                  CALL MPI_SEND(SEND, SCE%NPT*(SCE%NOPT+1), MPI_DOUBLE_PRECISION, &
                  IRANK, 1, MPI_COMM_WORLD, IERR)
               ENDDO
               WRITE(*,*) 'MASTER NODE DONE SENDING SHUTDOWN SIGNAL'

            ENDIF
            MESSAGE='MASTER NODE GOING TO SHUTDOWN_GRACEFULLY'
            CALL SHUTDOWN_GRACEFULLY(MESSAGE)
         ENDIF

      ENDDO !(SCE%ICALL .LT. SCE%MAXN) 
   
   ELSE !(.NOT. GA%MPI .OR. (GA%MPI .AND. GA%IRANK .EQ. 0) THEN

!Slave nodes

      DONE = .FALSE.

      DO WHILE (.NOT. DONE)

         CALL MPI_RECV(RECV, SCE%NOPT+SCE%NPT*(SCE%NOPT+1), MPI_DOUBLE_PRECISION, &
         0, 1, MPI_COMM_WORLD, ISTATUS, IERR)
         
         IF (RECV(1) .LT. -1D30) THEN
            DONE = .TRUE. 
            WRITE(THREE,'(I3.3)') GA%IRANK
            MESSAGE='IRANK ' // THREE // ' CALLING SHUTDOWN_GRACEFULLY'         
            CALL SHUTDOWN_GRACEFULLY(MESSAGE)
         ELSE
!Unpack XNSTD, X, and XF
            I = 1
            DO IGENE = 1, GA%NGENE
               SCE%XNSTD(IGENE) = RECV(I)
               I = I + 1
            ENDDO

            DO IINDIV = 1, SCE%NPT 
               DO IGENE = 1, GA%NGENE
                  SCE%X(IINDIV,IGENE) = RECV(I)
                  I = I + 1
               ENDDO
               SCE%XF(IINDIV) = RECV(I)
               I = I + 1 
            ENDDO

            SCE%ICALL = 0
            CALL SCE_OPTIMIZATION(GA%IRANK)

         !Pack up X and XF to send
            SEND(:) = 0.
            SEND(1) = REAL(SCE%ICALL)
             
            I = 1

            DO IINDIV = 1, SCE%NPT 
               DO IGENE = 1, SCE%NOPT
                  I = I + 1
                  SEND(I) = SCE%X(IINDIV,IGENE)
               ENDDO
               I = I + 1
               SEND(I) = SCE%XF(IINDIV)
            ENDDO
                
            CALL MPI_SEND(SEND, 1+SCE%NPT*(SCE%NOPT+1), MPI_DOUBLE_PRECISION, &
            0, 2, MPI_COMM_WORLD, IERR)
    
         ENDIF !IF (SEND(1) .LT. -8D30) 

      ENDDO ! .NOT. DONE

   ENDIF 

   WRITE(THREE,'(I3.3)') GA%IRANK
   MESSAGE='IRANK ' // THREE // ' shutting down'         
   CALL SHUTDOWN_GRACEFULLY(MESSAGE)
   
ENDIF

END PROGRAM