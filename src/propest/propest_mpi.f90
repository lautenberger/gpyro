MODULE GA_MPI

USE GPYRO_VARS
USE GA_VARS
USE MPI

IMPLICIT NONE

CONTAINS

! *****************************************************************************
SUBROUTINE SEND_POPULATION
! *****************************************************************************

INTEGER :: DEST, IINDIV, IERR

DO DEST = 1, GA%NPROC - 1
   DO IINDIV = 1,GA%NINDIV
      CALL MPI_SEND(GA%INDIV(IINDIV,1:GA%NGENE), GA%NGENE, &
      MPI_DOUBLE_PRECISION, DEST, IINDIV, MPI_COMM_WORLD, IERR)
   ENDDO
ENDDO
      
END SUBROUTINE SEND_POPULATION
! *****************************************************************************


! *****************************************************************************
SUBROUTINE RECV_POPULATION
! *****************************************************************************

INTEGER :: IINDIV, IERR
INTEGER :: ISTATUS(MPI_STATUS_SIZE)

DO IINDIV = 1,GA%NINDIV
   CALL MPI_RECV(GA%INDIV(IINDIV,1:GA%NGENE), GA%NGENE, &
   MPI_DOUBLE_PRECISION, 0, IINDIV, MPI_COMM_WORLD, ISTATUS, IERR)
ENDDO

END SUBROUTINE RECV_POPULATION
! *****************************************************************************

! *****************************************************************************
SUBROUTINE SEND_FITNESS_ANYNODE(IDEST,IINDIV)
! *****************************************************************************

INTEGER, INTENT(IN) :: IDEST, IINDIV
INTEGER :: IERR,IPHI
REAL(EB) :: SEND(1:GA%MPIFITSIZE)

SEND(1) = REAL(IDEST,EB)
SEND(2) = REAL(IINDIV,EB)
SEND(3) = GA%FIT(IINDIV,0) !overall

DO IPHI = 1, GA%NPHI
      SEND(3+IPHI) = GA%FIT(IINDIV,IPHI)
ENDDO

CALL MPI_SEND(SEND, GA%MPIFITSIZE, MPI_DOUBLE_PRECISION,0,IINDIV, MPI_COMM_WORLD, IERR)

END SUBROUTINE SEND_FITNESS_ANYNODE
! *****************************************************************************

! *****************************************************************************
SUBROUTINE RECV_FITNESS_ANYNODE(IDEST,IINDIV)
! *****************************************************************************

INTEGER, INTENT(OUT) :: IDEST, IINDIV
INTEGER :: ISTATUS(MPI_STATUS_SIZE), IERR,IPHI
REAL(EB) :: RECV(1:GA%MPIFITSIZE)

CALL MPI_RECV(RECV, GA%MPIFITSIZE, MPI_DOUBLE_PRECISION, &
     MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, ISTATUS, IERR)

IDEST            = INT(RECV(1))
IINDIV           = INT(RECV(2))
GA%FIT(IINDIV,0) = RECV(3)

DO IPHI = 1, GA%NPHI
   GA%FIT(IINDIV,IPHI) = RECV(3+IPHI)
ENDDO

END SUBROUTINE RECV_FITNESS_ANYNODE
! *****************************************************************************

! *****************************************************************************
SUBROUTINE SEND_IINDIV_TO_RUN(IDEST,IINDIV)
! *****************************************************************************

INTEGER, INTENT(IN) :: IDEST,IINDIV
INTEGER :: IERR

CALL MPI_SEND(IINDIV,1,MPI_INTEGER,IDEST,IINDIV,MPI_COMM_WORLD,IERR)

END SUBROUTINE SEND_IINDIV_TO_RUN
! *****************************************************************************

! *****************************************************************************
SUBROUTINE RECV_IINDIV_TO_RUN(IINDIV)
! *****************************************************************************

INTEGER, INTENT(OUT) :: IINDIV
INTEGER :: ISTATUS(MPI_STATUS_SIZE), IERR

CALL MPI_RECV(IINDIV,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,ISTATUS,IERR)

END SUBROUTINE RECV_IINDIV_TO_RUN
! *****************************************************************************
END MODULE GA_MPI
! *****************************************************************************