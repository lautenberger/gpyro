REM Set MPIRUN to full path to mpirun for MPI installation. 
SET MPIEXEC="C:\Program Files\Microsoft MPI\Bin\mpiexec"

REM Set GPYROBINDIR to directory containing gpyro binaries:
SET GPYROBINDIR=..\..\build\win

echo "Sample 01"
cd 01-dtg
%GPYROBINDIR%\gpyro sample_01.data > log.txt
cd ..

echo "Sample 02"
cd 02-thermoplastic
%GPYROBINDIR%\gpyro sample_02.data > log.txt
cd ..

echo "Sample 03"
cd 03-char
%GPYROBINDIR%\gpyro sample_03.data > log.txt
cd ..

echo "Sample 04"
cd 04-intumescent
%GPYROBINDIR%\gpyro sample_04.data > log.txt
cd ..

echo "Sample 05"
cd 05-smolder
%GPYROBINDIR%\gpyro sample_05.data > log.txt
cd ..

echo "Sample 06"
cd 06-oxidative_pyrolysis
%GPYROBINDIR%\gpyro sample_06.data > log.txt
cd ..

echo "Sample 07"
cd 07-2d
%GPYROBINDIR%\gpyro sample_07.data > log.txt
cd ..

echo "Sample 08"
cd 08-3d
%GPYROBINDIR%\gpyro ./sample_08.data > log.txt
cd ..

REM Uncomment EXIT below to prevent MPI samples from running:
REM EXIT

REM Below you will likely have to add a "-machinefile hosts" and then specify the
REM hosts to run in that file. 

echo "Sample 09"
cd 09-dtg_ga
%MPIEXEC% -n 9 %GPYROBINDIR%\gpyro_propest > log.txt
cd ..

echo "Sample 10"
cd 10-char_ga
%MPIEXEC% -n 17 %GPYROBINDIR%\gpyro_propest > log.txt
cd ..

echo "Sample 11"
cd 11-fds6_gpyro_coupled
%MPIEXEC% -np 16 %GPYROBINDIR%\fds6_gpyro crib.fds > log.txt
cd ..
