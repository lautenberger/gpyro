@echo off

set intelbin="%IFORT_COMPILER19%\bin"

IF "%SETUP_IFORT_COMPILER64%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER64=1

echo Setting up compiler environment
call %intelbin%\ifortvars intel64

:envexist

del *.exe

cd gpyro
del *.obj *.o *.mod *.exe
C:\MinGW\bin\mingw32-make.exe -f ../../Makefile_gpyro intel_win
move *.exe ..
del *.obj *.o *.mod

C:\MinGW\bin\mingw32-make.exe -f ../../Makefile_gpyro intel_win_openmp
move *.exe ..
del *.obj *.o *.mod

cd ..\gpyro_fds
del *.obj *.o *.mod *.exe
C:\MinGW\bin\mingw32-make.exe -f ../../Makefile_gpyro_fds intel_win_mpi
copy *.exe ..
del *.obj *.o *.mod

cd ..\gpyro_propest
del *.obj *.o *.mod *.exe
C:\MinGW\bin\mingw32-make.exe -f ../../Makefile_gpyro_propest intel_win_mpi
copy *.exe ..
del *.obj *.o *.mod
cd ..
