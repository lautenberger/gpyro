@echo off

cd gpyro_propest
del *.obj *.o *.mod *.exe
C:\MinGW\bin\mingw32-make.exe -f ../../Makefile_gpyro_propest gnu_win_mpi
del *.obj *.o *.mod
cd ..

cd gpyro
del *.obj *.o *.mod *.exe
C:\MinGW\bin\mingw32-make.exe -f ../../Makefile_gpyro gnu_win
del *.obj *.o *.mod
cd ..

cd gpyro_fds
del *.obj *.o *.mod *.exe
C:\MinGW\bin\mingw32-make.exe -f ../../Makefile_gpyro_fds gnu_win_mpi
del *.obj *.o *.mod
cd ..
