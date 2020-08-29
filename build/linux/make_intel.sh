#!/bin/bash

# Set gpyro version:
GPYRO_VER=0.8200

# Gpyro uses several environment variables for compilation. Since the default values
# specified on lines 17 - 27 below probably do not work on your system, add lines
# similar to the following to your ~/.bashrc file:
#
# export GPYRO_FCOMPL_SERIAL_INTEL=ifort
# export GPYRO_FCOMPL_MPI_INTEL=/usr/local/openmpi-1.10.3_intel/bin/mpifort
# export GPYRO_FCOMPL_SERIAL_GNU=gfortran
# export GPYRO_FCOMPL_MPI_GNU=/usr/local/openmpi-1.10.3_gnu/bin/mpifort
# export GPYRO_BINARY_DIRECTORY=/usr/local/bin

# Check if environment variables are set and export if not: 
if [ -z "$GPYRO_FCOMPL_SERIAL_INTEL" ]; then
   export GPYRO_FCOMPL_SERIAL_INTEL=ifort
fi

if [ -z "$GPYRO_FCOMPL_MPI_INTEL" ]; then
   export GPYRO_FCOMPL_MPI_INTEL=/usr/local/openmpi-2.1.0_intel/bin/mpifort
fi

if [ -z "$GPYRO_BINARY_DIRECTORY" ]; then
   GPYRO_BINARY_DIRECTORY=/usr/local/bin
fi

cd gpyro

rm -f *.o *.mod gpyro
make -f ../../Makefile_gpyro intel_linux
sudo cp -f gpyro $GPYRO_BINARY_DIRECTORY/gpyro_${GPYRO_VER}_intel
rm -f *.o *.mod

rm -f *.o *.mod gpyro
make -f ../../Makefile_gpyro intel_linux_openmp
sudo cp -f gpyro_openmp $GPYRO_BINARY_DIRECTORY/gpyro_openmp_${GPYRO_VER}_intel
rm -f *.o *.mod

rm -f *.o *.mod gpyro
make -f ../../Makefile_gpyro intel_linux_debug
sudo cp -f gpyro_debug $GPYRO_BINARY_DIRECTORY/gpyro_${GPYRO_VER}_debug_intel
rm -f *.o *.mod

cd ../gpyro_propest

rm -f *.o *.mod gpyro_propest
make -f ../../Makefile_gpyro_propest intel_linux_mpi
sudo cp -f gpyro_propest $GPYRO_BINARY_DIRECTORY/gpyro_propest_${GPYRO_VER}_intel
rm -f *.o *.mod

rm -f *.o *.mod gpyro_propest_debug
make -f ../../Makefile_gpyro_propest intel_linux_mpi_debug
sudo cp -f gpyro_propest $GPYRO_BINARY_DIRECTORY/gpyro_propest_${GPYRO_VER}_debug_intel
rm -f *.o *.mod

cd ../gpyro_fds

rm -f *.o *.mod
make -f ../../Makefile_gpyro_fds intel_linux_mpi
sudo cp -f fds6_mpi_gpyro $GPYRO_BINARY_DIRECTORY/fds6_mpi_22343_gpyro_${GPYRO_VER}_intel
rm -f *.o *.mod

rm -f *.o *.mod
make -f ../../Makefile_gpyro_fds intel_linux_openmp_mpi
sudo cp -f fds6_openmp_mpi_gpyro $GPYRO_BINARY_DIRECTORY/fds6_openmp_mpi_22343_gpyro_${GPYRO_VER}_intel
rm -f *.o *.mod

rm -f *.o *.mod
make -f ../../Makefile_gpyro_fds intel_linux_openmp_mpi_debug
sudo cp -f fds6_openmp_mpi_gpyro_debug $GPYRO_BINARY_DIRECTORY/fds6_openmp_mpi_22343_gpyro_${GPYRO_VER}_debug_intel
rm -f *.o *.mod

rm -f *.o *.mod
make -f ../../Makefile_gpyro_fds intel_linux_mpi_debug
sudo cp -f fds6_mpi_gpyro_debug $GPYRO_BINARY_DIRECTORY/fds6_mpi_22343_gpyro_${GPYRO_VER}_debug_intel
rm -f *.o *.mod

rm -f *.o *.mod fds6_gpyro
make -f ../../Makefile_gpyro_fds intel_linux
sudo cp -f fds6_gpyro $GPYRO_BINARY_DIRECTORY/fds6_22343_gpyro_${GPYRO_VER}_intel
rm -f *.o *.mod

exit 0
