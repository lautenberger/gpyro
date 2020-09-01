#!/bin/bash

# Set gpyro version:
GPYRO_VER=0.8200

# Gpyro uses several environment variables for compilation. Since the default values
# specified on lines 17 - 27 below probably do not work on your system, add lines
# similar to the following to your ~/.bashrc file:
#
# export GPYRO_FCOMPL_SERIAL_INTEL=ifort
# export GPYRO_FCOMPL_MPI_INTEL=/usr/local/openmpi-4.0.3_intel/bin/mpifort
# export GPYRO_FCOMPL_SERIAL_GNU=gfortran
# export GPYRO_FCOMPL_MPI_GNU=/usr/local/openmpi-4.0.3_gnu/bin/mpifort
# export GPYRO_BINARY_DIRECTORY=/usr/local/bin

# Check if environment variables are set and export if not:
if [ -z "$GPYRO_FCOMPL_SERIAL_GNU" ]; then
   export GPYRO_FCOMPL_SERIAL_GNU=gfortran
fi

if [ -z "$GPYRO_FCOMPL_MPI_GNU" ]; then
   export GPYRO_FCOMPL_MPI_GNU=/usr/local/openmpi-2.1.0_gnu/bin/mpifort
fi

if [ -z "$GPYRO_BINARY_DIRECTORY" ]; then
   GPYRO_BINARY_DIRECTORY=/usr/local/bin
fi

mkdir gpyro 2> /dev/null
cd gpyro

rm -f *.o *.mod gpyro_debug
make -f ../../Makefile_gpyro gnu_linux_debug
sudo cp -f gpyro_debug /usr/local/bin/gpyro_${GPYRO_VER}_debug_gnu
rm -f *.o *.mod

rm -f *.o *.mod gpyro
make -f ../../Makefile_gpyro gnu_linux
sudo cp -f gpyro /usr/local/bin/gpyro_${GPYRO_VER}_gnu
rm -f *.o *.mod

mkdir ../gpyro_propest 2> /dev/null
cd ../gpyro_propest

rm -f *.o *.mod gpyro_propest_debug
make -f ../../Makefile_gpyro_propest gnu_linux_mpi_debug
sudo cp -f gpyro_propest_debug /usr/local/bin/gpyro_propest_${GPYRO_VER}_debug_gnu
rm -f *.o *.mod

rm -f *.o *.mod gpyro_propest
make -f ../../Makefile_gpyro_propest gnu_linux_mpi
sudo cp -f gpyro_propest /usr/local/bin/gpyro_propest_${GPYRO_VER}_gnu
rm -f *.o *.mod

mkdir ../gpyro_fds 2> /dev/null
cd ../gpyro_fds

rm -f *.o *.mod
make -f ../../Makefile_gpyro_fds gnu_linux_mpi_debug
sudo cp -f fds6_mpi_gpyro_debug /usr/local/bin/fds6_mpi_22343_gpyro_${GPYRO_VER}_debug_gnu
rm -f *.o *.mod

rm -f *.o *.mod fds6_gpyro
make -f ../../Makefile_gpyro_fds gnu_linux
sudo cp -f fds6_gpyro /usr/local/bin/fds6_22343_gpyro_${GPYRO_VER}_gnu

rm -f *.o *.mod
make -f ../../Makefile_gpyro_fds gnu_linux_mpi
sudo cp -f fds6_mpi_gpyro /usr/local/bin/fds6_mpi_22343_gpyro_${GPYRO_VER}_gnu
rm -f *.o *.mod

exit 0
