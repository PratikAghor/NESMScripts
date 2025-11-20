#!/bin/bash

# Compiler and flags
FC=gfortran
FCFFLAGS="-fPIC -m64 -fno-second-underscore -fconvert=big-endian -O"

# Source files and libraries
SRC="nrl2nc.f zh_hp.f"
# on oncean8
# NCINCLUDE="-I/usr/include"
# NCLIBS="-L/lib/x86_64-linux-gnu/ -lnetcdf -lnetcdff"

# on PACE
# nf-config --all gives path
NCINCLUDE="-I/usr/local/pace-apps/spack/packages/linux-rhel7-x86_64/gcc-10.3.0/netcdf-fortran-4.5.4-yx5osuxluenmuvr3xnahmosfr3abeu2p/include"
NCLIBS="-L/usr/local/pace-apps/spack/packages/linux-rhel7-x86_64/gcc-10.3.0/netcdf-fortran-4.5.4-yx5osuxluenmuvr3xnahmosfr3abeu2p/lib -lnetcdff  -lnetcdf -lnetcdf -ldl -lm"
# Compile command
$FC $FCFFLAGS -o nrl2nc $SRC $NCINCLUDE $NCLIBS

