#!/bin/bash

#gfortran -g -O3 -fopenmp  -o xsect.exe \
gfortran -g -O3 -o xsect.exe \
types_module.f90 binomial_module.f90 \
abrasioncrs1.f90 abrasioncrs2.f90  density.f90  lgwt.f90 main.f90 \
-Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,invalid \
-ffree-line-length-none -finit-local-zero -fbounds-check -fcheck-array-temporaries



