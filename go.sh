#!/bin/bash
clear
rm *.dat
gfortran -Ofast *.f*
./a.out
