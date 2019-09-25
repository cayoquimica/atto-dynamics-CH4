gfortran-8 -c global_param.f90
gfortran-8 -c dyn.f90
gfortran-8 global_param.o dyn.o -o g-dyn.exe
./g-dyn.exe
