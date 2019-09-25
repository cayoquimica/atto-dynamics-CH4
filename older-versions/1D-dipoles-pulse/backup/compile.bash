ifort -c global_param.f90
ifort -c dyn.f90
ifort global_param.o dyn.o -o dyn.exe
./dyn.exe
