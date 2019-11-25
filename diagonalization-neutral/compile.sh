rm *.exe
ifort -c -mkl global_param.f90 ; ifort -c -mkl diag.f90 ; ifort -mkl diag.o global_param.o -o exe.exe;
