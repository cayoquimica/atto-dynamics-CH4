# atto-dynamics-CH4
This code is written to do the attoseconds dynamics of the Jahn-Teller effect in the sudden ionization of a CH4 molecule in a 2D scheme.

It has 3 fortran .f90 files and 15 .txt files. 

"dyn.f90" is the main code file which sets all initial variables and calls the time evolution of the wave-packet. It has the subroutine 'rkdumb', which treats each result of one time interaction (it has a commented part which enables to save the output files as .h5. The subroutines that does it - 'save_vector_h5' and 'save_matrix_h5' - are commented); 'rk4', which is the 4th order runge kutta;  'Ha_calc', which evaluates the operation H|Psi>; 'momentum_calc_q1' and 'momentum_calc_q2', which calculates (d/dq1)|Psi> and (d/dq2)|Psi> respectively - fisrt spatial derivatives, along each of our 2 dimensions (q1 and q2); and finally 'angular_momentum', which evaluates the operation L|Psi>. it also requires the file "csr_vectors", that is created executing the program "load_hamiltonian.f90".

"global_param.f90" loads all global variables and has the subroutine 'load_data', which loads all the data from the .txt files.

"load_hamiltonian.f90" is a main separated program that creates the Hamiltonian matrix and stores it into a set of vectors, following the CSR procedure for sparse matrix. These vectors are saved in the file "csr_vectors", so that this program needs to be executed only once if there is no modifications on the .txt files or in the "load_hamiltonian.f90" file itself. Once the file "csr_vectors" is created, you can use "dyn.f90"
