#ifort -qopenmp -c -fpic -mcmodel=large global_param.f90;ifort -qopenmp -c -fpic -mcmodel=large dyn.f90; #ifort -qopenmp -fpic -mcmodel=large global_param.o dyn.o -o dyn.exe -check bounds;
#export OMP_NUM_THREADS=4
#time ./dyn.exe
#h5pfc -c -assume bscc -qopenmp -mkl -fpic -mcmodel=large global_param.f90; h5pfc -c -assume bscc -qopenmp -mkl -fpic -mcmodel=large dyn.f90; h5pfc -assume bscc -qopenmp -mkl -fpic -mcmodel=large dyn.o global_param.o -o rundyn.exe
time gnuplot pl2.gpl &
time gnuplot pl3.gpl &
time gnuplot pl1.gpl
#sleep 150
#cp pop* figs/
#cp real* figs/
#cp imag* figs/
#rm pop*
#rm real*
#rm imag*
#cd figs
#ffmpeg -i pop%06d.png psi-square.mp4 &
#ffmpeg -i real%06d.png psi-real.mp4 &
#ffmpeg -i imag%06d.png psi-imag.mp4 
##ffmpeg -framerate 50 -i imag%06d.png psi-imag.mp4 
#sleep 15
#cd ..
#cp figs/psi-square.mp4 .
#cp figs/psi-real.mp4 .
#cp figs/psi-imag.mp4 .
#rm figs/*

