#ifort -qopenmp -c global_param.f90;ifort -qopenmp -c dyn.f90; ifort global_param.o dyn.o -o dyn.exe; 
#export OMP_NUM_THREADS=4
#time ./dyn.exe
ifort -c global_param.f90;ifort -c dyn.f90; ifort global_param.o dyn.o -o dyn.exe -check bounds; time ./dyn.exe
gnuplot pl2.gpl &
gnuplot pl3.gpl &
gnuplot pl1.gpl
sleep 10
cp pop* figs/
cp real* figs/
cp imag* figs/
rm pop*
rm real*
rm imag*
cd figs
ffmpeg -i pop%06d.png psi-square.mp4 &
ffmpeg -i real%06d.png psi-real.mp4 &
ffmpeg -i imag%06d.png psi-imag.mp4 
#ffmpeg -framerate 50 -i imag%06d.png psi-imag.mp4 
sleep 5
cd ..
cp figs/psi-square.mp4 .
cp figs/psi-real.mp4 .
cp figs/psi-imag.mp4 .
rm figs/*.*

