ifort -c global_param.f90;ifort -c dyn.f90; ifort global_param.o dyn.o -o dyn.exe;./dyn.exe
time gnuplot pl.gpl
rm figs/file*
cp file* figs/
cp real* figs/
cp imag* figs/
rm file*
rm real*
rm imag*
cd figs
ffmpeg -framerate 100 -i file%05d.png video.avi
ffmpeg -framerate 100 -i real%05d.png psi-real.avi
ffmpeg -framerate 100 -i imag%05d.png psi-imag.avi
cd ..
cp figs/video.avi .
cp figs/psi-real.avi .
cp figs/psi-imag.avi .
rm figs/video.avi
rm figs/psi-real.avi
rm figs/psi-imag.avi
