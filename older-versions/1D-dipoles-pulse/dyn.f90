program dyn
use global_param
implicit none
external rkdumb,derivs
integer n !m=dimension of the Hamiltonian (Nq1*Nq2*st x Nq1*Nq2*st)
real(kind=dp) :: a,e0,k0 !sig=width of gaussian for vec0; soma=variable for sums
complex(kind=dp) :: soma
!character(len=9) :: fmt
real(kind=dp) :: npoints !number of time steps to take in integration
real(kind=dp) :: t0,tf !t0=initial time for integration; tf=final time for integration
real(kind=dp) :: tt,ch,x,x0,kf,w,expo,c0,c1,tstep
complex(kind=dp) ::vec0(st*Nq1),vec1(st*Nq1)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
call hamiltonian
n=st*Nq1 ! m=Nq1*Nq2*st %size of the Hamiltonian

t0=0.d0
tf=800.d0 !40 = Approx. 1 femtosecond
npoints=8000.d0 !If 1000, the step will be approx. 1 attosecond
tstep=(tf-t0)/npoints ! define time step size in atomic units

e0=6.0d-3
k0= sqrt(2.d0*mtotal*e0) ! k0 = initial momentum (In atomic units)
k0=0.d0
tt=0.d0
!kf=4.d0*mtotal*e0**2.d0 ! force constant of harmonic potential
w=2.d0*e0 !w=sqrt(kf/mredu)
 c0 = ((mtotal*w)/(pi))**(1.d0/4.d0)
 c1 = (4.d0/pi*(mtotal*w)**3.d0)**(1.d0/4.d0)
vec0=0.d0
do i=1,Nq1
  ch=(i-(Nq1-1.d0)/2.d0-1.d0)*stepX !change in q1 plus a displacement for it not to pass through ch=0, where the NAC is infinity.
  expo = dexp(-(ch+1.40d0)**2.d0*mtotal*w/2.d0) !* exp(im*k0*(ch))
  vec0(i) = c0 * expo
!  vec0(i) = (c0 * expo + c1 * ch * expo) / sqrt(2.d0)
!  vec0(i) = c0 * (2.0d0*mtotal*w*ch**2.d0-1.d0) * expo
!  vec0(i+Nq1) = 0.0d0
!  vec0(i+2*Nq1) = 0.0d0
end do
!----------------------------------------!
!normalizing                             !
soma=0.0d0                               !
do i=1,Nq1                               !
  soma=soma+dconjg(vec0(i))*vec0(i)*stepX!
end do                                   !
vec0=vec0/sqrt(soma)                     !
!----------------------------------------!
write(*,'(a24,2(f23.15))')'normalization constant =',soma
!--------------------------------------------------------!
!Checking initial energy                                 !
call derivs(tt,vec0,n,vec1,tstep) !evaluating y'(t=0,y)!
soma=0.0d0                                               !
do i=1,Nq1                                               !
  soma=soma + dconjg(vec0(i)) * vec1(i)*im *stepX        ! 
end do                                                   !
!--------------------------------------------------------!
write(*,'(a24,2(f23.15),a8)')'inital energy =',soma,' hartree'
write(*,'(a24,(f23.15))')'initial k * dq =',sqrt(2*mtotal*e0)*stepX
write(*,'(a24,(f23.15),a4)')'vibrational frequency =',w,' au.'
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! nanoseconds = 1d-9 seconds
! picoseconds = 1d-12 seconds
! femtoseconds = 1d-15 seconds
! attoseconds = 1d-18 seconds
! 1 time atomic units = 24.18884326505 attoseconds
call rkdumb(vec0,n,t0,tf,npoints,tstep,derivs)
write(*,*) '************************************************************'
write(*,*) 'Program finished'
write(*,*) '************************************************************'
end program
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rkdumb(y0,n,t0,tf,npoints,h,derivs)
use global_param
integer n
real(kind=dp) :: t0,tf,npoints,tt(npoints),h,t,coordinate(n),soma,Te,aux1,aux2
complex(kind=dp) :: dydt(n),y(n),y0(n),ymatrix(n,npoints+2),ymatlab(n,npoints)
real(kind=dp) :: st1n(npoints),st1p(npoints),st2n(npoints),st2p(npoints),st3n(npoints),st3p(npoints)
real(kind=dp) :: st1ncount,st1pcount,st2ncount,st2pcount,st3ncount,st3pcount
real(kind=dp) :: pop1,pop2,pop3,enert(npoints+1),normt(npoints+1)
external derivs,rk4 !
!FIND A WAY TO CORRECTLY STORE tt AND ymatrix
do1: do i=1,n
  y(i)=y0(i) !loading initial conditions
  ymatrix(i,2)=y(i) !saving initial amplitudes
end do do1
tt(1)=t0 !saving initial time value
t=t0
write(*,'(a24,f23.15,a4)')'time step =',h,' au.'
write(*,'(a24,f23.15)')'step size in q1 =',stepX

call derivs(t,y,n,dydt,h) ! evaluating y'(t,y)
!-----------------------------------------!
!Checking norm and conservation of energy !
soma=0.0d0                                !
Te=0.d0                                   !
do i=1,n                                  !
  soma=soma+dconjg(y(i))*y(i)*stepX       !
  Te=Te + dconjg(y(i)) * dydt(i)*im *stepX!
end do                                    !
!-----------------------------------------!
aux1=soma
aux2=Te
normt(1)=soma
enert(1)=Te
write(*,*) '************************************************************'
write(*,*) 'Integrating amplitudes over time'
write(*,*) '************************************************************'
do3: do k=1,npoints
!write(*,'(a14,f23.15)')'current time =',t
  call derivs(t,y,n,dydt,h) ! evaluating y'(t,y)
  call rk4(y,dydt,n,t,h,y,derivs) !evaluating y(t+h)

  do2: do i=1,n
    ymatrix(i,k+2)=y(i) !saving the amplitude for time t
  end do do2
  t=t+h
  tt(k)=t !saving time value for each step in time

!-----------------------------------------!
!Saving norm and total energy in time     !
soma=0.0d0                                !
Te=0.d0                                   !
do i=1,n                                  !
  soma=soma+dconjg(y(i))*y(i)*stepX       !
  Te=Te + dconjg(y(i)) * dydt(i)*im *stepX!
end do                                    !
!-----------------------------------------!
normt(k+1)=soma
enert(k+1)=Te

!-----------------------------------------------------!
! saving the population at the borders for each state !
st1n(k)=dconjg(y(1))*y(1)*stepX                       !
st1p(k)=dconjg(y(Nq1))*y(Nq1)*stepX                   !
st2n(k)=dconjg(y(Nq1+1))*y(Nq1+1)*stepX               !
st2p(k)=dconjg(y(Nq1+Nq1))*y(Nq1+Nq1)*stepX           !
st3n(k)=dconjg(y(2*Nq1+1))*y(2*Nq1+1)*stepX           !
st3p(k)=dconjg(y(n))*y(n)*stepX                       !
!-----------------------------------------------------!
end do do3
!---------------------------------------------------------------------------------------------------------------------------!
! Determining the highest population at the simulation                                                                      !
! must be always small, > 10d-8. If not, consider changing grid borders and check the potential energy behaviour at borders !
st1ncount=1.d-100; st1pcount=1.d-100; st2ncount=1.d-100; st2pcount=1.d-100; st3ncount=1.d-100; st3pcount=1.d-100;           !
do k=1,npoints                                                                                                              !
  if (st1n(k) > st1ncount) then                                                                                             !
    st1ncount=st1n(k)                                                                                                       !
  end if                                                                                                                    !
  if (st1p(k) > st1pcount) then                                                                                             !
    st1pcount=st1p(k)                                                                                                       !
  end if                                                                                                                    !
  if (st2n(k) > st2ncount) then                                                                                             !
    st2ncount=st2n(k)                                                                                                       !
  end if                                                                                                                    !
  if (st2p(k) > st2pcount) then                                                                                             !
    st2pcount=st2p(k)                                                                                                       !
  end if                                                                                                                    !
  if (st3n(k) > st3ncount) then                                                                                             !
    st3ncount=st3n(k)                                                                                                       !
  end if                                                                                                                    !
  if (st3p(k) > st3pcount) then                                                                                             !
    st3pcount=st3p(k)                                                                                                       !
  end if                                                                                                                    !
end do                                                                                                                      !
!---------------------------------------------------------------------------------------------------------------------------!
write(*,'(a50,(es11.2e3))')'maximum population at negative border - state 1 =',st1ncount
write(*,'(a50,(es11.2e3))')'maximum population at positive border - state 1 =',st1pcount
write(*,'(a50,(es11.2e3))')'maximum population at negative border - state 2 =',st2ncount
write(*,'(a50,(es11.2e3))')'maximum population at positive border - state 2 =',st2pcount
write(*,'(a50,(es11.2e3))')'maximum population at negative border - state 3 =',st3ncount
write(*,'(a50,(es11.2e3))')'maximum population at positive border - state 3 =',st3pcount
!-----------------------------------------!
!Checking norm and conservation of energy !
soma=0.0d0                                !
Te=0.d0                                   !
call derivs(t,y,n,dydt,h)                 !
do i=1,n                                  !
  soma=soma+dconjg(y(i))*y(i)*stepX       !
  Te=Te + dconjg(y(i)) * dydt(i)*im *stepX!
end do                                    !
!-----------------------------------------!
write(*,'(a24,(e11.2e3))')'Norm conservation =',aux1-soma
write(*,'(a24,(e11.2e3),a8)')'Energy conservation =',aux2-Te,' hartree'
write(*,'(a24,(f24.15),a8)')'Final total energy =',Te,' hartree'
write(*,*) '************************************************************'
write(*,*) 'Saving amplitude over time'
write(*,*) '************************************************************'
do i=1,Nq1
  ymatrix(i,1)=(i-(Nq1-1.d0)/2.d0-1.d0)*stepX
!  ymatrix(i,1)=i
  coordinate(i)=ymatrix(i,1)
end do
!open(unit=10,file='amp-matrix.txt',status='unknown')
open(unit=20,file='time-amp.txt',status='unknown')
open(unit=30,file='amp-real.txt',status='unknown')
open(unit=40,file='amp-imag.txt',status='unknown')
open(unit=50,file='field-time.txt',status='unknown')
!do i=1,n
!  write(10,'(<2*(npoints+2)>(e24.15))')ymatrix(i,:)
!end do
t=0.d0
do j=1,npoints+1
Et = (E00/freq)*(-(t-t00)/sig**2.d0*sin(freq*(t-t00)+phase)+freq*cos(freq*(t-t00)+phase))*dexp(-(t-t00)**2.d0/(2.d0*sig**2.d0))
t=t+h
write(50,'(es26.16e3)') Et
  do i=1,Nq1
  write(20,'(4(es26.16e3))') coordinate(i),real(dconjg(ymatrix(i,j+1))*ymatrix(i,j+1)*stepX),&
real(dconjg(ymatrix(Nq1+i,j+1))*ymatrix(Nq1+i,j+1)*stepX),&
real(dconjg(ymatrix(2*Nq1+i,j+1)) * ymatrix(2*Nq1+i,j+1) *stepX)
    write(30,'(4(es26.16e3))')coordinate(i),(real(ymatrix(i,j+1))),(real(ymatrix(Nq1+i,j+1))),(real(ymatrix(2*Nq1+i,j+1)))
    write(40,'(4(es26.16e3))')coordinate(i),(aimag(ymatrix(i,j+1))),(aimag(ymatrix(Nq1+i,j+1))),(aimag(ymatrix(2*Nq1+i,j+1)))
  end do
  write(20,*)''
  write(20,*)''
  write(30,*)''
  write(30,*)''
  write(40,*)''
  write(40,*)''
end do
 close(unit=10)
 close(unit=20)
 close(unit=30)
 close(unit=40)
 close(unit=50)

open(unit=60,file='data.txt',status='unknown')
t=0.d0
do i=1,800/4+1 ! tf / step(4 = 0.1 femtoseconds)
  t=(i-1.d0)*4.d0
  Et = (E00/freq)*(-(t-t00)/sig**2.d0*sin(freq*(t-t00)+phase)+freq*cos(freq*(t-t00)+phase))*dexp(-(t-t00)**2.d0/(2.d0*sig**2.d0))
  pop1=0.d0
  pop2=0.d0
  pop3=0.d0
  do j=1,Nq1
    pop1=pop1 + dconjg(ymatrix(j,t*10+1)) * ymatrix(j,t*10+1) *stepX
    pop2=pop2 + dconjg(ymatrix(Nq1+j,t*10+1)) * ymatrix(Nq1+j,t*10+1) *stepX
    pop3=pop3 + dconjg(ymatrix(2*Nq1+j,t*10+1)) * ymatrix(2*Nq1+j,t*10+1) *stepX
  end do
 write(60,'(f5.1,6(es20.10e3))') t,Et,pop1,normt(t*10+1),enert(t*10+1),pop2,pop3
end do
 close(unit=60)

end subroutine rkdumb
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rk4(y0,dydt,n,t,h,Yout,derivs)
use global_param
implicit none
integer n
real(kind=dp) :: t,h,th
complex(kind=dp) :: dydt(n),y0(n),Yout(n),dym(n),dyt(n),yt(n)
external derivs 
! 4th order runge-kutta
!h = (tfinal - t0) / ntemp
!t+1= t + h
!k1 = h * y'(t , f(t))
!k2 = h *  y'(t + h/2.0d0 , y(t) + k1/2.0d0)
!k3 = h *  y'(t + h/2.0d0 , y(t) + k2/2.0d0)
!k4 = h *  y'(t + h , y(t) + k3)
!y(t+1) = y(t) + 1/6 (k1+ 2*k2 + 2*k3 + k4)
th=t+h/2.d0
!dydt = y'(t , y(t)) = k1/h ---- dydt MUST NOT BE PREVIOUSLY MULTIPLIED BY h ----
do1: do i=1,n
  yt(i) = y0(i) + (h/2.d0) * dydt(i) !Evaluation of y(t) + k1/2.0d0
end do do1
call derivs(th,yt,n,dyt,h/2.d0) !dyt = k2/h ---- Evaluation of y'(t + h/2.0d0 , y(t) + k1/2.0d0)
do2: do i=1,n
  yt(i) = y0(i) + (h/2.d0) * dyt(i)  !Evaluation of y(t) + k2/2.0d0
end do do2
call derivs(th,yt,n,dym,h/2.00) !dym = k3/h ---- Evaluation of y'(t + h/2.0d0 , y(t) + k2/2.0d0)
do3: do i=1,n
  yt(i) = y0(i) + h * dym(i) !Evaluation of y(t) + k3
  dym(i) = dyt(i) + dym(i)  ! making dym = (k2 + k3)/h ---- variable dyt free now
end do do3
call derivs(t+h,yt,n,dyt,h*2.d0) !dyt = k4/h ---- Evaluation of y'(t + h/2.0d0 , y(t) + k1/2.0d0)
do4: do i=1,n
  Yout(i) = y0(i) + h/6.d0 * (dydt(i) + 2.d0 * dym(i) + dyt(i))
end do do4
end subroutine rk4
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Subroutine computes the derivatives dPsi(t)/dt = H/i Psi
!This derivative routine calculates 
!Written as a simple matrix multiplication between the Hamiltonian matrix H and the amplitudes vector Psi
subroutine derivs(t,y,n,dydt,tstep)
use global_param
use omp_lib
implicit none
integer n,ii
real(kind=dp) :: t,tstep
complex(kind=dp) :: y(n),dydt(n)

Et = (E00/freq)*(-(t-t00)/sig**2.d0*sin(freq*(t-t00)+phase)+freq*cos(freq*(t-t00)+phase))*dexp(-(t-t00)**2.d0/(2.d0*sig**2.d0))

ham=Ha

do i=1,Nq1
  ham(i,i)=Ha(i,i)-(pdm1x(i)*Et)
  ham(Nq1+i,Nq1+i)=Ha(Nq1+i,Nq1+i)-(pdm2x(i)*Et)
  ham(2*Nq1+i,2*Nq1+i)=Ha(2*Nq1+i,2*Nq1+i)-(pdm3x(i)*Et)

  ham(i,Nq1+i)=-tdm21x(i) * Et
  ham(Nq1+i,i)=-tdm21x(i) * Et
  ham(i,2*Nq1+i)=-tdm31x(i) * Et
  ham(2*Nq1+i,i)=-tdm31x(i) * Et
  ham(Nq1+i,2*Nq1+i)=-tdm32x(i) * Et
  ham(2*Nq1+i,Nq1+i)=-tdm32x(i) * Et
end do

!$OMP parallel default(private) shared (ham,dydt,y)

!$OMP do
do i=1,n
  dydt(i)=dcmplx(0.d0,0.d0)
  do j=1,n
!    dydt(i) = dydt(i) +  (Ha(i,j)-(pdm1(i)*Et)) * y(j)/im
!    dydt(Nq1+i) = dydt(Nq1+i) +  (Ha(Nq1+i,Nq1+j)-(pdm2(i)*Et)) * y(Nq1+j)/im
!    dydt(2*Nq1+i) = dydt(2*Nq1+i) +  (Ha(2*Nq1+i,2*Nq1+j)-(pdm3(i)*Et)) * y(2*Nq1+j)/im

!    dydt(i) = dot_product( Ha(i,:),y(:) ) / im
!    dydt(i) = dydt(i) +  ham(i,j) * y(j)/im
!    dydt(Nq1+i) = dydt(Nq1+i) +  (ham(Nq1+i,Nq1+j)) * y(Nq1+j)/im
!    dydt(2*Nq1+i) = dydt(2*Nq1+i) +  (ham(2*Nq1+i,2*Nq1+j)) * y(2*Nq1+j)/im
    dydt(i) = dydt(i) +  ham(i,j) * y(j)/im
  end do
end do
!$OMP end do

!$OMP end parallel

end subroutine derivs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine momentum(mom)
use global_param
implicit none
real(kind=dp) :: mom(Nq1,Nq1),const

const = 1.d0/(2.d0*stepX)

mom(:,:)=0.0d0

mom(1,2) = const
mom(Nq1,Nq1-1) = - const

do i=2,-1 !from the third line to the antepenulmate one
  mom(i,i-1)=-const
  mom(i,i+1)=const
end do

end subroutine momentum


