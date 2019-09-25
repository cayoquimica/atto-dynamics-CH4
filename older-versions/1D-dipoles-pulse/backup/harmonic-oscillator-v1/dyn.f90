program dyn
use global_param
implicit none
external rkdumb,rk4,derivs !!!!! VER COMO CHAMAR !!!!!
integer, parameter :: dp = kind(1.d0) 
integer m !m=dimension of the Hamiltonian (Nq1*Nq2*st x Nq1*Nq2*st)
real(kind=dp) :: a,e0,k0,pi !sig=width of gaussian for vec0; soma=variable for sums
complex(kind=dp) :: soma,im
!character(len=9) :: fmt
integer npoints !number of time steps to take in integration
real(kind=dp) :: t0,tf !t0=initial time for integration; tf=final time for integration
real(kind=dp) :: tt,teta,phi,const,ch,x,x0,stepX,kf,w,expo,c0,c1,c2
complex(kind=dp) ::term1,term2,vec0(Nq1),vec1(Nq1),vec2(Nq1),vec3(Nq1)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im=cmplx(0.d0,1.d0)
pi = 3.141592653589793d0
!e0=0.00534042d0 ! E(t=0) = initial kinetic energy (In atomic units) --- Zero point energy
e0=0.005d0
!k0= sqrt(2.d0*mredu*e0) ! k0 = initial momentum (In atomic units)
k0=0.d0
x0=2.09970623d0 ! coordinate X of the minimum energy configuration in atomic units
stepX=0.01d0 !step in X
tt=0.d0
kf=4.d0*mredu*e0**2.d0 ! force constant of harmonic potential
!w=sqrt(kf/mredu)
w=2.d0*e0
 c0 = ((mredu*w)/pi)**(1.d0/4.d0)
 c1 = (4.d0/pi*(mredu*w)**3.d0)**(1.d0/4.d0)
 c2 = ((mredu*w)/(4.d0*pi))**(1.d0/4.d0)
!open(unit=3,file='n2harm.dat',status='unknown')
open(unit=4,file='n2.data',status='unknown')
do i=1,Nq1
  ch=(i-(Nq1-1.d0)/2.d0-1.d0)*stepX !change in x
  x=x0+ch ! x in atomic units
!write(3,'(3(f24.15))')x,0.5d0*kf*(x-x0)**2.0d0,e0
write(4,'(f24.15)') 0.5d0*kf*(x-x0)**2.0d0
  expo = dexp(-(x-x0)**2.d0*mredu*w/2.d0) !* exp(im*k0*(x-x0))
!  vec0(i) = c0 * expo
  vec0(i) = c0 * expo
  vec1(i) = c1 * (x-x0) * expo
!  vec3(i) = c2 * (2.d0*mredu*w*(x-x0)**2.d0 - 1.0d0) * expo
!  vec0(i) = (vec2(i) + vec1(i))/sqrt(2.d0)
!  vec0(i) = (vec2(i) + vec1(i) + vec3(i))/sqrt(3.d0)
end do
! close(unit=3)
 close(unit=4)
!vec0(:)=(vec1(:)+vec0(:))/sqrt(2.d0)
!vec0(:)=vec1(:)

!----------------------------------------!
!normalizing                             !
soma=0.0d0                               !
do i=1,Nq1                               !
  soma=soma+conjg(vec0(i))*vec0(i)*stepX !
end do                                   !
vec0=vec0/sqrt(soma)                     !
!----------------------------------------!

write(*,'(a24,2(f23.15))')'normalization constant =',soma



call derivs(tt,vec0,Nq1,vec1) ! evaluating y'(t,y)
!------------------------------------------------!
!Checking initial energy                         !
soma=0.0d0                                       !
do i=1,Nq1                                       !
  soma=soma + conjg(vec0(i)) * vec1(i)*im *stepX !
end do                                           !
!------------------------------------------------!
write(*,'(a24,2(f23.15))')'inital energy =',soma
write(*,'(a24,2(f23.15))')'initial velocity =',k0/mredu
!write(*,'(a24,2(f23.15))')'reduced mass =',mredu
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=Nq1 ! m=Nq1*Nq2*st
! nanoseconds = 1d-9 seconds
! picoseconds = 1d-12 seconds
! femtoseconds = 1d-15 seconds
! attoseconds = 1d-18 seconds
! 1 time atomic units = 24.18884326505 attoseconds
t0=0.d0
tf=1000.d0 !40 = Approx. 1 femtosecond
npoints=1000.d0 !If 1000, the step will be approx. 1 attosecond
call rkdumb(vec0,m,t0,tf,npoints,derivs,stepX)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end program
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rkdumb(y0,n,t0,tf,npoints,derivs,stepX)
integer, parameter :: dp = kind(1.d0) 
integer i,k,n
real(kind=dp) :: t0,tf,tt(npoints),h,t,coordinate(n),stepX
complex(kind=dp) :: soma,dydt(n),y(n),y0(n),ymatrix(n,npoints+2),ymatlab(n,npoints),Te,im,aux1,aux2
external derivs,rk4 !
im = dcmplx(0.d0,1.d0)
!FIND A WAY TO CORRECTLY STORE tt AND ymatrix
do1: do i=1,n
  y(i)=y0(i) !loading initial conditions
  ymatrix(i,2)=y(i) !saving initial amplitudes
end do do1
tt(1)=t0 !saving initial time value
t=t0
h=(tf-t0)/npoints ! define time step size in atomic units
write(*,'(a24,f23.15)')'time step =',h
write(*,'(a24,f23.15)')'step size in q1 =',stepX

call derivs(t,y,n,dydt) ! evaluating y'(t,y)
!-----------------------------------------!
!Checking norm and conservation of energy !
soma=0.0d0                                !
Te=0.d0                                   !
do i=1,n                                  !
  soma=soma+conjg(y(i))*y(i)*stepX        !
  Te=Te + conjg(y(i)) * dydt(i)*im *stepX !
end do                                    !
!-----------------------------------------!
aux1=soma
aux2=Te
do3: do k=1,npoints
!write(*,'(a14,f23.15)')'current time =',t
  call derivs(t,y,n,dydt) ! evaluating y'(t,y)
  call rk4(y,dydt,n,t,h,y0,derivs) !evaluating y(t+h)
y(:)=y0(:)

  do2: do i=1,n
    ymatrix(i,k+2)=y(i) !saving the amplitude for time t
  end do do2
  t=t+h
  tt(k)=t !saving time value for each step in time
end do do3

!-----------------------------------------!
!Checking norm and conservation of energy !
soma=0.0d0                                !
Te=0.d0
call derivs(t,y,n,dydt)                                   !
do i=1,n                                  !
  soma=soma+conjg(y(i))*y(i)*stepX        !
  Te=Te + conjg(y(i)) * dydt(i)*im *stepX !
end do                                    !
!-----------------------------------------!
write(*,'(a24,2(f23.15))')'Norm conservation =',aux1-soma
write(*,'(a24,2(f23.15),a8)')'Energy conservation =',aux2-Te,' hartree'
write(*,'(a24,2(f23.15),a8)')'Final total energy =',Te,' hartree'

do i=1,n
  ymatrix(i,1)=(i-(n-1)/2-1)*0.01d0 + 2.09970623d0
  coordinate(i)=ymatrix(i,1)
end do
open(unit=10,file='amp-matrix.txt',status='unknown')
open(unit=20,file='time-amp.txt',status='unknown')
open(unit=30,file='amp-real.txt',status='unknown')
open(unit=40,file='amp-imag.txt',status='unknown')
do i=1,n
  write(10,'(<2*(npoints+2)>(e24.15))')ymatrix(i,:)
end do
do j=1,npoints+1
  do i=1,n
    write(20,'(f24.14,(f24.15))')coordinate(i), (real(ymatrix(i,j+1))**2.0d0+aimag(ymatrix(i,j+1))**2.0d0)
    write(30,'(f24.14,(f24.15))')coordinate(i), (real(ymatrix(i,j+1)))
    write(40,'(f24.14,(f24.15))')coordinate(i), (aimag(ymatrix(i,j+1)))
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
end subroutine rkdumb
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rk4(y0,dydt,n,t,h,Yout,derivs)
implicit none
integer, parameter :: dp = kind(1.d0) 
integer :: i,n
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
th=t+h/2.0d0
!dydt = y'(t , y(t)) = k1/h ---- dydt MUST NOT BE PREVIOUSLY MULTIPLIED BY h ----
do1: do i=1,n
  yt(i) = y0(i) + (h/2.d0) * dydt(i) !Evaluation of y(t) + k1/2.0d0
end do do1
call derivs(th,yt,n,dyt) !dyt = k2/h ---- Evaluation of y'(t + h/2.0d0 , y(t) + k1/2.0d0)
do2: do i=1,n
  yt(i) = y0(i) + (h/2.0d0) * dyt(i)  !Evaluation of y(t) + k2/2.0d0
end do do2
call derivs(th,yt,n,dym) !dym = k3/h ---- Evaluation of y'(t + h/2.0d0 , y(t) + k2/2.0d0)
do3: do i=1,n
  yt(i) = y0(i) + h * dym(i) !Evaluation of y(t) + k3
  dym(i) = dyt(i) + dym(i)  ! making dym = (k2 + k3)/h ---- variable dyt free now
end do do3
call derivs(t+h,yt,n,dyt) !dyt = k4/h ---- Evaluation of y'(t + h/2.0d0 , y(t) + k1/2.0d0)
do4: do i=1,n
  Yout(i) = y0(i) + h/6.0d0 * (dydt(i) + 2.0d0 * dym(i) + dyt(i))
end do do4
end subroutine rk4
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Subroutine computes the derivatives dPsi(t)/dt = H/i Psi
!This derivative routine calculates 
!Written as a simple matrix multiplication between the Hamiltonian matrix H and the amplitudes vector Psi
subroutine derivs(t,y,n,dydt)
use global_param
implicit none
integer, parameter :: dp = kind(1.d0) 
integer n
real(kind=dp) :: Ha(n,n),t,totmass
complex(kind=dp) :: y(n),dydt(n),kine(n),im,cm(n)
external hamiltonian

call pot
call hamiltonian(Ha)
totmass=2.0d0*nitro

im = dcmplx(0.d0,1.d0)
do i=1,Nq1
  kine(i)=dcmplx(0.d0,0.d0)
  cm(i)=dcmplx(0.d0,0.d0)
  do j=1,Nq1
!    kine(i) = kine(i) + dcmplx(-1.d0/(2.d0 * mredu) * Ha(i,j) * real(y(j)),-1.d0/(2.d0 * mredu) * Ha(i,j) * aimag(y(j)))
    kine(i) = kine(i) - 1.d0/(2.0d0*mredu) * Ha(i,j) * y(j) ! kinetic energy of the internal motion
!    kine(i) = kine(i) - 1.d0/(mredu) * Ha(i,j) * y(j) ! kinetic energy of the internal motion
!    cm(i) = cm(i) - 1.d0/(2.d0 * totmass) * ((2.0d0*nitro**2.0d0)/(totmass**2.0d0)) * Ha(i,j) * y(j) !kinetic energy of the center of mass
  end do
end do
!dydt(:)=(kine(:)+cm(:)+pot1(:)*y(:))/im
dydt(:)=(kine(:)+pot1(:)*y(:))/im
!dydt(:)=(kine(:))/im

!write(*,'(a7,2(f23.15))')'dydt =', kine(33)

!write(*,'(a7,2(f23.15))')'dydt =', dot_product(Ha(31,:),y)

!read(*,*)

end subroutine derivs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine hamiltonian(Ha)
use global_param
implicit none
integer, parameter :: dp = kind(1.d0) 
real(kind=dp) :: Ha(Nq1,Nq1),ind1,ind2,ind3,step,const

step=0.01d0 !step in X in atomic units
const = 1.d0/(12.d0*step**2.d0)

Ha(:,:)=0.0d0
ind1=-30.d0*const
ind2= 16.d0*const
ind3=-1.d0*const
Ha(1,1)=ind1; Ha(1,2)=ind2; Ha(1,3)=ind3
Ha(Nq1,Nq1)=ind1; Ha(Nq1,Nq1-1)=ind2; Ha(Nq1,Nq1-2)=ind3
Ha(2,1)=ind2; Ha(2,2)=ind1; Ha(2,3)=ind2; Ha(2,4)=ind3
Ha(Nq1-1,Nq1)=ind2; Ha(Nq1-1,Nq1-1)=ind1; Ha(Nq1-1,Nq1-2)=ind2; Ha(Nq1-1,Nq1-3)=ind3
do i=3,Nq1-2 !from the third line to the antepenulmate one
  Ha(i,i-2)=ind3
  Ha(i,i-1)=ind2
  Ha(i,  i)=ind1
  Ha(i,i+1)=ind2
  Ha(i,i+2)=ind3
end do

end subroutine hamiltonian
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













 
