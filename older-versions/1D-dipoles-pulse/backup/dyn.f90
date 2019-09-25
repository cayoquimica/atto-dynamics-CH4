program dyn
use global_param
implicit none
external rkdumb,rk4,derivs !!!!! VER COMO CHAMAR !!!!!
integer m !m=dimension of the Hamiltonian (Nq1*Nq2*st x Nq1*Nq2*st)
real(kind=8) :: a,e0,k0,pi !sig=width of gaussian for vec0; soma=variable for sums
complex(kind=8), allocatable :: vec0(:) !vec0=
complex(kind=8) :: soma,im,aux
!character(len=9) :: fmt
integer npoints !number of time steps to take in integration
real(kind=8) :: t0,tf !t0=initial time for integration; tf=final time for integration
real(kind=8) :: Mtotal,tt,teta,phi,const
complex(kind=8) ::term1,term2,vec1(61),yy(61,1001)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!definition of the initial conditions
!vec0 is the vector with the initial amplitudes for the nuclear wave funtcion
!build a Nq1*Nq2*nst
! em 1 dimensão e 1 primeiro:
!The minimum happens in position q1=-1.36 (16 displacement) and q2=0
!Defining as a normalized gaussian centered in the minimum position
allocate(vec0(Nq1))
im=dcmplx(0.d0,1.d0)
pi = 3.141592653589793d0
Mtotal=29244.18598791975d-6
a=4.d0 !width of the gaussian - I choose 5 randomly.
e0=0.3674930360070d-1 ! E(t=0) = initial kinetic energy (In atomic units) --- 0.003674930360070 = approx  0.x1 ev
k0= sqrt(2.d0*Mtotal*e0) ! k0 = initial momentum (In atomic units) ! 150 = approx 0.385 hartree or 10.47 ev

tt=0.0d0
teta=( atan( (2*tt)/(a**2.d0*Mtotal) ) ) /2.d0
phi=-teta-(k0**2.d0/(2.d0*Mtotal))*tt
const = ( (2.d0*a**2.d0)/pi )**(1.d0/4.d0)
do i=1,Nq1
  term1 = exp(im*phi)/(a**4.d0+(4.d0*tt/Mtotal**2.d0))**(1.d0/4.d0) * exp(im*k0*(i-16))
  term2= exp( - ((i-16) - (k0*tt/Mtotal))**2.d0 / (a**2.d0 + (2*im*tt/Mtotal)) )
  vec0(i) = const * term1 * term2
!  vec0(i)=(2.d0/(a**2*pi))**(1.d0/4.d0) * exp( (im*k0*(i-16)) ) * exp(-(i-16)**2.d0/a**2.d0)  !16 is the position of the minimum in q1 coordinate
!  vec0(i)=(2.d0/(a**2*pi))**(1.d0/4.d0) * exp( (im*k0*i) - (i**2.d0/a**2.d0) ) !16 is the position of the minimum in q1 coordinate
!  vec0(i)=(dsqrt(2.d0)/a) * exp( (im*k0*(i-16)) - ((i-16)**2.d0/a**2.d0) ) !16 is the position of the minimum in q1 coordinate
end do
!write(*,'(a23,f23.16)')'normalization factor =',(a/dsqrt(2*pi))**(1.d0/2.d0)
!----------------------------------!
!normalizing                       !
soma=0.0d0                         !
do i=1,Nq1                         !
  soma=soma+conjg(vec0(i))*vec0(i) !
end do                             !
vec0=vec0/sqrt(soma)               !
!----------------------------------!
write(*,'(2(f23.16))')soma
write(*,'(a8,2(f23.16))')'vec 16 =',vec0(16)
write(*,'(a8,2(f23.16))')'vec 61 =',vec0(61)
write(*,'(2(f23.16))')vec0
read(*,*)


yy(:,1)=vec0(:)

do j=1,1000
tt=0.0d0+j
teta=( atan( (2*tt)/(a**2.d0*Mtotal) ) ) /2.d0
phi=-teta-(k0**2.d0/(2.d0*Mtotal))*tt
const = ( (2.d0*a**2.d0)/pi )**(1.d0/4.d0)
do i=1,Nq1
term1 = exp(im*phi)/(a**4.d0+(4.d0*tt/Mtotal**2.d0))**(1.d0/4.d0) * exp(im*k0*(i-16))
term2= exp( - ((i-16) - (k0*tt/Mtotal))**2.d0 / (a**2.d0 + (2*im*tt/Mtotal)) )
vec1(i) = const * term1 * term2
yy(i,j) = vec1(i)
end do

end do

open(unit=10,file='yy.txt',status='unknown')
do i=1,61
  write(10,'(<2*(npoints+1)>(f23.16))') yy(i,:)
end do

write(*,*)'done'
!read(*,*)
!----------------------------------!
!normalizing                       !
soma=0.0d0                         !
!write(*,'(2(f23.16))') soma
do i=1,Nq1                         !
  soma=soma+conjg(vec1(i))*vec1(i) !
end do                             !
vec1=vec1/sqrt(soma)               !
!----------------------------------!
!write(*,'(2(f23.16))') soma
!write(*,*)''
!write(*,*)''
!write(*,'(2(f23.16))') vec1










! for a gaussian wave packet Psi(x,0) = A e^(-a x^2)
!with A = (2a/pi)^(1/4)
! Psi (x,t) = (2a/pi)^(1/4) e^(-ax^2/(1 + 2 i hbar a t/m)) / sqrt( 1 + 2 i hbar a t/m )




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=Nq1 ! m=Nq1*Nq2*st
! nanoseconds = 1d-9 seconds
! picoseconds = 1d-12 seconds
! femtoseconds = 1d-15 seconds
! attoseconds = 1d-18 seconds
! 1 time atomic units = 24.18884326505 attoseconds
t0=0.d0
tf=4.d0 !40 = Approx. 1 femtosecond
npoints=1000.d0 !If 1000, the step will be approx. 1 attosecond
call rkdumb(vec0,m,t0,tf,npoints,derivs)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end program
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rkdumb(y0,n,t0,tf,npoints,derivs)
integer i,k,n
real(kind=8) :: t0,tf,tt(npoints),h,t,aux,coordinate(n)
complex(kind=8) :: soma,dydt(n),y(n),y0(n),ymatrix(n,npoints+1),ymatlab(n,npoints),Te
external derivs,rk4 !
!FIND A WAY TO CORRECTLY STORE tt AND ymatrix
do1: do i=1,n
  y(i)=y0(i) !loading initial conditions
  ymatrix(i,2)=y(i) !saving initial amplitudes
end do do1
tt(1)=t0 !saving initial time value for each step
t=t0
h=(tf-t0)/npoints ! define time step size






!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv CHECK POINT PART TO BE REMOVE LATER vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!----------------------------------!
!normalizing                       !
soma=0.0d0                         !
do i=1,n                           !
  soma=soma+conjg(y(i))*y(i)       !
end do                             !
!----------------------------------!
write(*,'(a14,2(f23.16))')'initial norm =',sqrt(soma)
write(*,'(a11,f23.16)')'time step =',h
aux=0.d0
open(unit=20,file='test.txt',status='unknown')
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ CHECK POINT PART TO BE REMOVE LATER ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


do3: do k=2,npoints
  call derivs(t,y,n,dydt) ! evaluating y'(t,y)
!Printting total energy for each time step
  Te=0.d0
Etotal:  do j=1,n
    Te=Te + conjg(y(j)) * dydt(j) ! IS IT CORRECT?  CHECK IF <Psi*|H|Psi> SHOULD BE WRITTEN LIKE THAT OR SHOULD I USE MODULUS OF COMPLEX NUMBERS <<<<<<<<<<<<<<<<<<<<------------------------
  end do Etotal
Te=dot_product(conjg(y),dydt)
write(*,'(a15,2(f23.16),a8)')'Total energy = ',Te,' hartree'
!write(*,'(a7,2(f23.16))')'dydt = ',dydt(j)
  call rk4(y,dydt,n,t,h,y,derivs) !evaluating y(t+h)



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv CHECK POINT PART TO BE REMOVE LATER vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
write(*,*)''
write(*,'(2(f23.16))')y
read(*,*)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ CHECK POINT PART TO BE REMOVE LATER ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



  if (t+h == t) write(*,*)'stepsize not significant in rkdumb subroutine'
  do2: do i=1,n
    ymatrix(i,k+1)=y(i) !saving the amplitude for time t
  end do do2
  t=t+h
  tt(k)=t !saving time value for each step in time


open(unit=30,file='hpsi.txt',status='unknown')
open(unit=40,file='psic.txt',status='unknown')

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv CHECK POINT PART TO BE REMOVE LATER vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!----------------------------------!
!normalizing                       !
soma=0.0d0                         !
do i=1,n                           !
  soma=soma+conjg(y(i))*y(i)       !
write(30,'(2(f23.16))')dydt(i)
write(40,'(2(f23.16))')conjg(y(i))
end do                             !
!----------------------------------!
!write(*,'(4(f23.16))')dydt,conjg(y)
!read(*,*)
!soma=dot_product(y,y)
!write(*,*)'aqui'
!write(*,'(a6,2(f25.16))')'norm =',sqrt(soma)
if (aux /= 0d0) write(20,'(e25.16)')real(aux-soma)
aux=soma
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ CHECK POINT PART TO BE REMOVE LATER ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


!write(*,*)'y'
!write(*,*)y(18)

!write(*,*)'dydt'
!write(*,*)dydt(18)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ CHECK POINT PART TO BE REMOVE LATER ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




end do do3


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv CHECK POINT PART TO BE REMOVE LATER vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!write(*,*)'y'
!write(*,*)y(18)

!write(*,*)'dydt'
!write(*,*)dydt(18)
!read(*,*)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ CHECK POINT PART TO BE REMOVE LATER ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^





do i=1,n
  ymatrix(i,1)=(i-(n-1)/2-1)*0.08d0
  coordinate(i)=ymatrix(i,1)
end do
open(unit=10,file='amp-matrix.txt',status='unknown')
do i=1,n
  write(10,'(<2*(npoints+1)>(e24.16))')ymatrix(i,:)
end do
!open(unit=20,file='test.txt',status='unknown')
!do i=1,npoints
!  do j=1,n
!    write(20,'(2(e24.15))')ymatrix(j,i+1)
!  end do
!  write(20,*)''
!  write(20,*)''
!end do
 close(unit=10)
 close(unit=20)
end subroutine rkdumb
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 4th order runge-kutta
!h = (tfinal - t0) / ntemp
!t+1= t + h
!k1 = h * y'(t , f(t))
!k2 = h *  y'(t + h/2.0d0 , y(t) + k1/2.0d0)
!k3 = h *  y'(t + h/2.0d0 , y(t) + k2/2.0d0)
!k4 = h *  y'(t + h , y(t) + k3)
!y(t+1) = y(t) + 1/6 (k1+ 2*k2 + 2*k3 + k4)
subroutine rk4(y0,dydt,n,t,h,Yout,derivs)
implicit none
integer :: i,n
real(kind=8) :: t,h,th
complex(kind=8) :: dydt(n),y0(n),Yout(n),dym(n),dyt(n),yt(n)
external derivs 
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
!Subroutine computes the derivatives dPsi(t)/dt = H Psi
!Written as a simple matrix multiplication between the Hamiltonian matrix H and the amplitudes vector Psi
subroutine derivs(t,y,n,dydt)
use global_param
implicit none
integer n
real(kind=8) :: step,sum1,t,pote(n),s(3*NA,1) !change for correct dimension - s(3*NA,2) if q1 and q2 are active
complex(kind=8) :: y(n),dydt(n),d2dq1(n),kine(n),im
call coord
call mass_vec

im = dcmplx(0.d0,1.d0)
!============================== Kinetic Energy ===========================================================!|
s(:,1)=q1                                                                                                 !|
!s(:,2)=q2                                                                                                !|
step=1.d0                                                                                                 !|
!------------------------------------------------------------------------------------------------------!  !|
! Secon derivative in finite difference                                                                !  !|
findiff: do i=1,n                                                                                      !  !|
  if (i == 1) then                                                                                     !  !|
    d2dq1(i) = 1.d0/(12.d0*step**2.d0) * (                       -30.d0*y(i) + 16.d0*y(i+1) - y(i+2))  !  !|
!    d2dq1(i) = 1.d0/(12.d0*step**2.d0) * ( 45.d0*y(i)-154.d0*y(i+1)+214.d0*y(i+2)-156.d0*y(i+3)+61.d0*y(i+4)-10.d0*y(i+5))  !  !|
  else if (i == 2) then                                                                                !  !|
    d2dq1(i) = 1.d0/(12.d0*step**2.d0) * (          16.d0*y(i-1) -30.d0*y(i) + 16.d0*y(i+1) - y(i+2))  !  !|
!    d2dq1(i) = 1.d0/(12.d0*step**2.d0) * ( 10.d0*y(i-1)-15.d0*y(i)-4.d0*y(i+1)+14.d0*y(i+2)-6.d0*y(i+3) + y(i+4))  !  !|
  else if (i == n-1) then                                                                              !  !|
    d2dq1(i) = 1.d0/(12.d0*step**2.d0) * (-y(i-2) + 16.d0*y(i-1) -30.d0*y(i) + 16.d0*y(i+1)         )  !  !|
!    d2dq1(i) = 1.d0/(12.d0*step**2.d0) * ( 10.d0*y(i+1)-15.d0*y(i)-4.d0*y(i-1)+14.d0*y(i-2)-6.d0*y(i-3) + y(i-4))  !  !|
  else if (i == n-2) then                                                                              !  !|
    d2dq1(i) = 1.d0/(12.d0*step**2.d0) * (-y(i-2) + 16.d0*y(i-1) -30.d0*y(i)                        )  !  !|
!    d2dq1(i) = 1.d0/(12.d0*step**2.d0) * ( 45.d0*y(i)-154.d0*y(i-1)+214.d0*y(i-2)-156.d0*y(i-3)+61.d0*y(i-4)-10.d0*y(i-5))  !  !|
  else                                                                                                 !  !|
    d2dq1(i) = 1.d0/(12.d0*step**2.d0) * (-y(i-2) + 16.d0*y(i-1) -30.d0*y(i) + 16.d0*y(i+1) - y(i+2))  !  !|
  end if                                                                                               !  !|
end do findiff                                                                                         !  !|
!------------------------------------------------------------------------------------------------------!  !|
do i=1,n                                                                                                  !|
  sum1=0.d0                                                                                               !|
  do j=1,3*NA                                                                                             !|
    !sum1 = sum1 + s(j,1)**2.d0/mass(j)                                                                    !|
! include NAC                                                                                             !|
  end do                                                                                                  !|
!  d2dq1(i) = sum1 * d2dq1(i)                                                                              !|
end do                                                                                                    !|
!kine(:) = -1.d0/2.d0 * d2dq1(:)                                                                           !|
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv CHECK POINT PART TO BE REMOVE LATER vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
kine(:) = -1.d0/(2.d0 * 29244.18598791975d-6) * d2dq1(:)                                                                  !|
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ CHECK POINT PART TO BE REMOVE LATER ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!=========================================================================================================!|
  
!============================== Potential Energy =======================================!|
call pot                                                                                !|
pote(:) = pot1(:,41)                                                                    !|
!pote(:) = pot1(33,:)                                                                    !|
pote=(pote-pot1(16,41))*y                                                               !|
!=======================================================================================!|

dydt(:)=(kine(:))/im



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv CHECK POINT PART TO BE REMOVE LATER vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!write(*,*)''
!write(*,*)kine(18)
!write(*,*)dydt(18)
!write(*,*)''
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ CHECK POINT PART TO BE REMOVE LATER ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!dydt(:)=(kine(:)+pote(:))/im



!dydt=H*y !In the end, it must be just a matrix multiplication
!deallocate(mass,pot1,pot2,pot3,q1,q2)
!deallocate(pot1,pot2,pot3)
end subroutine derivs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



















 
