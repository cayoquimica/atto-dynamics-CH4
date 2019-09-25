program dyna
implicit none
external rkdumb,rk4,derivs !!!!! VER COMO CHAMAR !!!!!
integer, parameter :: dp = kind(1.d0)
integer, parameter :: Nq1=61 !number of points of the grid along q1
integer m,i !m=dimension of the Hamiltonian (Nq1*Nq2*st x Nq1*Nq2*st)
real(kind=dp) :: a,e0,k0,pi !sig=width of gaussian for vec0; soma=variable for sums
complex(kind=dp), allocatable :: vec0(:) !vec0=
complex(kind=dp) :: soma,im,aux
!character(len=9) :: fmt
integer npoints !number of time steps to take in integration
real(kind=dp) :: t0,tf !t0=initial time for integration; tf=final time for integration
real(kind=dp) :: Mtotal,tt,teta,phi,const,w,stepX,ch
complex(kind=dp) ::term1,term2,vec1(61)
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
Mtotal=12000.d0
e0=0.005d0 ! E(t=0) = initial kinetic energy (In atomic units) --- 0.003674930360070 = approx  0.x1 ev
w=2.0d0*e0
a=sqrt(2.0d0/(Mtotal*w))
!a=4.0d0
!k0= sqrt(2.d0*Mtotal*e0) ! k0 = initial momentum (In atomic units) ! 150 = approx 0.385 hartree or 10.47 ev
k0=0.0d0
tt=0.0d0
teta=( atan( (2*tt)/(a**2.d0*Mtotal) ) ) /2.d0
phi=-teta-(k0**2.d0/(2.d0*Mtotal))*tt
const = ( (2.d0*a**2.d0)/pi )**(1.d0/4.d0) * exp(im*phi)/(a**4.d0+(4.d0*tt/Mtotal**2.d0))**(1.d0/4.d0)
stepX=0.01d0
do i=1,Nq1
  ch=(i-(Nq1-1.d0)/2.d0-1.d0)*stepX !change in x
  term1 = exp(im*k0*(ch))
  term2= exp( - ((ch) - (k0*tt/Mtotal))**2.d0 / (a**2.d0 + (2*im*tt/Mtotal)) )
  vec0(i) = const * term1 * term2
end do



!do i=1,Nq1
!  ch=(i-(Nq1-1.d0)/2.d0-1.d0)*stepX
!  x=x0+ch 
!  expo = dexp(-(x-x0)**2.d0*mredu*w/2.d0)
!  vec0(i) = c0 * expo
!end do
!----------------------------------------!
!normalizing                             !
soma=0.0d0                               !
do i=1,Nq1                               !
  soma=soma+conjg(vec0(i))*vec0(i)*stepX !
end do                                   !
vec0=vec0/sqrt(soma)                     !
!----------------------------------------!
write(*,'(a24,2(f23.15))')'normalization constant =',soma
end program
