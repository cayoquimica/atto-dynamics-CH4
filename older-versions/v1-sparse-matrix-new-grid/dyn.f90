program dyn
!maybe is useful: Function MERGE:
!test a condition, if true, assign a value X, if false assign value Y:
!the condition can be a boolean variable defined before. 
!result = merge( X, Y, i.eq.j)
!result, X and Y will be of same type
use global_param
use omp_lib
implicit none
external rkdumb,HA_calc,momentum_calc_q1,momentum_calc_q2,angular_momentum
integer cont,n,jj,ii !m=dimension of the Hamiltonian (Nq1*Nq2*Nst x Nq1*Nq2*Nst)
real(kind=dp) :: a,e0,k0 !sig=width of gaussian for vec0; soma=variable for sums
integer npoints !number of time steps to take in integration
real(kind=dp) :: t0,tf !t0=initial time for integration; tf=final time for integration
real(kind=dp) :: tt,ch1,ch2,x,x0,kf,w,expo,c0,c1,tstep
complex(kind=dp), allocatable :: vec0(:),vec1(:)
real(kind=dp) :: q10,q20 ! point along q1 and q2 where the minimum global C2v is
real(kind=dp) :: truni,trunf,start_time,stop_time,ompt0,ompt1,tt1,tt0 !time counters
complex(kind=dp) :: mom(0:Nst*Nq1*Nq2-1),am(0:Nst*Nq1*Nq2-1) !vectors to store operators acting on the wave function
real(kind=dp) :: soma,soma1,soma2,soma3,p1q1,p2q1,p3q1,p1q2,p2q2,p3q2,L1,L2,L3 !variables to store expectation values
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! I modified the Ha, moq1 and moq2 matrices building so that their indexes begin with 0 instead of 1. 
! It is not the ideal way, I use a auxiliar matrix to build with index 1 as I already had written in the code and then I put it
! in a matrix with index starting in 0.
! Later would be better and faster if I already build the matrices with index 0.     <<<---------------

!$ truni = omp_get_wtime()
open(unit=100,file='output',status='unknown') ! output for following the code running
write(100,*)'INITIATING SIMULATION'

n=Nst*Nq1*Nq2 ! m=Nq1*Nq2*st %size of the Hamiltonian

!$ ompt0 = omp_get_wtime()
call load_data
open(unit=20,file='csr_vectors',status='unknown')
read(20,'(i12)')k_moq1
allocate(moq1_val(0:k_moq1-1),moq1_rowc(0:n), moq1_row_col(0:k_moq1-1,0:1))
do i=0,k_moq1-1
  read(20,'(e50.40e3,2i12)')moq1_val(i),moq1_row_col(i,0),moq1_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')moq1_rowc(i)
end do

read(20,'(i12)')k_moq2
allocate(moq2_val(0:k_moq2-1),moq2_rowc(0:n), moq2_row_col(0:k_moq2-1,0:1))
do i=0,k_moq2-1
  read(20,'(e50.40e3,2i12)')moq2_val(i),moq2_row_col(i,0),moq2_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')moq2_rowc(i)
end do

read(20,'(i12)')k_Ha
allocate(Ha_val(0:k_Ha-1),Ha_rowc(0:n), Ha_row_col(0:k_Ha-1,0:1))
do i=0,k_Ha-1
  read(20,'(e50.40e3,2i12)')Ha_val(i),Ha_row_col(i,0),Ha_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')Ha_rowc(i)
end do

read(20,'(i12)')k_dip
allocate(dip_val(0:k_dip-1))
do i=0,k_dip-1
  read(20,'(e50.40e3)')dip_val(i)
end do

read(20,'(i12)')k_am
allocate(am_val(0:k_am-1),am_rowc(0:n), am_row_col(0:k_am-1,0:1))
do i=0,k_am-1
  read(20,'(e50.40e3,2i12)')am_val(i),am_row_col(i,0),am_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')am_rowc(i)
end do

!$ ompt1 = omp_get_wtime()
write(100,*) 'The time to load the CSR vectors=',ompt1 - ompt0, "seconds"

allocate ( vec0(0:Nst*Nq1*Nq2-1) )
allocate ( vec1(0:Nst*Nq1*Nq2-1) )
vec0=dcmplx(0.d0,0.d0)
vec1=dcmplx(0.d0,0.d0)
!! Follow mean value of momentum (k) through time and check if is always smaller than the maximum momentum from sampling
!############################################################################################################!#
! Buildin initial wave packet in 2D as the one of an Harmonic oscilator                                      !# 
t0=0.d0                                                                                                      !#      
tf=100.0d0               !40 = Approx. 1 femtosecond                                                        !#  
npoints=1000             !If 1000, the step will be approx. 1 attosecond                                    !#
tstep=(tf-t0)/npoints    ! define time step size in atomic units                                             !#
!                                                                                                            !#
e0=4.0d-4 !initial energy for a harmonic case                                                                !#
k0= sqrt(2.d0*mtotal*e0) ! k0 = initial momentum (In atomic units)                                           !#
k0=0.d0                                                                                                      !#
tt=0.d0                                                                                                      !#
!kf=4.d0*mtotal*e0**2.d0 ! force constant of harmonic potential                                              !#
w=2.d0*e0 !w=sqrt(kf/mredu) ---- 2.d0*e0 means the energy of the ground state of the harmonic oscilator      !#
c0 = ((mtotal*w)/(pi))**(1.d0/4.d0)                                                                          !#
c1 = (4.d0/pi*(mtotal*w)**3.d0)**(1.d0/4.d0)                                                                 !#
cont=0.d0                                                                                                    !#
q10=-0.4d0   ! mudei só pra ver o momento variando periodicamente                                            !#
q20= 0.d0                                                                                                    !#
do j=0,Nq2-1                                                                                                 !#
  do i=0,Nq1-1                                                                                               !#
    ch1=(i+1.d0-Nq1/2.d0-1.d0)*sq1+sq1/2.d0 !change in q1                                                    !#
    ch2=(j+1.d0-(114.d0)-1.d0)*sq2 !change in q2                                                             !#
!    expo = dexp(-(ch-q10)**2.d0*mtotal*w/2.d0) !* exp(im*k0*(ch))                                           !#
    expo = dexp((-(ch1-q10)**2.d0*mtotal*w/2.d0)+(-(ch2-q20)**2.d0*mtotal*w/2.d0)) !* exp(im*k0*(ch))        !#
    vec0(cont)= c0 * expo                                                                                    !#
!    vec0(cont+1*s)= c0 * expo                                                                                !#
!    vec0(cont+2*s)= c0 * expo                                                                                !#
    cont=cont+1;                                                                                             !#
  end do                                                                                                     !#
end do                                                                                                       !# 
!############################################################################################################!#
!--------------------------------------------!
!normalizing                                 !
soma=0.d0                                    !
do i=0,n-1                                   !
  soma=soma+dconjg(vec0(i))*vec0(i)          !
end do                                       !
vec0=vec0/sqrt(soma)                         !
!--------------------------------------------!
write(100,'(a24,(es15.7e3))')'normalization constant =',soma


!$ tt0 = omp_get_wtime()
call angular_momentum(vec0,n,am)
!$ tt1 = omp_get_wtime()


!--------------------------------------------------------!
!Checking initial energy and momentum                    !
!$ ompt0 = omp_get_wtime()                               !
call HA_calc(tt,vec0,n,vec1,tstep) !evaluating y'(t=0,y) !
!$ ompt1 = omp_get_wtime()                               !
!$ start_time = omp_get_wtime()                          !
call momentum_calc_q1(vec0,mom,n) ! evaluating dydq      !
!$ stop_time = omp_get_wtime()                           !
soma=0.d0                                                !
soma1=0.d0                                               !
soma2=0.d0                                               !
soma3=0.d0                                               !
p1q1=0.0d0                                               !
p2q1=0.0d0                                               !
p3q1=0.0d0                                               !
p1q2=0.0d0                                               !
p2q2=0.0d0                                               !
p3q2=0.0d0                                               !
L1=0.d0                                                  !
L2=0.d0                                                  !
L3=0.d0                                                  !
do i=0,s-1                                               !
  soma1=soma1 + dconjg(vec0(i)) * vec1(i)*im             !
  soma2=soma2 + dconjg(vec0(i+s)) * vec1(i+s)*im         !
  soma3=soma3 + dconjg(vec0(i+2*s)) * vec1(i+2*s)*im     !
  p1q1=p1q1 + dconjg(vec0(i)) * mom(i)                   !
  p2q1=p2q1 + dconjg(vec0(i+s)) * mom(i+s)               !
  p3q1=p3q1 + dconjg(vec0(i+2*s)) * mom(i+2*s)           !
  L1    = L1    + dconjg(vec0(i))     * am(i)            !
  L2    = L2    + dconjg(vec0(i+s))   * am(i+s)          !
  L3    = L3    + dconjg(vec0(i+2*s)) * am(i+2*s)        !
end do                                                   !
call momentum_calc_q2(vec0,mom,n) ! evaluating dydq      !
do i=0,s-1                                               !
  p1q2=p1q2 + dconjg(vec0(i)) * mom(i)                   !
  p2q2=p2q2 + dconjg(vec0(i+s)) * mom(i+s)               !
  p3q2=p3q2 + dconjg(vec0(i+2*s)) * mom(i+2*s)           !
end do                                                   !
write(100,'(a24,e23.15e3,a8)')'inital momentum in q1 =',p1q1+p2q1+p3q1,' au.'
write(100,'(a24,e23.15e3,a8)')'inital momentum in q2 =',p1q2+p2q2+p3q2,' au.'
!--------------------------------------------------------!
E_init=soma1+soma2+soma3 !Storing the initial energy 
write(100,*)'##############################################################################'
write(100,*)'Time the program took to do the operation H|Psi> is:',(ompt1-ompt0), "seconds"
write(100,*)'Be ready to wait probably ', npoints*(tt1-tt0+(stop_time-start_time)*2.d0+(ompt1-ompt0)*5.d0)/(3600.d0) , "hours"
write(100,*)'##############################################################################'
write(100,'(a24,3(f23.15))')'initial angular momentum = ',L1,L2,L3,L1+L2+L3
write(100,'(a24,(f23.15),a8)')'inital energy =',soma1+soma2+soma3,' hartree'
write(100,'(a24,(f23.15))')'initial 1 / dq1 =',1.d0/sq1
write(100,'(a24,(f23.15))')'initial 1 / dq2 =',1.d0/sq2
write(100,'(a24,(f23.15))')'initial k * dq =',dsqrt(2*mtotal*e0)*sq1*sq2 
write(100,'(a24,(f23.15),a4)')'vibrational frequency =',w,' au.'
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! nanoseconds = 1d-9 seconds
! picoseconds = 1d-12 seconds
! femtoseconds = 1d-15 seconds
! attoseconds = 1d-18 seconds
! 1 time atomic units = 24.18884326505 attoseconds
call rkdumb(vec0,n,t0,tf,npoints,tstep,HA_calc)
!$ trunf = omp_get_wtime()
write(100,*)'Time took to run totaly the program = ', (trunf-truni)/3600.d0, 'hours.'
write(100,*)'************************************************************'
write(100,*)'Program finished'
write(100,*)'************************************************************'
end program
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!subroutine rkdumb(y0,n,t0,tf,npoints,h,HA_calc)
subroutine rkdumb(y,n,t0,tf,npoints,h,HA_calc)
use global_param
external rk4,HA_calc,momentum_calc_q1,momentum_calc_q2,angular_momentum
integer           :: n,npoints,ii,gg
real(kind=dp)     :: t0,tf,tt(npoints),h,t!,coordinate1(s),coordinate2(s)
complex(kind=dp)  :: dydt(0:n-1),y(0:n-1)!,y0(n)!,ymatrix(n,npoints+1)!,ymatlab(n,npoints)
real(kind=dp)     :: soma,Te,aux1,aux2,Te1,Te2,Te3!,enert(npoints+1) ! This is fo saving total energy through time
real(kind=dp)     :: pq1_1,pq2_1,pq1_2,pq2_2,pq1_3,pq2_3 !This is for saving momentum through time
complex(kind=dp)  :: momq1(0:n-1),momq2(0:n-1),am(0:n-1) ! This is for saving momentum through time
real(kind=dp)     :: sum1,sum2,sum3,momq1t(npoints),momq2t(npoints) ! This is for saving norm through time
real(kind=dp)     :: e1,e2,e3,L1,L2,L3
character(len=23) :: fname
!real(kind=dp)     :: st1n(npoints),st1p(npoints),st2n(npoints),st2p(npoints),st3n(npoints),st3p(npoints)
!real(kind=dp)     :: st1ncount,st1pcount,st2ncount,st2pcount,st3ncount,st3pcount,pop1,pop2,pop3


!open(unit=10,file='norms-time.data',status='unknown')
!open(unit=11,file='momentum-time.data',status='unknown')
!open(unit=12,file='electric-field-time.data',status='unknown')
!open(unit=13,file='energy-time.data',status='unknown')
!open(unit=14,file='ang-mom.data',status='unknown')
open(unit=15,file='alldata.data',status='unknown')

t=t0
write(100,'(a27,f23.15,a4)')'time step =',h,' au.'
write(100,'(a27,f23.15)')'step size in coordinate 1 =',sq1
write(100,'(a27,f23.15)')'step size in coordinate 2 =',sq2

!--------------------------------------------------------!
!Checking momentum and saving norm                       !
call HA_calc(t,y,n,dydt,h) ! evaluating y'(t,y)
call momentum_calc_q1(y,momq1,n) ! evaluating dydq       !
call momentum_calc_q2(y,momq2,n) ! evaluating dydq       !
call angular_momentum(y,n,am)                            !
pq1_1=0.0d0 ! momentum in q1 state1                      !
pq2_1=0.0d0 ! momentum in q2 state1                      !
pq1_2=0.0d0 ! momentum in q1 state2                      !
pq2_2=0.0d0 ! momentum in q2 state2                      !
pq1_3=0.0d0 ! momentum in q1 state3                      !
pq2_3=0.0d0 ! momentum in q2 state3                      !
sum1=0.0d0  ! norm for state 1                           !
sum2=0.0d0  ! norm for state 2                           !
sum3=0.0d0  ! norm for state 3                           !
e1=0.d0     ! energy of state 1                          !
e2=0.d0     ! energy of state 2                          !
e3=0.d0     ! energy of state 3                          !
L1=0.d0     ! angular momentum of state 1                !
L2=0.d0     ! angular momentum of state 2                !
L3=0.d0     ! angular momentum of state 3                !
do i=0,s-1                                               !
  pq1_1 = pq1_1 + dconjg(y(i))     * momq1(i)            !
  pq2_1 = pq2_1 + dconjg(y(i))     * momq2(i)            !
  pq1_2 = pq1_2 + dconjg(y(i+s))   * momq1(i+s)          !
  pq2_2 = pq2_2 + dconjg(y(i+s))   * momq2(i+s)          !
  pq1_3 = pq1_3 + dconjg(y(i+2*s)) * momq1(i+2*s)        !
  pq2_3 = pq2_3 + dconjg(y(i+2*s)) * momq2(i+2*s)        !
  sum1  = sum1  + dconjg(y(i))     * y(i)                !
  sum2  = sum2  + dconjg(y(i+s))   * y(i+s)              !
  sum3  = sum3  + dconjg(y(i+2*s)) * y(i+2*s)            !
  e1    = e1    + dconjg(y(i))     * dydt(i)*im          !
  e2    = e2    + dconjg(y(i+s))   * dydt(i+s)*im        !
  e3    = e3    + dconjg(y(i+2*s)) * dydt(i+2*s)*im      !
  L1    = L1    + dconjg(y(i))     * am(i)               !
  L2    = L2    + dconjg(y(i+s))   * am(i+s)             !
  L3    = L3    + dconjg(y(i+2*s)) * am(i+2*s)           !
end do                                                   !
!--------------------------------------------------------!
soma=sum1+sum2+sum3
Te=e1+e2+e3
aux1=soma
aux2=Te
!nst1(1)=1.d0
!nst2(1)=0.d0
!nst3(1)=0.d0
!enert(1)=Te
write(15,'(a3,3(e26.14e3))')'# ',e1,e2,e3
write(15,*) '# time, Pulse, norm1, norm2, norm3, E1, E2, E3, angular momentum st1, angular momentum st2, angular momentum st3,&
linear momentum in q1 st1, linear momentum in q2 st1, linear momentum in q1 st2, linear momentum in q2 st2,&
linear momentum in q1 st3, linear momentum in q2 st3'
Et=(E00/freq)*(-(t-t00)/sig**2.d0*sin(freq*(t-t00)+phase)+freq*cos(freq*(t-t00)+phase))*dexp(-(t-t00)**2.d0/(2.d0*sig**2.d0))
write(15,'(17(es26.16e3))')t,Et,sum1,sum2,sum3,e1,e2,e3,L1,L2,L3,pq1_1,pq2_1,pq1_2,pq2_2,pq1_3,pq2_3


write(100,*) '************************************************************'
write(100,*) 'Integrating amplitudes over time'
write(100,*) '************************************************************'
write(100,*) 'time ,   Pulse  ,     E1    ,    E2     ,    E3     ,   ET-E0   ,   norm1   ,   norm2   ,   norm3   , NormT-1   &
,   Ltot    '
ii=0;gg=0
do ll=1,1000 ! for saving 1000 time samples
  gg=gg+1
  do k=1,(npoints/1000)
    ii=ii+1
    call HA_calc(t,y,n,dydt,h) ! evaluating y'(t,y)
    call rk4(y,dydt,n,t,h,y,HA_calc) !evaluating y(t+h)
    !---------------------------------------------!
    !Saving norm and total energy in time         !
    soma=0.0d0                                    !
    Te=0.d0                                       !
    do i=0,n-1                                    !
      soma=soma+dconjg(y(i))*y(i)                 !
      Te=Te + dconjg(y(i)) * dydt(i)*im           !
    end do                                        !
    !---------------------------------------------!
    !--------------------------------------------------------!
    !Checking momentum and saving norm                       !
    call momentum_calc_q1(y,momq1,n) ! evaluating dydq       !
    call momentum_calc_q2(y,momq2,n) ! evaluating dydq       !
    call angular_momentum(y,n,am)                            !
    pq1_1=0.0d0 ! momentum in q1 state1                      !
    pq2_1=0.0d0 ! momentum in q2 state1                      !
    pq1_2=0.0d0 ! momentum in q1 state2                      !
    pq2_2=0.0d0 ! momentum in q2 state2                      !
    pq1_3=0.0d0 ! momentum in q1 state3                      !
    pq2_3=0.0d0 ! momentum in q2 state3                      !
    sum1=0.0d0  ! norm for state 1                           !
    sum2=0.0d0  ! norm for state 2                           !
    sum3=0.0d0  ! norm for state 3                           !
    e1=0.d0     ! energy of state 1                          !
    e2=0.d0     ! energy of state 2                          !
    e3=0.d0     ! energy of state 3                          !
    L1=0.d0     ! angular momentum of state 1                !
    L2=0.d0     ! angular momentum of state 2                !
    L3=0.d0     ! angular momentum of state 3                !
    do i=0,s-1                                               !
      pq1_1 = pq1_1 + dconjg(y(i))     * momq1(i)            !
      pq2_1 = pq2_1 + dconjg(y(i))     * momq2(i)            !
      pq1_2 = pq1_2 + dconjg(y(i+s))   * momq1(i+s)          !
      pq2_2 = pq2_2 + dconjg(y(i+s))   * momq2(i+s)          !
      pq1_3 = pq1_3 + dconjg(y(i+2*s)) * momq1(i+2*s)        !
      pq2_3 = pq2_3 + dconjg(y(i+2*s)) * momq2(i+2*s)        !
      sum1  = sum1  + dconjg(y(i))     * y(i)                !
      sum2  = sum2  + dconjg(y(i+s))   * y(i+s)              !
      sum3  = sum3  + dconjg(y(i+2*s)) * y(i+2*s)            !
      e1    = e1    + dconjg(y(i))     * dydt(i)*im          !
      e2    = e2    + dconjg(y(i+s))   * dydt(i+s)*im        !
      e3    = e3    + dconjg(y(i+2*s)) * dydt(i+2*s)*im      !
      L1    = L1    + dconjg(y(i))     * am(i)               !
      L2    = L2    + dconjg(y(i+s))   * am(i+s)             !
      L3    = L3    + dconjg(y(i+2*s)) * am(i+2*s)           !
    end do                                                   !
    !--------------------------------------------------------!
!normt(k+1) =soma
!enert(ii+1) =Te
momq1t(ii)=pq1_1+pq1_2+pq1_3
momq2t(ii)=pq2_1+pq2_2+pq2_3
!nst1(ii+1)=sum1
!nst2(ii+1)=sum2
!nst3(ii+1)=sum3
!-----------------------------------------------------!
! saving the population at the borders for each state !
!st1n(k)=dconjg(y(1))*y(1)                             !
!st1p(k)=dconjg(y(s-54))*y(s-54)                         !
!st2n(k)=dconjg(y(s+1))*y(s+1)                     !
!st2p(k)=dconjg(y(s+s-54))*y(s+s-54)                 !
!st3n(k)=dconjg(y(2*s+1))*y(2*s+1)                 !
!st3p(k)=dconjg(y(n-54))*y(n-54)                             !
!-----------------------------------------------------!
    t=t+h
  end do
! ==================== STARTING SAVING PROCEDURE =====================

!fname='time-amp-real-000000.h5'
!write(fname(15:20),'(i0.6)') ii
!call save_vector_h5(real(y),n,fname,23)
!fname='time-amp-imag-000000.h5'
!write(fname(15:20),'(i0.6)') ii
!call save_vector_h5(aimag(y),n,fname,23)

  fname='amp-time-st1-000000.dat'
  write(fname(14:19),'(i0.6)') gg
  open(unit=20,file=fname,status='unknown')
  fname='amp-real-st1-000000.dat'
  write(fname(14:19),'(i0.6)') gg
  open(unit=30,file=fname,status='unknown')
  fname='amp-imag-st1-000000.dat'
  write(fname(14:19),'(i0.6)') gg
  open(unit=40,file=fname,status='unknown')
  fname='amp-time-st2-000000.dat'
  write(fname(14:19),'(i0.6)') gg
  open(unit=50,file=fname,status='unknown')
  fname='amp-real-st2-000000.dat'
  write(fname(14:19),'(i0.6)') gg
  open(unit=60,file=fname,status='unknown')
  fname='amp-imag-st2-000000.dat'
  write(fname(14:19),'(i0.6)') gg
  open(unit=70,file=fname,status='unknown')
  fname='amp-time-st3-000000.dat'
  write(fname(14:19),'(i0.6)') gg
  open(unit=80,file=fname,status='unknown')
  fname='amp-real-st3-000000.dat'
  write(fname(14:19),'(i0.6)') gg
  open(unit=90,file=fname,status='unknown')
  fname='amp-imag-st3-000000.dat'
  write(fname(14:19),'(i0.6)') gg
  open(unit=99,file=fname,status='unknown')
  do i=0,s-1
    write(20,'(3(es26.16e3))') real(dconjg(y(i))*y(i))
    write(30,'(3(es26.16e3))') (  real(y(i)) )
    write(40,'(3(es26.16e3))') ( aimag(y(i)) )
    write(50,'(3(es26.16e3))') real(dconjg(y(1*s+i))*y(1*s+i))
    write(60,'(3(es26.16e3))') (  real(y(1*s+i)) )
    write(70,'(3(es26.16e3))') ( aimag(y(1*s+i)) )
    write(80,'(3(es26.16e3))') real(dconjg(y(2*s+i)) * y(2*s+i))
    write(90,'(3(es26.16e3))') (  real(y(2*s+i)) )
    write(99,'(3(es26.16e3))') ( aimag(y(2*s+i)) )
  end do
  close(unit=20)
  close(unit=30)
  close(unit=40)
  close(unit=50)
  close(unit=60)
  close(unit=70)
  close(unit=80)
  close(unit=90)
  close(unit=99)

!  write(10,'(3(es26.16e3))') sum1,sum2,sum3
!  write(11,'(6(es26.16e3))') pq1_1,pq2_1,pq1_2,pq2_2,pq1_3,pq2_3
Et = (E00/freq)*(-(t-t00)/sig**2.d0*sin(freq*(t-t00)+phase)+freq*cos(freq*(t-t00)+phase))*dexp(-(t-t00)**2.d0/(2.d0*sig**2.d0))
!  write(12,'((es26.16e3))') Et
!  write(13,'(3(es26.16e3))') e1,e2,e3
!  write(14,'(3(es26.16e3))') L1,L2,L3
  write(15,'(17(es26.16e3))') t,Et,sum1,sum2,sum3,e1,e2,e3,L1,L2,L3,pq1_1,pq2_1,pq1_2,pq2_2,pq1_3,pq2_3
!  tt(ii)=t !saving time value for each step in time
  write(100,'(i9,11(e12.3e3))') ii,Et,e1,e2,e3,Te-E_init,sum1,sum2,sum3,sum1+sum2+sum3-1.d0,L1+L2+L3
  !WRITE(6,'(5(A))',ADVANCE="NO") "\b","\b","\b","\b","b"
  !write(6,'(f5.1,"%")',advance='no') ii/10.d0
  !flush(6)
end do

write(100,*)''
! close(unit=10)
! close(unit=11)
! close(unit=12)
! close(unit=13)
! close(unit=14)
 close(unit=15)
!---------------------------------------------------------------------------------------------------------------------------!
! Determining the highest population at the simulation                                                                      !
! must be always small, > 10d-8. If not, consider changing grid borders and check the potential energy behaviour at borders !
!st1ncount=1.d-100; st1pcount=1.d-100; st2ncount=1.d-100; st2pcount=1.d-100; st3ncount=1.d-100; st3pcount=1.d-100;           !
!do k=1,npoints                                                                                                              !
!  if (st1n(k) > st1ncount) then                                                                                             !
!    st1ncount=st1n(k)                                                                                                       !
!  end if                                                                                                                    !
!  if (st1p(k) > st1pcount) then                                                                                             !
!    st1pcount=st1p(k)                                                                                                       !
!  end if                                                                                                                    !
!  if (st2n(k) > st2ncount) then                                                                                             !
!    st2ncount=st2n(k)                                                                                                       !
!  end if                                                                                                                    !
!  if (st2p(k) > st2pcount) then                                                                                             !
!    st2pcount=st2p(k)                                                                                                       !
!  end if                                                                                                                    !
!  if (st3n(k) > st3ncount) then                                                                                             !
!    st3ncount=st3n(k)                                                                                                       !
!  end if                                                                                                                    !
!  if (st3p(k) > st3pcount) then                                                                                             !
!    st3pcount=st3p(k)                                                                                                       !
!  end if                                                                                                                    !
!end do                                                                                                                      !
!---------------------------------------------------------------------------------------------------------------------------!
!write(100,'(a50,(es11.2e3))')'maximum population at negative border - state 1 =',st1ncount
!write(100,'(a50,(es11.2e3))')'maximum population at positive border - state 1 =',st1pcount
!write(100,'(a50,(es11.2e3))')'maximum population at negative border - state 2 =',st2ncount
!write(100,'(a50,(es11.2e3))')'maximum population at positive border - state 2 =',st2pcount
!write(100,'(a50,(es11.2e3))')'maximum population at negative border - state 3 =',st3ncount
!write(100,'(a50,(es11.2e3))')'maximum population at positive border - state 3 =',st3pcount
!-----------------------------------------!
!Checking norm and conservation of energy !
soma=0.0d0                                !
sum1=0.0d0                                !
sum2=0.0d0                                !
sum3=0.0d0                                !
Te1=0.d0                                  !
Te2=0.d0                                  !
Te3=0.d0                                  !
call HA_calc(t,y,n,dydt,h)                !
do i=0,n-1                                !
  soma=soma+dconjg(y(i))*y(i)             !
end do                                    !
do i=0,s-1                                !
  sum1=sum1+dconjg(y(i))*y(i)             !
  sum2=sum2+dconjg(y(i+s))*y(i+s)         !
  sum3=sum3+dconjg(y(i+2*s))*y(i+2*s)     !
  Te1=Te1+dconjg(y(i))    *dydt(i)*im     !
  Te2=Te2+dconjg(y(i+s))  *dydt(i+s)*im   !
  Te3=Te3+dconjg(y(i+2*s))*dydt(i+2*s)*im !
end do                                    !
!-----------------------------------------!
write(100,'(a30,(e12.3e3))')'Norm conservation =',aux1-soma
write(100,'(a30,(e12.3e3))')'Final norm for state 1 =',sum1
write(100,'(a30,(e12.3e3))')'Final norm for state 2 =',sum2
write(100,'(a30,(e12.3e3))')'Final norm for state 3 =',sum3
write(100,'(a30,(e12.3e3),a8)')'Energy conservation =',aux2-Te,' hartree'
write(100,'(a30,(f20.15),a8)')'Final total energy state 1 =',Te1,' hartree'
write(100,'(a30,(f20.15),a8)')'Final total energy state 2 =',Te2,' hartree'
write(100,'(a30,(f20.15),a8)')'Final total energy state 3 =',Te3,' hartree'
!write(100,*)'************************************************************'
!write(100,*)'Saving amplitude over time'
!write(100,*)'************************************************************'
!call cpu_time(trun0)
!ii=0
!do j=1,Nq2
!  do i=1,Nq1
!    ii=ii+1
!    coordinate1(ii)=i
!    coordinate2(ii)=j
!  end do
!end do
!call save_vector_h5(coordinate1,Nq1*Nq2,'coordinate1.h5',14)
!call save_vector_h5(coordinate2,Nq1*Nq2,'coordinate2.h5',14)
!call save_vector_h5(coordinate2,Nq1*Nq2,'normst1.h5',14)
!call save_matrix_h5(ymatrix,n,npoints+1,'all-amp.h5',10)
!!call saveh5(coordinate1,coordinate2,ymatrix,npoints+1)

!open(unit=10,file='norm-time.txt',status='unknown')
!open(unit=20,file='amp-time.txt',status='unknown')
!open(unit=30,file='amp-real.txt',status='unknown')
!open(unit=40,file='amp-imag.txt',status='unknown')
!!open(unit=50,file='field-time.txt',status='unknown')
!!do i=1,n
!!  write(10,'(<2*(npoints+2)>(e24.15))')ymatrix(i,:)
!!end do
!!t=0.d0
!do j=1,npoints+1
!ett = (E00/freq)*(-(t-t00)/sig**2.d0*sin(freq*(t-t00)+phase)+freq*cos(freq*(t-t00)+phase))*dexp(-(t-t00)**2.d0/(2.d0*sig**2.d0))
!t=t+h
!  do i=1,s
!    write(20,'(3(es26.16e3))') real(dconjg(ymatrix(i,j))*ymatrix(i,j)),&
!real(dconjg(ymatrix(1*s+i,j))*ymatrix(1*s+i,j)),real(dconjg(ymatrix(2*s+i,j)) * ymatrix(2*s+i,j))
!    write(30,'(3(es26.16e3))') (real(ymatrix(i,j))),(real(ymatrix(1*s+i,j))),(real(ymatrix(2*s+i,j)))
!    write(40,'(3(es26.16e3))') (aimag(ymatrix(i,j))),(aimag(ymatrix(1*s+i,j))),(aimag(ymatrix(2*s+i,j)))
!  end do
!  write(20,*)''
!  write(20,*)''
!  write(30,*)''
!  write(30,*)''
!  write(40,*)''
!  write(40,*)''
!  write(10,'(3(es26.16e3))') nst1(j),nst2(j),nst3(j)
!end do

! close(unit=10)

!call cpu_time(trun1)
!write(100,*) 'Time to save final files= ', (trun1-trun0)/(16.d0*60.d0), 'minutes.'
write(100,'(a24,e12.3e3)') 'maximum momentum value =',maxval(momq1t)
write(100,'(a24,e12.3e3)') 'maximum momentum value =',maxval(momq2t)
!call save_vector_h5(ett,npoints+1,'electric-field.h5',17)
! close(unit=50)


!open(unit=60,file='data.txt',status='unknown')
!t=0.d0
!do i=1,100+1 ! tf / step(4 = 0.1 femtoseconds)
!  t=(i-1.d0)*1.d0
!  Et = (E00/freq)*(-(t-t00)/sig**2.d0*sin(freq*(t-t00)+phase)+freq*cos(freq*(t-t00)+phase))*dexp(-(t-t00)**2.d0/(2.d0*sig**2.d0))
!  pop1=0.d0
!  pop2=0.d0
!  pop3=0.d0
!  do j=1,s
!    pop1=pop1 + dconjg(ymatrix(j,t*10+1)) * ymatrix(j,t*10+1) *sq1*sq2
!    pop2=pop2 + dconjg(ymatrix(1*s+j,t*10+1)) * ymatrix(1*s+j,t*10+1) *sq1*sq2    ! NOT SURE OF THE TIME STEP HERE -----REVISE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    pop3=pop3 + dconjg(ymatrix(2*s+j,t*10+1)) * ymatrix(2*s+j,t*10+1) *sq1*sq2
!  end do
! write(60,'(f5.1,6(es20.10e3))') t,Et,pop1,normt(t*10+1),enert(t*10+1),pop2,pop3
!end do
! close(unit=60)

open(unit=222,file='final-wave-packet',status='unknown')
do i=0,n-1
  write(222,'(2(e26.14e3))') y(i)
end do
 close(unit=222)

end subroutine rkdumb
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rk4(y0,dydt,n,t,h,Yout,HA_calc)
use global_param
implicit none
integer n
real(kind=dp) :: t,h,th
complex(kind=dp) :: dydt(0:n-1),y0(0:n-1),Yout(0:n-1),dym(0:n-1),dyt(0:n-1),yt(0:n-1)
external HA_calc 
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
do1: do i=0,n-1
  yt(i) = y0(i) + (h/2.d0) * dydt(i) !Evaluation of y(t) + k1/2.0d0
end do do1
call HA_calc(th,yt,n,dyt,h/2.d0) !dyt = k2/h ---- Evaluation of y'(t + h/2.0d0 , y(t) + k1/2.0d0)
do2: do i=0,n-1
  yt(i) = y0(i) + (h/2.d0) * dyt(i)  !Evaluation of y(t) + k2/2.0d0
end do do2
call HA_calc(th,yt,n,dym,h/2.d0) !dym = k3/h ---- Evaluation of y'(t + h/2.0d0 , y(t) + k2/2.0d0)
do3: do i=0,n-1
  yt(i) = y0(i) + h * dym(i) !Evaluation of y(t) + k3
  dym(i) = dyt(i) + dym(i)  ! making dym = (k2 + k3)/h ---- variable dyt free now
end do do3
call HA_calc(t+h,yt,n,dyt,h*2.d0) !dyt = k4/h ---- Evaluation of y'(t + h/2.0d0 , y(t) + k1/2.0d0)
do4: do i=0,n-1
  Yout(i) = y0(i) + h/6.d0 * (dydt(i) + 2.d0 * dym(i) + dyt(i))
end do do4
end subroutine rk4
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Subroutine computes the derivatives dPsi(t)/dt = H/i Psi
!This derivative routine calculates 
!Written as a simple matrix multiplication between the Hamiltonian matrix H and the amplitudes vector Psi
subroutine HA_calc(t,y,n,dydt,tstep)
use global_param
use omp_lib
implicit none
integer n
real(kind=dp) :: t,tstep!,c1q1,c2q1,c3q1,c1q2,c2q2,c3q2,c1q3,c2q3,c3q3
complex(kind=dp) :: y(0:n-1),dydt(0:n-1)

!write(*,*) 'omp_get_max_threads= ', omp_get_max_threads ( )
!write(*,*) 'omp_get_num_procs = ', omp_get_num_procs ( )
!write(*,*) 'Time = ', omp_get_wtime ( )

Et = (E00/freq)*(-(t-t00)/sig**2.d0*dsin(freq*(t-t00)+phase)+freq*dcos(freq*(t-t00)+phase))*dexp(-(t-t00)**2.d0/(2.d0*sig**2.d0))
dydt=dcmplx(0.d0,0.d0)
!!$OMP PARALLEL DO shared(y,dydt,Ha_val,dip_val,Ha_row_col,Et)
do i=0,n-1
  do j=Ha_rowc(i),Ha_rowc(i+1)-1
    dydt(i) = dydt(i) + ( Ha_val(j) + (dip_val(j) * Et) ) * y(Ha_row_col(j,1)) / im
  end do
end do
!!$OMP end parallel do


!dydt=dcmplx(0.d0,0.d0)
!!$OMP PARALLEL DO shared(dydt,y,ham,n) 
!do i=0,n-1
!!dydt(i) = dot_product(ham(i,1:n),y(1:n))/im
!  do j=0,n-1
!    dydt(i) = dydt(i) +  ham(i,j) * y(j)/im
!  end do
!end do
!!$OMP end parallel do

end subroutine HA_calc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine momentum_calc_q1(y,mom,n)
use global_param
use omp_lib
implicit none
integer n
complex(kind=dp) :: mom(0:n-1),y(0:n-1)
mom=dcmplx(0.d0,0.d0)
!$OMP PARALLEL DO shared(y,mom,moq1_val,moq1_rowc,moq1_row_col)
do i=0,n-1
  do j=moq1_rowc(i),moq1_rowc(i+1)-1
    mom(i) = mom(i) + moq1_val(j)*ai*(-im) * y(moq1_row_col(j,1))
  end do
end do
!$OMP end parallel do
end subroutine momentum_calc_q1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine momentum_calc_q2(y,mom,n)
use global_param
use omp_lib
implicit none
integer n
complex(kind=dp) :: mom(0:n-1),y(0:n-1)
mom=dcmplx(0.d0,0.d0)
!$OMP PARALLEL DO shared(y,mom,moq2_val,moq2_rowc,moq2_row_col)
do i=0,n-1
  do j=moq2_rowc(i),moq2_rowc(i+1)-1
    mom(i) = mom(i) + moq2_val(j)*bi*(-im) * y(moq2_row_col(j,1))
  end do
end do
!$OMP end parallel do
end subroutine momentum_calc_q2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!subroutine gradient(y,mom,n)
!use global_param
!use omp_lib
!implicit none
!integer n
!complex(kind=dp) :: mom(n),y(n)
!mom=dcmplx(0.d0,0.d0)
!!$OMP PARALLEL DO shared(moq2,y,mom,n) 
!do i=1,n
!  do j=1,n
!mom(i) = dot_product( (moq1(i,j) + moq2(i,j) ) * (-im) * y(j))
!  do j=1,n
!    mom(i) = dydt(i) +  mo(i,j) * y(j) * (-im)
!  end do
!end do
!!$OMP end parallel do
!end subroutine gradient
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!subroutine save_vector_h5(x,n,fname,le)
!use HDF5
!use global_param
!implicit none
!integer(kind=4)    :: n,le
!real(kind=dp)      :: x(n)
!character(len=l)   :: fname         ! File name
!character(len=l-3) :: dsetname      ! dataset name
!integer(HID_T)     :: file_id       ! File identifier
!integer(HID_T)     :: dspace_id     ! Dataspace identifier
!integer(HID_T)     :: dset_id       ! Dataset identifier
!integer(HSIZE_T)   :: dims(1)       ! Dimensions for Dataset and Dataspace
!integer,parameter  :: rank = 1      ! Dataset rank = number of dimensions
!integer            :: error         ! Error flag
!dims=n
!write(dsetname,'(<le-3>a)') fname(1:le-3)
!! Initialize FORTRAN interface.
!call h5open_f(error)
!! Create a new file using default properties.
!call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)
!! Create the dataspace.
!call h5screate_simple_f(rank, dims, dspace_id, error)
!! Create the dataset with default properties.
!call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, dset_id,error)
!! Write the data to datset
!call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, dims, error)
!! End access to the dataset and release resources used by it.
!call h5dclose_f(dset_id, error)
!! Terminate access to the data space.
!call h5sclose_f(dspace_id, error)
!! Close the file.
!call h5fclose_f(file_id, error)
!! Close FORTRAN interface.
!call h5close_f(error)
!end subroutine save_vector_h5
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!subroutine save_matrix_h5(x,m,n,fname,le)
!use HDF5
!use global_param
!implicit none
!integer(kind=4)    :: m,n,le
!real(kind=dp)      :: x(m,n)
!character(len=l)   :: fname         ! File name
!character(len=l-3) :: dsetname      ! dataset name
!integer(HID_T)     :: file_id       ! File identifier
!integer(HID_T)     :: dspace_id     ! Dataspace identifier
!integer(HID_T)     :: dset_id       ! Dataset identifier
!integer(HSIZE_T)   :: dims(2)       ! Dimensions for Dataset and Dataspace
!integer,parameter  :: rank = 2      ! Dataset rank = number of dimensions
!integer            :: error         ! Error flag
!dims(1)=m
!dims(2)=n
!write(dsetname,'(<le-3>a)') fname(1:le-3)
!! Initialize FORTRAN interface.
!call h5open_f(error)
!! Create a new file using default properties.
!call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)
!! Create the dataspace.
!call h5screate_simple_f(rank, dims, dspace_id, error)
!! Create the dataset with default properties.
!call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, dset_id,error)
!! Write the data to datset
!call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, dims, error)
!! End access to the dataset and release resources used by it.
!call h5dclose_f(dset_id, error)
!! Terminate access to the data space.
!call h5sclose_f(dspace_id, error)
!! Close the file.
!call h5fclose_f(file_id, error)
!! Close FORTRAN interface.
!call h5close_f(error)
!end subroutine save_matrix_h5
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine angular_momentum(y,n,am)
use global_param
use omp_lib
implicit none
integer          :: n
complex(kind=dp) :: y(0:n-1),am(0:n-1)

am=dcmplx(0.d0,0.d0)

!$OMP PARALLEL DO shared(y,am,am_rowc,am_val,am_row_col)
do i=0,n-1
  do j=am_rowc(i),am_rowc(i+1)-1
    am(i) = am(i) + (-im) * am_val(j) * y(am_row_col(j,1))
  end do
end do
!$OMP end parallel do
end subroutine angular_momentum
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




