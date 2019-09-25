program dyn
use global_param
use omp_lib
implicit none
external rkdumb,HA_calc,first_derivative_matrix,momentum_calc
integer n,jj !m=dimension of the Hamiltonian (Nq1*Nq2*Nst x Nq1*Nq2*Nst)
real(kind=dp) :: a,e0,k0 !sig=width of gaussian for vec0; soma=variable for sums
complex(kind=dp) :: soma
!character(len=9) :: fmt
integer npoints !number of time steps to take in integration
real(kind=dp) :: t0,tf !t0=initial time for integration; tf=final time for integration
real(kind=dp) :: tt,ch1,ch2,x,x0,kf,w,expo,c0,c1,tstep,start_time,stop_time,ompt0,ompt1,tt1,tt0
complex(kind=dp), allocatable :: vec0(:),vec1(:,:),vec2(:),vec3(:)
real(kind=dp) :: q10,q20 ! point along q1 and q2 where the minimum global C2v is
integer cont
real(kind=dp) :: truni,trunf
complex(kind=dp) :: mom(Nst*Nq1*Nq2),am(Nst*Nq1*Nq2)
real(kind=dp) :: mom_expec_q1,mom_expec_q2,L1,L2,L3
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!$ truni = omp_get_wtime()
n=Nst*Nq1*Nq2 ! m=Nq1*Nq2*st %size of the Hamiltonian
!$ ompt0 = omp_get_wtime()
call hamiltonian_matrix
allocate (moq1(Nst*Nq1*Nq2,Nst*Nq1*Nq2))
allocate (moq2(Nst*Nq1*Nq2,Nst*Nq1*Nq2))
!$OMP parallel do shared(moq1,ham)
do i=1,n
  do j=1,n
    moq1(i,j)=0.d0
    ham(i,j) =0.d0
  end do 
end do
!$OMP end parallel do
!$OMP parallel do shared(moq2)
do i=1,n
  do j=1,n
    moq2(i,j)=0.d0
  end do 
end do
!$OMP end parallel do
call first_derivative_matrix_q1
call first_derivative_matrix_q2
call nac_ha_modify
!$OMP parallel do shared(ha,ham)
do i=1,n
  do j=1,n
    ham(i,j)=ha(i,j)
  end do
end do
!$OMP end parallel do
!$ ompt1 = omp_get_wtime()
write(*,*) 'The time to load the matrices=',ompt1 - ompt0, "seconds"

im = dcmplx(0.d0,1.d0)
n=Nst*Nq1*Nq2 ! m=Nq1*Nq2*st %size of the Hamiltonian
allocate ( vec0(Nst*Nq1*Nq2) )
allocate ( vec1(Nq1,Nq2)     )
allocate ( vec2(Nst*Nq1*Nq2) )
allocate ( vec3(Nq1*Nq2)     )
vec0=0.d0
vec1=0.d0
vec2=0.d0
vec3=0.d0

!! Follow mean value of momentum (k) through time and check if is always smaller than the maximum momentum from sampling


!############################################################################################################!#
! Buildin initial wave packet in 2D as the one of an Harmonic oscilator                                      !# 
t0=0.d0                                                                                                      !#      
tf=2000.0d0               !40 = Approx. 1 femtosecond                                                         !#  
npoints=10000             !If 1000, the step will be approx. 1 attosecond                                     !#
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
do j=1,Nq2                                                                                                   !#
  do i=1,Nq1                                                                                                 !#
    cont=cont+1;                                                                                             !#
    ch1=(i-Nq1/2.d0-1.d0)*sq1+sq1/2.d0 !change in q1                                                         !#
    ch2=(j-(119.d0)-1.d0)*sq2 !change in q2                                                                  !#
!    expo = dexp(-(ch-q10)**2.d0*mtotal*w/2.d0) !* exp(im*k0*(ch))                                           !#
    expo = dexp((-(ch1-q10)**2.d0*mtotal*w/2.d0)+(-(ch2-q20)**2.d0*mtotal*w/2.d0)) !* exp(im*k0*(ch))        !#
    vec1(i,j) = c0 * expo                                                                                    !#
    vec0(cont)=vec1(i,j)                                                                                     !#
!    vec0(cont+1*s)=vec1(i,j)                                                                                 !#
!    vec0(cont+2*s)=vec1(i,j)                                                                                 !#
    vec3(cont)=vec0(cont)                                                                                    !#
!    vec0(i) = (c0 * expo + c1 * ch * expo) / sqrt(2.d0)                                                     !#
!    vec0(i) = c0 * (2.0d0*mtotal*w*ch**2.d0-1.d0) * expo                                                    !#
!    vec0(i+Nq1) = 0.0d0                                                                                     !#
!    vec0(i+2*Nq1) = 0.0d0                                                                                   !#
  end do                                                                                                     !#
end do                                                                                                       !# 
!############################################################################################################!#
!--------------------------------------------!
!normalizing                                 !
soma=0.0d0                                   !
do i=1,n                                     !
  soma=soma+dconjg(vec0(i))*vec0(i)          !
end do                                       !
vec0=vec0/sqrt(soma)                         !
!--------------------------------------------!
write(*,'(a24,2(es15.7e3))')'normalization constant =',soma

!$ tt0 = omp_get_wtime()
call angular_momentum(vec0,n,am)
!$ tt1 = omp_get_wtime()

!--------------------------------------------------------!
!Checking initial energy and momentum                    !
!$ ompt0 = omp_get_wtime()
call HA_calc(tt,vec0,n,vec2,tstep) !evaluating y'(t=0,y) !
!$ ompt1 = omp_get_wtime()
!$ start_time = omp_get_wtime()
call momentum_calc_q1(vec0,mom,n,) ! evaluating dydq     !
!$ stop_time = omp_get_wtime()
soma=0.0d0                                               !
mom_expec_q1=0.0d0                                       !
mom_expec_q2=0.0d0                                       !
do i=1,s !n !!!! MUDEIIIIII                              !
  soma=soma + dconjg(vec0(i)) * vec2(i)*im               !
  mom_expec_q1=mom_expec_q1 + dconjg(vec0(i)) * mom(i)   !
  L1    = L1    + dconjg(vec0(i))     * am(i)
  L2    = L2    + dconjg(vec0(i+2))   * am(i+2)
  L3    = L3    + dconjg(vec0(i+2*s)) * am(i+2*s)
end do                                                   !
call momentum_calc_q2(vec0,mom,n) ! evaluating dydq      !
do i=1,s !n !!!! MUDEIIIIII                              !
  mom_expec_q2=mom_expec_q2 + dconjg(vec0(i)) * mom(i)   !
end do                                                   !
write(*,'(a24,e23.15e3,a8)')'inital momentum in q1 =',mom_expec_q1,' au.'
write(*,'(a24,e23.15e3,a8)')'inital momentum in q2 =',mom_expec_q2,' au.'
!--------------------------------------------------------!
write(*,*)'##############################################################################'
write(*,*)'The time the program took to do the operation H|Psi> is:',(ompt1-ompt0), "seconds"
write(*,*)'Be ready to wait probably ', npoints*(tt1-tt0+(stop_time-start_time)*2.d0+(ompt1-ompt0)*5.d0)/(3600.d0) , "hours"
write(*,*)'##############################################################################'
write(*,'(a24,3(f23.15))')'initial angular momentum = ',L1,L2,L3,L1+L2+L3
write(*,'(a24,2(f23.15),a8)')'inital energy =',soma,' hartree'
write(*,'(a24,(f23.15))')'initial 1 / dq1 =',1.d0/sq1
write(*,'(a24,(f23.15))')'initial 1 / dq2 =',1.d0/sq2
write(*,'(a24,(f23.15))')'initial k * dq =',dsqrt(2*mtotal*e0)*sq1*sq2 
write(*,'(a24,(f23.15),a4)')'vibrational frequency =',w,' au.'
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! nanoseconds = 1d-9 seconds
! picoseconds = 1d-12 seconds
! femtoseconds = 1d-15 seconds
! attoseconds = 1d-18 seconds
! 1 time atomic units = 24.18884326505 attoseconds
call rkdumb(vec0,n,t0,tf,npoints,tstep,HA_calc)
!$ trunf = omp_get_wtime()
write(*,*)'Time took to run totaly the program = ', (trunf-truni)/3600.d0, 'hours.'
write(*,*)'************************************************************'
write(*,*)'Program finished'
write(*,*)'************************************************************'
end program
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!subroutine rkdumb(y0,n,t0,tf,npoints,h,HA_calc)
subroutine rkdumb(y,n,t0,tf,npoints,h,HA_calc)
use global_param
external HA_calc,rk4,momentum_calc,angular_momentum !
integer           :: n,npoints,ii
real(kind=dp)     :: t0,tf,tt(npoints),h,t,coordinate1(s),coordinate2(s)
complex(kind=dp)  :: dydt(n),y(n)!,y0(n)!,ymatrix(n,npoints+1)!,ymatlab(n,npoints)
real(kind=dp)     :: soma,Te,aux1,aux2,Te1,Te2,Te3,enert(npoints+1) ! This is fo saving total energy through time
real(kind=dp)     :: pq1_1,pq2_1,pq1_2,pq2_2,pq1_3,pq2_3 !This is for saving momentum through time
complex(kind=dp)  :: momq1(n),momq2(n),am(n) ! This is for saving momentum through time
real(kind=dp)     :: sum1,sum2,sum3,momq1t(npoints),momq2t(npoints) ! This is for saving norm through time
real(kind=dp)     :: e1,e2,e3,L1,L2,L3
character(len=23) :: fname
real(kind=dp)     :: st1n(npoints),st1p(npoints),st2n(npoints),st2p(npoints),st3n(npoints),st3p(npoints)
real(kind=dp)     :: st1ncount,st1pcount,st2ncount,st2pcount,st3ncount,st3pcount,pop1,pop2,pop3


open(unit=10,file='norms-time.data',status='unknown')
open(unit=11,file='momentum-time.data',status='unknown')
open(unit=12,file='electric-field-time.data',status='unknown')
open(unit=13,file='energy-time.data',status='unknown')
open(unit=14,file='ang-mom.data',status='unknown')
open(unit=15,file='alldata.data',status='unknown')

!do1: do i=1,n
!  y(i)=y0(i) !loading initial conditions
!  ymatrix(i,1)=y(i) !saving initial amplitudes
!end do do1
tt(1)=t0 !saving initial time value
t=t0
write(*,'(a27,f23.15,a4)')'time step =',h,' au.'
write(*,'(a27,f23.15)')'step size in coordinate 1 =',sq1
write(*,'(a27,f23.15)')'step size in coordinate 2 =',sq2

call HA_calc(t,y,n,dydt,h) ! evaluating y'(t,y)
!-----------------------------------------!
!Checking norm and conservation of energy !
soma=0.0d0                                !
Te=0.d0                                   !
do i=1,n                                  !
  soma=soma+dconjg(y(i))*y(i)             !
  Te=Te + dconjg(y(i)) * dydt(i)*im       !
end do                                    !
!-----------------------------------------!
aux1=soma
aux2=Te
!nst1(1)=1.d0
!nst2(1)=0.d0
!nst3(1)=0.d0
enert(1)=Te

write(*,*) '************************************************************'
write(*,*) 'Integrating amplitudes over time'
write(*,*) '************************************************************'

do3: do k=1,npoints
!WRITE(6,'(5(A))',ADVANCE="NO") "\b","\b","\b","\b","b"
  call HA_calc(t,y,n,dydt,h) ! evaluating y'(t,y)
  call rk4(y,dydt,n,t,h,y,HA_calc) !evaluating y(t+h)
!  do2: do i=1,n
!    ymatrix(i,k+1)=y(i) !saving the amplitude for time t
!  end do do2
!---------------------------------------------!
!Saving norm and total energy in time         !
soma=0.0d0                                    !
Te=0.d0                                       !
do i=1,n                                      !
  soma=soma+dconjg(y(i))*y(i)                 !
  Te=Te + dconjg(y(i)) * dydt(i)*im           !
end do                                        !
!---------------------------------------------!
!--------------------------------------------------------!
!Checking momentum and saving norm                       !
call momentum_calc_q1(y,momq1,n) ! evaluating dydq       !
call momentum_calc_q2(y,momq2,n) ! evaluating dydq       !
call angular_momentum(y,n,am)
pq1_1=0.0d0 ! momentum in q1 state1                      !
pq2_1=0.0d0 ! momentum in q2 state1                      !
pq1_2=0.0d0 ! momentum in q1 state2                      !
pq2_2=0.0d0 ! momentum in q2 state2                      !
pq1_3=0.0d0 ! momentum in q1 state3                      !
pq2_3=0.0d0 ! momentum in q2 state3                      !
sum1=0.0d0  ! norm for state 1                           !
sum2=0.0d0  ! norm for state 2                           !
sum3=0.0d0  ! norm for state 3                           !
e1=0.d0
e2=0.d0
e3=0.d0
L1=0.d0
L2=0.d0
L3=0.d0
do i=1,s !n !!!! MUDEIIIIII                              !
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
  L1    = L1    + dconjg(y(i))     * am(i)
  L2    = L2    + dconjg(y(i+2))   * am(i+2)
  L3    = L3    + dconjg(y(i+2*s)) * am(i+2*s)
end do                                                   !
!--------------------------------------------------------!
!normt(k+1) =soma
enert(k+1) =Te
momq1t(k)=mom_expec_q1
momq2t(k)=mom_expec_q2
!nst1(k+1)=sum1
!nst2(k+1)=sum2
!nst3(k+1)=sum3

!-----------------------------------------------------!
! saving the population at the borders for each state !
st1n(k)=dconjg(y(1))*y(1)                             !
st1p(k)=dconjg(y(s-54))*y(s-54)                         !
st2n(k)=dconjg(y(s+1))*y(s+1)                     !
st2p(k)=dconjg(y(s+s-54))*y(s+s-54)                 !
st3n(k)=dconjg(y(2*s+1))*y(2*s+1)                 !
st3p(k)=dconjg(y(n-54))*y(n-54)                             !
!-----------------------------------------------------!

!fname='time-amp-real-000000.h5'
!write(fname(15:20),'(i0.6)') k
!call save_vector_h5(real(y),n,fname,23)
!fname='time-amp-imag-000000.h5'
!write(fname(15:20),'(i0.6)') k
!call save_vector_h5(aimag(y),n,fname,23)


  fname='amp-time-st1-000000.dat'
  write(fname(14:19),'(i0.6)') k
  open(unit=20,file=fname,status='unknown')
  fname='amp-real-st1-000000.dat'
  write(fname(14:19),'(i0.6)') k
  open(unit=30,file=fname,status='unknown')
  fname='amp-imag-st1-000000.dat'
  write(fname(14:19),'(i0.6)') k
  open(unit=40,file=fname,status='unknown')
  fname='amp-time-st2-000000.dat'
  write(fname(14:19),'(i0.6)') k
  open(unit=50,file=fname,status='unknown')
  fname='amp-real-st2-000000.dat'
  write(fname(14:19),'(i0.6)') k
  open(unit=60,file=fname,status='unknown')
  fname='amp-imag-st2-000000.dat'
  write(fname(14:19),'(i0.6)') k
  open(unit=70,file=fname,status='unknown')
  fname='amp-time-st3-000000.dat'
  write(fname(14:19),'(i0.6)') k
  open(unit=80,file=fname,status='unknown')
  fname='amp-real-st3-000000.dat'
  write(fname(14:19),'(i0.6)') k
  open(unit=90,file=fname,status='unknown')
  fname='amp-imag-st3-000000.dat'
  write(fname(14:19),'(i0.6)') k
  open(unit=99,file=fname,status='unknown')
  do i=1,s
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
  write(10,'(3(es26.16e3))') sum1,sum2,sum3
  write(11,'(6(es26.16e3))') pq1_1,pq2_1,pq1_2,pq2_2,pq1_3,pq2_3
Et = (E00/freq)*(-(t-t00)/sig**2.d0*sin(freq*(t-t00)+phase)+freq*cos(freq*(t-t00)+phase))*dexp(-(t-t00)**2.d0/(2.d0*sig**2.d0))
  write(12,'((es26.16e3))') Et
  write(13,'(3(es26.16e3))') e1,e2,e3
  write(14,'(3(es26.16e3))') L1,L2,L3
  write(15,'(7(es26.16e3))') t,Et,sum1,soma,Te,sum2,sum3,L1,L2,L3
  t=t+h
  tt(k)=t !saving time value for each step in time
  write(*,*) 'time, Pulse, E1, E2, E3, ET, norm1, norm2, norm3, NormT,Ltot   '
  write(*,'(i5,11(e12.3e3))') k,Et,e1,e2,e3,Te,sum1,sum2,sum3,sum1+sum2+sum3,L1+L2+L3
  !WRITE(6,'(5(A))',ADVANCE="NO") "\b","\b","\b","\b","b"
  !write(6,'(f5.1,"%")',advance='no') k/10.d0
  !flush(6)
end do do3

write(*,*)''
 close(unit=10)
 close(unit=11)
 close(unit=12)
 close(unit=13)
 close(unit=14)
 close(unit=15)
!---------------------------------------------------------------------------------------------------------------------------!
! Determining the highest population at the simulation                                                                      !
! must be always small, > 10d-8. If not, consider changing grid borders and check the potential energy behaviour at borders !
!st1ncount=1.d-100; st1pcount=1.d-100; st2ncount=1.d-100; st2pcount=1.d-100; st3ncount=1.d-100; st3pcount=1.d-100;           !
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
sum1=0.0d0                                !
sum2=0.0d0                                !
sum3=0.0d0                                !
Te1=0.d0                                  !
Te2=0.d0                                  !
Te3=0.d0                                  !
call HA_calc(t,y,n,dydt,h)                !
do i=1,n                                  !
  soma=soma+dconjg(y(i))*y(i)             !
end do                                    !
do i=1,s                                  !
  sum1=sum1+dconjg(y(i))*y(i)             !
  sum2=sum2+dconjg(y(i+s))*y(i+s)         !
  sum3=sum3+dconjg(y(i+2*s))*y(i+2*s)     !
  Te1=Te1+dconjg(y(i))    *dydt(i)*im     !
  Te2=Te2+dconjg(y(i+s))  *dydt(i+s)*im   !
  Te3=Te3+dconjg(y(i+2*s))*dydt(i+2*s)*im !
end do                                    !
!-----------------------------------------!
write(*,'(a30,(e12.3e3))')'Norm conservation =',aux1-soma
write(*,'(a30,(e12.3e3))')'Final norm for state 1 =',sum1
write(*,'(a30,(e12.3e3))')'Final norm for state 2 =',sum2
write(*,'(a30,(e12.3e3))')'Final norm for state 3 =',sum3
write(*,'(a30,(e12.3e3),a8)')'Energy conservation =',aux2-Te,' hartree'
write(*,'(a30,(f20.15),a8)')'Final total energy state 1 =',Te1,' hartree'
write(*,'(a30,(f20.15),a8)')'Final total energy state 2 =',Te2,' hartree'
write(*,'(a30,(f20.15),a8)')'Final total energy state 3 =',Te3,' hartree'
!write(*,*)'************************************************************'
!write(*,*)'Saving amplitude over time'
!write(*,*)'************************************************************'
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
!write(*,*) 'Time to save final files= ', (trun1-trun0)/(16.d0*60.d0), 'minutes.'
write(*,'(a24,e12.3e3)') 'maximum momentum value =',maxval(momq1t)
write(*,'(a24,e12.3e3)') 'maximum momentum value =',maxval(momq2t)
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

end subroutine rkdumb
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rk4(y0,dydt,n,t,h,Yout,HA_calc)
use global_param
implicit none
integer n
real(kind=dp) :: t,h,th
complex(kind=dp) :: dydt(n),y0(n),Yout(n),dym(n),dyt(n),yt(n)
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
do1: do i=1,n
  yt(i) = y0(i) + (h/2.d0) * dydt(i) !Evaluation of y(t) + k1/2.0d0
end do do1
call HA_calc(th,yt,n,dyt,h/2.d0) !dyt = k2/h ---- Evaluation of y'(t + h/2.0d0 , y(t) + k1/2.0d0)
do2: do i=1,n
  yt(i) = y0(i) + (h/2.d0) * dyt(i)  !Evaluation of y(t) + k2/2.0d0
end do do2
call HA_calc(th,yt,n,dym,h/2.d0) !dym = k3/h ---- Evaluation of y'(t + h/2.0d0 , y(t) + k2/2.0d0)
do3: do i=1,n
  yt(i) = y0(i) + h * dym(i) !Evaluation of y(t) + k3
  dym(i) = dyt(i) + dym(i)  ! making dym = (k2 + k3)/h ---- variable dyt free now
end do do3
call HA_calc(t+h,yt,n,dyt,h*2.d0) !dyt = k4/h ---- Evaluation of y'(t + h/2.0d0 , y(t) + k1/2.0d0)
do4: do i=1,n
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
integer n,ii
real(kind=dp) :: t,tstep!,c1q1,c2q1,c3q1,c1q2,c2q2,c3q2,c1q3,c2q3,c3q3
complex(kind=dp) :: y(n),dydt(n)

Et = (E00/freq)*(-(t-t00)/sig**2.d0*sin(freq*(t-t00)+phase)+freq*cos(freq*(t-t00)+phase))*dexp(-(t-t00)**2.d0/(2.d0*sig**2.d0))

ii=0
do j=1,Nq2
  do i=1,Nq1
    ii=ii+1
    ham(ii    ,ii    ) = Ha(ii    ,ii    ) - (pdm1x(i,j)*Et)
    ham(s+ii  ,s+ii  ) = Ha(s+ii  ,s+ii  ) - (pdm2x(i,j)*Et)
    ham(2*s+ii,2*s+ii) = Ha(2*s+ii,2*s+ii) - (pdm3x(i,j)*Et)

    ham(ii    ,s+ii  ) = Ha(ii    ,s+ii  ) - tdm21x(i,j) * Et
    ham(s+ii  ,ii    ) = ha(s+ii  ,ii    ) - tdm21x(i,j) * Et
    ham(ii    ,2*s+ii) = ha(ii    ,2*s+ii) - tdm31x(i,j) * Et
    ham(2*s+ii,ii    ) = ha(2*s+ii,ii    ) - tdm31x(i,j) * Et
    ham(s+ii  ,2*s+ii) = ha(s+ii  ,2*s+ii) - tdm32x(i,j) * Et
    ham(2*s+ii,s+ii  ) = ha(2*s+ii,s+ii  ) - tdm32x(i,j) * Et
  end do
end do

!write(*,*) 'omp_get_max_threads= ', omp_get_max_threads ( )
!write(*,*) 'omp_get_num_procs = ', omp_get_num_procs ( )
!write(*,*) 'Time = ', omp_get_wtime ( )
dydt=dcmplx(0.d0,0.d0)
!$OMP PARALLEL DO shared(dydt,y,ham,im,n) 
do i=1,n
!dydt(i) = dot_product(ham(i,1:n),y(1:n))/im
  do j=1,n
    dydt(i) = dydt(i) +  ham(i,j) * y(j)/im
  end do
end do
!$OMP end parallel do

end subroutine HA_calc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine first_derivative_matrix_q1
use global_param
implicit none
!=================================================================================================================================!
!Building the Momentun matrix                                                                                                     !
const1 = 1.d0/(60.d0*sq1)                                                                                                         !
!----------------------------------------------------------------------------------------------------------------------------|    !
!Doing the derivative along q1                                                                                               |    !
ind1= 45.d0*const1    !factor for 5 point first derivative finite diference in 1D along q1                                   |    !
ind2=  9.d0*const1    !factor for 5 point first derivative finite diference in 1D along q1                                   |    !
ind3=  1.d0*const1    !factor for 5 point first derivative finite diference in 1D along q1                                   |    !
!----------------------------------------------------------------------------------------------------------------------------|    !
do k=0,2*s,s ! running for all electronic states                                                                             |    !
  do j=0*Nq1,(Nq2-1)*Nq1,Nq1 ! Moving through the boxes, from the fisrt to the last one                                      |    !
!.....................................................................................................................!      |    !
! Doing the special cases of the 5 points finite difference that happens in the borders of each q1 sized box          !      |    !
! First line of each box                                                                                              !      |    !
moq1( k   + j       + 1     , k   + j       + 2     ) = moq1( k   + j       + 1     , k   + j       + 2     ) + ind1  !      |    !
moq1( k   + j       + 1     , k   + j       + 3     ) = moq1( k   + j       + 1     , k   + j       + 3     ) - ind2  !      |    !
moq1( k   + j       + 1     , k   + j       + 4     ) = moq1( k   + j       + 1     , k   + j       + 4     ) + ind3  !      |    !
! Second line of each box                                                                                             !      |    !
moq1( k   + j       + 2     , k   + j       + 1     ) = moq1( k   + j       + 2     , k   + j       + 1     ) - ind1  !      |    !
moq1( k   + j       + 2     , k   + j       + 3     ) = moq1( k   + j       + 2     , k   + j       + 3     ) + ind1  !      |    !
moq1( k   + j       + 2     , k   + j       + 4     ) = moq1( k   + j       + 2     , k   + j       + 4     ) - ind2  !      |    !
moq1( k   + j       + 2     , k   + j       + 5     ) = moq1( k   + j       + 2     , k   + j       + 5     ) + ind3  !      |    !
! Third line of each box                                                                                              !      |    !
moq1( k   + j       + 3     , k   + j       + 1     ) = moq1( k   + j       + 3     , k   + j       + 1     ) + ind2  !      |    !
moq1( k   + j       + 3     , k   + j       + 2     ) = moq1( k   + j       + 3     , k   + j       + 2     ) - ind1  !      |    !
moq1( k   + j       + 3     , k   + j       + 4     ) = moq1( k   + j       + 3     , k   + j       + 4     ) + ind1  !      |    !
moq1( k   + j       + 3     , k   + j       + 5     ) = moq1( k   + j       + 3     , k   + j       + 5     ) - ind2  !      |    !
moq1( k   + j       + 3     , k   + j       + 6     ) = moq1( k   + j       + 3     , k   + j       + 6     ) + ind3  !      |    !
! Antepenulmate line of each box                                                                                      !      |    !
moq1( k   + j       + Nq1-2 , k   + j       + Nq1-5 ) = moq1( k   + j       + Nq1-2 , k   + j       + Nq1-5 ) - ind3  !      |    !
moq1( k   + j       + Nq1-2 , k   + j       + Nq1-4 ) = moq1( k   + j       + Nq1-2 , k   + j       + Nq1-4 ) + ind2  !      |    !
moq1( k   + j       + Nq1-2 , k   + j       + Nq1-3 ) = moq1( k   + j       + Nq1-2 , k   + j       + Nq1-3 ) - ind1  !      |    !
moq1( k   + j       + Nq1-2 , k   + j       + Nq1-1 ) = moq1( k   + j       + Nq1-2 , k   + j       + Nq1-1 ) + ind1  !      |    !
moq1( k   + j       + Nq1-2 , k   + j       + Nq1   ) = moq1( k   + j       + Nq1-2 , k   + j       + Nq1   ) - ind2  !      |    !
! Penulmate line of each box                                                                                          !      |    !
moq1( k   + j       + Nq1-1 , k   + j       + Nq1-4 ) = moq1( k   + j       + Nq1-1 , k   + j       + Nq1-4 ) - ind3  !      |    !
moq1( k   + j       + Nq1-1 , k   + j       + Nq1-3 ) = moq1( k   + j       + Nq1-1 , k   + j       + Nq1-3 ) + ind2  !      |    !
moq1( k   + j       + Nq1-1 , k   + j       + Nq1-2 ) = moq1( k   + j       + Nq1-1 , k   + j       + Nq1-2 ) - ind1  !      |    !
moq1( k   + j       + Nq1-1 , k   + j       + Nq1   ) = moq1( k   + j       + Nq1-1 , k   + j       + Nq1   ) + ind1  !      |    !
! Last line of each box                                                                                               !      |    !
moq1( k   + j       + Nq1   , k   + j       + Nq1-3 ) = moq1( k   + j       + Nq1   , k   + j       + Nq1-3 ) - ind3  !      |    !
moq1( k   + j       + Nq1   , k   + j       + Nq1-2 ) = moq1( k   + j       + Nq1   , k   + j       + Nq1-2 ) + ind2  !      |    !
moq1( k   + j       + Nq1   , k   + j       + Nq1-1 ) = moq1( k   + j       + Nq1   , k   + j       + Nq1-1 ) - ind1  !      |    !
!.....................................................................................................................!      |    !
!________________________________________________________________________________________________________________________!   |    !
! From the fourth line to the ante-antepenulmate one of each box - non special cases                                     !   |    !
    do i=1+3,Nq1-3 ! do through single box                                                                               !   |    !
moq1( k   + j       + i     , k   + j       + i-3   ) = moq1( k   + j       + i     , k   + j       + i-3   ) - ind3     !   |    !
moq1( k   + j       + i     , k   + j       + i-2   ) = moq1( k   + j       + i     , k   + j       + i-2   ) + ind2     !   |    !
moq1( k   + j       + i     , k   + j       + i-1   ) = moq1( k   + j       + i     , k   + j       + i-1   ) - ind1     !   |    !
moq1( k   + j       + i     , k   + j       + i+1   ) = moq1( k   + j       + i     , k   + j       + i+1   ) + ind1     !   |    !
moq1( k   + j       + i     , k   + j       + i+2   ) = moq1( k   + j       + i     , k   + j       + i+2   ) - ind2     !   |    !
moq1( k   + j       + i     , k   + j       + i+3   ) = moq1( k   + j       + i     , k   + j       + i+3   ) + ind3     !   |    !
    end do ! end of do through single box                                                                                !   |    !
!________________________________________________________________________________________________________________________!   |    !
  end do ! end of do through boxes                                                                                           |    !
end do ! end of do through electronic states                                                                                 |    !
!----------------------------------------------------------------------------------------------------------------------------|    !
!=================================================================================================================================!
end subroutine first_derivative_matrix_q1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine first_derivative_matrix_q2
use global_param
implicit none
!=================================================================================================================================!
const2 = 1.d0/(60.d0*sq2)                                                                                                         !
!-------------------------------------------------------------------------------------------------------------------------|       !
! Doing derivative along q2                                                                                               |       !
ind1= 45.d0*const2    !factor for 5 point first derivative finite diference in 1D along q1                                |       !
ind2=  9.d0*const2    !factor for 5 point first derivative finite diference in 1D along q1                                |       !
ind3=  1.d0*const2    !factor for 5 point first derivative finite diference in 1D along q1                                |       !
!-------------------------------------------------------------------------------------------------------------------------|       !
do k=0,2*s,s ! Running through every state, the step size is Nq1*Nq2                                                      |       !
!....................................................................................................................!    |       !
! Doing the special cases of the 5 points finite difference that happens in first and last three boxes               !    |       !
  do i=1,Nq1 ! Moving inside the box                                                                                 !    |       !
! Upper and left part of the borders - first box of all                                                              !    |       !
moq2( k   +   0*Nq1 + i     , k   +   1*Nq1 + i     ) = moq2( k   +   0*Nq1 + i     , k   +   1*Nq1 + i     ) + ind1 !    |       !
moq2( k   +   0*Nq1 + i     , k   +   2*Nq1 + i     ) = moq2( k   +   0*Nq1 + i     , k   +   2*Nq1 + i     ) - ind2 !    |       !
moq2( k   +   0*Nq1 + i     , k   +   3*Nq1 + i     ) = moq2( k   +   0*Nq1 + i     , k   +   3*Nq1 + i     ) + ind3 !    |       !
! Upper and left part of the borders - second box of all                                                             !    |       !
moq2( k   +   1*Nq1 + i     , k   +   0*Nq1 + i     ) = moq2( k   +   1*Nq1 + i     , k   +   0*Nq1 + i     ) - ind1 !    |       !
moq2( k   +   1*Nq1 + i     , k   +   2*Nq1 + i     ) = moq2( k   +   1*Nq1 + i     , k   +   2*Nq1 + i     ) + ind1 !    |       !
moq2( k   +   1*Nq1 + i     , k   +   3*Nq1 + i     ) = moq2( k   +   1*Nq1 + i     , k   +   3*Nq1 + i     ) - ind2 !    |       !
moq2( k   +   1*Nq1 + i     , k   +   4*Nq1 + i     ) = moq2( k   +   1*Nq1 + i     , k   +   4*Nq1 + i     ) + ind3 !    |       !
! Upper and left part of the borders - third box of all                                                              !    |       !
moq2( k   +   2*Nq1 + i     , k   +   0*Nq1 + i     ) = moq2( k   +   2*Nq1 + i     , k   +   0*Nq1 + i     ) + ind2 !    |       !
moq2( k   +   2*Nq1 + i     , k   +   1*Nq1 + i     ) = moq2( k   +   2*Nq1 + i     , k   +   1*Nq1 + i     ) - ind1 !    |       !
moq2( k   +   2*Nq1 + i     , k   +   3*Nq1 + i     ) = moq2( k   +   2*Nq1 + i     , k   +   3*Nq1 + i     ) + ind1 !    |       !
moq2( k   +   2*Nq1 + i     , k   +   4*Nq1 + i     ) = moq2( k   +   2*Nq1 + i     , k   +   4*Nq1 + i     ) - ind2 !    |       !
moq2( k   +   2*Nq1 + i     , k   +   5*Nq1 + i     ) = moq2( k   +   2*Nq1 + i     , k   +   5*Nq1 + i     ) + ind3 !    |       !
! Botton and right part of the borders, last box of all                                                              !    |       !
moq2( k+s -   0*Nq1 - Nq1+i , k+s -   3*Nq1 - Nq1+i ) = moq2( k+s -   0*Nq1 - Nq1+i , k+s -   3*Nq1 - Nq1+i ) - ind3 !    |       !
moq2( k+s -   0*Nq1 - Nq1+i , k+s -   2*Nq1 - Nq1+i ) = moq2( k+s -   0*Nq1 - Nq1+i , k+s -   2*Nq1 - Nq1+i ) + ind2 !    |       !
moq2( k+s -   0*Nq1 - Nq1+i , k+s -   1*Nq1 - Nq1+i ) = moq2( k+s -   0*Nq1 - Nq1+i , k+s -   1*Nq1 - Nq1+i ) - ind1 !    |       !
! Botton and right part of the borders, penulmate box of all                                                         !    |       !
moq2( k+s -   1*Nq1 - Nq1+i , k+s -   4*Nq1 - Nq1+i ) = moq2( k+s -   1*Nq1 - Nq1+i , k+s -   4*Nq1 - Nq1+i ) - ind3 !    |       !
moq2( k+s -   1*Nq1 - Nq1+i , k+s -   3*Nq1 - Nq1+i ) = moq2( k+s -   1*Nq1 - Nq1+i , k+s -   3*Nq1 - Nq1+i ) + ind2 !    |       !
moq2( k+s -   1*Nq1 - Nq1+i , k+s -   2*Nq1 - Nq1+i ) = moq2( k+s -   1*Nq1 - Nq1+i , k+s -   2*Nq1 - Nq1+i ) - ind1 !    |       !
moq2( k+s -   1*Nq1 - Nq1+i , k+s -   0*Nq1 - Nq1+i ) = moq2( k+s -   1*Nq1 - Nq1+i , k+s -   0*Nq1 - Nq1+i ) + ind1 !    |       !
! Botton and right part of the borders, antepenulmate box of all                                                     !    |       !
moq2( k+s -   2*Nq1 - Nq1+i , k+s -   5*Nq1 - Nq1+i ) = moq2( k+s -   2*Nq1 - Nq1+i , k+s -   5*Nq1 - Nq1+i ) - ind3 !    |       !
moq2( k+s -   2*Nq1 - Nq1+i , k+s -   4*Nq1 - Nq1+i ) = moq2( k+s -   2*Nq1 - Nq1+i , k+s -   4*Nq1 - Nq1+i ) + ind2 !    |       !
moq2( k+s -   2*Nq1 - Nq1+i , k+s -   3*Nq1 - Nq1+i ) = moq2( k+s -   2*Nq1 - Nq1+i , k+s -   3*Nq1 - Nq1+i ) - ind1 !    |       !
moq2( k+s -   2*Nq1 - Nq1+i , k+s -   1*Nq1 - Nq1+i ) = moq2( k+s -   2*Nq1 - Nq1+i , k+s -   1*Nq1 - Nq1+i ) + ind1 !    |       !
moq2( k+s -   2*Nq1 - Nq1+i , k+s -   0*Nq1 - Nq1+i ) = moq2( k+s -   2*Nq1 - Nq1+i , k+s -   0*Nq1 - Nq1+i ) + ind1 !    |       !
  end do ! end of do through single box                                                                              !    |       !
!....................................................................................................................!    |       !
!_______________________________________________________________________________________________________________________! |       !
! Now doing the rest of the boxes                                                                                       ! |       !
  do j=3*Nq1,(Nq2-4)*Nq1,Nq1 ! Moving through the boxes, from the fourth to the ante-antepenulmate one                  ! |       !
    do i=1,Nq1 ! Moving inside the box                                                                                  ! |       !
moq2( k   + j       + i     , k   + j-3*Nq1 + i     ) = moq2( k   + j       + i     , k   + j-3*Nq1 + i     ) - ind3    ! |       !
moq2( k   + j       + i     , k   + j-2*Nq1 + i     ) = moq2( k   + j       + i     , k   + j-2*Nq1 + i     ) + ind2    ! |       !
moq2( k   + j       + i     , k   + j-1*Nq1 + i     ) = moq2( k   + j       + i     , k   + j-1*Nq1 + i     ) - ind1    ! |       !
moq2( k   + j       + i     , k   + j+1*Nq1 + i     ) = moq2( k   + j       + i     , k   + j+1*Nq1 + i     ) + ind1    ! |       !
moq2( k   + j       + i     , k   + j+2*Nq1 + i     ) = moq2( k   + j       + i     , k   + j+2*Nq1 + i     ) - ind2    ! |       !
moq2( k   + j       + i     , k   + j+3*Nq1 + i     ) = moq2( k   + j       + i     , k   + j+3*Nq1 + i     ) + ind3    ! |       !
    end do ! end of do through single box                                                                               ! |       !
  end do ! end of do through boxes                                                                                      ! |       !
!_______________________________________________________________________________________________________________________! |       !
end do ! end of do through electronic states                                                                              |       !
!-------------------------------------------------------------------------------------------------------------------------|       !
!=================================================================================================================================!
end subroutine first_derivative_matrix_q2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine momentum_calc_q1(y,mom,n)
use global_param
use omp_lib
implicit none
integer n
complex(kind=dp) :: mom(n),y(n)
mom=dcmplx(0.d0,0.d0)
!$OMP PARALLEL DO shared(moq1,y,mom,im,n) 
do i=1,n
mom(i) = dot_product(moq1(i,1:n)*ai*(-im),y(1:n))
!  do j=1,n
!    mom(i) = dydt(i) +  mo(i,j) * y(j) * (-im)
!  end do
end do
!$OMP end parallel do
end subroutine momentum_calc_q1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine momentum_calc_q2(y,mom,n)
use global_param
use omp_lib
implicit none
integer n
complex(kind=dp) :: mom(n),y(n)
mom=dcmplx(0.d0,0.d0)
!$OMP PARALLEL DO shared(moq2,y,mom,im,n) 
do i=1,n
mom(i) = dot_product(moq2(i,1:n)*bi*(-im),y(1:n))
!  do j=1,n
!    mom(i) = dydt(i) +  mo(i,j) * y(j) * (-im)
!  end do
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
!!$OMP PARALLEL DO shared(moq2,y,mom,im,n) 
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
subroutine save_vector_h5(x,n,fname,l)
use HDF5
use global_param
implicit none
integer(kind=4)    :: n,l
real(kind=dp)      :: x(n)
character(len=l)   :: fname         ! File name
character(len=l-3) :: dsetname      ! dataset name
integer(HID_T)     :: file_id       ! File identifier
integer(HID_T)     :: dspace_id     ! Dataspace identifier
integer(HID_T)     :: dset_id       ! Dataset identifier
integer(HSIZE_T)   :: dims(1)       ! Dimensions for Dataset and Dataspace
integer,parameter  :: rank = 1      ! Dataset rank = number of dimensions
integer            :: error         ! Error flag
dims=n
write(dsetname,'(<l-3>a)') fname(1:l-3)
! Initialize FORTRAN interface.
call h5open_f(error)
! Create a new file using default properties.
call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)
! Create the dataspace.
call h5screate_simple_f(rank, dims, dspace_id, error)
! Create the dataset with default properties.
call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, dset_id,error)
! Write the data to datset
call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, dims, error)
! End access to the dataset and release resources used by it.
call h5dclose_f(dset_id, error)
! Terminate access to the data space.
call h5sclose_f(dspace_id, error)
! Close the file.
call h5fclose_f(file_id, error)
! Close FORTRAN interface.
call h5close_f(error)
end subroutine save_vector_h5
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine save_matrix_h5(x,m,n,fname,l)
use HDF5
use global_param
implicit none
integer(kind=4)    :: m,n,l
real(kind=dp)      :: x(m,n)
character(len=l)   :: fname         ! File name
character(len=l-3) :: dsetname      ! dataset name
integer(HID_T)     :: file_id       ! File identifier
integer(HID_T)     :: dspace_id     ! Dataspace identifier
integer(HID_T)     :: dset_id       ! Dataset identifier
integer(HSIZE_T)   :: dims(2)       ! Dimensions for Dataset and Dataspace
integer,parameter  :: rank = 2      ! Dataset rank = number of dimensions
integer            :: error         ! Error flag
dims(1)=m
dims(2)=n
write(dsetname,'(<l-3>a)') fname(1:l-3)
! Initialize FORTRAN interface.
call h5open_f(error)
! Create a new file using default properties.
call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)
! Create the dataspace.
call h5screate_simple_f(rank, dims, dspace_id, error)
! Create the dataset with default properties.
call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, dset_id,error)
! Write the data to datset
call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, dims, error)
! End access to the dataset and release resources used by it.
call h5dclose_f(dset_id, error)
! Terminate access to the data space.
call h5sclose_f(dspace_id, error)
! Close the file.
call h5fclose_f(file_id, error)
! Close FORTRAN interface.
call h5close_f(error)
end subroutine save_matrix_h5
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine will include in the Hamiltonian variable "ham" the NAC
subroutine nac_ha_modify
use global_param
use omp_lib
integer :: cont
real(kind=dp), dimension(Nq1,Nq2) :: nac21q1,nac21q2,nac31q1,nac31q2,nac32q1,nac32q2
real(kind=dp), dimension(Nq1*Nq2) :: vec1,vec2,vec3,vec4,vec5,vec6

!-------------------------------------------------------------------!
!Loading NACS                                                       !
!THE NACS ARE ALREADY MASS SCALED                                   !
open(unit=1,file='nac21q1f.txt',status='old')                       !
open(unit=2,file='nac21q2f.txt',status='old')                       !
open(unit=3,file='nac31q1f.txt',status='old')                       !
open(unit=4,file='nac31q2f.txt',status='old')                       !
open(unit=5,file='nac32q1f.txt',status='old')                       !
open(unit=6,file='nac32q2f.txt',status='old')                       !
do i=1,Nq1                                                          !
  read(1,*) nac21q1(i,:)                                            !
  read(2,*) nac21q2(i,:)                                            !
  read(3,*) nac31q1(i,:)                                            !
  read(4,*) nac31q2(i,:)                                            !
  read(5,*) nac32q1(i,:)                                            !
  read(6,*) nac32q2(i,:)                                            !
end do                                                              !
 close(unit=1)                                                      !
 close(unit=2)                                                      !
 close(unit=3)                                                      !
 close(unit=4)                                                      !
 close(unit=5)                                                      !
 close(unit=6)                                                      !
!-------------------------------------------------------------------!
!vectorizing the grid nac matrix
cont=0
do j=1,Nq2
  do i=1,Nq1
    cont=cont+1;
    vec1(cont)=nac21q1(i,j)*ai*10.d0 ! The 'ai' and 'bi' terms is to transform back to cartesian
    vec2(cont)=nac21q2(i,j)*bi*10.d0 
    vec3(cont)=nac31q1(i,j)*ai*10.d0
    vec4(cont)=nac31q2(i,j)*bi*10.d0
    vec5(cont)=nac32q1(i,j)*ai*10.d0
    vec6(cont)=nac32q2(i,j)*bi*10.d0
  end do
end do
!vec1(1)=1.0d0;vec2(1)=1.0d0 ;vec3(1)=1.0d0  ;vec4(1)=1.0d0; vec5(1)=1.0d0  ;vec6(1)=1.0d0
!vec1(2)=1.0d0;vec2(2)=1.0d0 ;vec3(2)=1.0d0  ;vec4(2)=1.0d0; vec5(2)=1.0d0  ;vec6(2)=1.0d0   
!vec1(3)=1.0d0;vec2(3)=1.0d0 ;vec3(3)=1.0d0  ;vec4(3)=1.0d0; vec5(3)=1.0d0  ;vec6(3)=1.0d0   
!vec1(4)=1.0d0;vec2(4)=1.0d0 ;vec3(4)=1.0d0  ;vec4(4)=1.0d0; vec5(4)=1.0d0  ;vec6(4)=1.0d0      
!vec1(5)=2.0d0;vec2(5)=2.0d0 ;vec3(5)=2.0d0  ;vec4(5)=2.0d0; vec5(5)=2.0d0  ;vec6(5)=2.0d0      
!$OMP PARALLEL DO shared(moq1,moq2,ham,vec1,vec2,vec3,vec4,vec5,vec6)
do i=1,s
  do j=1,s
    ham(i+0*s,j+1*s)= vec1(i)*moq1(i,j) + vec2(i)*moq2(i,j)
    ham(j+1*s,i+0*s)= vec1(i)*moq1(i,j) + vec2(i)*moq2(i,j)
    ham(i+0*s,j+2*s)= vec3(i)*moq1(i,j) + vec4(i)*moq2(i,j)
    ham(j+2*s,i+0*s)= vec3(i)*moq1(i,j) + vec4(i)*moq2(i,j)
    ham(i+1*s,j+2*s)= vec5(i)*moq1(i,j) + vec6(i)*moq2(i,j)
    ham(j+2*s,i+1*s)= vec5(i)*moq1(i,j) + vec6(i)*moq2(i,j)
  end do
end do
!$OMP end parallel do

! Merging this matrix with the one that includes the second derivatives and the potential energy, Ha
!$OMP PARALLEL DO shared(ha,ham)
do i=1,Nst*Nq1*Nq2
  do j=1,Nst*Nq1*Nq2
    Ha(i,j)=Ha(i,j)-ham(i,j) ! Minus because of the sign of the kinetic energy term
  end do
end do
!$OMP END PARALLEL DO

end subroutine nac_ha_modify
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine angular_momentum(y,n,am)
use global_param
use omp_lib
implicit none
integer          :: n
complex(kind=dp) :: y(n),am(n)
real(kind=dp)    :: amxq1,amyq1,amzq1,amxq2,amyq2,amzq2

am=0.d0
!$OMP PARALLEL DO shared(moq1,moq2,am,im,y), private(j,k,amxq1,amxq2,amyq1,amyq2,amzq1,amzq2)
do i=1,n
  do j=1,n
    amxq1=0.d0
    amxq2=0.d0
    amyq1=0.d0
    amyq2=0.d0
    amzq1=0.d0
    amzq2=0.d0
    do k=1,NA
      amxq1=amxq1+ ( q1i(3*k-1)*co1(j)+q2i(3*k-1)*co2(j) ) * q1(3*k  ) - ( q1i(3*k  )*co1(j)+q2i(3*k  )*co2(j) ) * q1(3*k-1)
      amxq2=amxq2+ ( q1i(3*k-1)*co1(j)+q2i(3*k-1)*co2(j) ) * q2(3*k  ) - ( q1i(3*k  )*co1(j)+q2i(3*k  )*co2(j) ) * q2(3*k-1)
      amyq1=amyq1+ ( q1i(3*k  )*co1(j)+q2i(3*k  )*co2(j) ) * q1(3*k-2) - ( q1i(3*k-2)*co1(j)+q2i(3*k-2)*co2(j) ) * q1(3*k  )
      amyq2=amyq2+ ( q1i(3*k  )*co1(j)+q2i(3*k  )*co2(j) ) * q2(3*k-2) - ( q1i(3*k-2)*co1(j)+q2i(3*k-2)*co2(j) ) * q2(3*k  )
      amzq1=amzq1+ ( q1i(3*k-2)*co1(j)+q2i(3*k-2)*co2(j) ) * q1(3*k-1) - ( q1i(3*k-1)*co1(j)+q2i(3*k-1)*co2(j) ) * q1(3*k-2)
      amzq2=amzq2+ ( q1i(3*k-2)*co1(j)+q2i(3*k-2)*co2(j) ) * q2(3*k-1) - ( q1i(3*k-1)*co1(j)+q2i(3*k-1)*co2(j) ) * q2(3*k-2)
    end do
    am(i)=am(i) + (-im) * ( (amxq1+amyq1+amzq1)*moq1(i,j) + (amxq2+amyq2+amzq2)*moq2(i,j) ) * y(j)
  end do
end do
!$OMP END PARALLEL DO






end subroutine angular_momentum
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
