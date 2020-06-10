program dyn
!maybe is useful: Function MERGE:
!test a condition, if true, assign a value X, if false assign value Y:
!the condition can be a boolean variable defined before. 
!result = merge( X, Y, i.eq.j)
!result, X and Y will be of same type
use global_param
!use random_module
use omp_lib
implicit none
external                     :: rkdumb,HA_calc,momentum_calc_q1,momentum_calc_q2,angular_momentum
integer                      :: init_wf,ppp
integer                      :: cont,n,jj,ii !m=dimension of the Hamiltonian (Nq1*Nq2*Nst x Nq1*Nq2*Nst)
real(kind=dp),allocatable    :: wfinal(:,:)!pice(:),nwf0(:)
real(kind=dp)                :: q10,q20 ! point along q1 and q2 where the minimum global C2v is
real(kind=dp)                :: t,truni,trunf,start_time,stop_time,ompt0,ompt1,tt1,tt0 !time counters
complex(kind=dp)             :: mom(0:Nst*Nq1*Nq2-1),am(0:Nst*Nq1*Nq2-1) !vectors to store operators acting on the wave function
character(len=18)            :: fname
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! I modified the Ha, moq1 and moq2 matrices building so that their indexes begin with 0 instead of 1. 
! It is not the ideal way, I use a auxiliar matrix to build with index 1 as I already had written in the code and then I put it
! in a matrix with index starting in 0.
! Later would be better and faster if I already build the matrices with index 0.     <<<---------------
open(unit=33,file='simu_parameters',status='unknown') !File to save some parameters for the simulation
write(33,*)tf,npoints,nfiles,Nq1,Nq2

!$ truni = omp_get_wtime()
open(unit=100,file='output',status='unknown') ! output for following the code running
write(100,'(a)')'INITIATING SIMULATION'

n=Nst*Nq1*Nq2 ! m=Nq1*Nq2*st %size of the Hamiltonian

call load_data
open(unit=20,file='csr_vectors',status='unknown')
read(20,'(i12)')k_moq1
allocate(moq1_val(0:k_moq1-1),moq1_rowc(0:n), moq1_row_col(0:k_moq1-1,0:1))
do i=0,k_moq1-1
  read(20,'(e23.15e3,2i12)')moq1_val(i),moq1_row_col(i,0),moq1_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')moq1_rowc(i)
end do

read(20,'(i12)')k_moq2
allocate(moq2_val(0:k_moq2-1),moq2_rowc(0:n), moq2_row_col(0:k_moq2-1,0:1))
do i=0,k_moq2-1
  read(20,'(e23.15e3,2i12)')moq2_val(i),moq2_row_col(i,0),moq2_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')moq2_rowc(i)
end do

read(20,'(i12)')k_Ha
allocate(Ha_val(0:k_Ha-1),Ha_rowc(0:n), Ha_row_col(0:k_Ha-1,0:1))
do i=0,k_Ha-1
  read(20,'(e23.15e3,2i12)')Ha_val(i),Ha_row_col(i,0),Ha_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')Ha_rowc(i)
end do

read(20,'(i12)')k_dip
allocate(dip_val(0:k_dip-1,0:2))
do i=0,k_dip-1
  read(20,'(3e23.15e3)')dip_val(i,0),dip_val(i,1),dip_val(i,2)
end do

read(20,'(i12)')k_am
allocate(am_val(0:k_am-1),am_rowc(0:n), am_row_col(0:k_am-1,0:1))
do i=0,k_am-1
  read(20,'(e23.15e3,2i12)')am_val(i),am_row_col(i,0),am_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')am_rowc(i)
end do

close(unit=20)

open(unit=20,file='csr_vectors_only_NAC',status='unknown')
read(20,'(i12)')k_Ha2
allocate(Ha2_val(0:k_Ha2-1),Ha2_rowc(0:n), Ha2_row_col(0:k_Ha2-1,0:1))
do i=0,k_Ha2-1
  read(20,'(e23.15e3,2i12)')Ha2_val(i),Ha2_row_col(i,0),Ha2_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')Ha2_rowc(i)
end do
close(unit=20)

write(100,'(a)')'PARAMETERS OF THE DYNAMICS'
write(100,'(a50,i)')'Number of states =',Nst
write(100,'(a50,i)')'Number of points in coordinate 1 =',Nq1
write(100,'(a50,i)')'Number of points in coordinate 2 =',Nq2
write(100,'(a50,f11.4)')'Step size in coordinate 1 =',sq1
write(100,'(a50,f11.4)')'Step size in coordinate 2 =',sq2
write(100,'(a50,f11.4)')'Initial 1 / dq1 =',1.d0/sq1
write(100,'(a50,f11.4)')'Initial 1 / dq2 =',1.d0/sq2
write(100,'(a72,i5)')'Number of initial random orientations for the ionizing electric field = ',nsamples
write(100,'(a50,i7,a17,f6.1,a15)')'The time evolution goes until =',int(tf),' atomic units or ',tf*24.18884d0/1000.d0,' femtoseconds'
write(100,'(a50,i)')'Number of time evaluations =',npoints
write(100,'(a50,f6.4,a8,f10.2,a12)')'time step =',tstep,' au. or ',tstep*24.18884d0,' attoseconds'
write(100,'(a50,3f4.1)')'Orientation of the probing electric field =',ori
write(100,'(a50,(es15.7e3))')'Energy of the ionizing pulse =',e_ip
write(100,'(a50,(es15.7e3))')'Time where the electric field is centered =',t00
write(100,'(a50,(es15.7e3))')'Phase of the probing electric field =',phase
write(100,'(a50,(es15.7e3))')'Energy of the probing electric field =',freq
write(100,'(a50,(es15.7e3))')'Duration of the probing electric field (sigma) =',sig
write(100,'(a50,(es15.7e3))')'Intensity of the probing electric field =',E00
write(100,'(a)')'##############################################################################'
write(100,'(a)')'CSR vectors loaded'

allocate ( nwf0(0:Nst*Nq1*Nq2-1) )
nwf0=dcmplx(0.d0,0.d0)
open(newunit=init_wf,file='eigen-vec-neutral',status='unknown')
do i=0,s-1
  read(init_wf,*) nwf0(i) !reading the neutral initial eigen state and writting in the first third part of the vector
end do
allocate ( pice(0:Nst*Nq1*Nq2-1) )
allocate(pia(Nst,3))
allocate ( wf0(0:Nst*Nq1*Nq2-1) )
ij = 0
do j=q2_initial+20,q2_final+20  !The photoionization coefficiets were calculated only for q2=75+20:115+20 and q1=25+27:70+27
  do i=q1_initial+27,q1_final+27 !where the amplitudes of the eigen state of the neutral is non-zero (< 10^-8). This is the Frank-Condon region
    ij = ij + 1
    call p_i_a(i-27,j-20,pia) !Evaluating the photoionization coeficients for all electronic states
    phote(ij,:,:) = pia(:,:)
  end do
end do
write(100,'(a)')'Photoelecton coupling elements loaded'
write(100,'(a)')'Initiating dynamics'
allocate ( wfout(0:n-1,0:nfiles),wfinal(0:n-1,0:nfiles) )
allocate ( y_f(0:n-1,0:npoints) )
allocate( coh(0:n-1,0:nfiles),cohe(0:n-1,0:nfiles) )
wfinal = 0.d0
y_f = dcmplx(0.d0,0.d0)
cohe = dcmplx(0.d0,0.d0)
allocate( orie(nsamples,3) )
call generate_random_orientation
open(unit=44,file='orientations',status='unknown')
!CALL init_random_seed()         ! see example of RANDOM_SEED
fa(:)=0.d0;fb(:)=0.d0;fc(:)=0.d0;fd(:)=0.d0;fe(:)=0.d0;ff(:)=0.d0;fg(:)=0.d0;fh(:)=0.d0;fi(:)=0.d0;fj(:)=0.d0;fk(:)=0.d0
fl(:)=0.d0;fm(:)=0.d0;fn(:)=0.d0;fo(:)=0.d0
do ppp = 1,nsamples
!  if ( ppp <= nsamples/4 ) then
!    orie = [ 1.d0, 1.d0, 1.d0 ]
!    call generate_random_orientation
!    if ( nsamples == 4 ) then 
!      orientation = [ 1.d0, 1.d0, 1.d0]
!    end if
!  elseif ( ppp > nsamples/4 .and. ppp <= nsamples/4*2 ) then
!    orie = [-1.d0,-1.d0, 1.d0 ]
!    call generate_random_orientation
!    if ( nsamples == 4 ) then 
!      orientation = [-1.d0,-1.d0, 1.d0]
!    end if
!  elseif ( ppp > nsamples/4*2 .and. ppp <= nsamples/4*3 ) then
!    orie = [ 1.d0,-1.d0,-1.d0 ]
!    call generate_random_orientation
!    if ( nsamples == 4 ) then 
!      orientation = [ 1.d0,-1.d0,-1.d0]
!    end if
!  elseif (nsamples == 1 ) then
!    orientation = [ 1.d0, 1.d0, 1.d0]
!  else
!    orie = [-1.d0, 1.d0,-1.d0 ]
!    call generate_random_orientation
!    if ( nsamples == 4 ) then 
!      orientation = [-1.d0, 1.d0,-1.d0]
!    end if
!  end if
  if (nsamples == 4 ) then
    if ( ppp == 1 ) then
      orientation = [ 1.d0, 1.d0, 1.d0]
    elseif ( ppp == 2 ) then
      orientation = [-1.d0,-1.d0, 1.d0]
    elseif ( ppp == 3 ) then
      orientation = [ 1.d0,-1.d0,-1.d0]
    else
      orientation = [-1.d0, 1.d0,-1.d0]
    end if
  else
    orientation = orie(ppp,:)
  end if
  
  write(44,'(3f15.7)') orientation

  call generate_initial_wf
  wfout(:,0) = wf0(:)
  do i = 0,s-1
    coh(i,0)     = dconjg(wf0(i))   * wf0(i+s)
    coh(i+s,0)   = dconjg(wf0(i))   * wf0(i+2*s)
    coh(i+2*s,0) = dconjg(wf0(i+s)) * wf0(i+2*s)
  end do
  fname='iniWF-r-ori0000.h5'
  write(fname(12:15),'(i0.4)') ppp 
  call save_vector_h5(real(wf0),n,fname,18)
  fname='iniWF-i-ori0000.h5'
  write(fname(12:15),'(i0.4)') ppp 
  call save_vector_h5(aimag(wf0),n,fname,18)
  
  call rkdumb(wf0,n,HA_calc)
  maxmomq1(ppp) = maxval(momq1t)
  maxmomq2(ppp) = maxval(momq2t)
  !$OMP parallel do shared(wfout,wfinal)
  do j = 0,nfiles
    do i = 0,n-1
      wfinal(i,j) = wfinal(i,j) + dconjg( wfout(i,j) ) * wfout(i,j)
    end do
  end do
  !$OMP end parallel do
  cohe(0:s-1,:)   = cohe(0:s-1,:)   + coh(0:s-1,:)
  cohe(s:2*s-1,:) = cohe(s:2*s-1,:) + coh(s:2*s-1,:)
  cohe(2*s:n-1,:) = cohe(2*s:n-1,:) + coh(2*s:n-1,:)

  fa(:) = fa(:) + sum1(:)
  fb(:) = fb(:) + sum2(:)
  fc(:) = fc(:) + sum3(:)
  fd(:) = fd(:) + e1(:)
  fe(:) = fe(:) + e2(:)
  ff(:) = ff(:) + e3(:)
  fg(:) = fg(:) + L1(:)
  fh(:) = fh(:) + L2(:)
  fi(:) = fi(:) + L3(:)
  fj(:) = fj(:) + pq1_1(:)
  fk(:) = fk(:) + pq2_1(:)
  fl(:) = fl(:) + pq1_2(:)
  fm(:) = fm(:) + pq2_2(:)
  fn(:) = fn(:) + pq1_3(:)
  fo(:) = fo(:) + pq2_3(:)
write(100,'(a10,f7.2,a2)')'Progress: ',real(real(ppp)/real(nsamples)*100),' %'
end do
wfinal(:,:) = wfinal(:,:) / nsamples
cohe(:,:) = cohe(:,:) / nsamples
fa(:) = fa(:) / nsamples 
fb(:) = fb(:) / nsamples 
fc(:) = fc(:) / nsamples 
fd(:) = fd(:) / nsamples 
fe(:) = fe(:) / nsamples 
ff(:) = ff(:) / nsamples 
fg(:) = fg(:) / nsamples 
fh(:) = fh(:) / nsamples 
fi(:) = fi(:) / nsamples 
fj(:) = fj(:) / nsamples 
fk(:) = fk(:) / nsamples 
fl(:) = fl(:) / nsamples 
fm(:) = fm(:) / nsamples 
fn(:) = fn(:) / nsamples 
fo(:) = fo(:) / nsamples 

open(unit=15,file='alldata.data',status='unknown')

write(15,'(a3,3(e26.14e3))')'# ',fd(0),fe(0),ff(0)
write(15,'(a)') '# time, Pulse, norm1, norm2, norm3, E1, E2, E3, angular momentum st1, angular momentum st2, angular momentum st3,&
linear momentum in q1 st1, linear momentum in q2 st1, linear momentum in q1 st2, linear momentum in q2 st2,&
linear momentum in q1 st3, linear momentum in q2 st3'
write(100,'(a)') '   time  ,   Pulse   ,     E1    ,    E2     ,    E3     ,   ET-E0   ,   norm1   ,   norm2   ,   norm3   &
, NormT-1   ,   Ltot    '

t = t0
ii = 0; ll = 0!1
do ll=1,nfiles ! for saving 1000 time samples
  do k=1,(npoints/nfiles)
    ii = ii + 1
    t = t + tstep
  end do
  fname='diff-pop-000000.h5'
  write(fname(10:15),'(i0.6)') ll
  call save_vector_h5(y_f(:,ii-1),n,fname,18) ! saving the difference of the pop between 2 points in time to see pop transfer in the grid
  fname='time-pop-000000.h5'
  write(fname(10:15),'(i0.6)') ll
  call save_vector_h5(wfinal(:,ll-1),n,fname,18)
  fname='real-coh-000000.h5'
  write(fname(10:15),'(i0.6)') ll
  call save_vector_h5(real(cohe(:,ll-1)),n,fname,18)
  fname='imag-coh-000000.h5'
  write(fname(10:15),'(i0.6)') ll
  call save_vector_h5(aimag(cohe(:,ll-1)),n,fname,18)

  Et = (E00/freq)*(-(t-t00)/sig**2.d0*sin(freq*(t-t00)+phase)+freq*cos(freq*(t-t00)+phase))*dexp(-(t-t00)**2.d0/(2.d0*sig**2.d0))
  write(15,'(17(es26.16e3))') t-tstep,Et,fa(ll),fb(ll),fc(ll),fd(ll),fe(ll),ff(ll),fg(ll),fh(ll),fi(ll),fj(ll),fk(ll),&
fl(ll),fm(ll),fn(ll),fo(ll)
  write(100,'(f9.2,11(e12.3e3))') t,Et,fd(ll),fe(ll),ff(ll),fd(ll)+fe(ll)+ff(ll)-(fd(0)+fe(0)+ff(0)),fa(ll),fb(ll),fc(ll),&
fa(ll)+fb(ll)+fc(ll)-1.d0,fg(ll)+fh(ll)+fi(ll)
end do
write(100,'(a30,(e12.3e3))')'Norm conservation =',( fa(0)+fb(0)+fc(0) )-( fa(nfiles) + fb(nfiles) + fc(nfiles) )
write(100,'(a30,(e12.3e3))')'Final norm for state 1 =',fa(nfiles)
write(100,'(a30,(e12.3e3))')'Final norm for state 2 =',fb(nfiles)
write(100,'(a30,(e12.3e3))')'Final norm for state 3 =',fc(nfiles)
write(100,'(a30,(e12.3e3),a8)')'Energy conservation =',(fd(0)+fe(0)+ff(0))-(fd(nfiles)+fe(nfiles)+ff(nfiles)),' hartree'
write(100,'(a30,(f20.15),a8)')'Final total energy state 1 =',fd(nfiles),' hartree'
write(100,'(a30,(f20.15),a8)')'Final total energy state 2 =',fe(nfiles),' hartree'
write(100,'(a30,(f20.15),a8)')'Final total energy state 3 =',ff(nfiles),' hartree'

open(unit=222,file='final-wave-packet',status='unknown')
do i=0,n-1
  write(222,'(e26.14e3)') wfinal(i,nfiles)
end do
 close(unit=222)
write(100,'(a24,e12.3e3)') 'maximum momentum value =',maxval(maxmomq1)
write(100,'(a24,e12.3e3)') 'maximum momentum value =',maxval(maxmomq2)
!$ trunf = omp_get_wtime()
write(100,'(a38,f9.1,a7)')'Time took to run totaly the program = ', (trunf-truni)/3600.d0, ' hours.'
write(100,'(a)')'************************************************************'
write(100,'(a)')'Program finished'
write(100,'(a)')'************************************************************'
 close(unit=15)
 close(unit=100)
end program
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!subroutine rkdumb(y0,n,t0,tf,npoints,h,HA_calc)
subroutine rkdumb(y,n,HA_calc)
use omp_lib
use global_param
external rk4,HA_calc,momentum_calc_q1,momentum_calc_q2,angular_momentum
integer           :: n,ii
real(kind=dp)     :: tt(npoints),h,t!,coordinate1(s),coordinate2(s)
complex(kind=dp)  :: dydt(0:n-1),y(0:n-1)!,y0(n)!,ymatrix(n,npoints+1)!,ymatlab(n,npoints)
!real(kind=dp)     :: Te1,Te2,Te3!,enert(npoints+1) ! This is fo saving total energy through time
complex(kind=dp)  :: momq1(0:n-1),momq2(0:n-1),am(0:n-1) ! This is for saving momentum through time

h = tstep
t = t0
!-----------------------------------------------------------!
!Checking momentum and saving norm                          !
call HA_calc(t,y,n,dydt,h) ! evaluating y'(t,y)             !
call momentum_calc_q1(y,momq1,n) ! evaluating dydq          !
call momentum_calc_q2(y,momq2,n) ! evaluating dydq          !
call angular_momentum(y,n,am)                               !
pq1_1(0)=0.0d0 ! momentum in q1 state1                      !
pq2_1(0)=0.0d0 ! momentum in q2 state1                      !
pq1_2(0)=0.0d0 ! momentum in q1 state2                      !
pq2_2(0)=0.0d0 ! momentum in q2 state2                      !
pq1_3(0)=0.0d0 ! momentum in q1 state3                      !
pq2_3(0)=0.0d0 ! momentum in q2 state3                      !
sum1(0)=0.0d0  ! norm for state 1                           !
sum2(0)=0.0d0  ! norm for state 2                           !
sum3(0)=0.0d0  ! norm for state 3                           !
e1(0)=0.d0     ! energy of state 1                          !
e2(0)=0.d0     ! energy of state 2                          !
e3(0)=0.d0     ! energy of state 3                          !
L1(0)=0.d0     ! angular momentum of state 1                !
L2(0)=0.d0     ! angular momentum of state 2                !
L3(0)=0.d0     ! angular momentum of state 3                !
do i=0,s-1                                                  !
  pq1_1(0) = pq1_1(0) + dconjg(y(i))     * momq1(i)         !
  pq2_1(0) = pq2_1(0) + dconjg(y(i))     * momq2(i)         !
  pq1_2(0) = pq1_2(0) + dconjg(y(i+s))   * momq1(i+s)       !
  pq2_2(0) = pq2_2(0) + dconjg(y(i+s))   * momq2(i+s)       !
  pq1_3(0) = pq1_3(0) + dconjg(y(i+2*s)) * momq1(i+2*s)     !
  pq2_3(0) = pq2_3(0) + dconjg(y(i+2*s)) * momq2(i+2*s)     !
  sum1(0)  = sum1(0)  + dconjg(y(i))     * y(i)             !
  sum2(0)  = sum2(0)  + dconjg(y(i+s))   * y(i+s)           !
  sum3(0)  = sum3(0)  + dconjg(y(i+2*s)) * y(i+2*s)         !
  e1(0)    = e1(0)    + dconjg(y(i))     * dydt(i)*im       !
  e2(0)    = e2(0)    + dconjg(y(i+s))   * dydt(i+s)*im     !
  e3(0)    = e3(0)    + dconjg(y(i+2*s)) * dydt(i+2*s)*im   !
  L1(0)    = L1(0)    + dconjg(y(i))     * am(i)            !
  L2(0)    = L2(0)    + dconjg(y(i+s))   * am(i+s)          !
  L3(0)    = L3(0)    + dconjg(y(i+2*s)) * am(i+2*s)        !
end do                                                      !
!-----------------------------------------------------------!
!write(100,*)e1(0),e2(0),e3(0),sum1(0),sum2(0),sum3(0)
!read(*,*)
ii=0;
do ll=1,nfiles ! for saving 1000 time samples
  do k=1,(npoints/nfiles)
    ii = ii+1
    call HA_calc(t,y,n,dydt) ! evaluating y'(t,y)
    call rk4(y,dydt,n,t,h,y,HA_calc) !evaluating y(t+h)
    !$OMP parallel do shared(y_f,y_nac)
    do i = 0,Nst*s-1
      y_f(i,ii) = y_f(i,ii) + ( dconjg(y_nac(i)) * y_nac(i) )/nsamples
    end do
    !$OMP end parallel do
    t=t+h
  end do
  wfout(:,ll) = y(:)
  do i = 0,s-1
    coh(i,ll)     = dconjg(y(i))   * y(i+s)
    coh(i+s,ll)   = dconjg(y(i))   * y(i+2*s)
    coh(i+2*s,ll) = dconjg(y(i+s)) * y(i+2*s)
  end do
  !------------------------------------------------------------!
  !Checking momentum and saving norm                           !
  call momentum_calc_q1(y,momq1,n) ! evaluating dydq           !
  call momentum_calc_q2(y,momq2,n) ! evaluating dydq           !
  call angular_momentum(y,n,am)                                !
  pq1_1(ll)=0.0d0 ! momentum in q1 state1                      !
  pq2_1(ll)=0.0d0 ! momentum in q2 state1                      !
  pq1_2(ll)=0.0d0 ! momentum in q1 state2                      !
  pq2_2(ll)=0.0d0 ! momentum in q2 state2                      !
  pq1_3(ll)=0.0d0 ! momentum in q1 state3                      !
  pq2_3(ll)=0.0d0 ! momentum in q2 state3                      !
  sum1(ll)=0.0d0  ! norm for state 1                           !
  sum2(ll)=0.0d0  ! norm for state 2                           !
  sum3(ll)=0.0d0  ! norm for state 3                           !
  e1(ll)=0.d0     ! energy of state 1                          !
  e2(ll)=0.d0     ! energy of state 2                          !
  e3(ll)=0.d0     ! energy of state 3                          !
  L1(ll)=0.d0     ! angular momentum of state 1                !
  L2(ll)=0.d0     ! angular momentum of state 2                !
  L3(ll)=0.d0     ! angular momentum of state 3                !
  do i=0,s-1                                                   !
    pq1_1(ll) = pq1_1(ll) + dconjg(y(i))     * momq1(i)        !
    pq2_1(ll) = pq2_1(ll) + dconjg(y(i))     * momq2(i)        !
    pq1_2(ll) = pq1_2(ll) + dconjg(y(i+s))   * momq1(i+s)      !
    pq2_2(ll) = pq2_2(ll) + dconjg(y(i+s))   * momq2(i+s)      !
    pq1_3(ll) = pq1_3(ll) + dconjg(y(i+2*s)) * momq1(i+2*s)    !
    pq2_3(ll) = pq2_3(ll) + dconjg(y(i+2*s)) * momq2(i+2*s)    !
    sum1(ll)  = sum1(ll)  + dconjg(y(i))     * y(i)            !
    sum2(ll)  = sum2(ll)  + dconjg(y(i+s))   * y(i+s)          !
    sum3(ll)  = sum3(ll)  + dconjg(y(i+2*s)) * y(i+2*s)        !
    e1(ll)    = e1(ll)    + dconjg(y(i))     * dydt(i)*im      !
    e2(ll)    = e2(ll)    + dconjg(y(i+s))   * dydt(i+s)*im    !
    e3(ll)    = e3(ll)    + dconjg(y(i+2*s)) * dydt(i+2*s)*im  !
    L1(ll)    = L1(ll)    + dconjg(y(i))     * am(i)           !
    L2(ll)    = L2(ll)    + dconjg(y(i+s))   * am(i+s)         !
    L3(ll)    = L3(ll)    + dconjg(y(i+2*s)) * am(i+2*s)       !
  end do                                                       !
  !------------------------------------------------------------!
  momq1t(ll)=pq1_1(ll)+pq1_2(ll)+pq1_3(ll)
  momq2t(ll)=pq2_1(ll)+pq2_2(ll)+pq2_3(ll)
end do
end subroutine rkdumb
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
!Subroutine computes the derivatives dPsi(t)/dt = H/i Psi
!This derivative routine calculates 
!Written as a simple matrix multiplication between the Hamiltonian matrix H and the amplitudes vector Psi
subroutine HA_calc(t,y,n,dydt)
use global_param
use omp_lib
implicit none
integer n
real(kind=dp) :: t
complex(kind=dp) :: y(0:n-1),dydt(0:n-1)

!write(*,*) 'omp_get_max_threads= ', omp_get_max_threads ( )
!write(*,*) 'omp_get_num_procs = ', omp_get_num_procs ( )
!write(*,*) 'Time = ', omp_get_wtime ( )

Et = (E00/freq)*(-(t-t00)/sig**2.d0*dsin(freq*(t-t00)+phase)+freq*dcos(freq*(t-t00)+phase))*dexp(-(t-t00)**2.d0/(2.d0*sig**2.d0))
dydt=dcmplx(0.d0,0.d0)
y_nac=dcmplx(0.d0,0.d0)
!$OMP PARALLEL DO shared(y,dydt,Ha_val,dip_val,Ha_row_col,Et,y_nac,Ha2_val,Ha2_row_col)
do i=0,n-1
  do j=Ha_rowc(i),Ha_rowc(i+1)-1
    dydt(i) = dydt(i) + ( Ha_val(j) + (dot_product(orientation,dip_val(j,:)) * Et) ) * y(Ha_row_col(j,1)) / im
!    y_nac(i) = y_nac(i) + Ha2_val(j) * y(Ha2_row_col(j,1)) / im
  end do
end do
!$OMP end parallel do
!$OMP PARALLEL DO shared(y,dydt,Ha_val,dip_val,Ha_row_col,Et,y_nac,Ha2_val,Ha2_row_col)
do i=0,n-1
  do j=Ha2_rowc(i),Ha2_rowc(i+1)-1
!    dydt(i) = dydt(i) + ( Ha_val(j) + (dot_product(orientation,dip_val(j,:)) * Et) ) * y(Ha_row_col(j,1)) / im
    y_nac(i) = y_nac(i) + Ha2_val(j) * y(Ha2_row_col(j,1)) / im
  end do
end do
!$OMP end parallel do


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
subroutine save_vector_h5(x,n,fname,le)
use HDF5
use global_param
implicit none
integer(kind=4)     :: n,le
real(kind=dp)       :: x(n)
character(len=le)   :: fname         ! File name
character(len=le-3) :: dsetname      ! dataset name
integer(HID_T)      :: file_id       ! File identifier
integer(HID_T)      :: dspace_id     ! Dataspace identifier
integer(HID_T)      :: dset_id       ! Dataset identifier
integer(HSIZE_T)    :: dims(1)       ! Dimensions for Dataset and Dataspace
integer,parameter   :: rank = 1      ! Dataset rank = number of dimensions
integer             :: error         ! Error flag
dims=n
write(dsetname,'(<le-3>a)') fname(1:le-3)


!if ( all(x == 0.d0) ) then
!  write(100,*)'Error in wfout. It is all zeros '
!  read(*,*)
!endif
!write(*,*)'tudo certo até aqui,dentro do save'
!read(*,*)



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

subroutine progress_bar(jjjj,ffff)
implicit none
integer(kind=4)::jjjj,kkkk,ffff
character(len=30)::bar="?????% |                    | "
write(bar(1:5),'(f5.1)') 100.d0/real(ffff)*jjjj
do kkkk=1, int(real(jjjj)/real(ffff)*20.d0)
  bar(8+kkkk:8+kkkk)="*"
enddo
! print the progress bar.
write(6,'(a1,a30)',advance="no") char(13), bar
if (jjjj/=ffff) then
  flush(6)
else
  write(6,*)
endif
return
end subroutine progress_bar
