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
external                  :: rkdumb,HA_calc,momentum_calc_q1,momentum_calc_q2,angular_momentum
integer                   :: init_wf,ppp
integer,parameter         :: le = 18
integer                   :: cont,n,jj,ii !m=dimension of the Hamiltonian (Nq1*Nq2*Nst x Nq1*Nq2*Nst)
real(kind=dp),allocatable :: wfinal(:,:),temp_r(:),temp_i(:)!pice(:),nwf0(:)
real(kind=dp)             :: q10,q20,coord1(Nq1),coord2(Nq2),dnu1,dnu2 ! point along q1 and q2 where the minimum global C2v is
real(kind=dp)             :: t,truni,trunf,start_time,stop_time !time counters
complex(kind=dp)          :: mom(0:Nst*Nq1*Nq2-1),am(0:Nst*Nq1*Nq2-1) !vectors to store operators acting on the wave function
character(len=le)         :: fname
character(len=20)         :: fname1
character(len=20)         :: fname2
real(kind=dp)             :: norma_ori(nsamples)
real(kind=dp)             :: eig_val(nsamples)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! I modified the Ha, moq1 and moq2 matrices building so that their indexes begin with 0 instead of 1. 
! It is not the ideal way, I use a auxiliar matrix to build with index 1 as I already had written in the code and then I put it
! in a matrix with index starting in 0.
! Later would be better and faster if I already build the matrices with index 0.     <<<---------------

! The C2v point is in position 91,115 in the 146x184 sized matrix, and in 16734 in the vectorized wavefunction (starting at index=0)
! The D2d point is in position 73, 99 in the 146x184 sized matrix, and in 14380 in the vectorized wavefunction (starting at index=0)

!do i = 1,Nq1
!  coord1(i) = (i-Nq1/2.0_dp-1.0_dp)*sq1+sq1/2.0_dp!(Nq1-Nq11/2.0_dp-i)*sq1+sq1/2.0_dp
!end do
!do j=1,Nq2
!  coord2(j)=(j-(114.0_dp)-1.0_dp)*sq2
!end do
!write(*,*)coord1(1),coord1(ubound(coord1,1))
!write(*,*)coord2(1),coord2(ubound(coord2,1))
!dnu1 = 1 / ( coord1(ubound(coord1,1)) - coord1(1) )
!dnu2 = 1 / ( coord2(ubound(coord2,1)) - coord2(1) )
!write(*,*)dnu1,dnu2
!
!! ft_real_signal(b, size, a, fre, c )
!
!read(*,*)


n=Nst*Nq1*Nq2 ! m=Nq1*Nq2*st %size of the Hamiltonian

!$ truni = omp_get_wtime()
!$ start_time = omp_get_wtime()
open(unit=100,file='output',status='unknown') ! output for following the code running
write(100,'(a)')'INITIATING SIMULATION'

open(unit=33,file='simu_parameters',status='unknown') !File to save some parameters for the simulation
write(33,'(f16.3,2i12,2i8)')tf,npoints,nfiles,Nq1,Nq2
close(unit=33)

call load_data

open(unit=20,file='csr_vectors_NAC',status='unknown')
read(20,'(i12)')k_nac
allocate(nac_val(0:k_nac-1),nac_rowc(0:n), nac_row_col(0:k_nac-1,0:1))
do i=0,k_nac-1
  read(20,'(e23.15e3,2i12)')nac_val(i),nac_row_col(i,0),nac_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')nac_rowc(i)
end do
close(unit=20)

open(unit=20,file='csr_vectors_dip',status='unknown')
read(20,'(i12)')k_dipx
allocate(dipx_val(0:k_dipx-1),dipx_rowc(0:n), dipx_row_col(0:k_dipx-1,0:1))
do i=0,k_dipx-1
  read(20,'(e23.15e3,2i12)')dipx_val(i),dipx_row_col(i,0),dipx_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')dipx_rowc(i)
end do
read(20,'(i12)')k_dipy
allocate(dipy_val(0:k_dipy-1),dipy_rowc(0:n), dipy_row_col(0:k_dipy-1,0:1))
do i=0,k_dipy-1
  read(20,'(e23.15e3,2i12)')dipy_val(i),dipy_row_col(i,0),dipy_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')dipy_rowc(i)
end do
read(20,'(i12)')k_dipz
allocate(dipz_val(0:k_dipz-1),dipz_rowc(0:n), dipz_row_col(0:k_dipz-1,0:1))
do i=0,k_dipz-1
  read(20,'(e23.15e3,2i12)')dipz_val(i),dipz_row_col(i,0),dipz_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')dipz_rowc(i)
end do
close(unit=20)

open(unit=20,file='csr_vectors_kin',status='unknown')
read(20,'(i12)')k_kine
allocate(kine_val(0:k_kine-1),kine_rowc(0:n), kine_row_col(0:k_kine-1,0:1))
do i=0,k_kine-1
  read(20,'(e23.15e3,2i12)')kine_val(i),kine_row_col(i,0),kine_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')kine_rowc(i)
end do
close(unit=20)

open(unit=20,file='csr_vectors_pot',status='unknown')
read(20,'(i12)')k_pot
allocate(pot_val(0:k_pot-1),pot_rowc(0:n), pot_row_col(0:k_pot-1,0:1))
do i=0,k_pot-1
  read(20,'(e23.15e3,2i12)')pot_val(i),pot_row_col(i,0),pot_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')pot_rowc(i)
end do
close(unit=20)

open(unit=20,file='csr_vectors_anm',status='unknown')
read(20,'(i12)')k_am
allocate(am_val(0:k_am-1),am_rowc(0:n), am_row_col(0:k_am-1,0:1))
do i=0,k_am-1
  read(20,'(e23.15e3,2i12)')am_val(i),am_row_col(i,0),am_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')am_rowc(i)
end do
close(unit=20)

open(unit=20,file='csr_vectors_dq1',status='unknown')
read(20,'(i12)')k_moq1
allocate(moq1_val(0:k_kine-1),moq1_rowc(0:n), moq1_row_col(0:k_moq1-1,0:1))
do i=0,k_moq1-1
  read(20,'(e23.15e3,2i12)')moq1_val(i),moq1_row_col(i,0),moq1_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')moq1_rowc(i)
end do
close(unit=20)

open(unit=20,file='csr_vectors_dq2',status='unknown')
read(20,'(i12)')k_moq2
allocate(moq2_val(0:k_moq2-1),moq2_rowc(0:n), moq2_row_col(0:k_moq2-1,0:1))
do i=0,k_moq2-1
  read(20,'(e23.15e3,2i12)')moq2_val(i),moq2_row_col(i,0),moq2_row_col(i,1)
end do
do i=0,n
  read(20,'(i12)')moq2_rowc(i)
end do
close(unit=20)

write(100,'(a)')'PARAMETERS OF THE DYNAMICS'
write(100,'(a50,i)')'Number of states =',Nst
write(100,'(a50,i)')'Number of points in coordinate 1 =',Nq1
write(100,'(a50,i)')'Number of points in coordinate 2 =',Nq2
write(100,'(a50,f11.4)')'Step size in coordinate 1 =',sq1
write(100,'(a50,f11.4)')'Step size in coordinate 2 =',sq2
write(100,'(a50,f11.4)')'Initial 1 / dq1 =',1.0_dp/sq1
write(100,'(a50,f11.4)')'Initial 1 / dq2 =',1.0_dp/sq2
write(100,'(a72,i5)')'Number of initial random orientations for the ionizing electric field = ',nsamples
write(100,'(a50,i7,a17,f6.1,a15)')'The time evolution goes until =',int(tf),' atomic units or ',tf*24.18884_dp/1000.0_dp,' femtoseconds'
write(100,'(a50,i)')'Number of time evaluations =',npoints
write(100,'(a50,i)')'Number of snapshots to take of the wp =',nfiles
write(100,'(a50,f7.4,a8,f3.2,a12)')'time step = ',tstep,' au. or ',tstep*24.18884_dp,' attoseconds'
write(100,'(a,f7.4,a,f5.2,a)')'time step for saving files= ',tf/nfiles,' au. or ',tf/nfiles*24.18884_dp/1000.0_dp,' femtoseconds'
write(100,'(a50,3f4.1)')'Orientation of the probing electric field =',ori
write(100,'(a50,(es15.7e3))')'Energy of the ionizing pulse =',e_ip
write(100,'(a50,(es15.7e3))')'Time where the electric field is centered =',t00
write(100,'(a50,(es15.7e3))')'Phase of the probing electric field =',phase
write(100,'(a50,(es15.7e3))')'Energy of the probing electric field =',freq
write(100,'(a50,(es15.7e3))')'Duration of the probing electric field (sigma) =',sig
write(100,'(a50,(es15.7e3))')'Intensity of the probing electric field =',E00
write(100,'(a)')'##############################################################################'
write(100,'(a)')'CSR vectors loaded'

allocate ( wfout(0:n-1,0:nfiles),wfinal(0:n-1,0:nfiles) )
!allocate ( wfinal2(0:n-1,0:nfiles) )
allocate ( y_f (0:n-1,0:nfiles) )
allocate ( y_f2(0:n-1,0:nfiles) )
allocate( coh(0:n-1,0:nfiles),cohe(0:n-1,0:nfiles) )
!allocate( cohe2(0:n-1,0:nfiles) )
allocate ( pice(0:n-1) )
!allocate(pia(Nst,3))
allocate ( wf0(0:n-1) )
!allocate ( nwf0(0:n-1) )
!allocate( orie(nsamples,3) )
wfinal = 0.0_dp
y_f = 0.0_dp
cohe = 0.0_dp
!wfinal2 = 0.0_dp
y_f2 = 0.0_dp
!cohe2 = 0.0_dp
!phote(:,:,:,:) = dcmplx(0.0_dp,0.0_dp)
!nwf0=dcmplx(0.0_dp,0.0_dp)
!open(newunit=init_wf,file='eigen-vec-neutral',status='unknown')
!do i=0,s-1
!  read(init_wf,*) nwf0(i) !reading the neutral initial eigen state and writting in the first third part of the vector
!end do
!do j=q2_initial,q2_final !The photoionization coefficiets were calculated only for q2=75+20:115+20 and q1=25+27:70+27
!  do i=q1_initial,q1_final !where the amplitudes of the eigen state of the neutral is non-zero (< 10^-8). This is the Frank-Condon region
!    call p_i_a(i,j,pia) !Evaluating the photoionization coeficients for all electronic states
!    phote(i+q1i0,j+q2i0,:,:) = pia(:,:)
!  end do
!  call progress_bar(j-q2_initial+1,q2_final-q2_initial)
!end do
!write(100,'(a)')'Photoelecton coupling elements loaded'
!call generate_random_orientation
!open(unit=44,file='orientations',status='unknown')
!!$ trunf = omp_get_wtime()
!write(100,'(a,f9.1,a)')'Time to load PICEs = ', (trunf-truni)/60.0_dp, ' minutes.'

!CALL init_random_seed()         ! see example of RANDOM_SEED
fa(:)=0.0_dp;fb(:)=0.0_dp;fc(:)=0.0_dp;fd(:)=0.0_dp;fe(:)=0.0_dp;ff(:)=0.0_dp;fg(:)=0.0_dp;fh(:)=0.0_dp;fi(:)=0.0_dp;fj(:)=0.0_dp;fk(:)=0.0_dp
fl(:)=0.0_dp;fm(:)=0.0_dp;fn(:)=0.0_dp;fo(:)=0.0_dp
fp(:)=0.0_dp;fq(:)=0.0_dp;fr(:)=0.0_dp
fa2(:)=0.0_dp;fb2(:)=0.0_dp;fc2(:)=0.0_dp;fd2(:)=0.0_dp;fe2(:)=0.0_dp;ff2(:)=0.0_dp;fg2(:)=0.0_dp;fh2(:)=0.0_dp;fi2(:)=0.0_dp;fj2(:)=0.0_dp;fk2(:)=0.0_dp
fl2(:)=0.0_dp;fm2(:)=0.0_dp;fn2(:)=0.0_dp;fo2(:)=0.0_dp
fp2(:)=0.0_dp;fq2(:)=0.0_dp;fr2(:)=0.0_dp
!$ trunf = omp_get_wtime()
write(100,'(a14,f9.1,a9)')'Time so far = ', (trunf-truni)/60.0_dp, ' minutes.'
write(100,'(a)')'Initiating dynamics'
!open(unit=155,file='norma',status='unknown')
allocate(ax(n,2*nsamples))
!$ start_time = omp_get_wtime()

if ( e_ip == 13.0_dp/27.21138628_dp ) then
  fname1 = 'singular_values_13ev'
  fname2 = 'U_vec_13ev'
  write(*,'(a)')'Ionizing electric field corresponding to 13 eV'
else
  if     ( int(e_ip*27.21138628_dp) == 13 ) then
    fname1 = 'singular_values_09h'
    fname2 = 'U_vec_09h'
    write(*,'(a)')'Ionizing electric field corresponding to the 09th harmonic'
  elseif ( int(e_ip*27.21138628_dp) == 17 ) then
    fname1 = 'singular_values_11h'
    fname2 = 'U_vec_11h'
    write(*,'(a)')'Ionizing electric field corresponding to the 11th harmonic'
  elseif ( int(e_ip*27.21138628_dp) == 20 ) then
    fname1 = 'singular_values_13h'
    fname2 = 'U_vec_13h'
    write(*,'(a)')'Ionizing electric field corresponding to the 13th harmonic'
  elseif ( int(e_ip*27.21138628_dp) == 15 ) then
    fname1 = 'singular_values_15'
    fname2 = 'U_vec_15'
    write(*,'(a)')'Ionizing electric field of 15 eV'
  elseif ( int(e_ip*27.21138628_dp) == 21 ) then
    fname1 = 'singular_values_21'
    fname2 = 'U_vec_21'
    write(*,'(a)')'Ionizing electric field of 21.22 eV'
  elseif ( int(e_ip*27.21138628_dp) == 27 ) then
    fname1 = 'singular_values_27'
    fname2 = 'U_vec_27'
    write(*,'(a)')'Ionizing electric field of 27 eV'
  end if
end if

open(unit=136,file=fname1,status='unknown') ! <<<<<------------ Remember to edit this
do i=1,nsamples
  read(136,'(e25.16e2)')eig_val(i)
  write(*,*)eig_val(i)
end do
do ppp=1,3  !nsamples    <<<<<====================----------------==================----------------===============----------------================
  eig_val(ppp) = eig_val(ppp)**2.0_dp
  open(unit=131,file=fname2,status='unknown')
  do i=0,n-1
    read(131,'(<2*nsamples>e25.16e2)')ax(i+1,:)
    wf0(i) = dcmplx(ax(i+1,ppp*2-1),ax(i+1,ppp*2))
  end do
  close(unit=131)
  wfout(:,0) = wf0(:)
  do i = 0,s-1
    coh(i,0)     = 2.0_dp * dreal(dconjg(wf0(i))   * wf0(i+s)  ) 
    coh(i+s,0)   = 2.0_dp * dreal(dconjg(wf0(i))   * wf0(i+2*s)) 
    coh(i+2*s,0) = 2.0_dp * dreal(dconjg(wf0(i+s)) * wf0(i+2*s)) 
  end do
  call rkdumb(wf0,n,HA_calc)
!OMP parallel do shared(y_f2)
  do i=0,n-1
    y_f2(i,:) = y_f2(i,:) + y_f(i,:) * eig_val(ppp)
  end do
  !OMP end parallel do  
  maxmomq1(ppp) = maxval(momq1t)
  maxmomq2(ppp) = maxval(momq2t)
  !$OMP parallel do shared(wfout,wfinal)
  do j = 0,nfiles
    do i = 0,n-1
      wfinal(i,j) = wfinal(i,j) + dconjg( wfout(i,j) ) * wfout(i,j) * eig_val(ppp)
    end do
  end do
  !$OMP end parallel do
  cohe(0:s-1,:)   = cohe(0:s-1,:)   + coh(0:s-1,:)   * eig_val(ppp)  
  cohe(s:2*s-1,:) = cohe(s:2*s-1,:) + coh(s:2*s-1,:) * eig_val(ppp) 
  cohe(2*s:n-1,:) = cohe(2*s:n-1,:) + coh(2*s:n-1,:) * eig_val(ppp)

  fa(:) = fa(:) + sum1(:) * eig_val(ppp)   
  fb(:) = fb(:) + sum2(:) * eig_val(ppp)  
  fc(:) = fc(:) + sum3(:) * eig_val(ppp)  
  fd(:) = fd(:) + e1(:)   * eig_val(ppp)  
  fe(:) = fe(:) + e2(:)   * eig_val(ppp)  
  ff(:) = ff(:) + e3(:)   * eig_val(ppp)  
  fg(:) = fg(:) + L1(:)   * eig_val(ppp) 
  fh(:) = fh(:) + L2(:)   * eig_val(ppp)   
  fi(:) = fi(:) + L3(:)   * eig_val(ppp)  
  fj(:) = fj(:) + pq1_1(:)* eig_val(ppp)  
  fk(:) = fk(:) + pq2_1(:)* eig_val(ppp) 
  fl(:) = fl(:) + pq1_2(:)* eig_val(ppp)  
  fm(:) = fm(:) + pq2_2(:)* eig_val(ppp)  
  fn(:) = fn(:) + pq1_3(:)* eig_val(ppp) 
  fo(:) = fo(:) + pq2_3(:)* eig_val(ppp) 
  fp(:) = fp(:) + nac0(:) * eig_val(ppp)  
  fq(:) = fq(:) + nac1(:) * eig_val(ppp)  
  fr(:) = fr(:) + nac2(:) * eig_val(ppp)  
!!$ stop_time = omp_get_wtime()
!write(100,'(a10,f7.2,a6,f9.1,a8)')'Progress: ',real(ppp)/real(nsamples)*100,' % in ',(stop_time-start_time)/60.0_dp,' minutes'
  call progress_bar(ppp,nsamples)
end do
!$ stop_time = omp_get_wtime()
write(100,'(a,f9.1,a)')'Time for the dynamics to run: ',(stop_time-start_time)/60.0_dp,' minutes'

open(unit=15,file='alldata',status='unknown')

write(15,'(a3,3(e26.14e3))')'# ',fd(0),fe(0),ff(0)
write(15,'(a)') '# time, Pulse, norm1, norm2, norm3, E1, E2, E3, angular momentum st1, angular momentum st2, angular momentum st3,&
linear momentum in q1 st1, linear momentum in q2 st1, linear momentum in q1 st2, linear momentum in q2 st2,&
linear momentum in q1 st3, linear momentum in q2 st3'
write(100,'(a)') '   time  ,   Pulse   ,     E1    ,    E2     ,    E3     ,   ET-E0   ,   norm1   ,   norm2   ,   norm3   &
, NormT-1   ,   Ltot    '
write(15,'(17(e17.8e3))') t,Et,fa(0),fb(0),fc(0),fd(0),fe(0),ff(0),fg(0),fh(0),fi(0),fj(0),fk(0),&
fl(0),fm(0),fn(0),fo(0)

ll = 1
fname='diff-pop-000000.h5'
write(fname(10:15),'(i0.6)') ll
call save_vector_h5(y_f2(:,ll-1),n,fname,18) ! saving the difference of the pop between 2 points in time to see pop transfer in the grid
fname='time-pop-000000.h5'
write(fname(10:15),'(i0.6)') ll
call save_vector_h5(wfinal(:,ll-1),n,fname,18)
fname='real-coh-000000.h5'
write(fname(10:15),'(i0.6)') ll
call save_vector_h5(cohe(:,ll-1),n,fname,18)

!fname='dif2-pop-000000.h5'
!write(fname(10:15),'(i0.6)') ll
!call save_vector_h5(y_f2(:,ii-1),n,fname,18) ! saving the difference of the pop between 2 points in time to see pop transfer in the grid
!fname='tim2-pop-000000.h5'
!write(fname(10:15),'(i0.6)') ll
!call save_vector_h5(wfinal2(:,ll-1),n,fname,18)
!fname='rea2-coh-000000.h5'
!write(fname(10:15),'(i0.6)') ll
!call save_vector_h5(cohe2(:,ll-1),n,fname,18)

open(unit=74,file='norma',status='unknown')
do i=1,nsamples
  write(74,'(f20.12)') norma_ori(i)
end do

t = t0
ii = 0; 
do ll=1,nfiles ! for saving 1000 time samples
  do k=1,(npoints/nfiles)
    ii = ii + 1
    t = t + tstep
  end do
  fname='diff-pop-000000.h5'
  write(fname(10:15),'(i0.6)') ll+1
  call save_vector_h5(y_f(:,ll-1),n,fname,18) ! saving the difference of the pop between 2 points in time to see pop transfer in the grid
  fname='time-pop-000000.h5'
  write(fname(10:15),'(i0.6)') ll+1
  call save_vector_h5(wfinal(:,ll),n,fname,18)
  fname='real-coh-000000.h5'
  write(fname(10:15),'(i0.6)') ll+1
  call save_vector_h5(cohe(:,ll),n,fname,18)

!  fname='dif2-pop-000000.h5'
!  write(fname(10:15),'(i0.6)') ll+1
!  call save_vector_h5(y_f2(:,ii-1),n,fname,18) ! saving the difference of the pop between 2 points in time to see pop transfer in the grid
!  fname='tim2-pop-000000.h5'
!  write(fname(10:15),'(i0.6)') ll+1
!  call save_vector_h5(wfinal2(:,ll),n,fname,18)
!  fname='rea2-coh-000000.h5'
!  write(fname(10:15),'(i0.6)') ll+1
!  call save_vector_h5(cohe2(:,ll),n,fname,18)

  Et = (E00/freq)*(-(t-t00)/sig**2.0_dp*sin(freq*(t-t00)+phase)+freq*cos(freq*(t-t00)+phase))*dexp(-(t-t00)**2.0_dp/(2.0_dp*sig**2.0_dp))
  write(15,'(17(e17.8e3))') t,Et,fa(ll),fb(ll),fc(ll),fd(ll),fe(ll),ff(ll),fg(ll),fh(ll),fi(ll),fj(ll),fk(ll),&
fl(ll),fm(ll),fn(ll),fo(ll)
  write(100,'(f9.2,11(e12.3e3))') t,Et,fd(ll),fe(ll),ff(ll),fd(ll)+fe(ll)+ff(ll)-(fd(0)+fe(0)+ff(0)),fa(ll),fb(ll),fc(ll),&
fa(ll)+fb(ll)+fc(ll)-1.0_dp,fg(ll)+fh(ll)+fi(ll)
end do
write(100,'(a30,(e12.3e3))')'Norm conservation =',( fa(0)+fb(0)+fc(0) )-( fa(nfiles) + fb(nfiles) + fc(nfiles) )
write(100,'(a30,(e12.3e3))')'Norm conservation =',( fa2(0)+fb2(0)+fc2(0) )-( fa2(nfiles) + fb2(nfiles) + fc2(nfiles) )
write(100,'(a30,(e12.3e3))')'Final norm for state 1 =',fa(nfiles)
write(100,'(a30,(e12.3e3))')'Final norm for state 2 =',fb(nfiles)
write(100,'(a30,(e12.3e3))')'Final norm for state 3 =',fc(nfiles)
write(100,'(a30,(e12.3e3),a8)')'Energy conservation =',(fd(0)+fe(0)+ff(0))-(fd(nfiles)+fe(nfiles)+ff(nfiles)),' hartree'
write(100,'(a30,(e12.3e3),a8)')'Energy conservation =',(fd2(0)+fe2(0)+ff2(0))-(fd2(nfiles)+fe2(nfiles)+ff2(nfiles)),' hartree'
write(100,'(a30,(f20.15),a8)')'Final total energy state 1 =',fd(nfiles),' hartree'
write(100,'(a30,(f20.15),a8)')'Final total energy state 2 =',fe(nfiles),' hartree'
write(100,'(a30,(f20.15),a8)')'Final total energy state 3 =',ff(nfiles),' hartree'
write(100,'(a24,e12.3e3)') 'maximum momentum value =',maxval(maxmomq1)
write(100,'(a24,e12.3e3)') 'maximum momentum value =',maxval(maxmomq2)
!$ trunf = omp_get_wtime()
write(100,'(a38,f9.1,a7)')'Time took to run totaly the program = ', (trunf-truni)/60.0_dp, ' minutes.'
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
complex(kind=dp)  :: dydt(0:n-1),y(0:n-1) !,y0(n)!,ymatrix(n,npoints+1)!,ymatlab(n,npoints)
!real(kind=dp)     :: Te1,Te2,Te3!,enert(npoints+1) ! This is fo saving total energy through time
complex(kind=dp)  :: momq1(0:n-1),momq2(0:n-1),am(0:n-1) ! This is for saving momentum through time
complex(kind=dp)  :: dydt_nac(0:n-1),y_nac(0:n-1)
complex(kind=dp)  :: y_ant(0:n-1)
real(kind=dp)     :: tt1,tt0
h = tstep
t = t0
!-----------------------------------------------------------!
!Checking momentum and saving norm                          !
call HA_calc(t,y,n,dydt,y_nac,dydt_nac) ! evaluating y'(t,y)!
!call HA_calc(t,y,n,dydt) ! evaluating y'(t,y)               !
call momentum_calc_q1(y,momq1,n) ! evaluating dydq          !
call momentum_calc_q2(y,momq2,n) ! evaluating dydq          !
call angular_momentum(y,n,am)                               !
pq1_1(0)=0.0_dp ! momentum in q1 state1                      !
pq2_1(0)=0.0_dp ! momentum in q2 state1                      !
pq1_2(0)=0.0_dp ! momentum in q1 state2                      !
pq2_2(0)=0.0_dp ! momentum in q2 state2                      !
pq1_3(0)=0.0_dp ! momentum in q1 state3                      !
pq2_3(0)=0.0_dp ! momentum in q2 state3                      !
sum1(0)=0.0_dp  ! norm for state 1                           !
sum2(0)=0.0_dp  ! norm for state 2                           !
sum3(0)=0.0_dp  ! norm for state 3                           !
e1(0)=0.0_dp     ! energy of state 1                          !
e2(0)=0.0_dp     ! energy of state 2                          !
e3(0)=0.0_dp     ! energy of state 3                          !
L1(0)=0.0_dp     ! angular momentum of state 1                !
L2(0)=0.0_dp     ! angular momentum of state 2                !
L3(0)=0.0_dp     ! angular momentum of state 3                !
nac0(0)=0.0_dp   ! Change in the population of st 0           !
nac1(0)=0.0_dp   ! Change in the population of st 1           !
nac2(0)=0.0_dp   ! Change in the population of st 2           !
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
    !$omp parallel do shared(y_ant)
    do i=0,n-1
      y_ant(i) = y(i)
    end do
    !$omp end parallel do
    call HA_calc(t,y,n,dydt,y,dydt_nac) ! evaluating y'(t,y)
    call rk4(y,dydt,n,t,h,y,HA_calc,dydt_nac,y_nac) !evaluating y(t+h)
    t=t+h
  end do
  !$OMP parallel do shared(y_f,wfout)
  do i = 0,n-1
    y_f(i,ll) = y_f(i,ll) + 2.0_dp*real( dconjg(y_ant(i)) * y_nac(i) ,dp) / tstep
    wfout(i,ll) = y(i)
  end do
  !$OMP end parallel do
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
  pq1_1(ll)=0.0_dp ! momentum in q1 state1                      !
  pq2_1(ll)=0.0_dp ! momentum in q2 state1                      !
  pq1_2(ll)=0.0_dp ! momentum in q1 state2                      !
  pq2_2(ll)=0.0_dp ! momentum in q2 state2                      !
  pq1_3(ll)=0.0_dp ! momentum in q1 state3                      !
  pq2_3(ll)=0.0_dp ! momentum in q2 state3                      !
  sum1(ll)=0.0_dp  ! norm for state 1                           !
  sum2(ll)=0.0_dp  ! norm for state 2                           !
  sum3(ll)=0.0_dp  ! norm for state 3                           !
  e1(ll)=0.0_dp     ! energy of state 1                          !
  e2(ll)=0.0_dp     ! energy of state 2                          !
  e3(ll)=0.0_dp     ! energy of state 3                          !
  L1(ll)=0.0_dp     ! angular momentum of state 1                !
  L2(ll)=0.0_dp     ! angular momentum of state 2                !
  L3(ll)=0.0_dp     ! angular momentum of state 3                !
  nac0(ll)=0.0_dp   ! Change in the population of st 0           !
  nac1(ll)=0.0_dp   ! Change in the population of st 1           !
  nac2(ll)=0.0_dp   ! Change in the population of st 2           !
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
    nac0(ll)  = nac0(ll)  + (dconjg(y_nac(i))     * y_nac(i)) /tstep     ! EIIIIIIIIIIIIIIIIIIIIIIII
    nac1(ll)  = nac1(ll)  + (dconjg(y_nac(i+s))   * y_nac(i+s)) /tstep   ! EIIIIIIIIIIIIIIIIIIIIIIII
    nac2(ll)  = nac2(ll)  + (dconjg(y_nac(i+2*s)) * y_nac(i+2*s)) /tstep !  EIIIIIIIIIIIIIIIIIIIIIIII
  end do                                                       !
  !------------------------------------------------------------!
  momq1t(ll)=pq1_1(ll)+pq1_2(ll)+pq1_3(ll)
  momq2t(ll)=pq2_1(ll)+pq2_2(ll)+pq2_3(ll)
end do
end subroutine rkdumb
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rk4(y0,dydt,n,t,h,Yout,HA_calc,dydt_nac,Yout_nac)
!subroutine rk4(y0,dydt,n,t,h,Yout,HA_calc)
use global_param
use omp_lib
implicit none
integer n
real(kind=dp) :: t,h,th
complex(kind=dp) :: dydt(0:n-1),y0(0:n-1),Yout(0:n-1),dym(0:n-1),dyt(0:n-1),yt(0:n-1)
complex(kind=dp),dimension(0:n-1) :: dydt_nac,Yout_nac,dym_nac,dyt_nac,yt_nac
external HA_calc 
! 4th order runge-kutta
!h = (tfinal - t0) / ntemp
!t+1= t + h
!k1 = h * y'(t , f(t))
!k2 = h *  y'(t + h/2.0_dp , y(t) + k1/2.0_dp)
!k3 = h *  y'(t + h/2.0_dp , y(t) + k2/2.0_dp)
!k4 = h *  y'(t + h , y(t) + k3)
!y(t+1) = y(t) + 1/6 (k1+ 2*k2 + 2*k3 + k4)
th=t+h/2.0_dp
!dydt = y'(t , y(t)) = k1/h ---- dydt MUST NOT BE PREVIOUSLY MULTIPLIED BY h ----
!$OMP parallel do shared(yt,yt_nac)
!!$OMP parallel do shared(yt)
do1: do i=0,n-1
  yt(i) = y0(i) + (h/2.0_dp) * dydt(i) !Evaluation of y(t) + k1/2.0_dp
  yt_nac(i) = y0(i) + (h/2.0_dp) * dydt_nac(i) !Evaluation of y(t) + k1/2.0_dp
end do do1
!$OMP end parallel do
call HA_calc(th,yt,n,dyt,yt_nac,dyt_nac) !dyt = k2/h ---- Evaluation of y'(t + h/2.0_dp , y(t) + k1/2.0_dp)
!call HA_calc(th,yt,n,dyt) !dyt = k2/h ---- Evaluation of y'(t + h/2.0_dp , y(t) + k1/2.0_dp)
!$OMP parallel do shared(yt,yt_nac)
!!$OMP parallel do shared(yt)
do2: do i=0,n-1
  yt(i) = y0(i) + (h/2.0_dp) * dyt(i)  !Evaluation of y(t) + k2/2.0_dp
  yt_nac(i) = y0(i) + (h/2.0_dp) * dyt_nac(i)  !Evaluation of y(t) + k2/2.0_dp
end do do2
!$OMP end parallel do
call HA_calc(th,yt,n,dym,yt_nac,dym_nac) !dym = k3/h ---- Evaluation of y'(t + h/2.0_dp , y(t) + k2/2.0_dp)
!call HA_calc(th,yt,n,dym) !dym = k3/h ---- Evaluation of y'(t + h/2.0_dp , y(t) + k2/2.0_dp)
!$OMP parallel do shared(yt,dym,yt_nac,dym_nac)
!!$OMP parallel do shared(yt,dym)
do3: do i=0,n-1
  yt(i) = y0(i) + h * dym(i) !Evaluation of y(t) + k3
  yt_nac(i) = y0(i) + h * dym_nac(i) !Evaluation of y(t) + k3
  dym(i) = dyt(i) + dym(i)  ! making dym = (k2 + k3)/h ---- variable dyt free now
  dym_nac(i) = dyt_nac(i) + dym_nac(i)  ! making dym = (k2 + k3)/h ---- variable dyt free now
end do do3
!$OMP end parallel do
call HA_calc(t+h,yt,n,dyt,yt_nac,dyt_nac) !dyt = k4/h ---- Evaluation of y'(t + h/2.0_dp , y(t) + k1/2.0_dp)
!call HA_calc(t+h,yt,n,dyt) !dyt = k4/h ---- Evaluation of y'(t + h/2.0_dp , y(t) + k1/2.0_dp)
!$OMP parallel do shared(Yout,Yout_nac)
!!$OMP parallel do shared(Yout)
do4: do i=0,n-1
  Yout(i) = y0(i) + h/6.0_dp * (dydt(i) + 2.0_dp * dym(i) + dyt(i))
  Yout_nac(i) = h/6.0_dp * (dydt_nac(i) + 2.0_dp * dym_nac(i) + dyt_nac(i))
end do do4
!$OMP end parallel do

end subroutine rk4
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Subroutine computes the derivatives dPsi(t)/dt = H/i Psi
!This derivative routine calculates 
!Written as a simple matrix multiplication between the Hamiltonian matrix H and the amplitudes vector Psi
subroutine HA_calc(t,y,n,dydt,y_nac,dydt_nac)
!subroutine HA_calc(t,y,n,dydt)
use global_param
use omp_lib
implicit none
integer n
real(kind=dp) :: t
complex(kind=dp) :: y(0:n-1),dydt(0:n-1)
complex(kind=dp) :: dydt_nac(0:n-1),y_nac(0:n-1)

!write(*,*) 'omp_get_max_threads= ', omp_get_max_threads ( )
!write(*,*) 'omp_get_num_procs = ', omp_get_num_procs ( )
!write(*,*) 'Time = ', omp_get_wtime ( )

! Evaluating the external electric field at a given time
Et = (E00/freq)*(-(t-t00)/sig**2.0_dp*dsin(freq*(t-t00)+phase)+freq*dcos(freq*(t-t00)+phase))*dexp(-(t-t00)**2.0_dp/(2.0_dp*sig**2.0_dp))

dydt=dcmplx(0.0_dp,0.0_dp)
!$OMP PARALLEL DO shared(dydt)
do i=0,n-1
!  do j=Ha_rowc(i),Ha_rowc(i+1)-1 ! USING CSR VECTORS FOR THE SPARSE HMAILTONIAN MATRIX (CALCULATED IN load_hamiltonian.f90)
!    dydt(i) = dydt(i) + ( Ha_val(j) + (dot_product(orientation,di2_val(j,:)) * Et) ) * y(Ha_row_col(j,1)) / im ! dPsi(t)/dt = H/i Psi
!  end do
  do j=kine_rowc(i),kine_rowc(i+1)-1 ! USING CSR VECTORS FOR THE SPARSE HMAILTONIAN MATRIX (CALCULATED IN load_hamiltonian.f90)
    dydt(i) = dydt(i) + kine_val(j) * y(kine_row_col(j,1)) 
  end do
  do j=pot_rowc(i),pot_rowc(i+1)-1 ! USING CSR VECTORS FOR THE SPARSE HMAILTONIAN MATRIX (CALCULATED IN load_hamiltonian.f90)
    dydt(i) = dydt(i) + pot_val(j) * y(pot_row_col(j,1)) 
  end do
  do j=nac_rowc(i),nac_rowc(i+1)-1 ! USING CSR VECTORS FOR THE SPARSE HMAILTONIAN MATRIX (CALCULATED IN load_hamiltonian.f90)
    dydt(i) = dydt(i) + nac_val(j) * y(nac_row_col(j,1)) 
  end do
  do j=dipx_rowc(i),dipx_rowc(i+1)-1 ! USING CSR VECTORS FOR THE SPARSE HMAILTONIAN MATRIX (CALCULATED IN load_hamiltonian.f90)
  dydt(i) = dydt(i) + orientation(1) * dipx_val(j) * y(dipx_row_col(j,1)) * Et 
  end do
  do j=dipy_rowc(i),dipy_rowc(i+1)-1 ! USING CSR VECTORS FOR THE SPARSE HMAILTONIAN MATRIX (CALCULATED IN load_hamiltonian.f90)
  dydt(i) = dydt(i) + orientation(2) * dipy_val(j) * y(dipy_row_col(j,1)) * Et
  end do
  do j=dipz_rowc(i),dipz_rowc(i+1)-1 ! USING CSR VECTORS FOR THE SPARSE HMAILTONIAN MATRIX (CALCULATED IN load_hamiltonian.f90)
  dydt(i) = dydt(i) + orientation(3) * dipz_val(j) * y(dipz_row_col(j,1)) * Et
  end do
  dydt(i) = dydt(i) / im
end do
!$OMP end parallel do

!============================== EVALUATION OF THE NAC INFLUENCE ONLY ==============================!      
dydt_nac=dcmplx(0.0_dp,0.0_dp)
!$OMP PARALLEL DO shared(dydt_nac)
do i=0,n-1
  do j=nac_rowc(i),nac_rowc(i+1)-1
    dydt_nac(i) = dydt_nac(i) + nac_val(j) * y_nac(nac_row_col(j,1)) / im
  end do
end do
!$OMP end parallel do
!============================== EVALUATION OF THE NAC INFLUENCE ONLY ==============================!      

end subroutine HA_calc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!Subroutine computes the derivatives dPsi(t)/dt = H/i Psi
!!This derivative routine calculates 
!!Written as a simple matrix multiplication between the Hamiltonian matrix H and the amplitudes vector Psi
!subroutine HA_calc_nac_only(n,dydt_nac,y_nac)
!use global_param
!use omp_lib
!implicit none
!integer n
!complex(kind=dp) :: dydt_nac(0:n-1),y_nac(0:n-1)
!
!!============================== EVALUATION OF THE NAC INFLUENCE ONLY ==============================!      
!dydt_nac=dcmplx(0.0_dp,0.0_dp)
!!$OMP PARALLEL DO shared(dydt_nac)
!do i=0,n-1
!  do j=Ha2_rowc(i),Ha2_rowc(i+1)-1
!    dydt_nac(i) = dydt_nac(i) + Ha2_val(j) * y_nac(Ha2_row_col(j,1)) / im
!  end do
!end do
!!$OMP end parallel do
!!============================== EVALUATION OF THE NAC INFLUENCE ONLY ==============================!      
!
!end subroutine HA_calc_nac_only
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine momentum_calc_q1(y,mom,n)
use global_param
use omp_lib
implicit none
integer n
complex(kind=dp) :: mom(0:n-1),y(0:n-1)
mom=dcmplx(0.0_dp,0.0_dp)
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
mom=dcmplx(0.0_dp,0.0_dp)
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
!mom=dcmplx(0.0_dp,0.0_dp)
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
subroutine angular_momentum(y,n,am)
use global_param
use omp_lib
implicit none
integer          :: n
complex(kind=dp) :: y(0:n-1),am(0:n-1)

am=dcmplx(0.0_dp,0.0_dp)

!$OMP PARALLEL DO shared(y,am,am_rowc,am_val,am_row_col)
do i=0,n-1
  do j=am_rowc(i),am_rowc(i+1)-1
    am(i) = am(i) + (-im) * am_val(j) * y(am_row_col(j,1))
  end do
end do
!$OMP end parallel do
end subroutine angular_momentum
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function observable(x,y,n) 
use omp_lib
integer, parameter                           :: dp = kind(1.0d0)
real (kind=dp)                               :: observable
integer,intent(in)                           :: n
integer                                      :: i
complex(kind=dp),dimension(0:n-1),intent(in) :: x,y
observable = 0.0_dp
if (n < 184*146-1) then
  do i = 0,n-1
    observable = observable + dconjg(x(i)) * y(i)
  end do
else
  !$omp parallel do shared(observable)
  do i = 0,n-1
    observable = observable + dconjg(x(i)) * y(i)
  end do
  !$omp end parallel do
end if
!observable = dsqrt(observable)
return
end function observable
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ft_real_signal(b, n, a, fre, c )
integer,parameter :: dp = selected_real_kind(15, 307)
integer i,j,n,c,nu,t
complex(dp) :: ft_real_signal(n)
real(dp) :: a(n),b(n),da,dnu,fre(n)
real(dp),parameter :: pi = 3.141592653589793_dp

!define d(a)
da = a(2) - a(1)
ft = dcmplx(0_dp,0_dp)
!Do the integral
!OMP parallel do shared(ft)
do nu = 1,n
  do t = 1,n
    ft_real_signal(nu) = ft_real_signal(nu) + b(t) * exp(-(dcmplx(0_dp,1_dp)) * 2.0_dp* pi * fre(nu) * a(t)) * da
  end do
end do
!OMP end parallel do

return
end function ft_real_signal
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ft_complex_signal(b, n, a, fre, c )
integer,parameter :: dp = selected_real_kind(15, 307)
integer i,j,n,c,nu,t
complex(dp) :: ft_complex_signal(n),b(n)
real(dp) :: a(n),da,dnu,fre(n)
real(dp),parameter :: pi = 3.141592653589793_dp

!define d(a)
da = a(2) - a(1)
ft = dcmplx(0_dp,0_dp)
!Do the integral
!OMP parallel do shared(ft)
do nu = 1,n
  do t = 1,n
    ft_complex_signal(nu) = ft_complex_signal(nu) + b(t) * exp(-(dcmplx(0_dp,1_dp)) * 2.0_dp* pi * fre(nu) * a(t)) * da
  end do
end do
!OMP end parallel do

return
end function ft_complex_signal
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
