program photoelectron_spectrum
use omp_lib
use lib_array
implicit none
integer,parameter             :: dp = selected_real_kind(p=15,r=307)
integer                       :: i,j,ii,jj,kk,tt,a,ind
character(len=21)             :: fn1
character(len=18)             :: fn2
integer,parameter             :: nang = 512 !size of the angular grid for the photoelectron
integer,parameter             :: nk = 256 !sampling in kinetic energy for the photoelectron
integer,parameter             :: nq1 = 46 !Number of grid points along q1 in the franck-condon region
integer,parameter             :: nq2 = 41 !Number of grid points along q2 in the franck-condon region
integer,parameter             :: n1 = 146 !Number of grid points along q1 in the whole grid
integer,parameter             :: n2 = 184 !Number of grid points along q2 in the whole grid
real*8,parameter              :: ha2ev = 27.211386279352677_dp
real*8,parameter              :: pi = 3.141592653589793_dp
real*8,parameter              :: ener = 21.22_dp / ha2ev ! Frequency of the ionizing pulse in atomic units. 0.779820615611238 = 21.22eV
real*8,parameter              :: e00 = 0.05_dp ! Intensity of the ionizing pulse in atomic units
real*8,parameter              :: zpe_n = 0.00598956_dp ! Zero point energy of the ground electronic state of the neutral (above the neutral minimum energy)
real*8,parameter              :: c_au = 137.0359989995517_dp
real*8,parameter              :: d_omega = 4._dp * pi / Nang 
real*8,parameter              :: n = n1 * n2
character(len=2),parameter    :: orient = 'h4'  !Orientation of the ionizing electric field along one of the C-H bonds
real*8,dimension(3)           :: ori
real*8                        :: ang(nang,2),xmax,xmin,e_pe,k_pe,aux,ti,tf
real*8,dimension(nk,nang,18)  :: z !r0x,r0y,r0z,i0x,i0y,i0z,r1x,r1y,r1z,i1x,i1y,i1z,r2x,r2y,r2z,i2x,i2y,i2z
real*8,dimension(n1,n2)       :: v0neutral,v0cation,v1cation,v2cation,neutg,g
real*8,dimension(nq1,nq2)     :: v0n,v0c,v1c,v2c,nsg,sg
real*8,dimension(n)           :: neut
real*8,dimension(nk)          :: k,e,c0,c1,c2
real*8,dimension(nk*200)      :: k1,c00,c11,c22,p0,p1,p2
real*8,dimension(nq1,nk*200)  :: phot0,phot1,phot2
real*8,dimension(1000)        :: eip!,phot0,phot1,phot2 ! this is if I want to vary the energy of the ionizing pulse
real*8,dimension(n,300)       :: vec0,vec1,vec2
real*8,dimension(300)         :: val0,val1,val2
real*8,dimension(n1*n2)       :: vecst0,vecst1,vecst2 
real*8,dimension(:,:,:,:),allocatable :: x0,x1,x2

write(*,*) 'omp_get_max_threads= ', omp_get_max_threads ( )
write(*,*) 'omp_get_num_procs = ', omp_get_num_procs ( )

!$ ti=omp_get_wtime()
allocate(x0(nq1,nq2,nang*nk,6),x1(nq1,nq2,nang*nk,6),x2(nq1,nq2,nang*nk,6))


if ( orient == 'h1' ) then
  ori(1) = 1._dp
  ori(2) = 1._dp
  ori(3) = 1._dp
elseif (orient == 'h2' ) then
  ori(1) =-1._dp
  ori(2) =-1._dp
  ori(3) = 1._dp
elseif (orient == 'h3' ) then
  ori(1) = 1._dp
  ori(2) =-1._dp
  ori(3) =-1._dp
elseif (orient == 'h4' ) then
  ori(1) =-1._dp
  ori(2) = 1._dp
  ori(3) =-1._dp
endif

xmin = 0.005859375_dp
xmax = 1.505859375_dp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%
!do a linspace(xmin,xmax,nk)                                  !%
do i=1,nk                                                     !%
  k(i) = (xmax-xmin) * (i-1._dp) / (nk-1._dp) + xmin          !%
  e(i) = k(i)**2._dp / 2._dp * 27.211386279352677_dp          !%
end do                                                        !%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%
!do a linspace(xmin,xmax,nk)                                  !%
do i=1,nk*200                                                 !%
  k1(i) = (xmax-xmin) * (i-1._dp) / (nk*200-1_dp) + xmin      !% ! Linear momentum of the photoelectron
end do                                                        !%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%
xmin = 0.7198_dp
xmax = 1.30_dp/2._dp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!%
!do a linspace(xmin,xmax,nk)                                  !%
do i=1,1000                                                   !%
  eip(i) = (xmax-xmin) * (i-1._dp) / (1000._dp-1._dp) + xmin  !% ! Energy of the ionizing pulse.
end do                                                        !%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!% 
!------------------------ OPENING EVERY NEEDED FILE ----------------------
fn2 = 'test_sym_dist3.txt'
open(unit=100,file=fn2,status='old')
do i=1,nang
  read(100,*)ang(i,:)
enddo
close(unit=100)
open(unit=100,file='v1neutral.txt',status='old')
open(unit=101,file='v1f.txt',status='old')
open(unit=102,file='v2f.txt',status='old')
open(unit=103,file='v3f.txt',status='old')
do i=1,n1
  read(100,*)v0neutral(i,:)
  read(101,*)v0cation(i,:)
  read(102,*)v1cation(i,:)
  read(103,*)v2cation(i,:)
enddo
v0n(:,:) = v0neutral(52:97,95:135)
v0c(:,:) = v0cation(52:97,95:135)
v1c(:,:) = v1cation(52:97,95:135)
v2c(:,:) = v2cation(52:97,95:135)
close(unit=100)
close(unit=101)
close(unit=102)
close(unit=103)
open(unit=100,file='eigen-vec-neutral',status='old')
do i=1,n
  read(100,*)neut(i)
enddo
close(unit=100)
do j=1,n2  !PARALELIZÁVEL
  do i=1,n1
    neutg(i,j) = neut( (j-1)*n1+i )
  enddo
enddo
nsg(:,:)=neutg(52:97,95:135)

open(unit=100,file='eigen-vec-cation-st0',status='unknown')
open(unit=101,file='eigen-vec-cation-st1',status='unknown')
open(unit=102,file='eigen-vec-cation-st2',status='unknown')
do i=1,n
  read(100,*)vec0(i,:)
  read(101,*)vec1(i,:)
  read(102,*)vec2(i,:)
enddo
close(unit=100)
close(unit=101)
close(unit=102) 
open(unit=100,file='eigen-values-cation-st0',status='unknown')
open(unit=101,file='eigen-values-cation-st1',status='unknown')
open(unit=102,file='eigen-values-cation-st2',status='unknown')
read(100,*)val0(:)
read(101,*)val1(:)
read(102,*)val2(:)
close(unit=100)
close(unit=101)
close(unit=102)
!-------------------------------------------------------------------------

fn1 = 'pice_10000000_0_0.dat'
!do jj=1,nq2
!  do ii=1,nq1
!    fn1 = 'pice_10000000_0_0.dat'
!    write(fn1(7:9),'(i0.3)') jj+100+74
!    write(fn1(12:13),'(i0.2)') ii+24
!    write(fn1(17:17),'(i1)') 0
!    open(unit=(jj-1)*nq1+ii+100,file=fn1,status='old')
!    do i=1,nang*nk
!      read((jj-1)*nq1+ii+100,*) x1(ii,jj,i,:)
!    enddo
!    close(unit=(jj-1)*nq1+ii+100)
!  enddo
!enddo
!write(*,*)'fim do serial'
!!$OMP parallel do firstprivate(fn1) shared(x0) 
do jj=1,nq2
  do ii=1,nq1
    write(fn1(7:9),'(i0.3)') jj+100+74
    write(fn1(12:13),'(i0.2)') ii+24
    write(fn1(17:17),'(i1)') 0
    open(unit=(jj-1)*nq1+ii+100,file=fn1,status='old')
    write(fn1(17:17),'(i1)') 1
    open(unit=((jj-1)*nq1+ii)+2000,file=fn1,status='old')
    write(fn1(17:17),'(i1)') 2
    open(unit=((jj-1)*nq1+ii)+5000,file=fn1,status='old')
    do i=1,nang*nk
      read((jj-1)*nq1+ii+100,*) x0(ii,jj,i,:)
      read(((jj-1)*nq1+ii)+2000,*) x1(ii,jj,i,:)
      read(((jj-1)*nq1+ii)+5000,*) x2(ii,jj,i,:)
    enddo
    close(unit=(jj-1)*nq1+ii+100)
    close(unit=((jj-1)*nq1+ii)+2000)
    close(unit=((jj-1)*nq1+ii)+5000)
  enddo
enddo
!!$OMP end parallel do
!$ tf=omp_get_wtime()
write(*,*)'time to load everything = ',tf-ti
!aux=0._dp
!do jj=1,nq2
!  do ii=1,nq1
!    do i=1,nang*nk
!      aux = aux + x0(ii,jj,i,1) - x1(ii,jj,i,1)
!    enddo
!  enddo
!enddo
!write(*,*)aux
!stop

phot0(:,:) = 0._dp
phot1(:,:) = 0._dp
phot2(:,:) = 0._dp
!do tt = 1,1000 !sampling of possible energy of the ionizing pulse
!!$OMP parallel do default(shared) private(ii,i,j,kk,ind,aux,g,sg,z,c0,c1,c2,c00,c11,c22,e_pe,k_pe) reduction(+:phot0,phot1,phot2)
!!$OMP shared(nq1,nq2,nk,nang,d_omega,e00,freq,c_au,k,k1,ori,n1,n2,vec0,vec1,vec2,zpen,val0,val2,lav2,v0c.v1c,v2c,v0n,neutsg)
!!$ private(ii,i,j,kk,ind,aux,g,sg,z,c0,c1,c2,c00,c11,c22,e_pe,k_pe)
!!$OMP parallel do default(shared) private(ii,i,j,kk,ind,aux,g,sg,z,c0,c1,c2,c00,c11,c22,e_pe,k_pe)
  do jj=1,nq2
    do ii=1,nq1
      do i=1,nk
        do j=1,nang
          z(i,j, 1) = x0(ii,jj, (i-1)*nang+j, 1 )
          z(i,j, 2) = x0(ii,jj, (i-1)*nang+j, 2 )
          z(i,j, 3) = x0(ii,jj, (i-1)*nang+j, 3 )
          z(i,j, 4) = x0(ii,jj, (i-1)*nang+j, 4 )
          z(i,j, 5) = x0(ii,jj, (i-1)*nang+j, 5 )
          z(i,j, 6) = x0(ii,jj, (i-1)*nang+j, 6 )
          z(i,j, 7) = x1(ii,jj, (i-1)*nang+j, 1 )
          z(i,j, 8) = x1(ii,jj, (i-1)*nang+j, 2 )
          z(i,j, 9) = x1(ii,jj, (i-1)*nang+j, 3 )
          z(i,j,10) = x1(ii,jj, (i-1)*nang+j, 4 )
          z(i,j,11) = x1(ii,jj, (i-1)*nang+j, 5 )
          z(i,j,12) = x1(ii,jj, (i-1)*nang+j, 6 ) 
          z(i,j,13) = x2(ii,jj, (i-1)*nang+j, 1 )
          z(i,j,14) = x2(ii,jj, (i-1)*nang+j, 2 )
          z(i,j,15) = x2(ii,jj, (i-1)*nang+j, 3 )
          z(i,j,16) = x2(ii,jj, (i-1)*nang+j, 4 )
          z(i,j,17) = x2(ii,jj, (i-1)*nang+j, 5 )
          z(i,j,18) = x2(ii,jj, (i-1)*nang+j, 6 )
        enddo
      enddo
      c0 = 0._dp; c1 = 0._dp; c2 = 0._dp
      do i = 1,nk   
        do j = 1,nang
          c0(i) = c0(i) + ((z(i,j, 1)**2._dp+z(i,j, 2)**2._dp)*ori(1) + (z(i,j, 3)**2._dp+z(i,j, 4)**2._dp)*ori(2) + &
          (z(i,j, 5)**2._dp+z(i,j, 6)**2._dp)*ori(3) ) * d_omega / (-e00)**2._dp
          c1(i) = c1(i) + ((z(i,j, 7)**2._dp+z(i,j, 8)**2._dp)*ori(1) + (z(i,j, 9)**2._dp+z(i,j,10)**2._dp)*ori(2) + &
          (z(i,j,11)**2._dp+z(i,j,12)**2._dp)*ori(3) ) * d_omega / (-e00)**2._dp 
          c2(i) = c2(i) + ((z(i,j,13)**2._dp+z(i,j,14)**2._dp)*ori(1) + (z(i,j,15)**2._dp+z(i,j,16)**2._dp)*ori(2) + &
          (z(i,j,17)**2._dp+z(i,j,18)**2._dp)*ori(3) ) * d_omega / (-e00)**2._dp 
        end do
        c0(i) = c0(i) * k(i) / c_au * (2._dp * pi * ener ) ! total cross section for cation electronic ground state   ! CHANGE HERE IF WANNA VARY ENERGY OF THE PULSE
        c1(i) = c1(i) * k(i) / c_au * (2._dp * pi * ener ) ! total cross section for cation electronic state 1        ! CHANGE HERE IF WANNA VARY ENERGY OF THE PULSE
        c2(i) = c2(i) * k(i) / c_au * (2._dp * pi * ener ) ! total cross section fro cation electronic state 2        ! CHANGE HERE IF WANNA VARY ENERGY OF THE PULSE
      end do
      c00(:) = interp1d(k,c0,k1) !Interpolation for better resolution
      c11(:) = interp1d(k,c1,k1) !Interpolation for better resolution
      c22(:) = interp1d(k,c2,k1) !Interpolation for better resolution
  
!Selecting the best cross section element for the correct energy of the photoelectron of electronic state 0
      do kk = 1,300
!reshape each eigenvector
        do i = 1,n1         !paralelizável
          do j = 1,n2
           g(i,j) = vec0( (i-1)*n2+j, kk )
          enddo
        enddo
        sg(:,:) = g(52:97,95:135)
!get the index of momentum of the photoelectron
        if (ener - ( (v0c(ii,jj)+val0(kk)) - (v0n(ii,jj) + zpe_n) )  > 0._dp ) then                                              ! CHANGE HERE IF WANNA VARY ENERGY OF THE PULSE
          e_pe = ener - ( (v0c(ii,jj)+val0(kk)) - (v0n(ii,jj) + zpe_n) )                                                    ! CHANGE HERE IF WANNA VARY ENERGY OF THE PULSE
          k_pe = sqrt( 2._dp*e_pe ) ! value of the momentum of the photoelectron for that vibrational state of the cation
          aux = 20._dp
          do i = 1,nk*200
            if ( abs( k1(i) - k_pe ) < aux ) then
              aux = abs( k1(i) - k_pe )
              ind = i
            endif
          enddo
          phot0(jj,ind) = phot0(jj,ind) + c00(ind) * sqrt(nsg(ii,jj) * sg(ii,jj)**2._dp)                                         ! CHANGE HERE IF WANNA VARY ENERGY OF THE PULSE
        endif
      enddo
!Selecting the best cross section element for the correct energy of the photoelectron of electronic state 1
      do kk = 1,64
        do i = 1,n1         !paralelizável
          do j = 1,n2
           g(i,j) = vec1( (i-1)*n2+j, kk )
          enddo
        enddo
        sg(:,:) = g(52:97,95:135)
!get the index of momentum of the photoelectron
        if (ener - ( (v1c(ii,jj)+val1(kk)) - (v0n(ii,jj) + zpe_n) )  > 0._dp ) then                                              ! CHANGE HERE IF WANNA VARY ENERGY OF THE PULSE
          e_pe = ener - ( (v1c(ii,jj)+val1(kk)) - (v0n(ii,jj) + zpe_n) )                                                    ! CHANGE HERE IF WANNA VARY ENERGY OF THE PULSE
          k_pe = sqrt( 2._dp*e_pe )
          aux = 20._dp
          do i = 1,nk*200
            if ( abs( k1(i) - k_pe ) < aux ) then
              aux = abs( k1(i) - k_pe )
              ind = i
            endif
          enddo
          phot1(jj,ind) = phot1(jj,ind) + c11(ind) * sqrt(nsg(ii,jj) * sg(ii,jj)**2._dp)                                         ! CHANGE HERE IF WANNA VARY ENERGY OF THE PULSE
        endif 
      enddo 
!Selecting the best cross section element for the correct energy of the photoelectron of electronic state 2
      do kk = 1,17
        do i = 1,n1         !paralelizável
          do j = 1,n2
           g(i,j) = vec2( (i-1)*n2+j, kk )
          enddo
        enddo
        sg(:,:) = g(52:97,95:135)
!get the index of momentum of the photoelectron
        if (ener - ( (v0c(ii,jj)+val0(kk)) - (v0n(ii,jj) + zpe_n) ) > 0._dp ) then                                              ! CHANGE HERE IF WANNA VARY ENERGY OF THE PULSE
          e_pe = ener - ( (v2c(ii,jj)+val2(kk)) - (v0n(ii,jj) + zpe_n) )                                                    ! CHANGE HERE IF WANNA VARY ENERGY OF THE PULSE
          k_pe = sqrt( 2._dp*e_pe )
          aux = 20._dp
          do i = 1,nk*200
            if ( abs( k1(i) - k_pe ) < aux ) then
              aux = abs( k1(i) - k_pe )
              ind = i
            endif
          enddo
          phot2(jj,ind) = phot2(jj,ind) + c22(ind) * sqrt(nsg(ii,jj) * sg(ii,jj)**2._dp)                                                      ! CHANGE HERE IF WANNA VARY ENERGY OF THE PULSE
        endif
      enddo 
    enddo
!    write(*,*)real(jj)/real(nq2)*100._dp,'%'
  enddo 
!!$OMP end parallel do
!enddo 

do i=1,nk*200
  do j=1,nq2
    p0(i) = p0(i) + phot0(j,i)
    p1(i) = p1(i) + phot1(j,i)
    p2(i) = p2(i) + phot2(j,i)
  enddo
enddo



if (ori(1) == 1._dp .and. ori(2) == 1._dp .and. ori(3) == 1._dp ) then
  fn1='o-h1-spec_00.00eV.dat'
  write(fn1(11:15),'(f5.2)')ener * ha2ev
  open(unit=100,file=fn1,status='unknown')
  do i=1,nk*200
    write(100,'(4es24.15e3)')k1(i)**2._dp / 2._dp * ha2ev,p0(i),p1(i),p2(i)
  enddo
elseif(ori(1) ==-1._dp .and. ori(2) ==-1._dp .and. ori(3) == 1._dp ) then
  fn1='o-h2-spec_00.00eV.dat'
  write(fn1(11:15),'(f5.2)')ener * ha2ev
  open(unit=100,file=fn1,status='unknown')
  do i=1,nk*200
    write(100,'(4es24.15e3)')k1(i)**2._dp / 2._dp * ha2ev,p0(i),p1(i),p2(i)
  enddo
elseif(ori(1) == 1._dp .and. ori(2) ==-1._dp .and. ori(3) ==-1._dp ) then
  fn1='o-h3-spec_00.00eV.dat'
  write(fn1(11:15),'(f5.2)')ener * ha2ev
  open(unit=100,file=fn1,status='unknown')
  do i=1,nk*200
    write(100,'(4es24.15e3)')k1(i)**2._dp / 2._dp * ha2ev,p0(i),p1(i),p2(i)
  enddo
elseif(ori(1) ==-1._dp .and. ori(2) == 1._dp .and. ori(3) ==-1._dp ) then
  fn1='o-h4-spec_00.00eV.dat'
  write(fn1(11:15),'(f5.2)')ener * ha2ev
  open(unit=100,file=fn1,status='unknown')
  do i=1,nk*200
    write(100,'(4es24.15e3)')k1(i)**2._dp / 2._dp * ha2ev,p0(i),p1(i),p2(i)
  enddo
endif





end program photoelectron_spectrum

