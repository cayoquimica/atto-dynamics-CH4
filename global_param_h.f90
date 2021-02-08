module global_param
implicit none
integer, parameter                         :: dp = selected_real_kind(15, 307)
integer                                    :: i,j,k,l,ll,ij
integer, parameter                         :: NA=5 !Number of atoms
integer, parameter                         :: Nst=3 !number of electronic states
integer, parameter                         :: Nq1=146 !126 !number of points of the grid along q1
integer, parameter                         :: Nq2=184 !179 !number of points of the grid along q2
integer, parameter                         :: s=Nq1*Nq2 !counter to be used as index - means the size of a 1 state matrix, which is Nq1*Nq2
integer,parameter                          :: q1i0=27
integer,parameter                          :: q2i0=20
integer,parameter                          :: q1_initial = 25
integer,parameter                          :: q1_final   = 70
integer,parameter                          :: q2_initial = 75
integer,parameter                          :: q2_final   = 115
integer,parameter                          :: nsamples = 3 !Number of randon orientations to sample for the photoionization around each bond (the *4 is for each bond)
integer,parameter                          :: nfiles = 2000 !Number of snapshots in time to save
integer,parameter                          :: npoints = 200000 !Number of time steps to take in the simulation
complex(kind=dp), parameter                :: im=dcmplx(0.0_dp,1.0_dp) !imaginary unity
real(kind=dp),parameter                    :: t0 = 0.0_dp !Initial time
real(kind=dp),parameter                    :: tf = 2000.0_dp !Final time
real(kind=dp),parameter                    :: tstep = (tf-t0)/npoints !Time step
real(kind=dp),parameter                    :: pi = 3.141592653589793_dp
!real(kind=dp),parameter                    :: coneAng = 10.0_dp * pi / 180.0_dp !Angle around each bond to sample ramdonly for initial photoionizations
real(kind=dp),parameter                    :: sq1=0.08_dp !step in q1 in atomic units
real(kind=dp),parameter                    :: sq2=0.07_dp !step in q1 in atomic units
real(kind=dp),parameter                    :: saw2au=1822.889950851334_dp !mass conversion factor from Standard atomic weight to atomic units
real(kind=dp),parameter                    :: car=12.011_dp*saw2au !Carbon mass in atomic units
real(kind=dp),parameter                    :: hi=1.00794_dp*saw2au !Hidrogen mass in atomic units
real(kind=dp),parameter                    :: mtotal = car + 4.0_dp * hi ! total mass of CH4
real(kind=dp),parameter                    :: t00 = 800.0_dp !time where the pulse is centered
real(kind=dp),parameter                    :: phase = 0.0_dp !phase factor for the pulse related to the gaussian envelope
real(kind=dp),parameter                    :: freq =0.056937_dp !0.512285550500502 is the ionizating pulse !0.056937_dp !frequency of the pulse, in a.u. - 800 nm of wavelength
real(kind=dp),parameter                    :: sig = 100.0_dp  ! 50 approx 1200 attoseconds - width of the gaussian envelop
real(kind=dp),parameter                    :: E00 = 0.00_dp !0.05_dp !Electric field intensity
! The harmonics of 800nm pulse: 9 = 13.948222315188689 eV, 11 = 17.047827274119506, 13 = 20.147432233050328
! The fundamental frequency of 800nm pulse = 0.056954190556670 atomic units ------ 13.0_dp/27.21138628_dp if we want just to write in eV
real(kind=dp),parameter                    :: fre_800 = 0.056954190556670_dp
real(kind=dp),parameter                    :: e_ip =09.0_dp * fre_800 !Energy of the ionizing pulse in atomic units approx 21.22 ev - the Helium ressonance band
real(kind=dp),dimension(3), parameter      :: ori=[ 1.0_dp, 1.0_dp, 1.0_dp] !Orientation of the electric field of the probing pulse
real(kind=dp)                              :: mass(3*NA) !3N Vector with the masses of the atoms
real(kind=dp)                              :: q1(3*NA),q2(3*NA),q1i(3*NA),q2i(3*NA),ai,bi,aii,bii !Conversion vector from internal coordinates q1 and q2 to cartesian coordinates
real(kind=dp)                              :: co1(0:Nst*Nq1*Nq2-1),co2(0:Nst*Nq1*Nq2-1)
real(kind=dp),dimension(:,:),allocatable   :: pot1,pot2,pot3,e1neut !Matrices with each state potential energy
real(kind=dp),dimension(:,:),allocatable   :: pdm1x,pdm2x,pdm3x,pdm1y,pdm2y,pdm3y,pdm1z,pdm2z,pdm3z ! Permanent dipole moments
real(kind=dp),dimension(:,:),allocatable   :: tdm21x,tdm31x,tdm32x,tdm21y,tdm31y,tdm32y,tdm21z,tdm31z,tdm32z ! Transition dipole moments
real(kind=dp),dimension(:,:),allocatable   :: Ha! Hamiltonian matrix
real(kind=dp),dimension(:,:,:),allocatable :: ay! Hamiltonian matrix
real(kind=dp),dimension(:,:),allocatable   :: ham,ax ! Variable to store the Hamiltonian matrix temporary
real(kind=dp),dimension(:,:),allocatable   :: moq1,moq2! Momentum matrix
complex(kind=dp),allocatable               :: wfout(:,:)
real(kind=dp),allocatable                  :: y_f(:,:),coh(:,:),cohe(:,:)
real(kind=dp),allocatable       :: wfinal2(:,:),y_f2(:,:),cohe2(:,:)
real(kind=dp)                              :: mass1,mass2,mass3 !Reduced masses to be used in the second derivative
real(kind=dp)                              :: Et !Final electrical field of the pulse
real(kind=dp),dimension(3)                 :: orientation,u !Orientation of the electric field of the ionizing pulse
real(kind=dp),allocatable                  :: orie(:,:) !Vector with the random orientations of the electric field of the ionizing pulse
real(kind=dp)                              :: ind1,ind2,ind3,const1,const2,const3,E_init
complex(kind=dp),allocatable               :: pia(:,:),nwf0(:),pice(:),wf0(:)
complex(kind=dp)                           :: phote( Nq1,Nq2, Nst, 3 )
real(kind=dp),dimension(0:nfiles)          :: e1,e2,e3,L1,L2,L3,sum1,sum2,sum3,pq1_1,pq2_1,pq1_2,pq2_2,pq1_3,pq2_3,Te1,Te2,Te3
real(kind=dp),dimension(0:nfiles)          :: fa,fb,fc,fd,fe,ff,fg,fh,fi,fj,fk,fl,fm,fn,fo,fp,fq,fr
real(kind=dp),dimension(0:nfiles)          :: fa2,fb2,fc2,fd2,fe2,ff2,fg2,fh2,fi2,fj2,fk2,fl2,fm2,fn2,fo2,fp2,fq2,fr2
real(kind=dp)                              :: momq1t(nfiles),momq2t(nfiles),maxmomq1(nsamples),maxmomq2(nsamples) ! This is for saving norm through time
real(kind=dp),dimension(0:nfiles)          :: nac0,nac1,nac2
integer                                    :: k_Ha,k_di2,k_moq1,k_moq2,k_am,k_nac,k_pot,k_kine,k_dipx,k_dipy,k_dipz
real(kind=dp),allocatable,dimension(:)     :: Ha_val,dipx_val,dipy_val,dipz_val,am_val,moq1_val,moq2_val,pot_val,kine_val,nac_val !CSR vectors for sparse matrix multiplication
real(kind=dp),allocatable                  :: di2_val(:,:) !CSR vectors for sparse matrix multiplication
integer,allocatable                        :: Ha_rowc(:),Ha_row_col(:,:) !CSR vectors for sparse matrix multiplication
integer,allocatable                        :: moq1_rowc(:),moq1_row_col(:,:),moq2_rowc(:),moq2_row_col(:,:) !CSR vectors for sparse matrix multiplication
integer,allocatable                        :: kine_rowc(:),kine_row_col(:,:) !CSR vectors for sparse matrix multiplication
integer,allocatable :: dipx_rowc(:),dipx_row_col(:,:),dipy_rowc(:),dipy_row_col(:,:),dipz_rowc(:),dipz_row_col(:,:) !CSR vectors for sparse matrix multiplication
integer,allocatable                        :: pot_rowc(:),pot_row_col(:,:),nac_rowc(:),nac_row_col(:,:) !CSR vectors for sparse matrix multiplication
integer,allocatable                        :: am_rowc(:),am_row_col(:,:) !CSR vectors for sparse matrix multiplication
integer(kind=dp),dimension(0:s*Nst-1)      :: lol
real(kind=dp)                              :: norma
real(kind=dp),parameter                    :: const = 1_dp !sq1*sq2*0.0001_dp

contains

subroutine load_data
allocate (e1neut(Nq1,Nq2))
allocate (pot1(Nq1,Nq2))
allocate (pot2(Nq1,Nq2))
allocate (pot3(Nq1,Nq2))
allocate ( pdm1x(Nq1,Nq2),pdm2x(Nq1,Nq2),pdm3x(Nq1,Nq2),pdm1y(Nq1,Nq2),pdm2y(Nq1,Nq2),pdm3y(Nq1,Nq2) )
allocate ( pdm1z(Nq1,Nq2),pdm2z(Nq1,Nq2),pdm3z(Nq1,Nq2) )
!allocate (pdm1x(Nq1,Nq2))
!allocate (pdm2x(Nq1,Nq2))
!allocate (pdm3x(Nq1,Nq2))
allocate ( tdm21x(Nq1,Nq2),tdm31x(Nq1,Nq2),tdm32x(Nq1,Nq2),tdm21y(Nq1,Nq2),tdm31y(Nq1,Nq2),tdm32y(Nq1,Nq2) )
allocate ( tdm21z(Nq1,Nq2),tdm31z(Nq1,Nq2),tdm32z(Nq1,Nq2) )
!allocate (tdm21x(Nq1,Nq2))
!allocate (tdm31x(Nq1,Nq2))
!allocate (tdm32x(Nq1,Nq2))

ll=0
do k=1,Nst
  do j=1,Nq2
    do i=1,Nq1
      co1(ll)=(Nq1-Nq1/2.0_dp-i)*sq1+sq1/2.0_dp!(Nq1-Nq11/2.0_dp-i)*sq1+sq1/2.0_dp
      co2(ll)=(j-(114.0_dp)-1.0_dp)*sq2
      ll=ll+1
    end do
  end do
end do

!-------------------------------------------------------------------!
! Loading electronic structure data: Energy and dipole moments      !
open(unit=99,file='v1neutralf.txt',status='old')                     !
open(unit=1,file='v1f.txt',status='old')                            !
open(unit=2,file='v2f.txt',status='old')                            !
open(unit=3,file='v3f.txt',status='old')                            !
do i=1,Nq1                                                          !
  read(99,*) e1neut(i,:)                                            !
  read(1,*) pot1(i,:)                                               !
  read(2,*) pot2(i,:)                                               !
  read(3,*) pot3(i,:)                                               !
end do                                                              !
 close(unit=1)                                                      !
 close(unit=2)                                                      !
 close(unit=3)                                                      !
                                                                    !
open(unit=11,file='pdm1xf.txt',status='old')                        !
open(unit=12,file='pdm2xf.txt',status='old')                        !
open(unit=13,file='pdm3xf.txt',status='old')                        !
open(unit=14,file='pdm1yf.txt',status='old')                        !
open(unit=15,file='pdm2yf.txt',status='old')                        !
open(unit=16,file='pdm3yf.txt',status='old')                        !
open(unit=17,file='pdm1zf.txt',status='old')                        !
open(unit=18,file='pdm2zf.txt',status='old')                        !
open(unit=19,file='pdm3zf.txt',status='old')                        !
open(unit=21,file='tdm21xf.txt',status='old')                       !
open(unit=22,file='tdm31xf.txt',status='old')                       !
open(unit=23,file='tdm32xf.txt',status='old')                       !
open(unit=24,file='tdm21yf.txt',status='old')                       !
open(unit=25,file='tdm31yf.txt',status='old')                       !
open(unit=26,file='tdm32yf.txt',status='old')                       !
open(unit=27,file='tdm21zf.txt',status='old')                       !
open(unit=28,file='tdm31zf.txt',status='old')                       !
open(unit=29,file='tdm32zf.txt',status='old')                       !
do i=1,Nq1                                                          !
  read(11,*) pdm1x(i,:)                                             !
  read(12,*) pdm2x(i,:)                                             !
  read(13,*) pdm3x(i,:)                                             !
  read(14,*) pdm1y(i,:)                                             !
  read(15,*) pdm2y(i,:)                                             !
  read(16,*) pdm3y(i,:)                                             !
  read(17,*) pdm1z(i,:)                                             !
  read(18,*) pdm2z(i,:)                                             !
  read(19,*) pdm3z(i,:)                                             !
  read(21,*) tdm21x(i,:)                                            !
  read(22,*) tdm31x(i,:)                                            !
  read(23,*) tdm32x(i,:)                                            !
  read(24,*) tdm21y(i,:)                                            !
  read(25,*) tdm31y(i,:)                                            !
  read(26,*) tdm32y(i,:)                                            !
  read(27,*) tdm21z(i,:)                                            !
  read(28,*) tdm31z(i,:)                                            !
  read(29,*) tdm32z(i,:)                                            !
end do                                                              !
 close(unit=11)                                                     !
 close(unit=12)                                                     !
 close(unit=13)                                                     !
 close(unit=14)                                                     !
 close(unit=15)                                                     !
 close(unit=16)                                                     !
 close(unit=17)                                                     !
 close(unit=18)                                                     !
 close(unit=19)                                                     !
 close(unit=21)                                                     !
 close(unit=22)                                                     !
 close(unit=23)                                                     !
 close(unit=24)                                                     !
 close(unit=25)                                                     !
 close(unit=26)                                                     !
 close(unit=27)                                                     !
 close(unit=28)                                                     !
 close(unit=29)                                                     !
!-------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------------------------------------!
!Defining the coefficients for coordinates q1 and q2 in cartesian                                                                   !
q1(1)=0.005468383941194_dp;q1(2)=-0.005539667067822_dp;q1(3)=0.099940066929625_dp;q1(4)=0.352088100072972_dp;q1(5)=0.315861943182830_dp  !
q1(6)=-0.495856225504668_dp;q1(7)=-0.382798634963684_dp;q1(8)=-0.284702415219885_dp;q1(9)=-0.495881834018862_dp                         !
q1(10)=-0.056657307983618_dp;q1(11)=0.123960758550165_dp;q1(12)=-0.097711469008833_dp;q1(13)=0.022204480454270_dp                       !
q1(14)=-0.089107486989717_dp;q1(15)=-0.101474677166113_dp                                                                             !
                                                                                                                                    !
q2(1)=0.005503122218018_dp;q2(2)=-0.005574858175503_dp;q2(3)=0.000798651732665_dp;q2(4)=-0.173145514483078_dp;q2(5)=-0.209601800471665_dp!
q2(6)=0.414670087187852_dp;q2(7)=0.142239888853687_dp;q2(8)=0.240959271431981_dp;q2(9)=0.414644315993849_dp;q2(10)=-0.174783204984944_dp !
q2(11)=0.242514205248471_dp;q2(12)=-0.417522164824840_dp;q2(13)=0.140111513551201_dp;q2(14)=-0.207439525935984_dp                       !
q2(15)=-0.421309279015076_dp                                                                                                         !
!Defining the coefficients of the inverse transformation - from internal to cartesian                                               !
!This is just the inverse of the above (s^-1):                                                                                      !
!                                                                                                                                   !
!      | q1(1)  q1(2)  q1(3) ...  q1(i) |                                                                                           !
!  s = |                                |                                                                                           !
!      | q2(1)  q2(2)  q2(3) ...  q2(i) |                                                                                           !
!                                                                                                                                   !
q1i(1)=0.0_dp;q1i(2)=0.0_dp;q1i(3)=0.1404666820458888_dp;q1i(4)=0.327739015711321_dp;q1i(5)=0.327739015711321_dp                         !
q1i(6)=-0.382350835749710_dp;q1i(7)=-0.327739015711321_dp;q1i(8)=-0.327739015711321_dp;q1i(9)=-0.382350835749710_dp                     !
q1i(10)=-0.249049347411151_dp;q1i(11)=0.249049347411151_dp;q1i(12)=-0.454576619283899_dp;q1i(13)=0.249049347411151_dp                   !
q1i(14)=-0.249049347411151_dp;q1i(15)=-0.454576619283899_dp                                                                           !
                                                                                                                                    !
q2i(1)=0.0_dp;q2i(2)=0.0_dp;q2i(3)=0.075292919338856_dp;q2i(4)=-0.016803984583966_dp;q2i(5)=-0.016803984583966_dp                        !
q2i(6)=0.214461491870736_dp;q2i(7)=0.016803984583966_dp;q2i(8)=0.016803984583966_dp;q2i(9)=0.214461491870736_dp                         !
q2i(10)=-0.325973953345556_dp;q2i(11)=0.325973953345556_dp;q2i(12)=-0.663071158209507_dp;q2i(13)=0.325973953345556_dp                   !
q2i(14)=-0.325973953345556_dp;q2i(15)=-0.663071158209507_dp                                                                           !
!-----------------------------------------------------------------------------------------------------------------------------------!
!-----------------------------------!
!Defining a 3N mass vector          !
mass(1)=car;mass(2)=car;mass(3)=car !
do i=4,3*NA                         !
  mass(i)=hi                        !
end do                              !
!-----------------------------------!
!-----------------------------------------------------------------------------!
!Defining the reduced mass that is going to be used in the second derivatives !
!The coeficients ai and bi are going to be used in first derivatives          !
mass1=0.0_dp;mass2=0.0_dp;mass3=0.0_dp                                              !
ai=0.0_dp; bi=0.0_dp; aii=0.0_dp; bii=0.0_dp                                        !
do i=1,3*NA                                                                   !
  mass1=mass1 + q1(i)**2.0_dp / mass(i)                                         !
  mass2=mass2 + q2(i)**2.0_dp / mass(i)                                         !
  mass3=mass3 + q1(i)*q2(i)*2.0_dp / mass(i)                                    !
  ai=ai+q1(i) !NAC values are already divided by mass when loaded             !
  bi=bi+q2(i) !NAC values are already divided by mass when loaded             !
  aii=aii+q1i(i) !NAC values are already divided by mass when loaded          !
  bii=bii+q2i(i) !NAC values are already divided by mass when loaded          !
end do                                                                        !
!-----------------------------------------------------------------------------!

end subroutine load_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine generate_initial_wf
use omp_lib
real(kind=dp) :: soma,n
n = Nq1 * Nq2 * Nst 
!Now read the photoinization coefficients and project them into the initial neutral ground state wave function
pice=dcmplx(0.0_dp,0.0_dp)
!Projecting the PICEs on the ionizing electric field and vectorizing it
k = 0 
!$OMP parallel do shared(pice)
do j=1,Nq2
  do i=1,Nq1
    k = (j-1)*Nq1 + (i-1) ! We need to put -1 in 'i' because 'i' starts in 1 instead of 0 and we want k to start in 0. The same for 'j'
    pice(k)     = - dsqrt(2.0_dp) * 0.05_dp * dot_product( orientation , phote(i,j,1,:) )
    pice(k+s)   = - dsqrt(2.0_dp) * 0.05_dp * dot_product( orientation , phote(i,j,2,:) )
    pice(k+2*s) = - dsqrt(2.0_dp) * 0.05_dp * dot_product( orientation , phote(i,j,3,:) )
!    k = k + 1
  end do
end do
!$OMP end parallel do
!Projecting the photoionization coeficients into the neutral eigen state
do i=0,s-1
  wf0(i)     = nwf0(i) * pice(i)
  wf0(i+s)   = nwf0(i) * pice(i+s)
  wf0(i+2*s) = nwf0(i) * pice(i+2*s)
end do
!--------------------------------------------!
!normalizing                                 !
soma=0.0_dp                                    !
do i=0,n-1                                   !
  soma=soma+dconjg(wf0(i))*wf0(i) * const    !   
end do                                       ! 
!wf0=wf0/sqrt(soma)                           !
norma=soma                                   !
!--------------------------------------------!
end subroutine generate_initial_wf
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Subroutine that calculates the photoionization amplitudes for a given geometry and a given electric field
!Some parts of this function assumes that the number of electronic states is 3
subroutine p_i_a(i1,i2,pic)
implicit none
complex(kind=dp),dimension(Nst,3) :: pic 
integer                           :: i1,i2,ii,jj,kk
character(len=30)                 :: fname00
character(len=199)                :: fpath0,fpath1
integer                           :: file1,file2,file3,file4,file5
integer,parameter                 :: nang=512
integer,parameter                 :: nk=256
real(kind=dp),allocatable         :: vec1(:,:),vec2(:,:),vec3(:,:)
real(kind=dp)                     :: k0,k1,k2,e0,e1,e2,p_e(nk)
real(kind=dp)                     :: phi(nang),theta(nang),domega
complex(kind=dp),allocatable      :: r0(:,:,:),r1(:,:,:),r2(:,:,:),coef0(:,:),coef1(:,:),coef2(:,:)
real(kind=dp)                     :: ip0,ip1,ip2
integer                           :: aux1(Nst)
real(kind=dp)                     :: aux2(Nst),aux3(Nst)

allocate(vec1(nang*nk,Nst*2),vec2(nang*nk,Nst*2),vec3(nang*nk,Nst*2))
allocate(r0(nk,nang,Nst),r1(nk,nang,Nst),r2(nk,nang,Nst),coef0(nk,Nst),coef1(nk,Nst),coef2(nk,Nst))
call getcwd(fpath0) !getting the working directory path
fname00='pice_10000000_neut_cat_0_0.dat'
write(fname00(7:9),'(i0.3)') i2+100
write(fname00(12:13),'(i0.2)') i1
fpath1="~/pice_files/"//fname00
open(newunit=file1,file=fpath1,status='old')
write(fname00(26:26),'(i1)') 1
fpath1="~/pice_files/"//fname00
open(newunit=file2,file=fpath1,status='old')
write(fname00(26:26),'(i1)') 2
fpath1="~/pice_files/"//fname00
open(newunit=file3,file=fpath1,status='old')
!now read the x, y and z components of the photoionization coupling elements.
do ii=1,nk*nang
  read(file1,*) vec1(ii,:)!,vec1(i,2),vec1(i,3),vec1(i,4),vec1(i,5),vec1(i,6)
  read(file2,*) vec2(ii,:)!,vec2(i,2),vec2(i,3),vec2(i,4),vec2(i,5),vec2(i,6)
  read(file3,*) vec3(ii,:)!,vec3(i,2),vec3(i,3),vec3(i,4),vec3(i,5),vec3(i,6)
end do
do ii=1,nk
  do jj=1,nang
    r0(ii,jj,1)=dcmplx( vec1((ii-1)*nang+jj,1) , vec1((ii-1)*nang+jj,2) )
    r0(ii,jj,2)=dcmplx( vec1((ii-1)*nang+jj,3) , vec1((ii-1)*nang+jj,4) )
    r0(ii,jj,3)=dcmplx( vec1((ii-1)*nang+jj,5) , vec1((ii-1)*nang+jj,6) )
    r1(ii,jj,1)=dcmplx( vec2((ii-1)*nang+jj,1) , vec2((ii-1)*nang+jj,2) )
    r1(ii,jj,2)=dcmplx( vec2((ii-1)*nang+jj,3) , vec2((ii-1)*nang+jj,4) )
    r1(ii,jj,3)=dcmplx( vec2((ii-1)*nang+jj,5) , vec2((ii-1)*nang+jj,6) )
    r2(ii,jj,1)=dcmplx( vec3((ii-1)*nang+jj,1) , vec3((ii-1)*nang+jj,2) )
    r2(ii,jj,2)=dcmplx( vec3((ii-1)*nang+jj,3) , vec3((ii-1)*nang+jj,4) )
    r2(ii,jj,3)=dcmplx( vec3((ii-1)*nang+jj,5) , vec3((ii-1)*nang+jj,6) )
  end do
end do
open(newunit=file4,file='test_sym_dist3f.txt',status='old')
do ii=1,nang
  read(file4,*)theta(ii),phi(ii) !reading the angular distribution used to calculate the photoionization matrix elements 
end do

ip0 = pot1(i1+q1i0,i2+q2i0) - e1neut(i1+q1i0,i2+q2i0) !ionization potential for ground state of the cation
ip1 = pot2(i1+q1i0,i2+q2i0) - e1neut(i1+q1i0,i2+q2i0) !ionization potential for first excited state of the cation
ip2 = pot3(i1+q1i0,i2+q2i0) - e1neut(i1+q1i0,i2+q2i0) !ionization potential for second excited state of the cation
e0 = e_ip - ip0 !Energy of the ionized electron if the molecule goes for the cation ground state
e1 = e_ip - ip1 !Energy of the ionized electron if the molecule goes for the cation first excited state
e2 = e_ip - ip2 !Energy of the ionized electron if the molecule goes for the cation second excited state

domega=4*pi/nang
coef0=0.0_dp
coef1=0.0_dp
coef2=0.0_dp
do ii=1,nk
  do jj=1,nang
    coef0(ii,1) = coef0(ii,1) + r0(ii,jj,1) * domega !Doing the operation int( PICE(k) * domega )
    coef0(ii,2) = coef0(ii,2) + r0(ii,jj,2) * domega !Doing the operation int( PICE(k) * domega )
    coef0(ii,3) = coef0(ii,3) + r0(ii,jj,3) * domega !Doing the operation int( PICE(k) * domega )
    coef1(ii,1) = coef1(ii,1) + r1(ii,jj,1) * domega !Doing the operation int( PICE(k) * domega )
    coef1(ii,2) = coef1(ii,2) + r1(ii,jj,2) * domega !Doing the operation int( PICE(k) * domega )
    coef1(ii,3) = coef1(ii,3) + r1(ii,jj,3) * domega !Doing the operation int( PICE(k) * domega )
    coef2(ii,1) = coef2(ii,1) + r2(ii,jj,1) * domega !Doing the operation int( PICE(k) * domega )
    coef2(ii,2) = coef2(ii,2) + r2(ii,jj,2) * domega !Doing the operation int( PICE(k) * domega )
    coef2(ii,3) = coef2(ii,3) + r2(ii,jj,3) * domega !Doing the operation int( PICE(k) * domega )
  end do
  p_e(ii) = 0.005859375_dp + (ii-1) * 0.00588235294117647_dp !Momentum values in which the photoionization matrix elements are spanned
end do
!define the correct i, the correct momentum of the electron
aux2 = e_ip
do ii=1,nk
  if ( dsqrt( (p_e(ii)**2.0_dp/2.0_dp - (e_ip - ip0))**2.0_dp) < aux2(1) ) then
    aux1(1) = ii !Defining the value of the momentum of the leaving electron
    aux2(1) = dsqrt( (p_e(ii)**2.0_dp/2.0_dp - (e_ip - ip0))**2.0_dp)
    aux3(1) = p_e(ii)
  end if
  if ( dsqrt( (p_e(ii)**2.0_dp/2.0_dp - (e_ip - ip1))**2.0_dp) < aux2(2) ) then
    aux1(2) = ii !Defining the value of the momentum of the leaving electron
    aux2(2) = dsqrt( (p_e(ii)**2.0_dp/2.0_dp - (e_ip - ip1))**2.0_dp)
    aux3(2) = p_e(ii)
  end if
  if ( dsqrt( (p_e(ii)**2.0_dp/2.0_dp - (e_ip - ip2))**2.0_dp) < aux2(3) ) then
    aux1(3) = ii !Defining the value of the momentum of the leaving electron
    aux2(3) = dsqrt( (p_e(ii)**2.0_dp/2.0_dp - (e_ip - ip2))**2.0_dp)
    aux3(3) = p_e(ii)
  end if
end do

!It rests to do the operation -e * sqrt(2) * E * int(PICE(k) * domega)   --- 'e' is the electron charge, that in atomic units is 1
pic(1,:) = coef0(aux1(1),:)
pic(2,:) = coef1(aux1(2),:)
pic(3,:) = coef2(aux1(3),:)
pic(1,:) = pic(1,:)*dsqrt(aux3(1))!*mask1(i1,i2) !Here I include the multiplication by the electron density (that is just its momentum) and the mask for the sign correction of the wavefunction 
pic(2,:) = pic(2,:)*dsqrt(aux3(2))!*mask2(i1,i2) !Here I include the multiplication by the electron density (that is just its momentum) and the mask for the sign correction of the wavefunction 
pic(3,:) = pic(3,:)*dsqrt(aux3(3))!*mask3(i1,i2) !Here I include the multiplication by the electron density (that is just its momentum) and the mask for the sign correction of the wavefunction 
close(unit=file1)
close(unit=file2)
close(unit=file3)
close(unit=file4)
end subroutine p_i_a
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine generate_random_orientation
implicit none
integer :: n,i3
real(kind=dp),allocatable :: temp_m(:,:)
real(kind=dp) :: tmp(2)
n = nsamples/8
allocate(temp_m(n,3))
i3 = 0 
do while (i3 < n)
  call random_number(tmp)
  tmp = tmp*2-1
  if (tmp(1)**2.0_dp + tmp(2)**2.0_dp < 1.0_dp ) then
    i3 = i3 + 1 
    temp_m(i3,1) = 2.0_dp * tmp(1) * dsqrt(1.0_dp - tmp(1)**2.0_dp - tmp(2)**2.0_dp)
    temp_m(i3,2) = 2.0_dp * tmp(2) * dsqrt(1.0_dp - tmp(1)**2.0_dp - tmp(2)**2.0_dp)
    temp_m(i3,3) = 1.0_dp - 2.0_dp * (tmp(1)**2.0_dp + tmp(2)**2.0_dp)
  end if
end do
!call random_number(temp_m)
do i = 1,n 
  temp_m(i,:) = temp_m(i,:) / norm2(temp_m(i,:)) !* sqrt(3.0_dp)
end do
orie(1:n,:) = temp_m(:,:)
do i = 1,n 
  orie(i+n,:) = rotz( temp_m(i,:), 90.0_dp*pi/180.0_dp )
end do
do i = 1,2*n
  orie(i+n*2,:) = rotz( orie(i,:),180.0_dp*pi/180.0_dp )
end do
do i = 1,4*n
  orie(i+n*4,:) = - orie(i,:)
end do

end subroutine generate_random_orientation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotz(vec,ang)
!this subroutine takes a vector components x,y and z and rotates by an angle ang around the Z axis, returning the rotated vector
implicit none
real(kind=dp) :: vec(3),ang,matrix(3,3),rotz(3)

matrix(1,1) = cos(ang)
matrix(2,1) =-sin(ang)
matrix(3,1) = 0.0_dp
matrix(1,2) = sin(ang)
matrix(2,2) = cos(ang)
matrix(3,2) = 0.0_dp
matrix(1,3) = 0.0_dp
matrix(2,3) = 0.0_dp
matrix(3,3) = 1.0_dp

rotz(:) = matmul(matrix,vec)

return
end function rotz
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine save_matrix_h5(x,m,n,fname,le)
use HDF5
implicit none
integer             :: m,n,le
real(kind=dp)       :: x(m,n)
character(len=le)   :: fname         ! File name
character(len=le-3) :: dsetname      ! dataset name
integer(HID_T)      :: file_id       ! File identifier
integer(HID_T)      :: dspace_id     ! Dataspace identifier
integer(HID_T)      :: dset_id       ! Dataset identifier
integer(HSIZE_T)    :: dims(2)       ! Dimensions for Dataset and Dataspace
integer,parameter   :: rank = 2      ! Dataset rank = number of dimensions
integer*4           :: error         ! Error flag
dims(1)=m
dims(2)=n
write(dsetname,'(<le-3>a)') fname(1:le-3)
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
subroutine save_vector_h5(x,n,fname,le)
use HDF5
implicit none
integer             :: n,le
real(kind=dp)       :: x(n)
character(len=le)   :: fname         ! File name
character(len=le-3) :: dsetname      ! dataset name
integer(HID_T)      :: file_id       ! File identifier
integer(HID_T)      :: dspace_id     ! Dataspace identifier
integer(HID_T)      :: dset_id       ! Dataset identifier
integer(HSIZE_T)    :: dims(1)       ! Dimensions for Dataset and Dataspace
integer,parameter   :: rank = 1      ! Dataset rank = number of dimensions
integer*4           :: error         ! Error flag
dims=n
write(dsetname,'(<le-3>a)') fname(1:le-3)
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
subroutine open_vector_real_double_h5(x,n,fname,le)
use HDF5
implicit none
real(kind=dp)             :: x(n)
integer,parameter         :: rank = 1           ! Dataset rank = number of dimensions
integer                   :: n,le
character(len=le)         :: fname              ! File name
character(len=le-3)       :: dsetname           ! dataset name
integer(HID_T)            :: file_id            ! File identifier
integer(HID_T)            :: dspace_id          ! Dataspace identifier
integer(HID_T)            :: dset_id            ! Dataset identifier
integer(HSIZE_T)          :: dims(rank),maxdims(rank)       ! Dimensions for Dataset and Dataspace
integer*4                 :: error              ! Error flag
write(dsetname,'(<le-3>a)') fname(1:le-3) ! open a h5 file with a dataset that is the file name without the .h5
! Initialize FORTRAN interface.
call h5open_f(error)
! Open the file.
CALL h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
! Open an existing dataset.
CALL h5dopen_f(file_id, dsetname, dset_id, error)
! Get dataspace ID
CALL h5dget_space_f(dset_id, dspace_id, error)
! Get Dataspace dimensions
CALL h5sget_simple_extent_dims_f(dspace_id,dims, maxdims, error)
!n = dims(rank)
!allocate(x(n))
! Read the dataspace.
CALL h5dread_f(dset_id, h5kind_to_type(dp,H5_REAL_KIND), x, dims,  error)
! Close the file.
call h5fclose_f(file_id, error)
! Close FORTRAN interface.
call h5close_f(error)
end subroutine open_vector_real_double_h5
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine open_matrix_real_double_h5(x,m,n,fname,le)
use HDF5
implicit none
real(kind=dp)             :: x(m,n)
integer                   :: m,n,le
character(len=le)         :: fname              ! File name
character(len=le-3)       :: dsetname           ! dataset name
integer(HID_T)            :: file_id            ! File identifier
integer(HID_T)            :: dspace_id          ! Dataspace identifier
integer(HID_T)            :: dset_id            ! Dataset identifier
integer,parameter         :: rank = 2           ! Dataset rank = number of dimensions
integer(HSIZE_T)          :: dims(rank),maxdims(rank) ! Dimensions for Dataset and Dataspace
integer*4                 :: error              ! Error flag
dims(1)=m
dims(2)=n
write(dsetname,'(<le-3>a)') fname(1:le-3) ! open a h5 file with a dataset that is the file name without the .h5
! Initialize FORTRAN interface.
call h5open_f(error)
! Open the file.
CALL h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
! Open an existing dataset.
CALL h5dopen_f(file_id, dsetname, dset_id, error)
! Get dataspace ID
CALL h5dget_space_f(dset_id, dspace_id, error)
! Get Dataspace dimensions
CALL h5sget_simple_extent_dims_f(dspace_id,dims, maxdims, error)
!m = dims(1)
!n = dims(2)
!allocate(x(m,n))
! Read the dataspace.
CALL h5dread_f(dset_id, h5kind_to_type(dp,H5_REAL_KIND), x, dims,  error)
! Close the file.
call h5fclose_f(file_id, error)
! Close FORTRAN interface.
call h5close_f(error)
end subroutine open_matrix_real_double_h5
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine progress_bar(jjjj,ffff)
implicit none
integer        ::jjjj,kkkk,ffff
character(len=30)::bar="?????% |                    | "
write(bar(1:5),'(f5.1)') 100.0_dp/real(ffff)*jjjj
do kkkk=1, int(real(jjjj)/real(ffff)*20.0_dp)
  bar(8+kkkk:8+kkkk)="*"
enddo
! print the progress bar.
write(6,'(a1,a30)',advance="no") char(13), bar 
if (jjjj/=ffff) then
  flush(6)
else
  write(6,*)''
  write(*,*)''
endif
return
end subroutine progress_bar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end module global_param
