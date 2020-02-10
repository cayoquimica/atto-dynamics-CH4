module global_param
implicit none
integer, parameter                         :: dp = kind(1.d0) 
integer                                    :: i,j,k,l,ll,k_Ha,k_dip,k_moq1,k_moq2,k_am,ij
integer, parameter                         :: NA=5 !Number of atoms
integer, parameter                         :: Nst=3 !number of electronic states
integer, parameter                         :: Nq1=146 !126 !number of points of the grid along q1
integer, parameter                         :: Nq2=184 !179 !number of points of the grid along q2
integer, parameter                         :: s=Nq1*Nq2 !counter to be used as index - means the size of a 1 state matrix, which is Nq1*Nq2
integer,parameter                          :: q1_initial = 25
integer,parameter                          :: q1_final = 70
integer,parameter                          :: q2_initial = 75
integer,parameter                          :: q2_final = 115
integer,parameter                          :: nsamples = 25*8!25*4 !Number of randon orientations to sample for the photoionization around each bond (the *4 is for each bond)
integer,parameter                          :: nfiles = 2000 !Number of snapshots in time to save
integer,parameter                          :: npoints = 100000 !Number of time steps to take in the simulation
complex(kind=dp), parameter                :: im=dcmplx(0.d0,1.d0) !imaginary unity
real(kind=dp),parameter                    :: t0 = 0.d0 !Initial time
real(kind=dp),parameter                    :: tf = 2000.d0 !Final time
real(kind=dp),parameter                    :: tstep = (tf-t0)/npoints !Time step
real(kind=dp),parameter                    :: pi = 3.141592653589793d0
real(kind=dp),parameter                    :: coneAng = 10.d0 * pi / 180.d0 !Angle around each bond to sample ramdonly for initial photoionizations
real(kind=dp),parameter                    :: sq1=0.08d0 !step in q1 in atomic units
real(kind=dp),parameter                    :: sq2=0.07d0 !step in q1 in atomic units
real(kind=dp),parameter                    :: saw2au=1822.889950851334d0 !mass conversion factor from Standard atomic weight to atomic units
real(kind=dp),parameter                    :: car=12.011d0*saw2au !Carbon mass in atomic units
real(kind=dp),parameter                    :: hi=1.00794d0*saw2au !Hidrogen mass in atomic units
real(kind=dp),parameter                    :: mtotal = car + 4.d0 * hi ! total mass of CH4
real(kind=dp),parameter                    :: t00 = 800.d0 !time where the pulse is centered
real(kind=dp),parameter                    :: phase = 0.d0 !phase factor for the pulse related to the gaussian envelope
real(kind=dp),parameter                    :: freq =0.056937d0 !0.512285550500502 is the ionizating pulse !0.056937d0 !frequency of the pulse, in a.u. - 800 nm of wavelength
real(kind=dp),parameter                    :: sig = 100.d0  ! 50 approx 1200 attoseconds - width of the gaussian envelop
real(kind=dp),parameter                    :: E00 = 0.00d0 !0.05d0 !Electric field intensity
real(kind=dp),parameter                    :: e_ip =  1.d0 !Energy of the ionizing pulse in atomic units 1 = approx 27 ev
real(kind=dp),dimension(3), parameter      :: ori=[ 1.d0, 1.d0, 1.d0] !Orientation of the electric field of the probing pulse
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
complex(kind=dp),allocatable               :: wfout(:,:),coh(:,:),cohe(:,:)
real(kind=dp),allocatable                  :: Ha_val(:),dip_val(:,:),am_val(:),moq1_val(:),moq2_val(:) !CSR vectors for sparse matrix multiplication
integer,allocatable                        :: Ha_rowc(:),Ha_row_col(:,:) !CSR vectors for sparse matrix multiplication
integer,allocatable                        :: moq1_rowc(:),moq1_row_col(:,:),moq2_rowc(:),moq2_row_col(:,:) !CSR vectors for sparse matrix multiplication
integer,allocatable                        :: am_rowc(:),am_row_col(:,:) !CSR vectors for sparse matrix multiplication
real(kind=dp)                              :: mass1,mass2,mass3 !Reduced masses to be used in the second derivative
real(kind=dp)                              :: Et !Final electrical field of the pulse
real(kind=dp),dimension(3)                 :: orientation,u !Orientation of the electric field of the ionizing pulse
real(kind=dp),allocatable                  :: orie(:,:) !Vector with the random orientations of the electric field of the ionizing pulse
real(kind=dp)                              :: ind1,ind2,ind3,const1,const2,const3,E_init
complex(kind=dp),allocatable               :: pia(:,:),nwf0(:),pice(:),wf0(:)
complex(kind=dp)                           :: phote( (q1_final-q1_initial+1)*(q2_final-q2_initial+1), Nst, 3 )
real(kind=dp),dimension(0:nfiles)          :: e1,e2,e3,L1,L2,L3,sum1,sum2,sum3,pq1_1,pq2_1,pq1_2,pq2_2,pq1_3,pq2_3,Te1,Te2,Te3
real(kind=dp),dimension(0:nfiles)          :: fa,fb,fc,fd,fe,ff,fg,fh,fi,fj,fk,fl,fm,fn,fo 
real(kind=dp)                              :: momq1t(nfiles),momq2t(nfiles),maxmomq1(nsamples),maxmomq2(nsamples) ! This is for saving norm through time
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
      co1(ll)=(i-Nq1/2.d0-1.d0)*sq1+sq1/2.d0!(Nq1-Nq11/2.d0-i)*sq1+sq1/2.d0
      co2(ll)=(j-(114.d0)-1.d0)*sq2
      ll=ll+1
    end do
  end do
end do

!-------------------------------------------------------------------!
! Loading electronic structure data: Energy and dipole moments      !
open(unit=99,file='v1neutral.txt',status='old')                     !
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
q1(1)=0.005468383941194d0;q1(2)=-0.005539667067822d0;q1(3)=0.099940066929625d0;q1(4)=0.352088100072972d0;q1(5)=0.315861943182830d0  !
q1(6)=-0.495856225504668d0;q1(7)=-0.382798634963684d0;q1(8)=-0.284702415219885d0;q1(9)=-0.495881834018862d0                         !
q1(10)=-0.056657307983618d0;q1(11)=0.123960758550165d0;q1(12)=-0.097711469008833d0;q1(13)=0.022204480454270d0                       !
q1(14)=-0.089107486989717d0;q1(15)=-0.101474677166113d0                                                                             !
                                                                                                                                    !
q2(1)=0.005503122218018d0;q2(2)=-0.005574858175503d0;q2(3)=0.000798651732665d0;q2(4)=-0.173145514483078d0;q2(5)=-0.209601800471665d0!
q2(6)=0.414670087187852d0;q2(7)=0.142239888853687d0;q2(8)=0.240959271431981d0;q2(9)=0.414644315993849d0;q2(10)=-0.174783204984944d0 !
q2(11)=0.242514205248471d0;q2(12)=-0.417522164824840d0;q2(13)=0.140111513551201d0;q2(14)=-0.207439525935984d0                       !
q2(15)=-0.421309279015076d0                                                                                                         !
!Defining the coefficients of the inverse transformation - from internal to cartesian                                               !
!This is just the inverse of the above (s^-1):                                                                                      !
!                                                                                                                                   !
!      | q1(1)  q1(2)  q1(3) ...  q1(i) |                                                                                           !
!  s = |                                |                                                                                           !
!      | q2(1)  q2(2)  q2(3) ...  q2(i) |                                                                                           !
!                                                                                                                                   !
q1i(1)=0.0d0;q1i(2)=0.0d0;q1i(3)=0.1404666820458888d0;q1i(4)=0.327739015711321d0;q1i(5)=0.327739015711321d0                         !
q1i(6)=-0.382350835749710d0;q1i(7)=-0.327739015711321d0;q1i(8)=-0.327739015711321d0;q1i(9)=-0.382350835749710d0                     !
q1i(10)=-0.249049347411151d0;q1i(11)=0.249049347411151d0;q1i(12)=-0.454576619283899d0;q1i(13)=0.249049347411151d0                   !
q1i(14)=-0.249049347411151d0;q1i(15)=-0.454576619283899d0                                                                           !
                                                                                                                                    !
q2i(1)=0.0d0;q2i(2)=0.0d0;q2i(3)=0.075292919338856d0;q2i(4)=-0.016803984583966d0;q2i(5)=-0.016803984583966d0                        !
q2i(6)=0.214461491870736d0;q2i(7)=0.016803984583966d0;q2i(8)=0.016803984583966d0;q2i(9)=0.214461491870736d0                         !
q2i(10)=-0.325973953345556d0;q2i(11)=0.325973953345556d0;q2i(12)=-0.663071158209507d0;q2i(13)=0.325973953345556d0                   !
q2i(14)=-0.325973953345556d0;q2i(15)=-0.663071158209507d0                                                                           !
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
mass1=0.d0;mass2=0.d0;mass3=0.d0                                              !
ai=0.0d0; bi=0.0d0; aii=0.d0; bii=0.d0                                        !
do i=1,3*NA                                                                   !
  mass1=mass1 + q1(i)**2.d0 / mass(i)                                         !
  mass2=mass2 + q2(i)**2.d0 / mass(i)                                         !
  mass3=mass3 + q1(i)*q2(i)*2.d0 / mass(i)                                    !
  ai=ai+q1(i) !NAC values are already divided by mass when loaded             !
  bi=bi+q2(i) !NAC values are already divided by mass when loaded             !
  aii=aii+q1i(i) !NAC values are already divided by mass when loaded          !
  bii=bii+q2i(i) !NAC values are already divided by mass when loaded          !
end do                                                                        !
!-----------------------------------------------------------------------------!

end subroutine load_data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine generate_initial_wf
real(kind=dp) :: soma,n

n = Nq1 * Nq2 * Nst

!Now read the photoinization coefficients and project them into the initial neutral ground state wave function
pice=dcmplx(0.d0,0.d0)
k = (q2_initial-1+20) * Nq1 + (q1_initial+27) - 1 !The -1 is because the vector starts at index 0
ij = 0
!!$OMP parallel do default(private) shared(pice)
do j=q2_initial+20,q2_final+20  !The photoionization coefficiets were calculated only for q2=75+20:115+20 and q1=25+27:70+27
  do i=q1_initial+27,q1_final+27 !where the amplitudes of the eigen state of the neutral is non-zero (< 10^-8). This is the Frank-Condon region
    ij = ij +1 
    pice(k)     = - dsqrt(2.d0) * 0.05d0 * dot_product( orientation , phote(ij,1,:) )
    pice(k+s)   = - dsqrt(2.d0) * 0.05d0 * dot_product( orientation , phote(ij,2,:) )
    pice(k+2*s) = - dsqrt(2.d0) * 0.05d0 * dot_product( orientation , phote(ij,3,:) )
    k = k+1 
  end do
  k = k + (Nq1-(q1_final+27)) + (q1_initial+27) - 1 !Add the zeros values from 70+27 until Nq1 and from 1 to 25+27. The -1 here is to anulate the last -1 of the previous loop
end do
!!$OMP end parallel do
!Projecting the photoionization coeficients into the neutral eigen state
do i=0,s-1
  wf0(i)     = nwf0(i) * pice(i)
  wf0(i+s)   = nwf0(i) * pice(i+s)
  wf0(i+2*s) = nwf0(i) * pice(i+2*s)
end do
!--------------------------------------------!
!normalizing                                 !
soma=0.d0                                    !   
do i=0,n-1                                   !   
  soma=soma+dconjg(wf0(i))*wf0(i)            !   
end do                                       !   
wf0=wf0/sqrt(soma)                           !   
!--------------------------------------------!
end subroutine generate_initial_wf
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Subroutine that calculates the photoionization amplitudes for a given geometry and a given electric field
!Some parts of this function assumes that the number of electronic states is 3
subroutine p_i_a(i1,i2,pic)
implicit none
complex(kind=dp),dimension(Nst,3) :: pic
integer                           :: i1,i2,ii,jj,kk
character(len=21)                 :: fname00
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
real(kind=dp)                     :: aux2(Nst)

allocate(vec1(nang*nk,Nst*2),vec2(nang*nk,Nst*2),vec3(nang*nk,Nst*2))
allocate(r0(nk,nang,Nst),r1(nk,nang,Nst),r2(nk,nang,Nst),coef0(nk,Nst),coef1(nk,Nst),coef2(nk,Nst))
call getcwd(fpath0) !getting the working directory path
fname00='pice_10000000_0_0.dat'
write(fname00(7:9),'(i0.3)') i2+100
write(fname00(12:13),'(i0.2)') i1
fpath1="~/pice_files/"//fname00
open(newunit=file1,file=fpath1,status='old')
write(fname00(17:17),'(i1)') 1
fpath1="~/pice_files/"//fname00
open(newunit=file2,file=fpath1,status='old')
write(fname00(17:17),'(i1)') 2
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
open(newunit=file4,file='test_sym_dist3.txt',status='old')
do ii=1,nang
  read(file4,*)theta(ii),phi(ii) !reading the angular distribution used to calculate the photoionization matrix elements 
end do

ip0 = pot1(i1+27,i2+20) - e1neut(i1+27,i2+20) !ionization potential for ground state of the cation
ip1 = pot2(i1+27,i2+20) - e1neut(i1+27,i2+20) !ionization potential for first excited state of the cation
ip2 = pot3(i1+27,i2+20) - e1neut(i1+27,i2+20) !ionization potential for second excited state of the cation
e0 = e_ip - ip0 !Energy of the ionized electron if the molecule goes for the cation ground state
e1 = e_ip - ip1 !Energy of the ionized electron if the molecule goes for the cation first excited state
e2 = e_ip - ip2 !Energy of the ionized electron if the molecule goes for the cation second excited state

domega=4*pi/nang
coef0=0.d0
coef1=0.d0
coef2=0.d0
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
  p_e(ii) = 0.005859375d0 + (ii-1) * 0.00588235294117647d0 !Momentum values in which the photoionization matrix elements are spanned
end do

!define the correct i, the correct momentum of the electron
aux2 = e_ip
do ii=1,nk
  if ( dsqrt( (p_e(ii) - (e_ip - ip0))**2.d0) < aux2(1) ) then
    aux1(1) = ii !Defining the value of the momentum of the leaving electron
    aux2(1) = dsqrt( (p_e(ii) - (e_ip - ip0))**2.d0)
  end if
  if ( dsqrt( (p_e(ii) - (e_ip - ip1))**2.d0) < aux2(2) ) then
    aux1(2) = ii !Defining the value of the momentum of the leaving electron
    aux2(2) = dsqrt( (p_e(ii) - (e_ip - ip1))**2.d0)
  end if
  if ( dsqrt( (p_e(ii) - (e_ip - ip2))**2.d0) < aux2(3) ) then
    aux1(3) = ii !Defining the value of the momentum of the leaving electron
    aux2(3) = dsqrt( (p_e(ii) - (e_ip - ip2))**2.d0)
  end if
end do

!Do the operation -e * sqrt(2) * E * int(PICE(k) * domega)   --- 'e' is the electron charge, that in atomic units is 1
pic(1,:) = coef0(aux1(1),:)
pic(2,:) = coef1(aux1(2),:)
pic(3,:) = coef2(aux1(3),:)
close(unit=file1)
close(unit=file2)
close(unit=file3)
close(unit=file4)
end subroutine p_i_a
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!subroutine generate_random_orientation
!implicit none
!integer :: n
!real(kind=dp),allocatable :: temp_m(:,:)
!n = nsamples/8
!allocate(temp_m(n,3))
!call random_number(temp_m)
!do i = 1,n
!  temp_m(i,:) = temp_m(i,:) / norm2(temp_m(i,:)) * sqrt(3.d0)
!end do
!orie(1:n,:) = temp_m(:,:)
!do i = 1,n
!  orie(i+n,:) = rotz( temp_m(i,:), 90.d0*pi/180.d0 )
!end do
!do i = 1,2*n
!  orie(i+n*2,:) = rotz( orie(i,:),180.d0*pi/180.d0 )
!end do
!do i = 1,4*n
!  orie(i+n*4,:) = - orie(i,:)
!end do
!
!end subroutine generate_random_orientation
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotz(vec,ang)
!this subroutine takes a vector components x,y and z and rotates by an angle ang around the Z axis, returning the rotated vector
implicit none
real(kind=dp) :: vec(3),ang,matrix(3,3),rotz(3)

matrix(1,1) = cos(ang)
matrix(2,1) =-sin(ang)
matrix(3,1) = 0.d0
matrix(1,2) = sin(ang)
matrix(2,2) = cos(ang)
matrix(3,2) = 0.d0
matrix(1,3) = 0.d0
matrix(2,3) = 0.d0
matrix(3,3) = 1.d0

rotz(:) = matmul(matrix,vec)

return
end function rotz
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine generate_random_orientation
real(kind=dp) :: x,w,z,ang,u(3),rot,matrix(3,3),matrix2(3,3),eye(3,3),rotmatrix(3,3),normfactor
call random_number(z)
z = z * (1.d0 - dcos(coneAng)) + dcos(coneAng)
call random_number(ang)
ang = ang * 2.d0 * pi
x = dsqrt(1.d0 - z**2.d0) * dcos(ang)
w = dsqrt(1.d0 - z**2.d0) * dsin(ang) !The y component
!rotating x, w and z to be around each bond
u = [-1.d0 * orie(1,2), 1.d0 * orie(1,1), 0.d0 ]/norm2(orie(1,:)) !Cross product between [0,0,1] and orie
!normfactor = dsqrt( u(1)**2.d0 + u(2)**2.d0 + u(3)**2.d0 )
u = u / norm2(u)
if (u(1) /= u(1) .and. u(2) /= u(2) .and. u(3) /= u(3)) then
  u = [0.d0, 0.d0, 0.d0]
end if
rot = acos( dot_product(orie(1,:)/norm2(orie(1,:)),[0,0,1]) )
matrix(1,1) = 0.d0
matrix(1,2) =-u(3)
matrix(1,3) = u(2)
matrix(2,1) = u(3)
matrix(2,2) = 0.d0
matrix(2,3) =-u(1)
matrix(3,1) =-u(2)
matrix(3,2) = u(1)
matrix(3,3) = 0.d0

matrix2(1,1) = u(1)*u(1)
matrix2(1,2) = u(1)*u(2)
matrix2(1,3) = u(1)*u(3)
matrix2(2,1) = u(2)*u(1)
matrix2(2,2) = u(2)*u(2)
matrix2(2,3) = u(2)*u(3)
matrix2(3,1) = u(3)*u(1)
matrix2(3,2) = u(3)*u(2)
matrix2(3,3) = u(3)*u(3)

eye(1,1) = 1.d0
eye(1,2) = 0.d0
eye(1,3) = 0.d0
eye(2,1) = 0.d0
eye(2,2) = 1.d0
eye(2,3) = 0.d0
eye(3,1) = 0.d0
eye(3,2) = 0.d0
eye(3,3) = 1.d0

rotmatrix(:,:) = dcos(rot) * eye(:,:) + dsin(rot) * matrix + ( 1.d0 - dcos(rot) ) * matrix2

orientation(1) = rotmatrix(1,1)*x + rotmatrix(1,2)*w + rotmatrix(1,3)*z
orientation(2) = rotmatrix(2,1)*x + rotmatrix(2,2)*w + rotmatrix(2,3)*z
orientation(3) = rotmatrix(3,1)*x + rotmatrix(3,2)*w + rotmatrix(3,3)*z

!write(*,'(3f7.4)')orie(1,:)
!write(*,'(3f7.4)')u
!write(*,'(3f7.4)')orientation
!read(*,*)

end subroutine generate_random_orientation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module global_param
