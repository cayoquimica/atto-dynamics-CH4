module global_param
implicit none
integer, parameter                         :: dp = kind(1.d0) 
integer                                    :: i,j,k,l,ll,k_Ha,k_dip,k_moq1,k_moq2,k_am
integer, parameter                         :: NA=5 !Number of atoms
integer, parameter                         :: Nst=3 !number of electronic states
integer, parameter                         :: Nq1=146 !126 !number of points of the grid along q1
integer, parameter                         :: Nq2=184 !179 !number of points of the grid along q2
integer, parameter                         :: s=Nq1*Nq2 !counter to be used as index - means the size of a 1 state matrix, which is Nq1*Nq2
complex(kind=dp), parameter                :: im=dcmplx(0.d0,1.d0) !imaginary unity
real(kind=dp), parameter                   :: pi = 3.141592653589793d0
real(kind=dp), parameter                   :: sq1=0.08d0 !step in q1 in atomic units
real(kind=dp), parameter                   :: sq2=0.07d0 !step in q1 in atomic units
real(kind=dp), parameter                   :: saw2au=1822.889950851334d0 !mass conversion factor from Standard atomic weight to atomic units
real(kind=dp), parameter                   :: car=12.011d0*saw2au !Carbon mass in atomic units
real(kind=dp), parameter                   :: hi=1.00794d0*saw2au !Hidrogen mass in atomic units
real(kind=dp), parameter                   :: mtotal = car + 4.d0 * hi ! total mass of CH4
real(kind=dp)                              :: mass(3*NA) !3N Vector with the masses of the atoms
real(kind=dp)                              :: q1(3*NA),q2(3*NA),q1i(3*NA),q2i(3*NA),ai,bi,aii,bii !Conversion vector from internal coordinates q1 and q2 to cartesian coordinates
real(kind=dp)                              :: co1(0:Nst*Nq1*Nq2-1),co2(0:Nst*Nq1*Nq2-1)
real(kind=dp), dimension(:,:), allocatable :: pot1,pot2,pot3 !Matrices with each state potential energy
real(kind=dp), dimension(:,:), allocatable :: pdm1x,pdm2x,pdm3x,pdm1y,pdm2y,pdm3y ! Permanent dipole moments
real(kind=dp), dimension(:,:), allocatable :: tdm21x,tdm31x,tdm32x,tdm21y,tdm31y,tdm32y ! Transition dipole moments
real(kind=dp), dimension(:,:), allocatable :: Ha! Hamiltonian matrix
real(kind=dp), dimension(:,:), allocatable :: ham,ax ! Variable to store the Hamiltonian matrix temporary
real(kind=dp), dimension(:,:), allocatable :: moq1,moq2! Momentum matrix
real(kind=dp), allocatable                 :: Ha_val(:),dip_val(:),am_val(:),moq1_val(:),moq2_val(:) !CSR vectors for sparse matrix multiplication
integer, allocatable                       :: Ha_rowc(:),Ha_row_col(:,:) !CSR vectors for sparse matrix multiplication
integer, allocatable                       :: moq1_rowc(:),moq1_row_col(:,:),moq2_rowc(:),moq2_row_col(:,:) !CSR vectors for sparse matrix multiplication
integer, allocatable                       :: am_rowc(:),am_row_col(:,:) !CSR vectors for sparse matrix multiplication
real(kind=dp)                              :: mass1,mass2,mass3 !Reduced masses to be used in the second derivative
real(kind=dp)                              :: Et !Final electrical field of the pulse
real(kind=dp), parameter                   :: t00 = 0.d0 !time where the pulse is centered
real(kind=dp), parameter                   :: phase = 0.d0 !phase factor for the pulse related to the gaussian envelope
real(kind=dp), parameter                   :: freq =0.056937d0 !0.512285550500502 is the ionizating pulse !0.056937d0 !frequency of the pulse, in a.u. - 800 nm of wavelength
real(kind=dp), parameter                   :: sig = 125.d0  ! 50 approx 1200 attoseconds - width of the gaussian envelop
real(kind=dp), parameter                   :: E00 = 0.05d0 !0.05d0 !Electric field intensity
real(kind=dp)                              :: ind1,ind2,ind3,const1,const2,const3,E_init

contains

subroutine load_data
allocate (pot1(Nq1,Nq2))
allocate (pot2(Nq1,Nq2))
allocate (pot3(Nq1,Nq2))
!allocate ( pdm1x(Nq1,Nq2),pdm2x(Nq1,Nq2),pdm3x(Nq1,Nq2),pdm1y(Nq1,Nq2),pdm2y(Nq1,Nq2),pdm3y(Nq1,Nq2) )
allocate (pdm1x(Nq1,Nq2))
allocate (pdm2x(Nq1,Nq2))
allocate (pdm3x(Nq1,Nq2))
!allocate ( tdm21x(Nq1,Nq2),tdm31x(Nq1,Nq2),tdm32x(Nq1,Nq2),tdm21y(Nq1,Nq2),tdm31y(Nq1,Nq2),tdm32y(Nq1,Nq2) )
allocate (tdm21x(Nq1,Nq2))
allocate (tdm31x(Nq1,Nq2))
allocate (tdm32x(Nq1,Nq2))

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
open(unit=1,file='v1f.txt',status='old')                             !
open(unit=2,file='v2f.txt',status='old')                             !
open(unit=3,file='v3f.txt',status='old')                             !
do i=1,Nq1                                                          !
  read(1,*) pot1(i,:)                                               !
  read(2,*) pot2(i,:)                                               !
  read(3,*) pot3(i,:)                                               !
end do                                                              !
 close(unit=1)                                                      !
 close(unit=2)                                                      !
 close(unit=3)                                                      !
                                                                    !
open(unit=1,file='pdm1f.txt',status='old')                          !
open(unit=2,file='pdm2f.txt',status='old')                          !
open(unit=3,file='pdm3f.txt',status='old')                          !
open(unit=4,file='tdm21f.txt',status='old')                         !
open(unit=5,file='tdm31f.txt',status='old')                         !
open(unit=6,file='tdm32f.txt',status='old')                         !
do i=1,Nq1                                                          !
  read(1,*) pdm1x(i,:)                                              !
  read(2,*) pdm2x(i,:)                                              !
  read(3,*) pdm3x(i,:)                                              !
  read(4,*) tdm21x(i,:)                                             !
  read(5,*) tdm31x(i,:)                                             !
  read(6,*) tdm32x(i,:)                                             !
end do                                                              !
 close(unit=1)                                                      !
 close(unit=2)                                                      !
 close(unit=3)                                                      !
 close(unit=4)                                                      !
 close(unit=5)                                                      !
 close(unit=6)                                                      !
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

end module global_param
