module global_param
implicit none
integer, parameter :: dp = kind(1.d0) 
integer i,j,k,m,k_Ha,k_dip,k_moq1,k_moq2,k_am
integer, parameter :: NA=5 !Number of atoms
integer, parameter :: st=3 !number of electronic states
integer, parameter :: Nq1=252 !number of points of the grid along q1
integer, parameter :: Nq2=180 !number of points of the grid along q2
complex(kind=dp), parameter :: im = dcmplx(0.d0,1.d0)
real(kind=dp), parameter :: pi = 3.141592653589793d0
real(kind=dp), parameter :: stepX=0.04d0 !step in q1 in atomic units
real(kind=dp), parameter :: stepY=0.0675d0 !step in q1 in atomic units
real(kind=dp), parameter :: saw2au=1822.889950851334d0 !mass conversion factor from Standard atomic weight to atomic units
real(kind=dp), parameter :: car=12.011d0*saw2au !Carbon mass in atomic units
real(kind=dp), parameter :: hi=1.00794d0*saw2au !Hidrogen mass in atomic units
real(kind=dp), parameter :: mtotal = car + 4.d0 * hi ! total mass of CH4
real(kind=dp) :: mass(3*NA),q1(3*NA),q2(3*NA),pot1(Nq1),pot2(Nq1),pot3(Nq1) !Matrices with each state potential energy
real(kind=dp) :: pdm1x(Nq1),pdm2x(Nq1),pdm3x(Nq1),pdm1y(Nq1),pdm2y(Nq1),pdm3y(Nq1) ! Permanent dipole moment for state 1
real(kind=dp) :: tdm21x(Nq1),tdm31x(Nq1),tdm32x(Nq1),tdm21y(Nq1),tdm31y(Nq1),tdm32y(Nq1) ! Permanent dipole moment for state 1
real(kind=dp) :: Ha(st*Nq1,st*Nq1),ind1,ind2,ind3,const1,mass1
real(kind=dp) :: ham(st*Nq1,st*Nq1),Et,ax(0:st*Nq1-1,0:st*Nq1-1)
real(kind=dp), parameter :: t00 = 400.d0 !time where the pulse is centered
real(kind=dp), parameter :: phase = 0.d0 !phase factor for the pulse related to the gaussian envelope
real(kind=dp), parameter :: freq = 0.056937d0 !frequency of the pulse, in a.u.
real(kind=dp), parameter :: sig = 50 ! 50 approx 1200 attoseconds
real(kind=dp), parameter :: E00 = 0.05d0 !Electric field intensity


real(kind=dp), allocatable                 :: Ha_val(:),Ha_rowc(:),Ha_row_col(:,:),dip_val(:) !CSR vectors for sparse matrix multiplication


contains
!-------------------------------------------------------
subroutine hamiltonian

open(unit=1,file='e1int0.txt',status='old')
open(unit=2,file='e2int.txt',status='old')
open(unit=3,file='e3int.txt',status='old')
do i=1,Nq1
  read(1,*) pot1(i)
  read(2,*) pot2(i)
  read(3,*) pot3(i)
end do
 close(unit=1)
 close(unit=2)
 close(unit=3)

open(unit=1,file='pdm1x.txt',status='old')
open(unit=2,file='pdm2x.txt',status='old')
open(unit=3,file='pdm3x.txt',status='old')
open(unit=4,file='tdm21x.txt',status='old')
open(unit=5,file='tdm31x.txt',status='old')
open(unit=6,file='tdm32x.txt',status='old')
do i=1,Nq1
  read(1,*) pdm1x(i)
  read(2,*) pdm2x(i)
  read(3,*) pdm3x(i)
  read(4,*) tdm21x(i)
  read(5,*) tdm31x(i)
  read(6,*) tdm32x(i)
end do
 close(unit=1)
 close(unit=2)
 close(unit=3)
 close(unit=4)
 close(unit=5)
 close(unit=6)

open(unit=1,file='pdm1y.txt',status='old')
open(unit=2,file='pdm2y.txt',status='old')
open(unit=3,file='pdm3y.txt',status='old')
open(unit=4,file='tdm21y.txt',status='old')
open(unit=5,file='tdm31y.txt',status='old')
open(unit=6,file='tdm32y.txt',status='old')
do i=1,Nq1
  read(1,*) pdm1y(i)
  read(2,*) pdm2y(i)
  read(3,*) pdm3y(i)
  read(4,*) tdm21y(i)
  read(5,*) tdm31y(i)
  read(6,*) tdm32y(i)
end do
 close(unit=1)
 close(unit=2)
 close(unit=3)
 close(unit=4)
 close(unit=5)
 close(unit=6)

q1(1)=0.0d0;q1(2)=0.0d0;q1(3)=-0.100108188344758d0;q1(4)=-0.336746283751144d0;q1(5)=-0.336746283751144d0;q1(6)=0.497306436419762d0
q1(7)=0.336746283751144d0;q1(8)=0.336746283751144d0;q1(9)=0.497306436419762d0;q1(10)=0.074320862204978d0
q1(11)=-0.074320862204978d0;q1(12)=0.099157366092731d0;q1(13)=-0.074320862204978d0;q1(14)=0.074320862204978d0
q1(15)=0.099157366092731d0

q2(1)=0.0d0;q2(2)=0.0d0;q2(3)=0.0d0;q2(4)=0.192478578122720d0;q2(5)=0.192478578122720d0;q2(6)=-0.419409100911881d0
q2(7)=-0.192478578122720d0;q2(8)=-0.192478578122720d0;q2(9)=-0.419409100911881d0;q2(10)=0.192478578122720d0
q2(11)=-0.192478578122720d0;q2(12)=0.419409100911881d0;q2(13)=-0.192478578122720d0;q2(14)=0.192478578122720d0
q2(15)=0.419409100911881d0

mass(1)=car;mass(2)=car;mass(3)=car
do i=4,3*NA
  mass(i)=hi
end do

mass1=0.d0
do i=1,3*NA
  mass1=mass1 + q1(i)**2.d0 / mass(i)
end do

const1 = 1.d0/(12.d0*stepX**2.d0)

Ha(:,:)=0.0d0
ind1=-30.d0*const1
ind2= 16.d0*const1
ind3=-1.d0*const1
Ha(1,1)=ind1; Ha(1,2)=ind2; Ha(1,3)=ind3
Ha(Nq1,Nq1)=ind1; Ha(Nq1,Nq1-1)=ind2; Ha(Nq1,Nq1-2)=ind3
Ha(2,1)=ind2; Ha(2,2)=ind1; Ha(2,3)=ind2; Ha(2,4)=ind3
Ha(Nq1-1,Nq1)=ind2; Ha(Nq1-1,Nq1-1)=ind1; Ha(Nq1-1,Nq1-2)=ind2; Ha(Nq1-1,Nq1-3)=ind3

Ha(Nq1+1,Nq1+1)=ind1; Ha(Nq1+1,Nq1+2)=ind2; Ha(Nq1+1,Nq1+3)=ind3
Ha(2*Nq1,2*Nq1)=ind1; Ha(2*Nq1,2*Nq1-1)=ind2; Ha(2*Nq1,2*Nq1-2)=ind3
Ha(Nq1+2,Nq1+1)=ind2; Ha(Nq1+2,Nq1+2)=ind1; Ha(Nq1+2,Nq1+3)=ind2; Ha(Nq1+2,Nq1+4)=ind3
Ha(2*Nq1-1,2*Nq1)=ind2; Ha(2*Nq1-1,2*Nq1-1)=ind1; Ha(2*Nq1-1,2*Nq1-2)=ind2; Ha(2*Nq1-1,2*Nq1-3)=ind3

Ha(2*Nq1+1,2*Nq1+1)=ind1; Ha(2*Nq1+1,2*Nq1+2)=ind2; Ha(2*Nq1+1,2*Nq1+3)=ind3
Ha(3*Nq1,3*Nq1)=ind1; Ha(3*Nq1,3*Nq1-1)=ind2; Ha(3*Nq1,3*Nq1-2)=ind3
Ha(2*Nq1+2,2*Nq1+1)=ind2; Ha(2*Nq1+2,2*Nq1+2)=ind1; Ha(2*Nq1+2,2*Nq1+3)=ind2; Ha(2*Nq1+2,2*Nq1+4)=ind3
Ha(3*Nq1-1,3*Nq1)=ind2; Ha(3*Nq1-1,3*Nq1-1)=ind1; Ha(3*Nq1-1,3*Nq1-2)=ind2; Ha(3*Nq1-1,3*Nq1-3)=ind3

do i=1+2,Nq1-2 !from the third line to the antepenulmate one
  Ha(i,i-2)=ind3
  Ha(i,i-1)=ind2
  Ha(i,  i)=ind1
  Ha(i,i+1)=ind2
  Ha(i,i+2)=ind3

  Ha(Nq1+i,Nq1+i-2)=ind3
  Ha(Nq1+i,Nq1+i-1)=ind2
  Ha(Nq1+i,  Nq1+i)=ind1
  Ha(Nq1+i,Nq1+i+1)=ind2
  Ha(Nq1+i,Nq1+i+2)=ind3

  Ha(2*Nq1+i,2*Nq1+i-2)=ind3
  Ha(2*Nq1+i,2*Nq1+i-1)=ind2
  Ha(2*Nq1+i,  2*Nq1+i)=ind1
  Ha(2*Nq1+i,2*Nq1+i+1)=ind2
  Ha(2*Nq1+i,2*Nq1+i+2)=ind3
end do

Ha(:,:)= - 1.d0/2.d0 * mass1 * Ha(:,:)

do i=1,Nq1
!  Ha(i,i)=Ha(i,i)+pot1(i)
!  Ha(Nq1+i,Nq1+i)=Ha(Nq1+i,Nq1+i)+pot2(i)
!  Ha(2*Nq1+i,2*Nq1+i)=Ha(2*Nq1+i,2*Nq1+i)+pot3(i)

!  Ha(i,Nq1+i)=-tdm21(i) * Et
!  Ha(Nq1+i,i)=-tdm21(i) * Et
!  Ha(i,2*Nq1+i)=-tdm31(i) * Et
!  Ha(2*Nq1+i,i)=-tdm31(i) * Et
!  Ha(Nq1+i,2*Nq1+i)=-tdm32(i) * Et
!  Ha(2*Nq1+i,Nq1+i)=-tdm32(i) * Et
end do

m=st*Nq1
do i=1,m
  do j=1,m
    ax(i-1,j-1)=ha(i,j)
  end do
end do

!Creating the CSR vectors --------------------------------------------------------------------------------------------------!
!This is for the the hamiltonian - dipoles will be in a separate vector                                                     !
k=0                                                                                                                         !
do i=0,m-1 !running through rows                                                                                            !
  do j=0,m-1 !running through columns                                                                                       !
    if (ax(i,j) /= 0.d0) then                                                                                              !
      k=k+1                                                                                                                 !
    end if                                                                                                                  !
  end do                                                                                                                    !
end do                                                                                                                      !
k_Ha=k                                                                                                                      !
k_dip=k                                                                                                                     !
write(*,*)'k = ',k
write(*,'(99f13.6)')ax(0,0),ax(0,1),ax(0,2),ax(0,3)
write(*,'(99f13.6)')ax(1,0),ax(1,1),ax(1,2),ax(1,3),ax(1,4)
write(*,'(99f13.6)')ax(2,0),ax(2,1),ax(2,2),ax(2,3),ax(2,4),ax(2,5)
write(*,'(99f13.6)')ax(3,0),ax(3,1),ax(3,2),ax(3,3),ax(3,4),ax(3,5),ax(3,6)
write(*,'(4f13.6)')ax(Nq1+0,Nq1+0),ax(Nq1+0,Nq1+1),ax(Nq1+0,Nq1+2),ax(Nq1+0,Nq1+3)
write(*,'(5f13.6)')ax(Nq1+1,Nq1+0),ax(Nq1+1,Nq1+1),ax(Nq1+1,Nq1+2),ax(Nq1+1,Nq1+3),ax(Nq1+1,Nq1+4)
write(*,'(6f13.6)')ax(Nq1+2,Nq1+0),ax(Nq1+2,Nq1+1),ax(Nq1+2,Nq1+2),ax(Nq1+2,Nq1+3),ax(Nq1+2,Nq1+4),ax(Nq1+2,Nq1+5)
write(*,'(7f13.6)')ax(Nq1+3,Nq1+0),ax(Nq1+3,Nq1+1),ax(Nq1+3,Nq1+2),ax(Nq1+3,Nq1+3),ax(Nq1+3,Nq1+4),ax(Nq1+3,Nq1+5),ax(Nq1+3,Nq1+6)
allocate(Ha_val(0:k-1),Ha_rowc(0:m),Ha_row_col(0:k-1,0:1) )                                                !
!                                                                                                                           !
Ha_rowc=0.d0                                                                                                                !
k=0                                                                                                                         !
do i=0,m-1 !running through rows                                                                                            !
  do j=0,m-1 !running through columns                                                                                       !
    if (ax(i,j) /= 0.d0) then                                                                                              !
      Ha_val(k)=ax(i,j)  !storing each non-zero element                                                                     !
      Ha_row_col(k,0)=i !storing the row index of each non-zero element                                                     !
      Ha_row_col(k,1)=j !storing the column index of each non-zero element                                                  !
      k=k+1                                                                                                                 !
!write(*,'(2i6,7f13.6)')i,j,ax(i,j),Ha_val(k-1),Ha_row_col(k-1,:)
!read(*,*)
    end if                                                                                                                  !
  end do                                                                                                                    !
  Ha_rowc(i+1)=k !storing the counting of non-zero elements in each row                                                     ! 
!write(*,'(3i6,f13.6)')i,j,k,Ha_rowc(i+1)
end do                                                                                                                      !
!---------------------------------------------------------------------------------------------------------------------------!

!do j=0,k-1
!write(*,*)j
!write(*,*)'Ha e col = ',Ha_val(j),Ha_row_col(j,1)

!end do

!do i=0,m
!write(*,*)'lines = ',Ha_rowc(i)
!end do
!read(*,*)

!write(*,'(<Nq1>(f6.3))') Ha(:,:)
!read(*,*)

end subroutine hamiltonian

end module global_param
