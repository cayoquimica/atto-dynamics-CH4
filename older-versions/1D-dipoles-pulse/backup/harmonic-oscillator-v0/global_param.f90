module global_param
implicit none 
integer i,j,k
integer, parameter :: NA=2 !Number of atoms
integer, parameter :: st=1 !number of electronic states
integer, parameter :: Nq1=61 !number of points of the grid along q1
integer, parameter :: Nq2=61 !number of points of the grid along q2
real(kind=8), parameter :: saw2au=1822.889950851334d0 !mass conversion factor from Standard atomic weight to atomic units
real(kind=8), parameter :: car=12.011d0*saw2au !Carbon mass in atomic units
real(kind=8), parameter :: nitro=14.0067d0*saw2au !Carbon mass in atomic units
real(kind=8), parameter :: mredu=nitro**2.0d0/(2.0d0*nitro) ! reduced mass of N2
real(kind=8), parameter :: hi=1.00794d0*saw2au !Hidrogen mass in atomic units
real(kind=8) :: mass(3*NA),q1(3*NA),q2(3*NA)
real(kind=8) :: pot1(Nq1) !Matrices with each state potential energy

contains
!-------------------------------------------------------
subroutine mass_vec
implicit none
mass(1)=car;mass(2)=car;mass(3)=car
do i=4,3*NA
  mass(i)=hi
end do
end subroutine mass_vec
!-------------------------------------------------------
subroutine pot
implicit none
open(unit=1,file='n2.data',status='old')
do i=1,Nq1
  read(1,*) pot1(i)
end do
 close(unit=1)
end subroutine pot
!-------------------------------------------------------
subroutine coord
implicit none
q1(1)=0.0d0;q1(2)=0.0d0;q1(3)=-0.100108188344758d0;q1(4)=-0.336746283751144d0;q1(5)=-0.336746283751144d0;q1(6)=0.497306436419762d0
q1(7)=0.336746283751144d0;q1(8)=0.336746283751144d0;q1(9)=0.497306436419762d0;q1(10)=0.074320862204978d0
q1(11)=-0.074320862204978d0;q1(12)=0.099157366092731d0;q1(13)=-0.074320862204978d0;q1(14)=0.074320862204978d0
q1(15)=0.099157366092731d0

q2(1)=0.0d0;q2(2)=0.0d0;q2(3)=0.0d0;q2(4)=0.192478578122720d0;q2(5)=0.192478578122720d0;q2(6)=-0.419409100911881d0
q2(7)=-0.192478578122720d0;q2(8)=-0.192478578122720d0;q2(9)=-0.419409100911881d0;q2(10)=0.192478578122720d0
q2(11)=-0.192478578122720d0;q2(12)=0.419409100911881d0;q2(13)=-0.192478578122720d0;q2(14)=0.192478578122720d0
q2(15)=0.419409100911881d0
end subroutine coord
!-------------------------------------------------------
end module global_param
