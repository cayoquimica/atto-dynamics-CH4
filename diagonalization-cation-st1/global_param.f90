module global_param
implicit none
integer, parameter                         :: dp = kind(1.d0) 
integer                                    :: i,j,k,ll
integer, parameter                         :: NA=5 !Number of atoms
integer, parameter                         :: Nst=1 !number of electronic states
integer, parameter                         :: Nq1=146 !126 !number of points of the grid along q1
integer, parameter                         :: Nq2=184 !179 !number of points of the grid along q2
integer, parameter                         :: s=Nq1*Nq2 !counter to be used as index - means the size of a 1 state matrix, which is Nq1*Nq2
complex(kind=dp)                           :: im
real(kind=dp), parameter                   :: pi = 3.141592653589793d0
real(kind=dp), parameter                   :: sq1=0.08d0 !step in q1 in atomic units
real(kind=dp), parameter                   :: sq2=0.07d0 !step in q1 in atomic units
real(kind=dp), parameter                   :: saw2au=1822.889950851334d0 !mass conversion factor from Standard atomic weight to atomic units
real(kind=dp), parameter                   :: car=12.011d0*saw2au !Carbon mass in atomic units
real(kind=dp), parameter                   :: hi=1.00794d0*saw2au !Hidrogen mass in atomic units
real(kind=dp), parameter                   :: mtotal = car + 4.d0 * hi ! total mass of CH4
real(kind=dp)                              :: mass(3*NA) !3N Vector with the masses of the atoms
real(kind=dp)                              :: q1(3*NA),q2(3*NA),q1i(3*NA),q2i(3*NA),ai,bi,aii,bii !Conversion vector from internal coordinates q1 and q2 to cartesian coordinates
real(kind=dp)                              :: co1(Nst*Nq1*Nq2),co2(Nst*Nq1*Nq2)
real(kind=dp), dimension(:,:), allocatable :: pot1
real(kind=dp), dimension(:,:), allocatable :: Ha! Hamiltonian matrix
!real(kind=dp), dimension(:,:), allocatable :: ham ! Variable to store the Hamiltonian matrix temporary
!real(kind=dp), dimension(:,:), allocatable :: moq1,moq2! Momentum matrix
real(kind=dp)                              :: mass1,mass2,mass3 !Reduced masses to be used in the second derivative
!real(kind=dp)                              :: Et !Final electrical field of the pulse
!real(kind=dp), parameter                   :: t00 = 30.d0 !time where the pulse is centered
!real(kind=dp), parameter                   :: phase = 0.d0 !phase factor for the pulse related to the gaussian envelope
!real(kind=dp), parameter                   :: freq = 0.056937d0 !frequency of the pulse, in a.u.
!real(kind=dp), parameter                   :: sig = 5 ! 50 approx 1200 attoseconds - width of the gaussian envelop
!real(kind=dp), parameter                   :: E00 = 0.0d0 !0.05d0 !Electric field intensity
real(kind=dp)                              :: ind1,ind2,ind3,const1,const2,const3

contains

subroutine hamiltonian_matrix
! This subroutine builds the hamiltonian matrrix and stores it in variable Ha. 
! It also creates a temporary copy of Ha in the variable ham, that can be modified to include the time dependent interaction of the dipoles with the pulse.
! So, ham is the real matrix that will be used in the code and both ham and Ha will be global variables.
allocate (pot1(Nq1,Nq2))
allocate (Ha(Nst*Nq1*Nq2,Nst*Nq1*Nq2))
!allocate (ham(Nst*Nq1*Nq2,Nst*Nq1*Nq2))

ll=0
do k=1,Nst
  do j=1,Nq2
    do i=1,Nq1
      ll=ll+1
      co1(ll)=(i-Nq1/2.d0-1.d0)*sq1+sq1/2.d0!(Nq1-Nq11/2.d0-i)*sq1+sq1/2.d0
      co2(ll)=(j-(Nq2/2.d0)-1.d0)*sq2
    end do
  end do
end do

!-------------------------------------------------------------------!
! Loading electronic structure data: Energy and dipole moments      !
open(unit=1,file='v1f.txt',status='old')                            !
do i=1,Nq1                                                          !
  read(1,*) pot1(i,:)                                               !
end do                                                              !
 close(unit=1)                                                      !
!-------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------------------------------------!
!Defining the coefficients for coordinates q1 and q2 in cartesian                                                                   !
q1(1)=0.0d0;q1(2)=0.0d0;q1(3)=-0.100108188344758d0;q1(4)=-0.336746283751144d0;q1(5)=-0.336746283751144d0;q1(6)=0.497306436419762d0  !
q1(7)=0.336746283751144d0;q1(8)=0.336746283751144d0;q1(9)=0.497306436419762d0;q1(10)=0.074320862204978d0                            !
q1(11)=-0.074320862204978d0;q1(12)=0.099157366092731d0;q1(13)=-0.074320862204978d0;q1(14)=0.074320862204978d0                       !
q1(15)=0.099157366092731d0                                                                                                          !
                                                                                                                                    !
q2(1)=0.0d0;q2(2)=0.0d0;q2(3)=0.0d0;q2(4)=0.192478578122720d0;q2(5)=0.192478578122720d0;q2(6)=-0.419409100911881d0                  !
q2(7)=-0.192478578122720d0;q2(8)=-0.192478578122720d0;q2(9)=-0.419409100911881d0;q2(10)=0.192478578122720d0                         !
q2(11)=-0.192478578122720d0;q2(12)=0.419409100911881d0;q2(13)=-0.192478578122720d0;q2(14)=0.192478578122720d0                       !
q2(15)=0.419409100911881d0                                                                                                          !
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
!=================================================================================================================================!
!Building the Hamiltonian matrix                                                                                                  !
!Here, the Hamiltonian matrix will include kinetic and potential energy only.                                                     !
!NAC and dipoles are inserted in the main program.                                                                                !
!It is build with blocks as:                                                                                                      !
! | q1,q2=1,st1       0               0               0               0               0      |                                    !
! |      0       q1,q2=2,st1          0               0               0               0      |                                    !
! |     ...          ...             ...             ...             ...             ...     |                                    !
! |      0            0          q1,q2=1,st2          0               0               0      |                                    !
! |      0            0               0          q1,q2=2,st2          0               0      |                                    !
! |     ...          ...             ...             ...             ...             ...     |                                    !
! |      0            0               0               0          q1,q2=1,st3          0      |                                    !
! |      0            0               0               0               0          q1,q2=2,st3 |                                    !
! |     ...          ...             ...             ...             ...             ...     |                                    !
!-----------------------------------------------------------------------------!                                                   !
!Constants for the finite difference second derivatives, mass scalled already !                                                   ! 
const1 = 1.d0/(12.d0*sq1**2.d0)*mass1 ! for q1 in 1 dimension                 !                                                   !
const2 = 1.d0/(12.d0*sq2**2.d0)*mass2 ! for q1 in 1 dimension                 !                                                   !
const3 =-1.d0/( 2.d0*sq1*sq2  )*mass3 ! for derivative in q1 and q2           !                                                   !
!-----------------------------------------------------------------------------!                                                   !
Ha=0.d0 ! Creating the whole matrix and defining every element as zero                                                            !
!----------------------------------------------------------------------------------------------------------------------------|    !
! 5 step finite difference second derivative in 1D along q1 (d^2/dq1^2)                                                      |    !
! Each element of the above matrix is going to be a square q1 sized box in which we have q1 varying for a fixed q2           |    !
! | x x x 0 0 0 0 ..... ....... 0 0 0 .....|                                                                                 |    !
! | x x x x 0 0 0 ..... ....... 0 0 0 .....|                                                                                 |    !
! | x x x x x 0 0 ..... ....... 0 0 0 .....|                                                                                 |    !
! | 0 x x x x x 0 ..... ....... 0 0 0 .....|                                                                                 |    !
! | . . . . . . . ..... . . . . . . . .....|                                                                                 |    !
! |..... 0 0 0 ....... ..... 0 x x x x x 0 |                                                                                 |    !
! |..... 0 0 0 ....... ..... 0 0 x x x x x |                                                                                 |    !
! |..... 0 0 0 ....... ..... 0 0 0 x x x x |                                                                                 |    !
! |..... 0 0 0 ....... ..... 0 0 0 0 x x x |                                                                                 |    !
ind1=-30.d0*const1 !factor for 5 point second derivative finite diference in 1D along q1                                     |    !
ind2= 16.d0*const1 !factor for 5 point second derivative finite diference in 1D along q1                                     |    !
ind3=- 1.d0*const1 !factor for 5 point second derivative finite diference in 1D along q1                                     |    !
do k=0,0 ! running for all electronic states                                                                             |    !
  do j=0*Nq1,(Nq2-1)*Nq1,Nq1 ! Moving through the boxes, from the fisrt to the last one                                      |    !
!.....................................................................................................................!      |    !
! Doing the special cases of the 5 points finite difference that happens in the borders of each q1 sized box          !      |    !
! Upper and left part of the borders                                                                                  !      |    !
Ha( k   + j       + 1     , k   + j       + 1     ) = Ha( k   + j       + 1     , k   + j       + 1     ) + ind1      !      |    !
Ha( k   + j       + 1     , k   + j       + 2     ) = Ha( k   + j       + 1     , k   + j       + 2     ) + ind2      !      |    !
Ha( k   + j       + 1     , k   + j       + 3     ) = Ha( k   + j       + 1     , k   + j       + 3     ) + ind3      !      |    !
! Upper and left part of the borders, second row                                                                      !      |    !
Ha( k   + j       + 2     , k   + j       + 1     ) = Ha( k   + j       + 2     , k   + j       + 1     ) + ind2      !      |    !
Ha( k   + j       + 2     , k   + j       + 2     ) = Ha( k   + j       + 2     , k   + j       + 2     ) + ind1      !      |    !
Ha( k   + j       + 2     , k   + j       + 3     ) = Ha( k   + j       + 2     , k   + j       + 3     ) + ind2      !      |    !
Ha( k   + j       + 2     , k   + j       + 4     ) = Ha( k   + j       + 2     , k   + j       + 4     ) + ind3      !      |    !
! Botton and right part of the borders, antepenulmate row                                                             !      |    !
Ha( k   + j       + Nq1-1 , k   + j       + Nq1-3 ) = Ha( k   + j       + Nq1-1 , k   + j       + Nq1-3 ) + ind3      !      |    !
Ha( k   + j       + Nq1-1 , k   + j       + Nq1-2 ) = Ha( k   + j       + Nq1-1 , k   + j       + Nq1-2 ) + ind2      !      |    !
Ha( k   + j       + Nq1-1 , k   + j       + Nq1-1 ) = Ha( k   + j       + Nq1-1 , k   + j       + Nq1-1 ) + ind1      !      |    !
Ha( k   + j       + Nq1-1 , k   + j       + Nq1   ) = Ha( k   + j       + Nq1-1 , k   + j       + Nq1   ) + ind2      !      |    !
! Botton and right part of the borders, last row                                                                      !      |    !
Ha( k   + j       + Nq1   , k   + j       + Nq1-2 ) = Ha( k   + j       + Nq1   , k   + j       + Nq1-2 ) + ind3      !      |    !
Ha( k   + j       + Nq1   , k   + j       + Nq1-1 ) = Ha( k   + j       + Nq1   , k   + j       + Nq1-1 ) + ind2      !      |    !
Ha( k   + j       + Nq1   , k   + j       + Nq1   ) = Ha( k   + j       + Nq1   , k   + j       + Nq1   ) + ind1      !      |    !
!.....................................................................................................................!      |    !
!________________________________________________________________________________________________________________________!   |    !
! From the third line to the antepenulmate one of each box - non special cases                                           !   |    !
    do i=1+2,Nq1-2 ! do through single box                                                                               !   |    !
Ha( k   + j       + i     , k   + j       + i-2   ) = Ha( k   + j       + i     , k   + j       + i-2   ) + ind3         !   |    !
Ha( k   + j       + i     , k   + j       + i-1   ) = Ha( k   + j       + i     , k   + j       + i-1   ) + ind2         !   |    !
Ha( k   + j       + i     , k   + j       + i     ) = Ha( k   + j       + i     , k   + j       + i     ) + ind1         !   |    !
Ha( k   + j       + i     , k   + j       + i+1   ) = Ha( k   + j       + i     , k   + j       + i+1   ) + ind2         !   |    !
Ha( k   + j       + i     , k   + j       + i+2   ) = Ha( k   + j       + i     , k   + j       + i+2   ) + ind3         !   |    !
    end do ! end of do through single box                                                                                !   |    !
!________________________________________________________________________________________________________________________!   |    !
  end do ! end of do through boxes                                                                                           |    !
end do ! end of do through electronic states                                                                                 |    !
!----------------------------------------------------------------------------------------------------------------------------|    !
!---------------------------------------------------------------------------------------------------------------------|           !
!                                                                                                                     |           !
! 5 step finite difference second derivative in 1D along q2 (d^2/dq2^2)                                               |           !
! The simbol 'o' is the center of the 5 points finite difference                                                      |           !
! | |o        ||x        ||x        ||         ||         ||         ||         ||         ||         |               |           !
! | | x       || o       || x       || x       ||         ||         ||         ||         ||         |               |           !
! | |  x      ||  x      ||  o      ||  x      ||  x      ||         ||         ||         ||         | q2=1          |           !
! | |         ||   x     ||   x     ||   o     ||   x     ||   x     ||         ||         ||         |               |           !
! | |         ||         ||    x    ||    x    ||    o    ||    x    ||    x    ||         ||         |               |           !
! |                                                                                                   |               |           !
! | |x        ||o        ||x        ||x        ||         ||         ||         ||         ||         |               |           !
! | | x       || x       || o       || x       ||         ||         ||         ||         ||         |               |           !
! | |         ||  x      ||  x      ||  o      ||  x      ||         ||         ||         ||         | q2=2          |           !
! | |         ||         ||   x     ||   x     ||   o     ||   x     ||   x     ||         ||         |               |           !
! | |         ||         ||         ||    x    ||    x    ||    o    ||    x    ||         ||         |               |           !
! |                                                                                                   |               |           !
! | |x        ||x        ||o        ||x        ||x        ||         ||         ||         ||         |               |           !
! | |         || x       || x       || o       || x       || x       ||         ||         ||         |               |           !
! | |         ||         ||  x      ||  x      ||  o      ||  x      ||  x      ||         ||         | q2=3          |           !
! | | q1 size ||         ||         ||   x     ||   x     ||   o     ||   x     ||   x     ||         |               |           !
! | |   box   ||         ||         ||         ||    x    ||    x    ||    o    ||    x    ||    x    |               |           !
! |                                                                                                   |               |           !
! | |         ||x        ||x        ||o        ||x        ||x        ||         ||         ||         |               |           !
! | |         ||         || x       || x       || o       || x       || x       ||         ||         |               |           !
! | |         ||         ||         ||  x      ||  x      ||  o      ||  x      ||  x      ||         | q2=4          |           !
! | | q1 size || q1 size ||         ||         ||   x     ||   x     ||   o     ||   x     ||   x     |               |           !
! | |   box   ||   box   ||         ||         ||         ||    x    ||    x    ||    o    ||    x    |               |           !
!                                                                                                                     |           !
ind1=-30.d0*const2 !factor for 5 point second derivative finite diference in 1D along q2                              |           !
ind2= 16.d0*const2 !factor for 5 point second derivative finite diference in 1D along q2                              |           !
ind3=- 1.d0*const2 !factor for 5 point second derivative finite diference in 1D along q2                              |           !
!                                                                                                                     |           !
! Each element will be defined by the following indexes:                                                              |           !
!Ha( k   + j       + i     , k   + j       + i     ) = Ha( k   + j       + i     , k   + j       + i     ) + ind1     |           !
! where 'k' is running through the different electronic states, so it is going to be:                                 |           !
!                    It has to start on zero, so thats why the -1                                                     |           !
!k=(st-1)*Nq1*Nq2    And when st changes it means that we moved Nq1*Nq2 along the matrix                              |           !
!                                                                                                                     |           !
! 'j' is going to count the q1 sized boxes along the matrix, so it is going to be:                                    |           !
!                    n is the real counter for the boxes, so it goes from 0 to Nq2-1 if one wants the                 |           !
!                    beginning of the box or from 1 to Nq2 if one wants the end of the box                            |           !
!j=(n-1)*Nq1         j is actually pointing the row or column at the beginning or the end of the box                  |           !
!                    Each box has Nq1 x Nq1 size, so changing box means move Nq1 rows or columns                      |           !
!                                                                                                                     |           !
! 'i' is counting rows or columns inside a q1 sized box, so it simply goes from 1 to Nq1                              |           !
!                                                                                                                     |           !
! The '+s' will mean that we are going to the last row of the last box                                                |           !
!                                                                                                                     |           !
do k=0,0  !2*s,s ! Running through every state, the step size is Nq1*Nq2                                                  |           !
!................................................................................................................!    |           !
! Doing the special cases of the 5 points finite difference that happens in first and last two boxes             !    |           !
  do i=1,Nq1 ! Moving inside the box                                                                             !    |           !
! Upper and left part of the borders - first box of all                                                          !    |           !
Ha( k   +   0*Nq1 + i     , k   +   0*Nq1 + i     ) = Ha( k   +   0*Nq1 + i     , k   +   0*Nq1 + i     ) + ind1 !    |           !
Ha( k   +   0*Nq1 + i     , k   +   1*Nq1 + i     ) = Ha( k   +   0*Nq1 + i     , k   +   1*Nq1 + i     ) + ind2 !    |           !
Ha( k   +   0*Nq1 + i     , k   +   2*Nq1 + i     ) = Ha( k   +   0*Nq1 + i     , k   +   2*Nq1 + i     ) + ind3 !    |           !
! Botton and right part of the borders, last box of all                                                          !    |           !
Ha( k+s -   0*Nq1 - Nq1+i , k+s -   2*Nq1 - Nq1+i ) = Ha( k+s -   0*Nq1 - Nq1+i , k+s -   2*Nq1 - Nq1+i ) + ind3 !    |           !
Ha( k+s -   0*Nq1 - Nq1+i , k+s -   1*Nq1 - Nq1+i ) = Ha( k+s -   0*Nq1 - Nq1+i , k+s -   1*Nq1 - Nq1+i ) + ind2 !    |           !
Ha( k+s -   0*Nq1 - Nq1+i , k+s -   0*Nq1 - Nq1+i ) = Ha( k+s -   0*Nq1 - Nq1+i , k+s -   0*Nq1 - Nq1+i ) + ind1 !    |           !
! Upper and left part of the borders - second box of all                                                         !    |           !
Ha( k   +   1*Nq1 + i     , k   +   0*Nq1 + i     ) = Ha( k   +   1*Nq1 + i     , k   +   0*Nq1 + i     ) + ind2 !    |           !
Ha( k   +   1*Nq1 + i     , k   +   1*Nq1 + i     ) = Ha( k   +   1*Nq1 + i     , k   +   1*Nq1 + i     ) + ind1 !    |           !
Ha( k   +   1*Nq1 + i     , k   +   2*Nq1 + i     ) = Ha( k   +   1*Nq1 + i     , k   +   2*Nq1 + i     ) + ind2 !    |           !
Ha( k   +   1*Nq1 + i     , k   +   3*Nq1 + i     ) = Ha( k   +   1*Nq1 + i     , k   +   3*Nq1 + i     ) + ind3 !    |           !
! Botton and right part of the borders, penulmate box of all                                                     !    |           !
Ha( k+s -   1*Nq1 - Nq1+i , k+s -   3*Nq1 - Nq1+i ) = Ha( k+s -   1*Nq1 - Nq1+i , k+s -   3*Nq1 - Nq1+i ) + ind3 !    |           !
Ha( k+s -   1*Nq1 - Nq1+i , k+s -   2*Nq1 - Nq1+i ) = Ha( k+s -   1*Nq1 - Nq1+i , k+s -   2*Nq1 - Nq1+i ) + ind2 !    |           !
Ha( k+s -   1*Nq1 - Nq1+i , k+s -   1*Nq1 - Nq1+i ) = Ha( k+s -   1*Nq1 - Nq1+i , k+s -   1*Nq1 - Nq1+i ) + ind1 !    |           !
Ha( k+s -   1*Nq1 - Nq1+i , k+s -   0*Nq1 - Nq1+i ) = Ha( k+s -   1*Nq1 - Nq1+i , k+s -   0*Nq1 - Nq1+i ) + ind2 !    |           !
  end do ! end of do through single box                                                                          !    |           !
!................................................................................................................!    |           !
!___________________________________________________________________________________________________________________! |           !
! Now doing the rest of the boxes                                                                                   ! |           !
  do j=2*Nq1,(Nq2-3)*Nq1,Nq1 ! Moving through the boxes, from the third to the antipenulmate                        ! |           !
    do i=1,Nq1 ! Moving inside the box                                                                              ! |           !
Ha( k   + j       + i     , k   + j-2*Nq1 + i     ) = Ha( k   + j       + i     , k   + j-2*Nq1 + i     ) + ind3    ! |           !
Ha( k   + j       + i     , k   + j-1*Nq1 + i     ) = Ha( k   + j       + i     , k   + j-1*Nq1 + i     ) + ind2    ! |           !
Ha( k   + j       + i     , k   + j       + i     ) = Ha( k   + j       + i     , k   + j       + i     ) + ind1    ! |           !
Ha( k   + j       + i     , k   + j+1*Nq1 + i     ) = Ha( k   + j       + i     , k   + j+1*Nq1 + i     ) + ind2    ! |           !
Ha( k   + j       + i     , k   + j+2*Nq1 + i     ) = Ha( k   + j       + i     , k   + j+2*Nq1 + i     ) + ind3    ! |           !
    end do ! end of do through single box                                                                           ! |           !
  end do ! end of do through boxes                                                                                  ! |           !
!___________________________________________________________________________________________________________________! |           !
end do ! end of do through electronic states                                                                          |           !
!---------------------------------------------------------------------------------------------------------------------|           !
!-------------------------------------------------------------------------------------------------------------------------|       !
! 7 steps finite difference second derivative in 2D                                                                       |       !
!              |q1                                                                                                        |       !
!              |                                                                                                          |       !
!              |                                                                                                          |       !
!            7 2             q2                                                                                           |       !
!----------- 5 1 4 ------------>                                                                                          |       !
!              3 6                                                                                                        |       !
!              |                                                                                                          |       !
!              |                                                                                                          |       !
!              |                                                                                                          |       !
!              V                                                                                                          |       !
!                                                                                                                         |       !
! 1 3 _ _ _ _    4 6 _ _ _ _                                                                                              |       !
! 2 1 3 _ _ _    _ 4 6 _ _ _                                                                                              |       !
! _ 2 1 3 _ _    _ _ 4 6 _ _                                                                                              |       !
! _ _ 2 1 3 _    _ _ _ 4 6 _                                                                                              |       !
! _ _ _ 2 1 3    _ _ _ _ 4 6                                                                                              |       !
! _ _ _ _ 2 1    _ _ _ _ _ 4                                                                                              |       !
!                                                                                                                         |       !
! 5 _ _ _ _ _    1 3 _ _ _ _    4 6 _ _ _ _                                                                               |       !
! 7 5 _ _ _ _    2 1 3 _ _ _    _ 4 6 _ _ _                                                                               |       !
! _ 7 5 _ _ _    _ 2 1 3 _ _    _ _ 4 6 _ _                                                                               |       !
! _ _ 7 5 _ _    _ _ 2 1 3 _    _ _ _ 4 6 _                                                                               |       !
! _ _ _ 7 5 _    _ _ _ 2 1 3    _ _ _ _ 4 6                                                                               |       !
! _ _ _ _ 7 5    _ _ _ _ 2 1    _ _ _ _ _ 4                                                                               |       !
!                                                                                                                         |       !
ind1=- 2.d0*const3 !factor for 7 point second derivative finite diference in 2D                                           |       !
ind2=  1.d0*const3 !factor for 7 point second derivative finite diference in 2D                                           |       !
ind3=- 1.d0*const3 !factor for 7 point second derivative finite diference in 2D                                           |       !
! Second derivative along q1 and q2 for state 1                                                                           |       !
! Upper and left part of the borders                                                                                      |       !
! Running along q1                                                                                                        |       !
!                                                                                                                         |       !
! Doing the special case of the 7 points finite difference                                                                |       !
! that happens when moving along q1 in first and last row and collumn of each box                                         |       !
do k=0,0 !2*s,s ! running for all electronic states                                                                          |       !
!..................................................................................................................!      |       !
! Doing the special cases that happens in the borders of each q1 sized box when moving along q1                    !      |       !
! Upper and left part of the borders - first row of each box                                                       !      |       !
  do j=0*Nq1,(Nq2-1)*Nq1,Nq1 ! Moving through the boxes, from the fisrt to the last one                            !      |       !
Ha( k   + j       + 1     , k   + j       + 2     ) = Ha( k   + j       + 1     , k   + j       + 2     ) + ind2   !      |       !
! Botton and right part of the borders, last row of each box                                                       !      |       !
Ha( k   + j       + Nq1   , k   + j       + Nq1-1 ) = Ha( k   + j       + Nq1   , k   + j       + Nq1-1 ) + ind2   !      |       !
  end do ! end of do through boxes                                                                                 !      |       !
! Now doing the rest of the lines of the first and last box                                                        !      |       ! 
  do i=2,Nq1-1 ! Moving inside the first and last box, but avoiding the first and last row                         !      |       !
Ha( k   + 0       + i     , k   + 0       + i-1   ) = Ha( k   + 0       + i     , k   + 0       + i-1   ) + ind2   !      |       !
Ha( k   + 0       + i     , k   + 0       + i+1   ) = Ha( k   + 0       + i     , k   + 0       + i+1   ) + ind2   !      |       !
Ha( k+s -Nq1      + i     , k+s -Nq1      + i-1   ) = Ha( k+s -Nq1      + i     , k+s -Nq1      + i-1   ) + ind2   !      |       !
Ha( k+s -Nq1      + i     , k+s -Nq1      + i+1   ) = Ha( k+s -Nq1      + i     , k+s -Nq1      + i+1   ) + ind2   !      |       !
  end do ! end of do through single box                                                                            !      |       !
!..................................................................................................................!      |       !
! Doing the special cases that happens on the first and last box for each state when moving along q2               !      |       !
! Upper and left part of the borders - first box                                                                   !      |       !
  do i=1,Nq1 ! Moving inside the box                                                                               !      |       !
Ha( k   + 0       + i     , k   + Nq1     + i     ) = Ha( k   + 0       + i     , k   + Nq1     + i     ) + ind2   !      |       !
! Botton and right part of the borders, last box                                                                   !      |       !
Ha( k+s - Nq1     + i     , k+s - 2*Nq1   + i     ) = Ha( k+s - Nq1     + i     , k+s - 2*Nq1   + i     ) + ind2   !      |       !
  end do ! end of do through box                                                                                   !      |       !
! Now doing the first and last row of each of the others boxes                                                     !      |       !
  do j=1*Nq1,(Nq2-2)*Nq1,Nq1 ! Moving through the boxes, from the second to penulmate one                          !      |       !
Ha( k   + j       + 1     , k   + j-Nq1   + 1     ) = Ha( k   + j       + 1     , k   + j-Nq1   + 1     ) + ind2   !      |       !
Ha( k   + j       + 1     , k   + j+Nq1   + 1     ) = Ha( k   + j       + 1     , k   + j+Nq1   + 1     ) + ind2   !      |       !
Ha( k   + j       + Nq1   , k   + j-Nq1   + Nq1   ) = Ha( k   + j       + Nq1   , k   + j-Nq1   + Nq1   ) + ind2   !      |       !
Ha( k   + j       + Nq1   , k   + j+Nq1   + Nq1   ) = Ha( k   + j       + Nq1   , k   + j+Nq1   + Nq1   ) + ind2   !      |       !
  end do ! end of do through boxes                                                                                 !      |       !
!..................................................................................................................!      |       !
! Now the special cases when moving along q1 and q2 at the same time                                               !      |       !
! Cross terms in the first box and first row                                                                       !      |       !
Ha( k   + 0       + 1     , k   + Nq1     + 2     ) = Ha( k   + 0       + 1     , k   + Nq1     + 2     ) + ind3   !      |       !
! Cross terms in the last box and last row                                                                         !      |       !
Ha( k+s + 0       + 0     , k+s - Nq1     - 1     ) = Ha( k+s + 0       + 0     , k+s - Nq1     - 1     ) + ind3   !      |       !
! Now for the rest of the rows of the first and last box                                                           !      |       !
  do i=1,Nq1-2 ! Moving inside the first and last box, but avoiding the first and last row                         !      |       !
Ha( k   + 0       + i+1   , k   + Nq1     + i+2   ) = Ha( k   + 0       + i+1   , k   + Nq1     + i+2   ) + ind3   !      |       !
Ha( k+s + 0       - i     , k+s - 1*Nq1   - (i+1) ) = Ha( k+s + 0       - i     , k+s - 1*Nq1   - (i+1) ) + ind3   !      |       !
  end do ! end of do through single box                                                                            !      |       !
! Now for the first and last rows of each of the others boxes                                                      !      |       !
  do j=1*Nq1,(Nq2-2)*Nq1,Nq1 ! Moving through the boxes, from the second to penulmate one                          !      |       !
Ha( k   + j       + 1     , k   + j + Nq1 + 1+1   ) = Ha( k   + j       + 1     , k   + j + Nq1 + 1+1   ) + ind3   !      |       !
Ha( k   + j       + Nq1   , k   + j - Nq1 + Nq1-1 ) = Ha( k   + j       + Nq1   , k   + j - Nq1 + Nq1-1 ) + ind3   !      |       !
  end do ! end of do through boxes                                                                                 !      |       !
!..................................................................................................................!      |       !
! Finally, doing for the rest of the boxes, always avoinding the first and last row                                !      |       !
  do j=1*Nq1,(Nq2-2)*Nq1,Nq1 ! Moving through the boxes, from the second to penulmate one                          !      |       !
    do i=2,Nq1-1 ! Moving inside the first and last box, but avoiding the first and last row                       !      |       !
Ha( k   + j       + i     , k   + j + Nq1 + i+1   ) = Ha( k   + j       + i     , k   + j + Nq1 + i+1   ) + ind3   !      |       !
Ha( k   + j       + i     , k   + j       + i+1   ) = Ha( k   + j       + i     , k   + j       + i+1   ) + ind2   !      |       !
Ha( k   + j       + i     , k   + j       + i-1   ) = Ha( k   + j       + i     , k   + j       + i-1   ) + ind2   !      |       !
Ha( k   + j       + i     , k   + j + Nq1 + i     ) = Ha( k   + j       + i     , k   + j + Nq1 + i     ) + ind2   !      |       !
Ha( k   + j       + i     , k   + j - Nq1 + i     ) = Ha( k   + j       + i     , k   + j - Nq1 + i     ) + ind2   !      |       !
Ha( k   + j       + i     , k   + j - Nq1 + i-1   ) = Ha( k   + j       + i     , k   + j - Nq1 + i-1   ) + ind3   !      |       !
    end do ! end of do through single box                                                                          !      |       !
  end do ! end of do through boxes                                                                                 !      |       !
!..................................................................................................................!      |       !
! Now doing only the diagonal of each box                                                                          !      |       !
  do j=0*Nq1,(Nq2-1)*Nq1,Nq1 ! Moving through the boxes, from the fisrt to the last one                            !      |       !
    do i=1,Nq1 ! Moving inside the box                                                                             !      |       !
Ha( k   + j       + i     , k   + j       + i     ) = Ha( k   + j       + i     , k   + j       + i     ) + ind1   !      |       !
    end do ! end of do through single box                                                                          !      |       !
  end do ! end of do through boxes                                                                                 !      |       !
end do ! end of do through electronic states                                                                       !      |       !
!-------------------------------------------------------------------------------------------------------------------------|       !
! Finishing the kinetic energy part                                                                                               !
Ha(: , :)= - 1.d0/2.d0 * Ha(: , :)                                                                                                !
! Adding potential energy - moving only along the diagonal of the Hamiltonian                                                     !
do j=0*Nq1,(Nq2-1)*Nq1,Nq1 ! Moving through the boxes, from the fisrt to the last one                                             !
  do i=1,Nq1 ! Moving inside the box                                                                                              !
    Ha( 0*s + j + i , 0*s + j + i ) = Ha( 0*s + j + i , 0*s + j + i ) + pot1 (i,(j+Nq1)/Nq1)                                      !
  end do ! end of do through single box                                                                                           !
end do ! end of do through boxes                                                                                                  !
!=================================================================================================================================!                       

!ham=ha
!ham=0.d0

end subroutine hamiltonian_matrix

end module global_param
