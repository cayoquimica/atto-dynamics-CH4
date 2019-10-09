program load_hamiltonian
!maybe is useful: Function MERGE:
!test a condition, if true, assign a value X, if false assign value Y:
!the condition can be a boolean variable defined before. 
!result = merge( X, Y, i.eq.j)
!result, X and Y will be of same type
use global_param
use omp_lib
implicit none
!external rkdumb,HA_calc,first_derivative_matrix,momentum_calc
integer n,jj,ii,nana !m=dimension of the Hamiltonian (Nq1*Nq2*Nst x Nq1*Nq2*Nst)
real(kind=dp) :: ompt0,ompt1,x1,x2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! I modified the Ha, moq1 and moq2 matrices building so that their indexes begin with 0 instead of 1. 
! It is not the ideal way, I use a auxiliar matrix to build with index 1 as I already had written in the code and then I put it
! in a matrix with index starting in 0.
! Later would be better and faster if I already build the matrices with index 0.     <<<---------------

n=Nst*Nq1*Nq2 ! m=Nq1*Nq2*st %size of the Hamiltonian

!$ ompt0 = omp_get_wtime()
call load_data

call hamiltonian_matrix
write(*,*)'Kinetic and potential energy loaded'
call first_derivative_matrix_q1
!$ x1 = omp_get_wtime()
!Creating the CSR vectors --------------------------------------------------------------------------------------------------!
!This is for the first derivative with respect of q1                                                                        !
k=0                                                                                                                         !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (moq1(i,j) /= 0.d0) then                                                                                             !
      k=k+1                                                                                                                 !
    end if                                                                                                                  !
  end do                                                                                                                    !
end do                                                                                                                      !
k_moq1=k                                                                                                                    !
write(*,*)'k moq1= ',k                                                                                                      !
allocate(moq1_val(0:k-1),moq1_rowc(0:n), moq1_row_col(0:k-1,0:1))                                                           !
!                                                                                                                           !
moq1_rowc=0.d0                                                                                                              !
k=0                                                                                                                         !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (moq1(i,j) /= 0.d0) then                                                                                             !
      moq1_val(k)=moq1(i,j)  !storing each non-zero element                                                                 !
      moq1_row_col(k,0)=i !storing the row index of each non-zero element                                                   !
      moq1_row_col(k,1)=j !storing the column index of each non-zero element                                                !
      k=k+1                                                                                                                 !
    end if                                                                                                                  !
  end do                                                                                                                    !
  moq1_rowc(i+1)=k !storing the counting of non-zero elements in each row                                                   !
end do                                                                                                                      !
!---------------------------------------------------------------------------------------------------------------------------!
!$ x2 = omp_get_wtime()
write(*,*) 'Time to create CSR vectors = ',x2 - x1, "seconds"
write(*,*)'Csr vectors for first derivative along q1 created'
call first_derivative_matrix_q2
!Creating the CSR vectors --------------------------------------------------------------------------------------------------!
!This is for the first derivative with respect of q1                                                                        !
k=0                                                                                                                         !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (moq2(i,j) /= 0.d0) then                                                                                             !
      k=k+1                                                                                                                 !
    end if                                                                                                                  !
  end do                                                                                                                    !
end do                                                                                                                      !
k_moq2=k                                                                                                                    !
write(*,*)'k moq2= ',k                                                                                                      !
allocate(moq2_val(0:k-1),moq2_rowc(0:n), moq2_row_col(0:k-1,0:1))                                                           !
!                                                                                                                           !
moq2_rowc=0.d0                                                                                                              !
k=0                                                                                                                         !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (moq2(i,j) /= 0.d0) then                                                                                             !
      moq2_val(k)=moq2(i,j)  !storing each non-zero element                                                                 !
      moq2_row_col(k,0)=i !storing the row index of each non-zero element                                                   !
      moq2_row_col(k,1)=j !storing the column index of each non-zero element                                                !
      k=k+1                                                                                                                 !
    end if                                                                                                                  !
  end do                                                                                                                    !
  moq2_rowc(i+1)=k !storing the counting of non-zero elements in each row                                                   !
end do                                                                                                                      !
!---------------------------------------------------------------------------------------------------------------------------!

write(*,*)'CSR vectors for first derivative along q2 created'
call nac_ha_modify

deallocate(moq1,moq2)                                                                                                                             !

!I still calculate angular momentum without using the sparse matrix CSR method.                                                                   !
!Deallocate these matrices to salve RAM at this point. They will be reallocated later again just to be used in the angular momentum subroutine    !
!Not perfect way to do it, but these matrices will be loaded two times now in order to not use to much RAM memory.                                !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!==============================================!
! Loading dipoles                              !
allocate (ax(0:n-1,0:n-1))                     !
!$OMP PARALLEL DO shared(ax)                   !
do i=0,n-1                                     !
  do j=0,n-1                                   !
    ax(i,j)=0.d0                               !
  end do                                       !
end do                                         !
!$OMP end parallel do                          !
ii=0                                           !
do i=0,Nq1-1                                   !
  do j=0,Nq2-1                                 !
    ax(ii    ,ii    ) = - dot_product( orientation, [ pdm1x(i+1,j+1),pdm1y(i+1,j+1),pdm1z(i+1,j+1) ] )  !
    ax(s+ii  ,s+ii  ) = - dot_product( orientation, [ pdm2x(i+1,j+1),pdm2y(i+1,j+1),pdm2z(i+1,j+1) ] )  !
    ax(2*s+ii,2*s+ii) = - dot_product( orientation, [ pdm3x(i+1,j+1),pdm3y(i+1,j+1),pdm3z(i+1,j+1) ] )  !
!                                              !
    ax(ii    ,s+ii  ) = - dot_product( orientation, [ tdm21x(i+1,j+1),tdm21y(i+1,j+1),tdm21z(i+1,j+1) ] )  !
    ax(s+ii  ,ii    ) = - dot_product( orientation, [ tdm21x(i+1,j+1),tdm21y(i+1,j+1),tdm21z(i+1,j+1) ] )  !
    ax(ii    ,2*s+ii) = - dot_product( orientation, [ tdm31x(i+1,j+1),tdm31y(i+1,j+1),tdm31z(i+1,j+1) ] )  !
    ax(2*s+ii,ii    ) = - dot_product( orientation, [ tdm31x(i+1,j+1),tdm31y(i+1,j+1),tdm31z(i+1,j+1) ] )  !
    ax(s+ii  ,2*s+ii) = - dot_product( orientation, [ tdm32x(i+1,j+1),tdm32y(i+1,j+1),tdm32z(i+1,j+1) ] )  !
    ax(2*s+ii,s+ii  ) = - dot_product( orientation, [ tdm32x(i+1,j+1),tdm32y(i+1,j+1),tdm32z(i+1,j+1) ] )  !
    ii=ii+1                                    !
  end do                                       !
end do                                         !
!==============================================!

allocate (ham(0:n-1,0:n-1))
!!$OMP PARALLEL DO shared(ax)
do i=0,n-1
  do j=0,n-1
    ham(i,j)=Ha(i,j)+ax(i,j)
  end do
end do

!!$OMP end parallel do
!Creating the CSR vectors --------------------------------------------------------------------------------------------------!
!This is for the the hamiltonian - dipoles will be in a separate vector                                                     !
k=0                                                                                                                         !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (ham(i,j) /= 0.d0) then                                                                                              !
      k=k+1                                                                                                                 !
    end if                                                                                                                  !
  end do                                                                                                                    !
end do                                                                                                                      !
k_Ha=k                                                                                                                      !
k_dip=k                                                                                                                     !
write(*,*)'k_Ha = ',k_Ha                                                                                                    !
allocate(Ha_val(0:k-1),Ha_rowc(0:n),Ha_row_col(0:k-1,0:1),dip_val(0:k-1))                                                   !
!                                                                                                                           !
Ha_rowc=0.d0                                                                                                                !
k=0                                                                                                                         !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (ham(i,j) /= 0.d0) then                                                                                              !
      Ha_val(k)=Ha(i,j)  !storing each non-zero element                                                                     !
      Ha_row_col(k,0)=i !storing the row index of each non-zero element                                                     !
      Ha_row_col(k,1)=j !storing the column index of each non-zero element                                                  !
      dip_val(k)=ax(i,j)  !storing each non-zero element                                                                    !
      k=k+1                                                                                                                 !
    end if                                                                                                                  !
  end do                                                                                                                    !
  Ha_rowc(i+1)=k !storing the counting of non-zero elements in each row                                                     ! 
end do                                                                                                                      !
!---------------------------------------------------------------------------------------------------------------------------!
deallocate(Ha,ham,ax)
!deallocate(ax)
call first_derivative_matrix_q1
call first_derivative_matrix_q2
!$ ompt1 = omp_get_wtime()
write(*,*) 'The time to load the matrices=',ompt1 - ompt0, "seconds"

open(unit=20,file='csr_vectors',status='unknown')
write(20,'(i12)')k_moq1
do i=0,k_moq1-1
  write(20,'(e23.15e3,2i12)')moq1_val(i),moq1_row_col(i,0),moq1_row_col(i,1)
end do
do i=0,n
  write(20,'(i12)')moq1_rowc(i)
end do

write(20,'(i12)')k_moq2
do i=0,k_moq2-1
  write(20,'(e23.15e3,2i12)')moq2_val(i),moq2_row_col(i,0),moq2_row_col(i,1)
end do
do i=0,n
  write(20,'(i12)')moq2_rowc(i)
end do

write(20,'(i12)')k_Ha
do i=0,k_Ha-1
  write(20,'(e23.15e3,2i12)')Ha_val(i),Ha_row_col(i,0),Ha_row_col(i,1)
end do
do i=0,n
  write(20,'(i12)')Ha_rowc(i)
end do

write(20,'(i12)')k_dip
do i=0,k_dip-1
  write(20,'(e23.15e3)')dip_val(i)
end do


call angular_momentum(n)

write(20,'(i12)')k_am
do i=0,k_am-1
  write(20,'(e23.15e3,2i12)')am_val(i),am_row_col(i,0),am_row_col(i,1)
end do
do i=0,n
  write(20,'(i12)')am_rowc(i)
end do

write(*,*)'FINISHED'

 close(unit=20)

end program load_hamiltonian



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine first_derivative_matrix_q1
use global_param
use omp_lib
implicit none
integer n
n=Nst*Nq1*Nq2
allocate (ax(n,n))
!$OMP parallel do shared(ax)
do i=1,n
  do j=1,n
    ax(i,j)=0.d0
  end do 
end do
!$OMP end parallel do
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
ax( k   + j       + 1     , k   + j       + 2     ) = ax( k   + j       + 1     , k   + j       + 2     ) + ind1      !      |    !
ax( k   + j       + 1     , k   + j       + 3     ) = ax( k   + j       + 1     , k   + j       + 3     ) - ind2      !      |    !
ax( k   + j       + 1     , k   + j       + 4     ) = ax( k   + j       + 1     , k   + j       + 4     ) + ind3      !      |    !
! Second line of each box                                                                                             !      |    !
ax( k   + j       + 2     , k   + j       + 1     ) = ax( k   + j       + 2     , k   + j       + 1     ) - ind1      !      |    !
ax( k   + j       + 2     , k   + j       + 3     ) = ax( k   + j       + 2     , k   + j       + 3     ) + ind1      !      |    !
ax( k   + j       + 2     , k   + j       + 4     ) = ax( k   + j       + 2     , k   + j       + 4     ) - ind2      !      |    !
ax( k   + j       + 2     , k   + j       + 5     ) = ax( k   + j       + 2     , k   + j       + 5     ) + ind3      !      |    !
! Third line of each box                                                                                              !      |    !
ax( k   + j       + 3     , k   + j       + 1     ) = ax( k   + j       + 3     , k   + j       + 1     ) + ind2      !      |    !
ax( k   + j       + 3     , k   + j       + 2     ) = ax( k   + j       + 3     , k   + j       + 2     ) - ind1      !      |    !
ax( k   + j       + 3     , k   + j       + 4     ) = ax( k   + j       + 3     , k   + j       + 4     ) + ind1      !      |    !
ax( k   + j       + 3     , k   + j       + 5     ) = ax( k   + j       + 3     , k   + j       + 5     ) - ind2      !      |    !
ax( k   + j       + 3     , k   + j       + 6     ) = ax( k   + j       + 3     , k   + j       + 6     ) + ind3      !      |    !
! Antepenulmate line of each box                                                                                      !      |    !
ax( k   + j       + Nq1-2 , k   + j       + Nq1-5 ) = ax( k   + j       + Nq1-2 , k   + j       + Nq1-5 ) - ind3      !      |    !
ax( k   + j       + Nq1-2 , k   + j       + Nq1-4 ) = ax( k   + j       + Nq1-2 , k   + j       + Nq1-4 ) + ind2      !      |    !
ax( k   + j       + Nq1-2 , k   + j       + Nq1-3 ) = ax( k   + j       + Nq1-2 , k   + j       + Nq1-3 ) - ind1      !      |    !
ax( k   + j       + Nq1-2 , k   + j       + Nq1-1 ) = ax( k   + j       + Nq1-2 , k   + j       + Nq1-1 ) + ind1      !      |    !
ax( k   + j       + Nq1-2 , k   + j       + Nq1   ) = ax( k   + j       + Nq1-2 , k   + j       + Nq1   ) - ind2      !      |    !
! Penulmate line of each box                                                                                          !      |    !
ax( k   + j       + Nq1-1 , k   + j       + Nq1-4 ) = ax( k   + j       + Nq1-1 , k   + j       + Nq1-4 ) - ind3      !      |    !
ax( k   + j       + Nq1-1 , k   + j       + Nq1-3 ) = ax( k   + j       + Nq1-1 , k   + j       + Nq1-3 ) + ind2      !      |    !
ax( k   + j       + Nq1-1 , k   + j       + Nq1-2 ) = ax( k   + j       + Nq1-1 , k   + j       + Nq1-2 ) - ind1      !      |    !
ax( k   + j       + Nq1-1 , k   + j       + Nq1   ) = ax( k   + j       + Nq1-1 , k   + j       + Nq1   ) + ind1      !      |    !
! Last line of each box                                                                                               !      |    !
ax( k   + j       + Nq1   , k   + j       + Nq1-3 ) = ax( k   + j       + Nq1   , k   + j       + Nq1-3 ) - ind3      !      |    !
ax( k   + j       + Nq1   , k   + j       + Nq1-2 ) = ax( k   + j       + Nq1   , k   + j       + Nq1-2 ) + ind2      !      |    !
ax( k   + j       + Nq1   , k   + j       + Nq1-1 ) = ax( k   + j       + Nq1   , k   + j       + Nq1-1 ) - ind1      !      |    !
!.....................................................................................................................!      |    !
!________________________________________________________________________________________________________________________!   |    !
! From the fourth line to the ante-antepenulmate one of each box - non special cases                                     !   |    !
    do i=1+3,Nq1-3 ! do through single box                                                                               !   |    !
ax( k   + j       + i     , k   + j       + i-3   ) = ax( k   + j       + i     , k   + j       + i-3   ) - ind3         !   |    !
ax( k   + j       + i     , k   + j       + i-2   ) = ax( k   + j       + i     , k   + j       + i-2   ) + ind2         !   |    !
ax( k   + j       + i     , k   + j       + i-1   ) = ax( k   + j       + i     , k   + j       + i-1   ) - ind1         !   |    !
ax( k   + j       + i     , k   + j       + i+1   ) = ax( k   + j       + i     , k   + j       + i+1   ) + ind1         !   |    !
ax( k   + j       + i     , k   + j       + i+2   ) = ax( k   + j       + i     , k   + j       + i+2   ) - ind2         !   |    !
ax( k   + j       + i     , k   + j       + i+3   ) = ax( k   + j       + i     , k   + j       + i+3   ) + ind3         !   |    !
    end do ! end of do through single box                                                                                !   |    !
!________________________________________________________________________________________________________________________!   |    !
  end do ! end of do through boxes                                                                                           |    !
end do ! end of do through electronic states                                                                                 |    !
!----------------------------------------------------------------------------------------------------------------------------|    !
!=================================================================================================================================!
allocate (moq1(0:Nst*Nq1*Nq2-1,0:Nst*Nq1*Nq2-1))
!$OMP PARALLEL DO shared(ha,ax) 
do i=1,n
  do j=1,n
    moq1(i-1,j-1)=ax(i,j)
  end do
end do
!$OMP end parallel do
deallocate(ax)
end subroutine first_derivative_matrix_q1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine first_derivative_matrix_q2
use global_param
use omp_lib
implicit none
integer n
n=Nst*Nq1*Nq2
allocate (ax(n,n))
!$OMP parallel do shared(ax)
do i=1,n
  do j=1,n
    ax(i,j)=0.d0
  end do 
end do
!$OMP end parallel do
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
ax( k   +   0*Nq1 + i     , k   +   1*Nq1 + i     ) = ax( k   +   0*Nq1 + i     , k   +   1*Nq1 + i     ) + ind1     !    |       !
ax( k   +   0*Nq1 + i     , k   +   2*Nq1 + i     ) = ax( k   +   0*Nq1 + i     , k   +   2*Nq1 + i     ) - ind2     !    |       !
ax( k   +   0*Nq1 + i     , k   +   3*Nq1 + i     ) = ax( k   +   0*Nq1 + i     , k   +   3*Nq1 + i     ) + ind3     !    |       !
! Upper and left part of the borders - second box of all                                                             !    |       !
ax( k   +   1*Nq1 + i     , k   +   0*Nq1 + i     ) = ax( k   +   1*Nq1 + i     , k   +   0*Nq1 + i     ) - ind1     !    |       !
ax( k   +   1*Nq1 + i     , k   +   2*Nq1 + i     ) = ax( k   +   1*Nq1 + i     , k   +   2*Nq1 + i     ) + ind1     !    |       !
ax( k   +   1*Nq1 + i     , k   +   3*Nq1 + i     ) = ax( k   +   1*Nq1 + i     , k   +   3*Nq1 + i     ) - ind2     !    |       !
ax( k   +   1*Nq1 + i     , k   +   4*Nq1 + i     ) = ax( k   +   1*Nq1 + i     , k   +   4*Nq1 + i     ) + ind3     !    |       !
! Upper and left part of the borders - third box of all                                                              !    |       !
ax( k   +   2*Nq1 + i     , k   +   0*Nq1 + i     ) = ax( k   +   2*Nq1 + i     , k   +   0*Nq1 + i     ) + ind2     !    |       !
ax( k   +   2*Nq1 + i     , k   +   1*Nq1 + i     ) = ax( k   +   2*Nq1 + i     , k   +   1*Nq1 + i     ) - ind1     !    |       !
ax( k   +   2*Nq1 + i     , k   +   3*Nq1 + i     ) = ax( k   +   2*Nq1 + i     , k   +   3*Nq1 + i     ) + ind1     !    |       !
ax( k   +   2*Nq1 + i     , k   +   4*Nq1 + i     ) = ax( k   +   2*Nq1 + i     , k   +   4*Nq1 + i     ) - ind2     !    |       !
ax( k   +   2*Nq1 + i     , k   +   5*Nq1 + i     ) = ax( k   +   2*Nq1 + i     , k   +   5*Nq1 + i     ) + ind3     !    |       !
! Botton and right part of the borders, last box of all                                                              !    |       !
ax( k+s -   0*Nq1 - Nq1+i , k+s -   3*Nq1 - Nq1+i ) = ax( k+s -   0*Nq1 - Nq1+i , k+s -   3*Nq1 - Nq1+i ) - ind3     !    |       !
ax( k+s -   0*Nq1 - Nq1+i , k+s -   2*Nq1 - Nq1+i ) = ax( k+s -   0*Nq1 - Nq1+i , k+s -   2*Nq1 - Nq1+i ) + ind2     !    |       !
ax( k+s -   0*Nq1 - Nq1+i , k+s -   1*Nq1 - Nq1+i ) = ax( k+s -   0*Nq1 - Nq1+i , k+s -   1*Nq1 - Nq1+i ) - ind1     !    |       !
! Botton and right part of the borders, penulmate box of all                                                         !    |       !
ax( k+s -   1*Nq1 - Nq1+i , k+s -   4*Nq1 - Nq1+i ) = ax( k+s -   1*Nq1 - Nq1+i , k+s -   4*Nq1 - Nq1+i ) - ind3     !    |       !
ax( k+s -   1*Nq1 - Nq1+i , k+s -   3*Nq1 - Nq1+i ) = ax( k+s -   1*Nq1 - Nq1+i , k+s -   3*Nq1 - Nq1+i ) + ind2     !    |       !
ax( k+s -   1*Nq1 - Nq1+i , k+s -   2*Nq1 - Nq1+i ) = ax( k+s -   1*Nq1 - Nq1+i , k+s -   2*Nq1 - Nq1+i ) - ind1     !    |       !
ax( k+s -   1*Nq1 - Nq1+i , k+s -   0*Nq1 - Nq1+i ) = ax( k+s -   1*Nq1 - Nq1+i , k+s -   0*Nq1 - Nq1+i ) + ind1     !    |       !
! Botton and right part of the borders, antepenulmate box of all                                                     !    |       !
ax( k+s -   2*Nq1 - Nq1+i , k+s -   5*Nq1 - Nq1+i ) = ax( k+s -   2*Nq1 - Nq1+i , k+s -   5*Nq1 - Nq1+i ) - ind3     !    |       !
ax( k+s -   2*Nq1 - Nq1+i , k+s -   4*Nq1 - Nq1+i ) = ax( k+s -   2*Nq1 - Nq1+i , k+s -   4*Nq1 - Nq1+i ) + ind2     !    |       !
ax( k+s -   2*Nq1 - Nq1+i , k+s -   3*Nq1 - Nq1+i ) = ax( k+s -   2*Nq1 - Nq1+i , k+s -   3*Nq1 - Nq1+i ) - ind1     !    |       !
ax( k+s -   2*Nq1 - Nq1+i , k+s -   1*Nq1 - Nq1+i ) = ax( k+s -   2*Nq1 - Nq1+i , k+s -   1*Nq1 - Nq1+i ) + ind1     !    |       !
ax( k+s -   2*Nq1 - Nq1+i , k+s -   0*Nq1 - Nq1+i ) = ax( k+s -   2*Nq1 - Nq1+i , k+s -   0*Nq1 - Nq1+i ) - ind2     !    |       !
  end do ! end of do through single box                                                                              !    |       !
!....................................................................................................................!    |       !
!_______________________________________________________________________________________________________________________! |       !
! Now doing the rest of the boxes                                                                                       ! |       !
  do j=3*Nq1,(Nq2-4)*Nq1,Nq1 ! Moving through the boxes, from the fourth to the ante-antepenulmate one                  ! |       !
    do i=1,Nq1 ! Moving inside the box                                                                                  ! |       !
ax( k   + j       + i     , k   + j-3*Nq1 + i     ) = ax( k   + j       + i     , k   + j-3*Nq1 + i     ) - ind3        ! |       !
ax( k   + j       + i     , k   + j-2*Nq1 + i     ) = ax( k   + j       + i     , k   + j-2*Nq1 + i     ) + ind2        ! |       !
ax( k   + j       + i     , k   + j-1*Nq1 + i     ) = ax( k   + j       + i     , k   + j-1*Nq1 + i     ) - ind1        ! |       !
ax( k   + j       + i     , k   + j+1*Nq1 + i     ) = ax( k   + j       + i     , k   + j+1*Nq1 + i     ) + ind1        ! |       !
ax( k   + j       + i     , k   + j+2*Nq1 + i     ) = ax( k   + j       + i     , k   + j+2*Nq1 + i     ) - ind2        ! |       !
ax( k   + j       + i     , k   + j+3*Nq1 + i     ) = ax( k   + j       + i     , k   + j+3*Nq1 + i     ) + ind3        ! |       !
    end do ! end of do through single box                                                                               ! |       !
  end do ! end of do through boxes                                                                                      ! |       !
!_______________________________________________________________________________________________________________________! |       !
end do ! end of do through electronic states                                                                              |       !
!-------------------------------------------------------------------------------------------------------------------------|       !
!=================================================================================================================================!
allocate (moq2(0:Nst*Nq1*Nq2-1,0:Nst*Nq1*Nq2-1))
!$OMP PARALLEL DO shared(ha,ax) 
do i=1,n
  do j=1,n
    moq2(i-1,j-1)=ax(i,j)
  end do
end do
!$OMP end parallel do
deallocate(ax)
end subroutine first_derivative_matrix_q2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine will include in the Hamiltonian variable "Ha" the NAC
subroutine nac_ha_modify
use global_param
use omp_lib
integer                               :: cont,n
real(kind=dp), dimension(Nq1,Nq2)     :: nac21q1,nac21q2,nac31q1,nac31q2,nac32q1,nac32q2
real(kind=dp), dimension(0:Nq1*Nq2-1) :: vec1,vec2,vec3,vec4,vec5,vec6

n=Nst*Nq1*Nq2
allocate (ham(0:Nst*Nq1*Nq2-1,0:Nst*Nq1*Nq2-1))
!$OMP parallel do shared(ham)
do i=0,n-1
  do j=0,n-1
    ham(i,j) =0.d0
  end do 
end do
!$OMP end parallel do
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
do j=0,Nq2-1
  do i=0,Nq1-1
    vec1(cont)=nac21q1(i+1,j+1)*ai       ! The 'ai' and 'bi' terms is to transform back to cartesian
    vec2(cont)=nac21q2(i+1,j+1)*bi 
    vec3(cont)=nac31q1(i+1,j+1)*ai
    vec4(cont)=nac31q2(i+1,j+1)*bi
    vec5(cont)=nac32q1(i+1,j+1)*ai
    vec6(cont)=nac32q2(i+1,j+1)*bi
    cont=cont+1
  end do
end do
!$OMP PARALLEL DO shared(moq1,moq2,ham,vec1,vec2,vec3,vec4,vec5,vec6)
do i=0,s-1
  do j=0,s-1
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
do i=0,Nst*Nq1*Nq2-1
  do j=0,Nst*Nq1*Nq2-1
    Ha(i,j)=Ha(i,j)-ham(i,j) ! Minus because of the sign of the kinetic energy term
  end do
end do
!$OMP END PARALLEL DO
deallocate(ham)
write(*,*)'NAC included'
end subroutine nac_ha_modify
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine angular_momentum(n)
use global_param
use omp_lib
implicit none
integer          :: n
real(kind=dp)    :: amxq1,amyq1,amzq1,amxq2,amyq2,amzq2

allocate (ax(0:n-1,0:n-1))
!$OMP PARALLEL DO shared(ax,Ha)
do i=0,n-1
  do j=0,n-1
    ax(i,j)=0.d0
  end do
end do
!$OMP end parallel do

!!$OMP PARALLEL DO shared(moq1,moq2,ax), private(j,k,amxq1,amxq2,amyq1,amyq2,amzq1,amzq2)
do i=0,n-1
  do j=0,n-1
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
!    am(i)=am(i) + (-im) * ( (amxq1+amyq1+amzq1)*moq1(i,j) + (amxq2+amyq2+amzq2)*moq2(i,j) ) * y(j)
    ax(i,j)=ax(i,j) + ( (amxq1+amyq1+amzq1)*moq1(i,j) + (amxq2+amyq2+amzq2)*moq2(i,j) )
  end do
end do
!!$OMP END PARALLEL DO

!Creating the CSR vectors --------------------------------------------------------------------------------------------------!
!This is for the the hamiltonian - dipoles will be in a separate vector                                                     !
k=0                                                                                                                         !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (ax(i,j) /= 0.d0) then                                                                                               !
      k=k+1                                                                                                                 !
    end if                                                                                                                  !
  end do                                                                                                                    !
end do                                                                                                                      !
k_am=k                                                                                                                      !
write(*,*)'k_am = ',k_am                                                                                                    !
allocate(am_val(0:k-1),am_rowc(0:n),am_row_col(0:k-1,0:1))                                                                  !
am_rowc=0.d0                                                                                                                !
k=0                                                                                                                         !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (ax(i,j) /= 0.d0) then                                                                                               !
      am_val(k)=ax(i,j)  !storing each non-zero element                                                                     !
      am_row_col(k,0)=i !storing the row index of each non-zero element                                                     !
      am_row_col(k,1)=j !storing the column index of each non-zero element                                                  !
      k=k+1                                                                                                                 !
    end if                                                                                                                  !
  end do                                                                                                                    !
  am_rowc(i+1)=k !storing the counting of non-zero elements in each row                                                     !
end do                                                                                                                      !
!---------------------------------------------------------------------------------------------------------------------------!

end subroutine angular_momentum
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine hamiltonian_matrix
use global_param

! This subroutine builds the hamiltonian matrrix and stores it in variable Ha. 
! It also creates a temporary copy of Ha in the variable ham, that can be modified to include the time dependent interaction of the dipoles with the pulse.
! So, ham is the real matrix that will be used in the code and both ham and Ha will be global variables.
allocate (ax(Nst*Nq1*Nq2,Nst*Nq1*Nq2))
!$OMP parallel do shared(ax)
do i=1,Nst*Nq1*Nq2
  do j=1,Nst*Nq1*Nq2
    ax(i,j)=0.d0
  end do 
end do
!$OMP end parallel do
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
ax=0.d0 ! Creating the whole matrix and defining every element as zero                                                            !
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
do k=0,2*s,s ! running for all electronic states                                                                             |    !
  do j=0*Nq1,(Nq2-1)*Nq1,Nq1 ! Moving through the boxes, from the fisrt to the last one                                      |    !
!.....................................................................................................................!      |    !
! Doing the special cases of the 5 points finite difference that happens in the borders of each q1 sized box          !      |    !
! Upper and left part of the borders                                                                                  !      |    !
ax( k   + j       + 1     , k   + j       + 1     ) = ax( k   + j       + 1     , k   + j       + 1     ) + ind1      !      |    !
ax( k   + j       + 1     , k   + j       + 2     ) = ax( k   + j       + 1     , k   + j       + 2     ) + ind2      !      |    !
ax( k   + j       + 1     , k   + j       + 3     ) = ax( k   + j       + 1     , k   + j       + 3     ) + ind3      !      |    !
! Upper and left part of the borders, second row                                                                      !      |    !
ax( k   + j       + 2     , k   + j       + 1     ) = ax( k   + j       + 2     , k   + j       + 1     ) + ind2      !      |    !
ax( k   + j       + 2     , k   + j       + 2     ) = ax( k   + j       + 2     , k   + j       + 2     ) + ind1      !      |    !
ax( k   + j       + 2     , k   + j       + 3     ) = ax( k   + j       + 2     , k   + j       + 3     ) + ind2      !      |    !
ax( k   + j       + 2     , k   + j       + 4     ) = ax( k   + j       + 2     , k   + j       + 4     ) + ind3      !      |    !
! Botton and right part of the borders, antepenulmate row                                                             !      |    !
ax( k   + j       + Nq1-1 , k   + j       + Nq1-3 ) = ax( k   + j       + Nq1-1 , k   + j       + Nq1-3 ) + ind3      !      |    !
ax( k   + j       + Nq1-1 , k   + j       + Nq1-2 ) = ax( k   + j       + Nq1-1 , k   + j       + Nq1-2 ) + ind2      !      |    !
ax( k   + j       + Nq1-1 , k   + j       + Nq1-1 ) = ax( k   + j       + Nq1-1 , k   + j       + Nq1-1 ) + ind1      !      |    !
ax( k   + j       + Nq1-1 , k   + j       + Nq1   ) = ax( k   + j       + Nq1-1 , k   + j       + Nq1   ) + ind2      !      |    !
! Botton and right part of the borders, last row                                                                      !      |    !
ax( k   + j       + Nq1   , k   + j       + Nq1-2 ) = ax( k   + j       + Nq1   , k   + j       + Nq1-2 ) + ind3      !      |    !
ax( k   + j       + Nq1   , k   + j       + Nq1-1 ) = ax( k   + j       + Nq1   , k   + j       + Nq1-1 ) + ind2      !      |    !
ax( k   + j       + Nq1   , k   + j       + Nq1   ) = ax( k   + j       + Nq1   , k   + j       + Nq1   ) + ind1      !      |    !
!.....................................................................................................................!      |    !
!________________________________________________________________________________________________________________________!   |    !
! From the third line to the antepenulmate one of each box - non special cases                                           !   |    !
    do i=1+2,Nq1-2 ! do through single box                                                                               !   |    !
ax( k   + j       + i     , k   + j       + i-2   ) = ax( k   + j       + i     , k   + j       + i-2   ) + ind3         !   |    !
ax( k   + j       + i     , k   + j       + i-1   ) = ax( k   + j       + i     , k   + j       + i-1   ) + ind2         !   |    !
ax( k   + j       + i     , k   + j       + i     ) = ax( k   + j       + i     , k   + j       + i     ) + ind1         !   |    !
ax( k   + j       + i     , k   + j       + i+1   ) = ax( k   + j       + i     , k   + j       + i+1   ) + ind2         !   |    !
ax( k   + j       + i     , k   + j       + i+2   ) = ax( k   + j       + i     , k   + j       + i+2   ) + ind3         !   |    !
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
!ax( k   + j       + i     , k   + j       + i     ) = ax( k   + j       + i     , k   + j       + i     ) + ind1     |           !
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
do k=0,2*s,s ! Running through every state, the step size is Nq1*Nq2                                                  |           !
!................................................................................................................!    |           !
! Doing the special cases of the 5 points finite difference that happens in first and last two boxes             !    |           !
  do i=1,Nq1 ! Moving inside the box                                                                             !    |           !
! Upper and left part of the borders - first box of all                                                          !    |           !
ax( k   +   0*Nq1 + i     , k   +   0*Nq1 + i     ) = ax( k   +   0*Nq1 + i     , k   +   0*Nq1 + i     ) + ind1 !    |           !
ax( k   +   0*Nq1 + i     , k   +   1*Nq1 + i     ) = ax( k   +   0*Nq1 + i     , k   +   1*Nq1 + i     ) + ind2 !    |           !
ax( k   +   0*Nq1 + i     , k   +   2*Nq1 + i     ) = ax( k   +   0*Nq1 + i     , k   +   2*Nq1 + i     ) + ind3 !    |           !
! Botton and right part of the borders, last box of all                                                          !    |           !
ax( k+s -   0*Nq1 - Nq1+i , k+s -   2*Nq1 - Nq1+i ) = ax( k+s -   0*Nq1 - Nq1+i , k+s -   2*Nq1 - Nq1+i ) + ind3 !    |           !
ax( k+s -   0*Nq1 - Nq1+i , k+s -   1*Nq1 - Nq1+i ) = ax( k+s -   0*Nq1 - Nq1+i , k+s -   1*Nq1 - Nq1+i ) + ind2 !    |           !
ax( k+s -   0*Nq1 - Nq1+i , k+s -   0*Nq1 - Nq1+i ) = ax( k+s -   0*Nq1 - Nq1+i , k+s -   0*Nq1 - Nq1+i ) + ind1 !    |           !
! Upper and left part of the borders - second box of all                                                         !    |           !
ax( k   +   1*Nq1 + i     , k   +   0*Nq1 + i     ) = ax( k   +   1*Nq1 + i     , k   +   0*Nq1 + i     ) + ind2 !    |           !
ax( k   +   1*Nq1 + i     , k   +   1*Nq1 + i     ) = ax( k   +   1*Nq1 + i     , k   +   1*Nq1 + i     ) + ind1 !    |           !
ax( k   +   1*Nq1 + i     , k   +   2*Nq1 + i     ) = ax( k   +   1*Nq1 + i     , k   +   2*Nq1 + i     ) + ind2 !    |           !
ax( k   +   1*Nq1 + i     , k   +   3*Nq1 + i     ) = ax( k   +   1*Nq1 + i     , k   +   3*Nq1 + i     ) + ind3 !    |           !
! Botton and right part of the borders, penulmate box of all                                                     !    |           !
ax( k+s -   1*Nq1 - Nq1+i , k+s -   3*Nq1 - Nq1+i ) = ax( k+s -   1*Nq1 - Nq1+i , k+s -   3*Nq1 - Nq1+i ) + ind3 !    |           !
ax( k+s -   1*Nq1 - Nq1+i , k+s -   2*Nq1 - Nq1+i ) = ax( k+s -   1*Nq1 - Nq1+i , k+s -   2*Nq1 - Nq1+i ) + ind2 !    |           !
ax( k+s -   1*Nq1 - Nq1+i , k+s -   1*Nq1 - Nq1+i ) = ax( k+s -   1*Nq1 - Nq1+i , k+s -   1*Nq1 - Nq1+i ) + ind1 !    |           !
ax( k+s -   1*Nq1 - Nq1+i , k+s -   0*Nq1 - Nq1+i ) = ax( k+s -   1*Nq1 - Nq1+i , k+s -   0*Nq1 - Nq1+i ) + ind2 !    |           !
  end do ! end of do through single box                                                                          !    |           !
!................................................................................................................!    |           !
!___________________________________________________________________________________________________________________! |           !
! Now doing the rest of the boxes                                                                                   ! |           !
  do j=2*Nq1,(Nq2-3)*Nq1,Nq1 ! Moving through the boxes, from the third to the antipenulmate                        ! |           !
    do i=1,Nq1 ! Moving inside the box                                                                              ! |           !
ax( k   + j       + i     , k   + j-2*Nq1 + i     ) = ax( k   + j       + i     , k   + j-2*Nq1 + i     ) + ind3    ! |           !
ax( k   + j       + i     , k   + j-1*Nq1 + i     ) = ax( k   + j       + i     , k   + j-1*Nq1 + i     ) + ind2    ! |           !
ax( k   + j       + i     , k   + j       + i     ) = ax( k   + j       + i     , k   + j       + i     ) + ind1    ! |           !
ax( k   + j       + i     , k   + j+1*Nq1 + i     ) = ax( k   + j       + i     , k   + j+1*Nq1 + i     ) + ind2    ! |           !
ax( k   + j       + i     , k   + j+2*Nq1 + i     ) = ax( k   + j       + i     , k   + j+2*Nq1 + i     ) + ind3    ! |           !
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
do k=0,2*s,s ! running for all electronic states                                                                          |       !
!..................................................................................................................!      |       !
! Doing the special cases that happens in the borders of each q1 sized box when moving along q1                    !      |       !
! Upper and left part of the borders - first row of each box                                                       !      |       !
  do j=0*Nq1,(Nq2-1)*Nq1,Nq1 ! Moving through the boxes, from the fisrt to the last one                            !      |       !
ax( k   + j       + 1     , k   + j       + 2     ) = ax( k   + j       + 1     , k   + j       + 2     ) + ind2   !      |       !
! Botton and right part of the borders, last row of each box                                                       !      |       !
ax( k   + j       + Nq1   , k   + j       + Nq1-1 ) = ax( k   + j       + Nq1   , k   + j       + Nq1-1 ) + ind2   !      |       !
  end do ! end of do through boxes                                                                                 !      |       !
! Now doing the rest of the lines of the first and last box                                                        !      |       ! 
  do i=2,Nq1-1 ! Moving inside the first and last box, but avoiding the first and last row                         !      |       !
ax( k   + 0       + i     , k   + 0       + i-1   ) = ax( k   + 0       + i     , k   + 0       + i-1   ) + ind2   !      |       !
ax( k   + 0       + i     , k   + 0       + i+1   ) = ax( k   + 0       + i     , k   + 0       + i+1   ) + ind2   !      |       !
ax( k+s -Nq1      + i     , k+s -Nq1      + i-1   ) = ax( k+s -Nq1      + i     , k+s -Nq1      + i-1   ) + ind2   !      |       !
ax( k+s -Nq1      + i     , k+s -Nq1      + i+1   ) = ax( k+s -Nq1      + i     , k+s -Nq1      + i+1   ) + ind2   !      |       !
  end do ! end of do through single box                                                                            !      |       !
!..................................................................................................................!      |       !
! Doing the special cases that happens on the first and last box for each state when moving along q2               !      |       !
! Upper and left part of the borders - first box                                                                   !      |       !
  do i=1,Nq1 ! Moving inside the box                                                                               !      |       !
ax( k   + 0       + i     , k   + Nq1     + i     ) = ax( k   + 0       + i     , k   + Nq1     + i     ) + ind2   !      |       !
! Botton and right part of the borders, last box                                                                   !      |       !
ax( k+s - Nq1     + i     , k+s - 2*Nq1   + i     ) = ax( k+s - Nq1     + i     , k+s - 2*Nq1   + i     ) + ind2   !      |       !
  end do ! end of do through box                                                                                   !      |       !
! Now doing the first and last row of each of the others boxes                                                     !      |       !
  do j=1*Nq1,(Nq2-2)*Nq1,Nq1 ! Moving through the boxes, from the second to penulmate one                          !      |       !
ax( k   + j       + 1     , k   + j-Nq1   + 1     ) = ax( k   + j       + 1     , k   + j-Nq1   + 1     ) + ind2   !      |       !
ax( k   + j       + 1     , k   + j+Nq1   + 1     ) = ax( k   + j       + 1     , k   + j+Nq1   + 1     ) + ind2   !      |       !
ax( k   + j       + Nq1   , k   + j-Nq1   + Nq1   ) = ax( k   + j       + Nq1   , k   + j-Nq1   + Nq1   ) + ind2   !      |       !
ax( k   + j       + Nq1   , k   + j+Nq1   + Nq1   ) = ax( k   + j       + Nq1   , k   + j+Nq1   + Nq1   ) + ind2   !      |       !
  end do ! end of do through boxes                                                                                 !      |       !
!..................................................................................................................!      |       !
! Now the special cases when moving along q1 and q2 at the same time                                               !      |       !
! Cross terms in the first box and first row                                                                       !      |       !
ax( k   + 0       + 1     , k   + Nq1     + 2     ) = ax( k   + 0       + 1     , k   + Nq1     + 2     ) + ind3   !      |       !
! Cross terms in the last box and last row                                                                         !      |       !
ax( k+s + 0       + 0     , k+s - Nq1     - 1     ) = ax( k+s + 0       + 0     , k+s - Nq1     - 1     ) + ind3   !      |       !
! Now for the rest of the rows of the first and last box                                                           !      |       !
  do i=1,Nq1-2 ! Moving inside the first and last box, but avoiding the first and last row                         !      |       !
ax( k   + 0       + i+1   , k   + Nq1     + i+2   ) = ax( k   + 0       + i+1   , k   + Nq1     + i+2   ) + ind3   !      |       !
ax( k+s + 0       - i     , k+s - 1*Nq1   - (i+1) ) = ax( k+s + 0       - i     , k+s - 1*Nq1   - (i+1) ) + ind3   !      |       !
  end do ! end of do through single box                                                                            !      |       !
! Now for the first and last rows of each of the others boxes                                                      !      |       !
  do j=1*Nq1,(Nq2-2)*Nq1,Nq1 ! Moving through the boxes, from the second to penulmate one                          !      |       !
ax( k   + j       + 1     , k   + j + Nq1 + 1+1   ) = ax( k   + j       + 1     , k   + j + Nq1 + 1+1   ) + ind3   !      |       !
ax( k   + j       + Nq1   , k   + j - Nq1 + Nq1-1 ) = ax( k   + j       + Nq1   , k   + j - Nq1 + Nq1-1 ) + ind3   !      |       !
  end do ! end of do through boxes                                                                                 !      |       !
!..................................................................................................................!      |       !
! Finally, doing for the rest of the boxes, always avoinding the first and last row                                !      |       !
  do j=1*Nq1,(Nq2-2)*Nq1,Nq1 ! Moving through the boxes, from the second to penulmate one                          !      |       !
    do i=2,Nq1-1 ! Moving inside the first and last box, but avoiding the first and last row                       !      |       !
ax( k   + j       + i     , k   + j + Nq1 + i+1   ) = ax( k   + j       + i     , k   + j + Nq1 + i+1   ) + ind3   !      |       !
ax( k   + j       + i     , k   + j       + i+1   ) = ax( k   + j       + i     , k   + j       + i+1   ) + ind2   !      |       !
ax( k   + j       + i     , k   + j       + i-1   ) = ax( k   + j       + i     , k   + j       + i-1   ) + ind2   !      |       !
ax( k   + j       + i     , k   + j + Nq1 + i     ) = ax( k   + j       + i     , k   + j + Nq1 + i     ) + ind2   !      |       !
ax( k   + j       + i     , k   + j - Nq1 + i     ) = ax( k   + j       + i     , k   + j - Nq1 + i     ) + ind2   !      |       !
ax( k   + j       + i     , k   + j - Nq1 + i-1   ) = ax( k   + j       + i     , k   + j - Nq1 + i-1   ) + ind3   !      |       !
    end do ! end of do through single box                                                                          !      |       !
  end do ! end of do through boxes                                                                                 !      |       !
!..................................................................................................................!      |       !
! Now doing only the diagonal of each box                                                                          !      |       !
  do j=0*Nq1,(Nq2-1)*Nq1,Nq1 ! Moving through the boxes, from the fisrt to the last one                            !      |       !
    do i=1,Nq1 ! Moving inside the box                                                                             !      |       !
ax( k   + j       + i     , k   + j       + i     ) = ax( k   + j       + i     , k   + j       + i     ) + ind1   !      |       !
    end do ! end of do through single box                                                                          !      |       !
  end do ! end of do through boxes                                                                                 !      |       !
end do ! end of do through electronic states                                                                       !      |       !
!-------------------------------------------------------------------------------------------------------------------------|       !
! Finishing the kinetic energy part                                                                                               !
ax(: , :)= - 1.d0/2.d0 * ax(: , :)                                                                                                !
! Adding potential energy - moving only along the diagonal of the axmiltonian                                                     !
do j=0*Nq1,(Nq2-1)*Nq1,Nq1 ! Moving through the boxes, from the fisrt to the last one                                             !
  do i=1,Nq1 ! Moving inside the box                                                                                              !
    ax( 0*s + j + i , 0*s + j + i ) = ax( 0*s + j + i , 0*s + j + i ) + pot1 (i,(j+Nq1)/Nq1)                                      !
    ax( 1*s + j + i , 1*s + j + i ) = ax( 1*s + j + i , 1*s + j + i ) + pot2 (i,(j+Nq1)/Nq1)                                      !
    ax( 2*s + j + i , 2*s + j + i ) = ax( 2*s + j + i , 2*s + j + i ) + pot3 (i,(j+Nq1)/Nq1)                                      !
  end do ! end of do through single box                                                                                           !
end do ! end of do through boxes                                                                                                  !
!=================================================================================================================================!

allocate (Ha(0:Nst*Nq1*Nq2-1,0:Nst*Nq1*Nq2-1))

!$OMP PARALLEL DO shared(ha,ax) 
do i=1,Nst*Nq1*Nq2
  do j=1,Nst*Nq1*Nq2
    Ha(i-1,j-1)=ax(i,j)
  end do
end do
!$OMP end parallel do

deallocate(ax)

end subroutine hamiltonian_matrix









