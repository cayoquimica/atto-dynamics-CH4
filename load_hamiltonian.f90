program load_hamiltonian
!maybe is useful: Function MERGE:
!test a condition, if true, assign a value X, if false assign value Y:
!the condition can be a boolean variable defined before. 
!result = merge( X, Y, i.eq.j)
!result, X and Y will be of same type
use global_param
use omp_lib
implicit none
integer                      :: n,jj,ii,nana !m=dimension of the Hamiltonian (Nq1*Nq2*Nst x Nq1*Nq2*Nst)
real(kind=dp)                :: ompt0,ompt1,x1,x2,soma
integer(kind=dp)             :: lolx(0:Nq1*Nq2*Nst-1),loly(0:Nq1*Nq2*Nst-1),lolz(0:Nq1*Nq2*Nst-1)
integer(kind=dp)             :: kx,ky,kz
!integer,parameter            :: q1_initial = 25
!integer,parameter            :: q1_final = 70
!integer,parameter            :: q2_initial = 75
!integer,parameter            :: q2_final = 115
!complex(kind=dp),allocatable :: pia(:),nwf0(:),pice(:),wf0(:)
integer                      :: init_wf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! I modified the Ha, moq1 and moq2 matrices building so that their indexes begin with 0 instead of 1. 
! It is not the ideal way, I use a auxiliar matrix to build with index 1 as I already had written in the code and then I put it
! in a matrix with index starting in 0.
! Later would be better and faster if I already build the matrices with index 0.     <<<---------------

n=Nst*Nq1*Nq2 ! m=Nq1*Nq2*st %size of the Hamiltonian

!$ ompt0 = omp_get_wtime()
call load_data

!!$ ompt0 = omp_get_wtime()
!allocate ( nwf0(0:Nst*Nq1*Nq2-1) )
!nwf0=dcmplx(0.d0,0.d0)
!open(newunit=init_wf,file='eigen-vec-neutral',status='unknown')
!do i=0,s-1
!  read(init_wf,*) nwf0(i) !reading the neutral initial eigen state and writting in the first third part of the vector
!end do
!!Now read the photoinization coefficients and project them into the initial neutral ground state wave function
!allocate ( pice(0:Nst*Nq1*Nq2-1) )
!pice=dcmplx(0.d0,0.d0)
!allocate(pia(Nst))
!k = (q2_initial-1+20) * Nq1 + (q1_initial+27) - 1 !The -1 is because the vector starts at index 0
!!!$OMP parallel do default(private) shared(pice)
!do j=q2_initial+20,q2_final+20  !The photoionization coefficiets were calculated only for q2=75+20:115+20 and q1=25+27:70+27
!  do i=q1_initial+27,q1_final+27 !where the amplitudes of the eigen state of the neutral is non-zero (< 10^-8). This is the Frank-Condon region
!    call p_i_a(i-27,j-20,pia) !Evaluating the photoionization coeficients for all electronic states
!    pice(k)     = pia(1)
!    pice(k+s)   = pia(2)
!    pice(k+2*s) = pia(3)
!    k = k+1
!  end do
!  k = k + (Nq1-(q1_final+27)) + (q1_initial+27) - 1 !Add the zeros values from 70+27 until Nq1 and from 1 to 25+27. The -1 here is to anulate the last -1 of the previous loop
!end do
!!!$OMP end parallel do
!allocate ( wf0(0:Nst*Nq1*Nq2-1) )
!!wf0=dcmplx(0.d0,0.d0)
!!Projecting the photoionization coeficients into the neutral eigen state
!do i=0,s-1
!  wf0(i)     = nwf0(i) * pice(i)
!  wf0(i+s)   = nwf0(i) * pice(i+s)
!  wf0(i+2*s) = nwf0(i) * pice(i+2*s)
!end do
!!--------------------------------------------!
!!normalizing                                 !
!soma=0.d0                                    !
!do i=0,n-1                                   !
!  soma=soma+dconjg(wf0(i))*wf0(i)            !
!end do                                       !
!wf0=wf0/sqrt(soma)                           !
!!--------------------------------------------!
!open(unit=88,file='wfINIT',status='unknown')
!do i=0,n-1
!  write(88,*)'(',dreal(wf0(i)),',',dimag(wf0(i)),')'
!end do
!close(unit=88)
!!$ ompt1 = omp_get_wtime()
!write(*,*)'The time to prepare the initial wavepacket=',ompt1 - ompt0, "seconds"

call hamiltonian_matrix
write(*,*)'Kinetic and potential energy loaded. Their CSR vectors are created.'
call first_derivative_matrix_q1
!$ x1 = omp_get_wtime()
!Creating the CSR vectors --------------------------------------------------------------------------------------------------!
!This is for the first derivative with respect of q1                                                                        !
lol=0.d0                                                                                                                    !
k=0                                                                                                                         !
!$OMP parallel do shared(lol)                                                                                               !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (moq1(i,j) /= 0.d0) then                                                                                             !
      lol(i) = lol(i) + 1                                                                                                   !
    end if                                                                                                                  !
  end do                                                                                                                    !
end do                                                                                                                      !
!$OMP end parallel do                                                                                                       !
k = sum(lol)                                                                                                                !
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
lol=0.d0                                                                                                                    !
k=0                                                                                                                         !
!$OMP parallel do shared(lol)                                                                                               !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (moq2(i,j) /= 0.d0) then                                                                                             !
      lol(i) = lol(i) + 1                                                                                                   !
    end if                                                                                                                  !
  end do                                                                                                                    !
end do                                                                                                                      !
!$OMP end parallel do                                                                                                       !
k = sum(lol)                                                                                                                !
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
deallocate(moq1,moq2)

!I still calculate angular momentum without using the sparse matrix CSR method.                                                                   !
!Deallocate these matrices to salve RAM at this point. They will be reallocated later again just to be used in the angular momentum subroutine    !
!Not perfect way to do it, but these matrices will be loaded two times now in order to not use to much RAM memory.                                !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!==============================================!
! Loading dipoles                              !
allocate (ay(0:n-1,0:n-1,0:2))                 !
!$OMP PARALLEL DO shared(ay)                   !
do i=0,n-1                                     !
  do j=0,n-1                                   !
    ay(i,j,:)=0.d0                             !
  end do                                       !
end do                                         !
!$OMP end parallel do                          !
!$OMP PARALLEL DO shared(ay)                   !
do i=0,Nq1-1                                   !
  do j=0,Nq2-1                                 !
    ii = (i*(Nq2-1)+j)
    ay(ii    ,ii    ,0) = - pdm1x(i+1,j+1)     !
    ay(s+ii  ,s+ii  ,0) = - pdm2x(i+1,j+1)     !
    ay(2*s+ii,2*s+ii,0) = - pdm3x(i+1,j+1)     !
!
    ay(ii    ,s+ii  ,0) = - tdm21x(i+1,j+1)    !
    ay(s+ii  ,ii    ,0) = - tdm21x(i+1,j+1)    !
    ay(ii    ,2*s+ii,0) = - tdm31x(i+1,j+1)    !
    ay(2*s+ii,ii    ,0) = - tdm31x(i+1,j+1)    !
    ay(s+ii  ,2*s+ii,0) = - tdm32x(i+1,j+1)    !
    ay(2*s+ii,s+ii  ,0) = - tdm32x(i+1,j+1)    !


    ay(ii    ,ii    ,1) = - pdm1y(i+1,j+1)     !
    ay(s+ii  ,s+ii  ,1) = - pdm2y(i+1,j+1)     !
    ay(2*s+ii,2*s+ii,1) = - pdm3y(i+1,j+1)     !
!
    ay(ii    ,s+ii  ,1) = - tdm21y(i+1,j+1)    !
    ay(s+ii  ,ii    ,1) = - tdm21y(i+1,j+1)    !
    ay(ii    ,2*s+ii,1) = - tdm31y(i+1,j+1)    !
    ay(2*s+ii,ii    ,1) = - tdm31y(i+1,j+1)    !
    ay(s+ii  ,2*s+ii,1) = - tdm32y(i+1,j+1)    !
    ay(2*s+ii,s+ii  ,1) = - tdm32y(i+1,j+1)    !


    ay(ii    ,ii    ,2) = - pdm1z(i+1,j+1)     !
    ay(s+ii  ,s+ii  ,2) = - pdm2z(i+1,j+1)     !
    ay(2*s+ii,2*s+ii,2) = - pdm3z(i+1,j+1)     !
!
    ay(ii    ,s+ii  ,2) = - tdm21z(i+1,j+1)    !
    ay(s+ii  ,ii    ,2) = - tdm21z(i+1,j+1)    !
    ay(ii    ,2*s+ii,2) = - tdm31z(i+1,j+1)    !
    ay(2*s+ii,ii    ,2) = - tdm31z(i+1,j+1)    !
    ay(s+ii  ,2*s+ii,2) = - tdm32z(i+1,j+1)    !
    ay(2*s+ii,s+ii  ,2) = - tdm32z(i+1,j+1)    !
  end do                                       !
end do                                         !
!$OMP end parallel do                          !
!==============================================!
allocate (ham(0:n-1,0:n-1))
!$OMP PARALLEL DO shared(ay,ham,Ha)
do i=0,n-1
  do j=0,n-1
    ham(i,j)=Ha(i,j)+ay(i,j,0)+ay(i,j,1)+ay(i,j,2)
  end do
end do
!$OMP end parallel do
!Creating the CSR vectors --------------------------------------------------------------------------------------------------!
lolx=0.d0                                                                                                                   !
loly=0.d0                                                                                                                   !
lolz=0.d0                                                                                                                   !
!$OMP parallel do shared(lol)                                                                                               !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (ay(i,j,0) /= 0.d0) then                                                                                             !
      lolx(i) = lolx(i) + 1                                                                                                 !
    end if                                                                                                                  !
    if (ay(i,j,1) /= 0.d0) then                                                                                             !
      loly(i) = loly(i) + 1                                                                                                 !
    end if                                                                                                                  !
    if (ay(i,j,2) /= 0.d0) then                                                                                             !
      lolz(i) = lolz(i) + 1                                                                                                 !
    end if                                                                                                                  !
  end do                                                                                                                    !
end do                                                                                                                      !
!$OMP end parallel do                                                                                                       !
kx = sum(lolx)                                                                                                              !
ky = sum(loly)                                                                                                              !
kz = sum(lolz)                                                                                                              !
k_dipx=kx                                                                                                                   !
k_dipy=ky                                                                                                                   !
k_dipz=kz                                                                                                                   !
allocate(dipx_val(0:kx-1),dipx_rowc(0:n),dipx_row_col(0:kx-1,0:1))                                                          !
allocate(dipy_val(0:ky-1),dipy_rowc(0:n),dipy_row_col(0:ky-1,0:1))                                                          !
allocate(dipz_val(0:kz-1),dipz_rowc(0:n),dipz_row_col(0:kz-1,0:1))                                                          !
!                                                                                                                           !
dipx_rowc=0.d0                                                                                                              !
k=0                                                                                                                         !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (ay(i,j,0) /= 0.d0) then                                                                                             !
      dipx_val(k)=ay(i,j,0)  !storing each non-zero element                                                                 !
      dipx_row_col(k,0)=i !storing the row index of each non-zero element                                                   !
      dipx_row_col(k,1)=j !storing the column index of each non-zero element                                                !
      k=k+1                                                                                                                 !
    end if                                                                                                                  !
  end do                                                                                                                    !
  dipx_rowc(i+1)=k !storing the counting of non-zero elements in each row                                                   !
end do                                                                                                                      !
dipy_rowc=0.d0                                                                                                              !
k=0                                                                                                                         !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (ay(i,j,1) /= 0.d0) then                                                                                             !
      dipy_val(k)=ay(i,j,1)  !storing each non-zero element                                                                 !
      dipy_row_col(k,0)=i !storing the row index of each non-zero element                                                   !
      dipy_row_col(k,1)=j !storing the column index of each non-zero element                                                !
      k=k+1                                                                                                                 !
    end if                                                                                                                  !
  end do                                                                                                                    !
  dipy_rowc(i+1)=k !storing the counting of non-zero elements in each row                                                   !
end do                                                                                                                      !
dipz_rowc=0.d0                                                                                                              !
k=0                                                                                                                         !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (ay(i,j,2) /= 0.d0) then                                                                                             !
      dipz_val(k)=ay(i,j,2)  !storing each non-zero element                                                                 !
      dipz_row_col(k,0)=i !storing the row index of each non-zero element                                                   !
      dipz_row_col(k,1)=j !storing the column index of each non-zero element                                                !
      k=k+1                                                                                                                 !
    end if                                                                                                                  !
  end do                                                                                                                    !
  dipz_rowc(i+1)=k !storing the counting of non-zero elements in each row                                                   !
end do                                                                                                                      !
!---------------------------------------------------------------------------------------------------------------------------!
!Creating the CSR vectors --------------------------------------------------------------------------------------------------!
!This is for the the hamiltonian - dipoles will be in a separate vector                                                     !
lol=0.d0                                                                                                                    !
k=0                                                                                                                         !
!$OMP parallel do shared(lol)                                                                                               !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (ham(i,j) /= 0.d0) then                                                                                              !
      lol(i) = lol(i) + 1                                                                                                   !
    end if                                                                                                                  !
  end do                                                                                                                    !
end do                                                                                                                      !
!$OMP end parallel do                                                                                                       !
k = sum(lol)                                                                                                                !
k_Ha=k                                                                                                                      !
k_di2=k                                                                                                                     !
write(*,*)'k_Ha = ',k_Ha                                                                                                    !
allocate(Ha_val(0:k-1),Ha_rowc(0:n),Ha_row_col(0:k-1,0:1),di2_val(0:k-1,0:2))                                               !
!                                                                                                                           !
Ha_rowc=0.d0                                                                                                                !
k=0                                                                                                                         !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (ham(i,j) /= 0.d0) then                                                                                              !
      Ha_val(k)=Ha(i,j)  !storing each non-zero element                                                                     !
      Ha_row_col(k,0)=i !storing the row index of each non-zero element                                                     !
      Ha_row_col(k,1)=j !storing the column index of each non-zero element                                                  !
      di2_val(k,0)=ay(i,j,0)  !storing each non-zero element                                                                !
      di2_val(k,1)=ay(i,j,1)  !storing each non-zero element                                                                !
      di2_val(k,2)=ay(i,j,2)  !storing each non-zero element                                                                !
      k=k+1                                                                                                                 !
    end if                                                                                                                  !
  end do                                                                                                                    !
  Ha_rowc(i+1)=k !storing the counting of non-zero elements in each row                                                     ! 
end do                                                                                                                      !
!---------------------------------------------------------------------------------------------------------------------------!
deallocate(Ha,ham,ay)
call first_derivative_matrix_q1
call first_derivative_matrix_q2
call angular_momentum(n)
!$ ompt1 = omp_get_wtime()
write(*,*) 'The time to load the matrices=',ompt1 - ompt0, "seconds"

write(*,*) 'Initiating saving dq1'

open(unit=20,file='csr_vectors_dq1',status='unknown')
write(20,'(i12)')k_moq1
do i=0,k_moq1-1
  write(20,'(e23.15e3,2i12)')moq1_val(i),moq1_row_col(i,0),moq1_row_col(i,1)
end do
do i=0,n
  write(20,'(i12)')moq1_rowc(i)
end do
close(unit=20)

write(*,*) 'Initiating saving dq2'

open(unit=20,file='csr_vectors_dq2',status='unknown')
write(20,'(i12)')k_moq2
do i=0,k_moq2-1
  write(20,'(e23.15e3,2i12)')moq2_val(i),moq2_row_col(i,0),moq2_row_col(i,1)
end do
do i=0,n
  write(20,'(i12)')moq2_rowc(i)
end do
close(unit=20)

write(*,*) 'Initiating saving kin'

open(unit=20,file='csr_vectors_kin',status='unknown')
write(20,'(i12)')k_kine
do i=0,k_kine-1
  write(20,'(e23.15e3,2i12)')kine_val(i),kine_row_col(i,0),kine_row_col(i,1)
end do

do i=0,n
  write(20,'(i12)')kine_rowc(i)
end do
close(unit=20)

write(*,*) 'Initiating saving pot'

open(unit=20,file='csr_vectors_pot',status='unknown')
write(20,'(i12)')k_pot
do i=0,k_pot-1
  write(20,'(e23.15e3,2i12)')pot_val(i),pot_row_col(i,0),pot_row_col(i,1)
end do
do i=0,n
  write(20,'(i12)')pot_rowc(i)
end do
close(unit=20)

write(*,*) 'Initiating saving dip'

open(unit=20,file='csr_vectors_dip',status='unknown')
write(20,'(i12)')k_dipx
do i=0,k_dipx-1
  write(20,'(e23.15e3,2i12)')dipx_val(i),dipx_row_col(i,0),dipx_row_col(i,1)
end do
do i=0,n
  write(20,'(i12)')dipx_rowc(i)
end do
write(20,'(i12)')k_dipy
do i=0,k_dipy-1
  write(20,'(e23.15e3,2i12)')dipy_val(i),dipy_row_col(i,0),dipy_row_col(i,1)
end do
do i=0,n
  write(20,'(i12)')dipy_rowc(i)
end do
write(20,'(i12)')k_dipz
do i=0,k_dipz-1
  write(20,'(e23.15e3,2i12)')dipz_val(i),dipz_row_col(i,0),dipz_row_col(i,1)
end do
do i=0,n
  write(20,'(i12)')dipz_rowc(i)
end do
close(unit=20)

write(*,*) 'Initiating saving nac'

open(unit=20,file='csr_vectors_NAC',status='unknown')
write(20,'(i12)')k_nac
do i=0,k_nac-1
  write(20,'(e23.15e3,2i12)')nac_val(i),nac_row_col(i,0),nac_row_col(i,1)
end do
do i=0,n
  write(20,'(i12)')nac_rowc(i)
end do
close(unit=20)

write(*,*) 'Initiating saving anm'

open(unit=20,file='csr_vectors_anm',status='unknown')
write(20,'(i12)')k_am
do i=0,k_am-1
  write(20,'(e23.15e3,2i12)')am_val(i),am_row_col(i,0),am_row_col(i,1)
end do
do i=0,n
  write(20,'(i12)')am_rowc(i)
end do
close(unit=20)


write(*,*) 'Initiating saving normal csr vectors'

!DOING THE ORIGINAL WITH THE WHOLE HAMILTONIAN, THE DIPOLES APART, FIRST DERIVATIVES AND ANGULAR MOMENTUM
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

write(20,'(i12)')k_di2
do i=0,k_di2-1
  write(20,'(3e23.15e3)')di2_val(i,:)
end do

write(20,'(i12)')k_am
do i=0,k_am-1
  write(20,'(e23.15e3,2i12)')am_val(i),am_row_col(i,0),am_row_col(i,1)
end do
do i=0,n
  write(20,'(i12)')am_rowc(i)
end do

close(unit=20)

write(*,*)'FINISHED'


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
real(kind=dp), dimension(Nq1,Nq2)     :: nac10c0x, nac10h1x, nac10h2x, nac10h3x, nac10h4x 
real(kind=dp), dimension(Nq1,Nq2)     :: nac10c0y, nac10h1y, nac10h2y, nac10h3y, nac10h4y 
real(kind=dp), dimension(Nq1,Nq2)     :: nac10c0z, nac10h1z, nac10h2z, nac10h3z, nac10h4z 
real(kind=dp), dimension(Nq1,Nq2)     :: nac20c0x, nac20h1x, nac20h2x, nac20h3x, nac20h4x 
real(kind=dp), dimension(Nq1,Nq2)     :: nac20c0y, nac20h1y, nac20h2y, nac20h3y, nac20h4y 
real(kind=dp), dimension(Nq1,Nq2)     :: nac20c0z, nac20h1z, nac20h2z, nac20h3z, nac20h4z 
real(kind=dp), dimension(Nq1,Nq2)     :: nac21c0x, nac21h1x, nac21h2x, nac21h3x, nac21h4x 
real(kind=dp), dimension(Nq1,Nq2)     :: nac21c0y, nac21h1y, nac21h2y, nac21h3y, nac21h4y 
real(kind=dp), dimension(Nq1,Nq2)     :: nac21c0z, nac21h1z, nac21h2z, nac21h3z, nac21h4z 
real(kind=dp), dimension(0:Nq1*Nq2-1) :: vec1,vec2,vec3,vec4,vec5,vec6

real(kind=dp) :: truni,trunf

!$ truni = omp_get_wtime()


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
open(unit=21,file='nac10c0xf.txt',status='old')                      !
open(unit=22,file='nac10c0yf.txt',status='old')                      !
open(unit=23,file='nac10c0zf.txt',status='old')                      !
open(unit=24,file='nac10h1xf.txt',status='old')                      !
open(unit=25,file='nac10h1yf.txt',status='old')                      !
open(unit=26,file='nac10h1zf.txt',status='old')                      !
open(unit=27,file='nac10h2xf.txt',status='old')                      !
open(unit=28,file='nac10h2yf.txt',status='old')                      !
open(unit=29,file='nac10h2zf.txt',status='old')                      !
open(unit=30,file='nac10h3xf.txt',status='old')                      !
open(unit=31,file='nac10h3yf.txt',status='old')                      !
open(unit=32,file='nac10h3zf.txt',status='old')                      !
open(unit=33,file='nac10h4xf.txt',status='old')                      !
open(unit=34,file='nac10h4yf.txt',status='old')                      !
open(unit=35,file='nac10h4zf.txt',status='old')                      !
open(unit=36,file='nac20c0xf.txt',status='old')                      !
open(unit=37,file='nac20c0yf.txt',status='old')                      !
open(unit=38,file='nac20c0zf.txt',status='old')                      !
open(unit=39,file='nac20h1xf.txt',status='old')                      !
open(unit=40,file='nac20h1yf.txt',status='old')                      !
open(unit=41,file='nac20h1zf.txt',status='old')                      !
open(unit=42,file='nac20h2xf.txt',status='old')                      !
open(unit=43,file='nac20h2yf.txt',status='old')                      !
open(unit=44,file='nac20h2zf.txt',status='old')                      !
open(unit=45,file='nac20h3xf.txt',status='old')                      !
open(unit=46,file='nac20h3yf.txt',status='old')                      !
open(unit=47,file='nac20h3zf.txt',status='old')                      !
open(unit=48,file='nac20h4xf.txt',status='old')                      !
open(unit=49,file='nac20h4yf.txt',status='old')                      !
open(unit=50,file='nac20h4zf.txt',status='old')                      !
open(unit=51,file='nac21c0xf.txt',status='old')                      !
open(unit=52,file='nac21c0yf.txt',status='old')                      !
open(unit=53,file='nac21c0zf.txt',status='old')                      !
open(unit=54,file='nac21h1xf.txt',status='old')                      !
open(unit=55,file='nac21h1yf.txt',status='old')                      !
open(unit=56,file='nac21h1zf.txt',status='old')                      !
open(unit=57,file='nac21h2xf.txt',status='old')                      !
open(unit=58,file='nac21h2yf.txt',status='old')                      !
open(unit=59,file='nac21h2zf.txt',status='old')                      !
open(unit=60,file='nac21h3xf.txt',status='old')                      !
open(unit=61,file='nac21h3yf.txt',status='old')                      !
open(unit=62,file='nac21h3zf.txt',status='old')                      !
open(unit=63,file='nac21h4xf.txt',status='old')                      !
open(unit=64,file='nac21h4yf.txt',status='old')                      !
open(unit=65,file='nac21h4zf.txt',status='old')                      !
do i=1,Nq1                                                          !
  read(21,*) nac10c0x(i,:)                                          !
  read(22,*) nac10c0y(i,:)                                          !
  read(23,*) nac10c0z(i,:)                                          !
  read(24,*) nac10h1x(i,:)                                          !
  read(25,*) nac10h1y(i,:)                                          !
  read(26,*) nac10h1z(i,:)                                          !
  read(27,*) nac10h2x(i,:)                                          !
  read(28,*) nac10h2y(i,:)                                          !
  read(29,*) nac10h2z(i,:)                                          !
  read(30,*) nac10h3x(i,:)                                          !
  read(31,*) nac10h3y(i,:)                                          !
  read(32,*) nac10h3z(i,:)                                          !
  read(33,*) nac10h4x(i,:)                                          !
  read(34,*) nac10h4y(i,:)                                          !
  read(35,*) nac20h4z(i,:)                                          !
  read(36,*) nac20c0x(i,:)                                          !
  read(37,*) nac20c0y(i,:)                                          !
  read(38,*) nac20c0z(i,:)                                          !
  read(39,*) nac20h1x(i,:)                                          !
  read(40,*) nac20h1y(i,:)                                          !
  read(41,*) nac20h1z(i,:)                                          !
  read(42,*) nac20h2x(i,:)                                          !
  read(43,*) nac20h2y(i,:)                                          !
  read(44,*) nac20h2z(i,:)                                          !
  read(45,*) nac20h3x(i,:)                                          !
  read(46,*) nac20h3y(i,:)                                          !
  read(47,*) nac20h3z(i,:)                                          !
  read(48,*) nac20h4x(i,:)                                          !
  read(49,*) nac20h4y(i,:)                                          !
  read(50,*) nac20h4z(i,:)                                          !
  read(51,*) nac21c0x(i,:)                                          !
  read(52,*) nac21c0y(i,:)                                          !
  read(53,*) nac21c0z(i,:)                                          !
  read(54,*) nac21h1x(i,:)                                          !
  read(55,*) nac21h1y(i,:)                                          !
  read(56,*) nac21h1z(i,:)                                          !
  read(57,*) nac21h2x(i,:)                                          !
  read(58,*) nac21h2y(i,:)                                          !
  read(59,*) nac21h2z(i,:)                                          !
  read(60,*) nac21h3x(i,:)                                          !
  read(61,*) nac21h3y(i,:)                                          !
  read(62,*) nac21h3z(i,:)                                          !
  read(63,*) nac21h4x(i,:)                                          !
  read(64,*) nac21h4y(i,:)                                          !
  read(65,*) nac21h4z(i,:)                                          !
end do                                                              !
do i=21,65                                                          !
  close(unit=i)                                                     !
end do                                                              !
!-------------------------------------------------------------------!
!vectorizing the grid nac matrix
cont=0
vec1=0_dp;vec2=0_dp;vec3=0_dp;vec4=0_dp;vec5=0_dp;vec6=0_dp;
do j=0,Nq2-1
  do i=0,Nq1-1
vec1(cont)=nac10c0x(i+1,j+1)*q1(1)/mass(1)+nac10c0y(i+1,j+1)*q1(2)/mass(2)+nac10c0z(i+1,j+1)*q1(3)/mass(3)+nac10h1x(i+1,j+1)*q1(4)/&
mass(4)+nac10h1y(i+1,j+1)*q1(5)/mass(5)+nac10h1z(i+1,j+1)*q1(6)/mass(6)+nac10h2x(i+1,j+1)*q1(7)/mass(7)+nac10h2y(i+1,j+1)*q1(8)/&
mass(8)+nac10h2z(i+1,j+1)*q1(9)/mass(9)+nac10h3x(i+1,j+1)*q1(10)/mass(10)+nac10h3y(i+1,j+1)*q1(11)/mass(11)+nac10h3z(i+1,j+1)*&
q1(12)/mass(12)+nac10h4x(i+1,j+1)*q1(13)/mass(13)+nac10h4y(i+1,j+1)*q1(14)/mass(14)+nac10h4z(i+1,j+1)*q1(15)/mass(15)
vec2(cont)=nac10c0x(i+1,j+1)*q2(1)/mass(1)+nac10c0y(i+1,j+1)*q2(2)/mass(2)+nac10c0z(i+1,j+1)*q2(3)/mass(3)+nac10h1x(i+1,j+1)*q2(4)/&
mass(4)+nac10h1y(i+1,j+1)*q2(5)/mass(5)+nac10h1z(i+1,j+1)*q2(6)/mass(6)+nac10h2x(i+1,j+1)*q2(7)/mass(7)+nac10h2y(i+1,j+1)*q2(8)/&
mass(8)+nac10h2z(i+1,j+1)*q2(9)/mass(9)+nac10h3x(i+1,j+1)*q2(10)/mass(10)+nac10h3y(i+1,j+1)*q2(11)/mass(11)+nac10h3z(i+1,j+1)*&
q2(12)/mass(12)+nac10h4x(i+1,j+1)*q2(13)/mass(13)+nac10h4y(i+1,j+1)*q2(14)/mass(14)+nac10h4z(i+1,j+1)*q2(15)/mass(15)
vec3(cont)=nac20c0x(i+1,j+1)*q1(1)/mass(1)+nac20c0y(i+1,j+1)*q1(2)/mass(2)+nac20c0z(i+1,j+1)*q1(3)/mass(3)+nac20h1x(i+1,j+1)*q1(4)/&
mass(4)+nac20h1y(i+1,j+1)*q1(5)/mass(5)+nac20h1z(i+1,j+1)*q1(6)/mass(6)+nac20h2x(i+1,j+1)*q1(7)/mass(7)+nac20h2y(i+1,j+1)*q1(8)/&
mass(8)+nac20h2z(i+1,j+1)*q1(9)/mass(9)+nac20h3x(i+1,j+1)*q1(10)/mass(10)+nac20h3y(i+1,j+1)*q1(11)/mass(11)+nac20h3z(i+1,j+1)*&
q1(12)/mass(12)+nac20h4x(i+1,j+1)*q1(13)/mass(13)+nac20h4y(i+1,j+1)*q1(14)/mass(14)+nac20h4z(i+1,j+1)*q1(15)/mass(15)
vec4(cont)=nac20c0x(i+1,j+1)*q2(1)/mass(1)+nac20c0y(i+1,j+1)*q2(2)/mass(2)+nac20c0z(i+1,j+1)*q2(3)/mass(3)+nac20h1x(i+1,j+1)*q2(4)/&
mass(4)+nac20h1y(i+1,j+1)*q2(5)/mass(5)+nac20h1z(i+1,j+1)*q2(6)/mass(6)+nac20h2x(i+1,j+1)*q2(7)/mass(7)+nac20h2y(i+1,j+1)*q2(8)/&
mass(8)+nac20h2z(i+1,j+1)*q2(9)/mass(9)+nac20h3x(i+1,j+1)*q2(10)/mass(10)+nac20h3y(i+1,j+1)*q2(11)/mass(11)+nac20h3z(i+1,j+1)*&
q2(12)/mass(12)+nac20h4x(i+1,j+1)*q2(13)/mass(13)+nac20h4y(i+1,j+1)*q2(14)/mass(14)+nac20h4z(i+1,j+1)*q2(15)/mass(15)
vec5(cont)=nac21c0x(i+1,j+1)*q1(1)/mass(1)+nac21c0y(i+1,j+1)*q1(2)/mass(2)+nac21c0z(i+1,j+1)*q1(3)/mass(3)+nac21h1x(i+1,j+1)*q1(4)/&
mass(4)+nac21h1y(i+1,j+1)*q1(5)/mass(5)+nac21h1z(i+1,j+1)*q1(6)/mass(6)+nac21h2x(i+1,j+1)*q1(7)/mass(7)+nac21h2y(i+1,j+1)*q1(8)/&
mass(8)+nac21h2z(i+1,j+1)*q1(9)/mass(9)+nac21h3x(i+1,j+1)*q1(10)/mass(10)+nac21h3y(i+1,j+1)*q1(11)/mass(11)+nac21h3z(i+1,j+1)*&
q1(12)/mass(12)+nac21h4x(i+1,j+1)*q1(13)/mass(13)+nac21h4y(i+1,j+1)*q1(14)/mass(14)+nac21h4z(i+1,j+1)*q1(15)/mass(15)
vec6(cont)=nac21c0x(i+1,j+1)*q2(1)/mass(1)+nac21c0y(i+1,j+1)*q2(2)/mass(2)+nac21c0z(i+1,j+1)*q2(3)/mass(3)+nac21h1x(i+1,j+1)*q2(4)/&
mass(4)+nac21h1y(i+1,j+1)*q2(5)/mass(5)+nac21h1z(i+1,j+1)*q2(6)/mass(6)+nac21h2x(i+1,j+1)*q2(7)/mass(7)+nac21h2y(i+1,j+1)*q2(8)/&
mass(8)+nac21h2z(i+1,j+1)*q2(9)/mass(9)+nac21h3x(i+1,j+1)*q2(10)/mass(10)+nac21h3y(i+1,j+1)*q2(11)/mass(11)+nac21h3z(i+1,j+1)*&
q2(12)/mass(12)+nac21h4x(i+1,j+1)*q2(13)/mass(13)+nac21h4y(i+1,j+1)*q2(14)/mass(14)+nac21h4z(i+1,j+1)*q2(15)/mass(15)
    cont=cont+1
  end do
end do
!$OMP PARALLEL DO shared(moq1,moq2,ham)
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
!!$OMP PARALLEL DO shared(ha,ham)
!do i=0,Nst*Nq1*Nq2-1
!  do j=0,Nst*Nq1*Nq2-1
!    Ha(i,j)=Ha(i,j)-ham(i,j) ! Minus because of the sign of the kinetic energy term
!  end do
!end do
!!$OMP END PARALLEL DO

!$ trunf = omp_get_wtime()
write(*,'(a38,f9.1,a9)')'Time = ', (trunf-truni)/60.d0, ' minutes.'
!$ truni = omp_get_wtime()

write(*,*)'NAC EVALUATED, STARTING TO CREATE CSR VECTORS'

!Creating the CSR vectors --------------------------------------------------------------------------------------------------!
k=0                                                                                                                         !
lol=0.d0                                                                                                                    !
!$OMP parallel do shared(lol)                                                                                               !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (ham(i,j) /= 0.d0) then                                                                                              !
      lol(i) = lol(i)+1                                                                                                     !
    end if                                                                                                                  !
  end do                                                                                                                    !
end do                                                                                                                      !
!$OMP end parallel do                                                                                                       !
k = sum(lol)                                                                                                                !
k_nac=k                                                                                                                     !
write(*,*)'k_Ha (only NAC) = ',k_nac                                                                                        !
!$ trunf = omp_get_wtime()                                                                                                  !
write(*,'(a38,f9.1,a9)')'Time = ', (trunf-truni)/60.d0, ' minutes.'                                                         !
!$ truni = omp_get_wtime()                                                                                                  !
                                                                                                                            !
allocate(nac_val(0:k-1),nac_rowc(0:n),nac_row_col(0:k-1,0:1))                                                               !
nac_val = 1.23456789d0                                                                                                      !
!                                                                                                                           !
nac_rowc=0.d0                                                                                                               !
!$ trunf = omp_get_wtime()                                                                                                  !
write(*,'(a38,f9.1,a9)')'Time = ', (trunf-truni)/60.d0, ' minutes.'                                                         !
!$ truni = omp_get_wtime()                                                                                                  !
k=0                                                                                                                         !
do i=0,n-1 !running through rows                                                                                            !
  do j=0,n-1 !running through columns                                                                                       !
    if (ham(i,j) /= 0.d0) then                                                                                              !
      nac_val(k) = - ham(i,j)  !storing each non-zero element                                                               !
      nac_row_col(k,0)=i !storing the row index of each non-zero element                                                    !
      nac_row_col(k,1)=j !storing the column index of each non-zero element                                                 !
      k=k+1                                                                                                                 !
    end if                                                                                                                  !
  end do                                                                                                                    !
  nac_rowc(i+1)=k !storing the counting of non-zero elements in each row                                                    !
end do                                                                                                                      !
!---------------------------------------------------------------------------------------------------------------------------!
!$ trunf = omp_get_wtime()
write(*,'(a38,f9.1,a9)')'Time = ', (trunf-truni)/60.d0, ' minutes.'
!$ truni = omp_get_wtime()
write(*,*)'Parte só do NAC feita'

! Merging this matrix with the one that includes the second derivatives and the potential energy, Ha
!$OMP PARALLEL DO shared(ha,ham)
do i=0,Nst*Nq1*Nq2-1
  do j=0,Nst*Nq1*Nq2-1
    Ha(i,j)=Ha(i,j)-ham(i,j) ! Minus because of the sign of the kinetic energy term
  end do
end do
!$OMP END PARALLEL DO
!$ trunf = omp_get_wtime()
write(*,'(a38,f9.1,a9)')'Time = ', (trunf-truni)/60.d0, ' minutes.'
!$ truni = omp_get_wtime()
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
integer n

n = s * Nst
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
!Creating the CSR vectors --------------------------------------------------------------------------------------------------!     !
lol=0.d0                                                                                                                    !     !
k=0                                                                                                                         !     !
!$OMP parallel do shared(lol)                                                                                               !     !
do i=0,n-1 !running through rows                                                                                            !     !
  do j=0,n-1 !running through columns                                                                                       !     !
    if (ax(i+1,j+1) /= 0.d0) then                                                                                           !     !
      lol(i) = lol(i) + 1                                                                                                   !     !
    end if                                                                                                                  !     !
  end do                                                                                                                    !     !
end do                                                                                                                      !     !
!$OMP end parallel do                                                                                                       !     !
k = sum(lol)                                                                                                                !     !
k_kine = k                                                                                                                  !     !
write(*,*)'k kine= ',k                                                                                                      !     !
allocate(kine_val(0:k-1),kine_rowc(0:n),kine_row_col(0:k-1,0:1))                                                            !     !
!                                                                                                                           !     !
kine_rowc=0.d0                                                                                                              !     !
k=0                                                                                                                         !     !
do i=0,n-1 !running through rows                                                                                            !     !
  do j=0,n-1 !running through columns                                                                                       !     !
    if (ax(i+1,j+1) /= 0.d0) then                                                                                           !     !
      kine_val(k)=ax(i+1,j+1)  !storing each non-zero element                                                               !     !
      kine_row_col(k,0)=i !storing the row index of each non-zero element                                                   !     !
      kine_row_col(k,1)=j !storing the column index of each non-zero element                                                !     !
      k=k+1                                                                                                                 !     !
    end if                                                                                                                  !     !
  end do                                                                                                                    !     !
  kine_rowc(i+1)=k !storing the counting of non-zero elements in each row                                                   !     !
end do                                                                                                                      !     !
!---------------------------------------------------------------------------------------------------------------------------!     !
allocate (Ha(0:Nst*Nq1*Nq2-1,0:Nst*Nq1*Nq2-1))                                                                                    !
!$OMP PARALLEL DO shared(Ha)                                                                                                      !
do i=0,Nst*Nq1*Nq2-1                                                                                                              !
  do j=0,Nst*Nq1*Nq2-1                                                                                                            !
    Ha(i,j)=0.0d0                                                                                                                 !
  end do                                                                                                                          !
end do                                                                                                                            !
!$OMP end parallel do                                                                                                             !
do j=0*Nq1,(Nq2-1)*Nq1,Nq1 ! Moving through the boxes, from the fisrt to the last one                                             !
  do i=1,Nq1 ! Moving inside the box                                                                                              !
    Ha( 0*s + j + i -1, 0*s + j + i -1) = Ha( 0*s + j + i -1, 0*s + j + i -1) + pot1 (i,(j+Nq1)/Nq1)                              !
    Ha( 1*s + j + i -1, 1*s + j + i -1) = Ha( 1*s + j + i -1, 1*s + j + i -1) + pot2 (i,(j+Nq1)/Nq1)                              !
    Ha( 2*s + j + i -1, 2*s + j + i -1) = Ha( 2*s + j + i -1, 2*s + j + i -1) + pot3 (i,(j+Nq1)/Nq1)                              !
  end do ! end of do through single box                                                                                           !
end do ! end of do through boxes                                                                                                  !
!Creating the CSR vectors --------------------------------------------------------------------------------------------------!     !
lol=0.d0                                                                                                                    !     !
k=0                                                                                                                         !     !
!$OMP parallel do shared(lol)                                                                                               !     !
do i=0,n-1 !running through rows                                                                                            !     !
  do j=0,n-1 !running through columns                                                                                       !     !
    if (Ha(i,j) /= 0.d0) then                                                                                               !     !
      lol(i) = lol(i) + 1                                                                                                   !     !
    end if                                                                                                                  !     !
  end do                                                                                                                    !     !
end do                                                                                                                      !     !
!$OMP end parallel do                                                                                                       !     !
k = sum(lol)                                                                                                                !     !
k_pot=k                                                                                                                     !     !
write(*,*)'k pote= ',k                                                                                                      !     !
allocate(pot_val(0:k-1),pot_rowc(0:n),pot_row_col(0:k-1,0:1))                                                               !     !
!                                                                                                                           !     !
pot_rowc=0.d0                                                                                                               !     !
k=0                                                                                                                         !     !
do i=0,n-1 !running through rows                                                                                            !     !
  do j=0,n-1 !running through columns                                                                                       !     !
    if (Ha(i,j) /= 0.d0) then                                                                                               !     !
      pot_val(k)=Ha(i,j)  !storing each non-zero element                                                                    !     !
      pot_row_col(k,0)=i !storing the row index of each non-zero element                                                    !     !
      pot_row_col(k,1)=j !storing the column index of each non-zero element                                                 !     !
      k=k+1                                                                                                                 !     !
    end if                                                                                                                  !     !
  end do                                                                                                                    !     !
  pot_rowc(i+1)=k !storing the counting of non-zero elements in each row                                                    !     !
end do                                                                                                                      !     !
!---------------------------------------------------------------------------------------------------------------------------!     !
! Adding potential energy - moving only along the diagonal of the axmiltonian                                                     !
do j=0*Nq1,(Nq2-1)*Nq1,Nq1 ! Moving through the boxes, from the fisrt to the last one                                             !
  do i=1,Nq1 ! Moving inside the box                                                                                              !
    ax( 0*s + j + i , 0*s + j + i ) = ax( 0*s + j + i , 0*s + j + i ) + pot1 (i,(j+Nq1)/Nq1)                                      !
    ax( 1*s + j + i , 1*s + j + i ) = ax( 1*s + j + i , 1*s + j + i ) + pot2 (i,(j+Nq1)/Nq1)                                      !
    ax( 2*s + j + i , 2*s + j + i ) = ax( 2*s + j + i , 2*s + j + i ) + pot3 (i,(j+Nq1)/Nq1)                                      !
  end do ! end of do through single box                                                                                           !
end do ! end of do through boxes                                                                                                  !
!=================================================================================================================================!

!$OMP PARALLEL DO shared(Ha,ax) 
do i=1,Nst*Nq1*Nq2
  do j=1,Nst*Nq1*Nq2
    Ha(i-1,j-1)=ax(i,j) ! Reseting Ha to include only the value of the kinetic energy and potential energy in this point.
  end do
end do
!$OMP end parallel do

deallocate(ax)

end subroutine hamiltonian_matrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!Subroutine that calculates the photoionization amplitudes for a given geometry and a given electric field
!!Some parts of this function assumes that the number of electronic states is 3
!subroutine p_i_a(i1,i2,pic)
!use global_param
!implicit none 
!complex(kind=dp),dimension(Nst) :: pic
!integer                         :: i1,i2,ii,jj,kk
!character(len=21)               :: fname00
!character(len=199)              :: fpath0,fpath1
!integer                         :: file1,file2,file3,file4,file5
!integer,parameter               :: nang=512
!integer,parameter               :: nk=256
!real(kind=dp),allocatable       :: vec1(:,:),vec2(:,:),vec3(:,:)
!real(kind=dp)                   :: k0,k1,k2,e0,e1,e2,p_e(nk)
!real(kind=dp)                   :: phi(nang),theta(nang),domega
!complex(kind=dp),allocatable    :: r0(:,:,:),r1(:,:,:),r2(:,:,:),coef0(:,:),coef1(:,:),coef2(:,:)
!real(kind=dp)                   :: ip0,ip1,ip2
!integer                         :: aux1(Nst)
!real(kind=dp)                   :: aux2(Nst)
!
!allocate(vec1(nang*nk,Nst*2),vec2(nang*nk,Nst*2),vec3(nang*nk,Nst*2))
!allocate(r0(nk,nang,Nst),r1(nk,nang,Nst),r2(nk,nang,Nst),coef0(nk,Nst),coef1(nk,Nst),coef2(nk,Nst))
!call getcwd(fpath0) !getting the working directory path
!fname00='pice_10000000_0_0.dat'
!write(fname00(7:9),'(i0.3)') i2+100
!write(fname00(12:13),'(i0.2)') i1
!fpath1="~/pice_files/"//fname00
!open(newunit=file1,file=fpath1,status='old')
!write(fname00(17:17),'(i1)') 1
!fpath1="~/pice_files/"//fname00
!open(newunit=file2,file=fpath1,status='old')
!write(fname00(17:17),'(i1)') 2
!fpath1="~/pice_files/"//fname00
!open(newunit=file3,file=fpath1,status='old')
!!now read the x, y and z components of the photoionization coupling elements.
!do ii=1,nk*nang
!  read(file1,*) vec1(ii,:)!,vec1(i,2),vec1(i,3),vec1(i,4),vec1(i,5),vec1(i,6)
!  read(file2,*) vec2(ii,:)!,vec2(i,2),vec2(i,3),vec2(i,4),vec2(i,5),vec2(i,6)
!  read(file3,*) vec3(ii,:)!,vec3(i,2),vec3(i,3),vec3(i,4),vec3(i,5),vec3(i,6)
!end do
!kk=1
!do ii=1,nk
!  do jj=1,nang
!    r0(ii,jj,1)=dcmplx( vec1(kk,1) , vec1(kk,2) )
!    r0(ii,jj,2)=dcmplx( vec1(kk,3) , vec1(kk,4) )
!    r0(ii,jj,3)=dcmplx( vec1(kk,5) , vec1(kk,6) )
!    r1(ii,jj,1)=dcmplx( vec2(kk,1) , vec2(kk,2) )
!    r1(ii,jj,2)=dcmplx( vec2(kk,3) , vec2(kk,4) )
!    r1(ii,jj,3)=dcmplx( vec2(kk,5) , vec2(kk,6) )
!    r2(ii,jj,1)=dcmplx( vec3(kk,1) , vec3(kk,2) )
!    r2(ii,jj,2)=dcmplx( vec3(kk,3) , vec3(kk,4) )
!    r2(ii,jj,3)=dcmplx( vec3(kk,5) , vec3(kk,6) )
!    kk=kk+1
!  end do
!end do
!open(newunit=file4,file='test_sym_dist3.txt',status='old')
!do ii=1,nang
!  read(file4,*)theta(ii),phi(ii) !reading the angular distribution used to calculate the photoionization matrix elements 
!end do
!
!ip0 = pot1(i1+27,i2+20) - e1neut(i1+27,i2+20) !ionization potential for ground state of the cation
!ip1 = pot2(i1+27,i2+20) - e1neut(i1+27,i2+20) !ionization potential for first excited state of the cation
!ip2 = pot3(i1+27,i2+20) - e1neut(i1+27,i2+20) !ionization potential for second excited state of the cation
!e0 = e_ip - ip0 !Energy of the ionized electron if the molecule goes for the cation ground state
!e1 = e_ip - ip1 !Energy of the ionized electron if the molecule goes for the cation first excited state
!e2 = e_ip - ip2 !Energy of the ionized electron if the molecule goes for the cation second excited state
!
!domega=4*pi/nang
!coef0=0.d0
!coef1=0.d0
!coef2=0.d0
!do ii=1,nk
!  do jj=1,nang
!    coef0(ii,1) = coef0(ii,1) + r0(ii,jj,1) * domega !Doing the operation int( PICE(k) * domega )
!    coef0(ii,2) = coef0(ii,2) + r0(ii,jj,2) * domega !Doing the operation int( PICE(k) * domega )
!    coef0(ii,3) = coef0(ii,3) + r0(ii,jj,3) * domega !Doing the operation int( PICE(k) * domega )
!    coef1(ii,1) = coef1(ii,1) + r1(ii,jj,1) * domega !Doing the operation int( PICE(k) * domega )
!    coef1(ii,2) = coef1(ii,2) + r1(ii,jj,2) * domega !Doing the operation int( PICE(k) * domega )
!    coef1(ii,3) = coef1(ii,3) + r1(ii,jj,3) * domega !Doing the operation int( PICE(k) * domega )
!    coef2(ii,1) = coef2(ii,1) + r2(ii,jj,1) * domega !Doing the operation int( PICE(k) * domega )
!    coef2(ii,2) = coef2(ii,2) + r2(ii,jj,2) * domega !Doing the operation int( PICE(k) * domega )
!    coef2(ii,3) = coef2(ii,3) + r2(ii,jj,3) * domega !Doing the operation int( PICE(k) * domega )
!  end do
!  p_e(ii) = 0.005859375d0 + (ii-1) * 0.00588235294117647d0 !Momentum values in which the photoionization matrix elements are spanned
!end do
!
!!define the correct i, the correct momentum of the electron
!aux2 = e_ip
!do ii=1,nk
!  if ( dsqrt( (p_e(ii) - (e_ip - ip0))**2.d0) < aux2(1) ) then
!    aux1(1) = ii !Defining the value of the momentum of the leaving electron
!    aux2(1) = dsqrt( (p_e(ii) - (e_ip - ip0))**2.d0)
!  end if
!  if ( dsqrt( (p_e(ii) - (e_ip - ip1))**2.d0) < aux2(2) ) then
!    aux1(2) = ii !Defining the value of the momentum of the leaving electron
!    aux2(2) = dsqrt( (p_e(ii) - (e_ip - ip1))**2.d0)
!  end if
!  if ( dsqrt( (p_e(ii) - (e_ip - ip2))**2.d0) < aux2(3) ) then
!    aux1(3) = ii !Defining the value of the momentum of the leaving electron
!    aux2(3) = dsqrt( (p_e(ii) - (e_ip - ip2))**2.d0)
!  end if
!end do
!
!!Do the operation -e * sqrt(2) * E * int(PICE(k) * domega)   --- 'e' is the electron charge, that in atomic units is 1
!pic(1) = - dsqrt(2.d0) * E00 * dot_product( orientation , coef0(aux1(1),:) )
!pic(2) = - dsqrt(2.d0) * E00 * dot_product( orientation , coef1(aux1(2),:) )
!pic(3) = - dsqrt(2.d0) * E00 * dot_product( orientation , coef2(aux1(3),:) )
!close(unit=file1)
!close(unit=file2)
!close(unit=file3)
!close(unit=file4)
!end subroutine p_i_a






