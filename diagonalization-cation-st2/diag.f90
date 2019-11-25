program dyn
use global_param
use omp_lib
implicit none
integer n,jj !m=dimension of the Hamiltonian (Nq1*Nq2*Nst x Nq1*Nq2*Nst)
integer :: LDA,LWORK,LIWORK,INFO
integer,allocatable :: IWORK(:)
real(kind=dp), allocatable :: EIG(:),WORK(:)
character(len=1) :: JOBZ,UPLO
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=Nst*Nq1*Nq2 ! m=Nq1*Nq2*st %size of the Hamiltonian
call hamiltonian_matrix


JOBZ='V'
UPLO='U'
N=Nq1*Nq2*Nst
allocate(EIG(N), source = 0.d0 )
LDA=N
LWORK= 1 + 6*N + 2*N**2    + 10*s
!LWORK = 351230976
allocate(WORK(LWORK))
WORK(:)=0.d0
LIWORK=3 + 5*N             + 10*s
!LIWORK = 18608
allocate(IWORK(LIWORK))
IWORK(:)=0
INFO=0

write(*,*)'INFO = ',INFO
write(*,*)'UPLO = ',UPLO
write(*,*)'LWORK = ',LWORK
write(*,*)'LIWORK = ',LIWORK


call dsyevd(JOBZ,UPLO,N,ha,LDA,EIG,WORK,LWORK,IWORK,LIWORK,INFO)

write(*,*)'Finished'
write(*,*)'INFO = ',INFO,' if = 0 it is fine'
write(*,*)'UPLO = ',UPLO
write(*,*)'LWORK = ',LWORK
write(*,*)'LIWORK = ',LIWORK
write(*,*)'WORK(1) = ',int(WORK(1))
write(*,*)'IWORK(1) = ',IWORK(1)

write(*,'(20f12.8)')EIG(1:10)
open(unit=100,file='eigen-vec-cation-st2',status='unknown')
do i=1,s
  write(100,*)  ha(i,1:100)
end do

end program
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
