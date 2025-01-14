module kroenecker

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use stdlib_sparse


    implicit none
    private 
    public :: kron, kron_iden, get_state

    interface kron
        procedure kron_dense
        procedure kron_sp
    end interface kron

    interface  kron_iden
        procedure kron_iden_right_dense
        procedure kron_iden_left_dense
        procedure kron_iden_right_sp
        procedure kron_iden_left_sp
    end interface kron_iden

    contains

    function kron_dense(A,B) result(C)
    ! Returns Kroenecker product of two dense matrices A and B
       complex(dp), intent(in)  :: A(:,:), B(:,:)
       complex(dp), allocatable :: C(:,:)
       
       integer :: i = 0, j = 0, k = 0, l = 0
       integer :: m = 0, n = 0, p = 0, q = 0

        allocate(C(size(A,1)*size(B,1),size(A,2)*size (B,2)))
        C = cmplx(0.0_dp, 0.0_dp, kind=dp)

        do i = 1,size(A,1)
            do j = 1,size(A,2)
                n=(i-1)*size(B,1) + 1
                m=n+size(B,1) - 1
                p=(j-1)*size(B,2) + 1
                q=p+size(B,2) - 1
                C(n:m,p:q) = A(i,j)*B
            end do
        end do
       
    end function kron_dense

    function kron_sp(A, B) result(C)
    ! Returns Kroenecker product of sparse matrix A with sparse matrix B
        type(COO_cdp_type), intent(in)  :: A, B
        type(COO_cdp_type)              :: C

        integer :: i = 0, j = 0, k

        integer :: nrows, ncols
        integer :: nnz

        nrows = A%nrows * B%nrows 
        ncols = A%ncols * B%ncols
        nnz = A%nnz * B%nnz

        call C%malloc(nrows, ncols, nnz)

        k = 1

        do i = 1, size(A%index, dim=2)
            do j = 1, size(B%index, dim=2)
                C%index(1,k) = B%nrows*(A%index(1,i)-1) + B%index(1,j)
                C%index(2,k) = B%ncols*(A%index(2,i)-1) + B%index(2,j)
                C%data(k) = A%data(i) * B%data(j)
                k=k+1
            end do
        end do

    end function kron_sp

    
    function kron_iden_right_dense(A,N) result(B)
    ! Returns Kroenecker product of matrix A with identity of size ZxZ
       complex(dp), intent(in)  :: A(:,:)
       integer, intent(in)      :: N
       complex(dp), allocatable :: B(:,:)
       
       complex(dp), allocatable :: C(:,:,:,:)
       integer                  :: i                     

        allocate(B(size(A,1)*N,size(A,2)*N))
        allocate(C(N,size(A,1),N,size(A,2)))
            
        C = cmplx(0.0_dp, 0.0_dp, kind=dp)  

        do i = 1,N
            C(i,:,i,:) = A
        end do 

        B = reshape(C, [size(A,1)*N,size(A,2)*N])
       
    end function kron_iden_right_dense


    function kron_iden_left_dense(A,N) result(B)
    ! Returns Kroenecker product of matrix N with identity of size AxA (labels are switched because of way interfacing works)
       integer, intent(in)      :: A
       complex(dp), intent(in)  :: N(:,:)
       complex(dp), allocatable :: B(:,:)
       
       integer :: i = 0, j = 0, k = 0, l = 0
       integer :: m = 0, z = 0, p = 0, q = 0


        allocate(B(size(N,1)*A,size(N,2)*A))
       
        B = cmplx(0.0_dp, 0.0_dp, kind=dp)  
       
        do i = 1,A
                z=A*(i-1) + 1
                m=z+size(N,dim=1) - 1
                p=z+size(N,dim=2) - 1
                B(z:m,z:p) = N
        end do
       
    end function kron_iden_left_dense

    function kron_iden_right_sp(A,N) result(B)
    ! Returns Kroenecker product of matrix A with identity of size ZxZ
        type(COO_cdp_type), intent(in)  :: A
        integer, intent(in)             :: N
        type(COO_cdp_type)              :: B

        integer :: i = 0, j = 0, k = 1

        integer :: nrows, ncols
        integer :: nnz

        nrows = A%nrows * N 
        ncols = A%ncols * N
        nnz = A%nnz * N

        call B%malloc(nrows, ncols, nnz)
       
        do i=1,size(A%data)
            B%index(1, 1+(i-1)*N:i*N) = [(j, j=1+(A%index(1, i)-1)*N, A%index(1, i)*N)]  
            B%index(2, 1+(i-1)*N:i*N) = [(j, j=1+(A%index(2, i)-1)*N, A%index(2, i)*N)]
            B%data(1+(i-1)*N:i*N) = A%data(i)
        end do

    end function kron_iden_right_sp

    function kron_iden_left_sp(A,N) result(B)
    ! Returns Kroenecker product of matrix A with identity of size ZxZ
        integer, intent(in)       :: A
        type(COO_cdp_type), intent(in) :: N
        type(COO_cdp_type)             :: B

        integer :: i = 0, j = 0, k = 1

        integer :: nrows, ncols
        integer :: nnz

        nrows = A * N%nrows 
        ncols = A * N%ncols
        nnz = A * N%nnz

        call B%malloc(nrows, ncols, nnz)
        ! print*, B%index
        ! print*, size(B%data)

        do i=1,A
            B%index(:, 1+(i-1)*size(N%data):i*size(N%data)) = N%index + (i-1)*size(N%data)
            B%data(1+(i-1)*size(N%data):i*size(N%data)) = N%data
        end do
       
    end function kron_iden_left_sp

    subroutine get_state(el_state, sample, Z)

        complex(dp), intent(in) :: el_state(4)
        complex(dp), intent(in) :: sample(:,:)

    end subroutine get_state

end module kroenecker 