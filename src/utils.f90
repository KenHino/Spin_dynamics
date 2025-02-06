module utils
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use stdlib_sparse

    implicit none
    public :: eye_cp, print_matrix, coo_to_csr
    private :: scale_sparse_cdp, scale_sparse_dp, add_sparse


    real(dp) :: pi = 4.0_dp*atan(1.0_dp)


    interface operator(*)
        procedure scale_sparse_cdp
        procedure scale_sparse_dp
    end interface

    interface operator(+)
        procedure add_sparse
    end interface

    contains


    subroutine coo_to_csr(H_coo, H_csr)
    ! Converts unordered COO sparse matrix to CSR
        type(COO_cdp_type), intent(inout)  :: H_coo
        type(CSR_cdp_type), intent(out)    :: H_csr

        integer, allocatable :: row_elem(:) ! Array containing the number of non-zero elements in each row of H 
        integer              :: indx        
        
        integer :: i

        call H_csr%malloc(H_coo%nrows, H_coo%ncols, H_coo%nnz)
        allocate(row_elem(H_coo%nrows), source=0)

        ! Count number of elements in each row
        do i=1,H_coo%nnz
            row_elem(H_coo%index(1,i)) = row_elem(H_coo%index(1,i)) + 1
        end do

        ! Convert number of elements to the rowpointer format for a CSR matrix
        H_csr%rowptr(1) = 1
        do i=1,H_coo%nrows
            H_csr%rowptr(i+1) = H_csr%rowptr(i) + row_elem(i)  
        end do

        ! Start moving elements from COO to CSR matrix
        do i=1,H_coo%nnz
            ! Find index where the row of COO element i starts
            indx = H_csr%rowptr(H_coo%index(1,i)) 
            indx = indx + row_elem(H_coo%index(1,i)) - 1
            row_elem(H_coo%index(1,i)) = row_elem(H_coo%index(1,i)) - 1 
            H_csr%col(indx) = H_coo%index(2,i) 
            H_csr%data(indx) = H_coo%data(i) 
        end do

        deallocate(H_coo%data)
        deallocate(H_coo%index)

    end subroutine coo_to_csr

    function scale_sparse_cdp(c, A) result(B)
    ! Multiplies sparse matrix A by constant c
        complex(dp), intent(in)        :: c
        type(COO_cdp_type), intent(in) :: A
        type(COO_cdp_type)             :: B

        call B%malloc(num_rows=A%nrows, num_cols=A%ncols, nnz=A%nnz)
        B%index = A%index
        B%data = c*A%data

    end function scale_sparse_cdp
    
    function scale_sparse_dp(c, A) result(B)
    ! Multiplies sparse matrix A by constant c
        real(dp), intent(in)           :: c
        type(COO_cdp_type), intent(in) :: A
        type(COO_cdp_type)             :: B

        call B%malloc(num_rows=A%nrows, num_cols=A%ncols, nnz=A%nnz)
        B%index = A%index
        B%data = c*A%data

    end function scale_sparse_dp
    

    function add_sparse(A, B) result(C)
    ! Adds two sparse matrices with the same index structure
        type(COO_cdp_type), intent(in) :: A, B
        type(COO_cdp_type)             :: C

        if (all(A%index == B%index)) then
            call C%malloc(num_rows=A%nrows, num_cols=A%ncols, nnz=A%nnz)
            C%index = A%index
            C%data  = A%data + B%data
        else
            print*, 'Cannot add two sparse matrices with different sets of arrays'
        end if

    end function add_sparse

    function eye_cp(N) result(id)
    ! Generates an indentity matrix of Nension N
        integer, intent(in)      :: N
        complex(dp), allocatable :: id(:, :)
        
        integer :: i

        allocate(id(N, N))
        
        id = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, N
            id(i, i) = cmplx(1.0_dp, 0.0_dp, kind=dp)
        end do

    end function eye_cp

    subroutine print_matrix(a)
    ! Prints 2D array as a matrix
        complex(dp), intent(in) :: a(:,:)
        ! integer, intent(in) :: a(:,:)

        integer :: i

        do i = 1, size(a, 1)
             write(*, "(*('('sf6.4xspf6.4x'i)':x))") a(i,:)
            ! print*, a(i, :) 
        end do
        
    end subroutine 

end module utils