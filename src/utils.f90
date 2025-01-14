module utils
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use stdlib_sparse

    implicit none
    public :: eye, print_matrix
    private :: scale_sparse, add_sparse

    real(dp) :: pi = 4.0_dp*atan(1.0_dp)

    interface operator(*)
        procedure scale_sparse
    end interface

    interface operator(+)
        procedure add_sparse
    end interface

    contains

    function scale_sparse(c, A) result(B)
    ! Multiplies sparse matrix A by constant c
        complex(dp), intent(in)        :: c
        type(COO_cdp_type), intent(in) :: A
        type(COO_cdp_type)             :: B

        call B%malloc(num_rows=A%nrows, num_cols=A%ncols, nnz=A%nnz)
        B%index = A%index
        B%data = c*A%data

    end function scale_sparse
    
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


    function eye(N) result(id)
    ! Generates an indentity matrix of Nension N
        integer, intent(in)   :: N
        complex(dp), allocatable :: id(:, :)
        
        integer :: i

        allocate(id(N, N))
        
        id = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, N
            id(i, i) = cmplx(1.0_dp, 0.0_dp, kind=dp)
        end do

    end function eye

    subroutine print_matrix(a)
    ! Prints 2D array as a matrix
        ! complex(dp), intent(in) :: a(:,:)
        integer, intent(in) :: a(:,:)

        integer :: i

        do i = 1, size(a, 1)
            print*, a(i, :)
        end do
        
    end subroutine 

end module utils