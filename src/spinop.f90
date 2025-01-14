module spinop
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use utils
    use stdlib_sparse

    implicit none
    private 
    public :: getSpinop

    interface getSpinop
        ! Construct Sx,Sy and Sz operators for particle with spin multiplicity g_I 
        procedure getSpinop_dense
        procedure getSpinop_sp
    end interface getSpinop

    contains

    subroutine getSpinop_dense(g_I, Sx, Sy, Sz)
        integer, intent(in)                     :: g_I
        complex(dp), allocatable, intent(inout) :: Sx(:, :), Sy(:, :), Sz(:, :) 

        real(dp)              :: S 
        real(dp), allocatable :: Ms(:)
        integer               :: i

        S = (real(g_I, kind=dp)-1.0_dp)/2.0_dp

        allocate(Sx(g_I, g_I))
        allocate(Sy(g_I, g_I))
        allocate(Sz(g_I, g_I))
        allocate(Ms(g_I))

        Sx = cmplx(0.0_dp, 0.0_dp)
        Sy = cmplx(0.0_dp, 0.0_dp)
        Sz = cmplx(0.0_dp, 0.0_dp)

        Ms(1) = S
        do i = 2,g_I
            Ms(i) = Ms(i-1) - 1.0_dp
        end do

        do i = 1,g_I-1
            Sx(i, i+1) = cmplx(0.5_dp*sqrt(S*(S+1) - Ms(i+1)*Ms(i)), 0.0_dp)  
            Sx(i+1, i) = Sx(i, i+1) 
        end do

        do i = 1,g_I-1
            Sy(i, i+1) = cmplx(0.0_dp, -0.5_dp*sqrt(S*(S+1) - Ms(i+1)*Ms(i)))  
            Sy(i+1, i) = conjg(Sy(i, i+1)) 
        end do

        do i = 1,g_I
            Sz(i, i) = Ms(i)
        end do

    end subroutine getSpinop_dense

 
    subroutine getSpinop_sp(g_I, Sx, Sy, Sz)
        integer, intent(in)               :: g_I
        type(COO_cdp_type), intent(inout) :: Sx, Sy, Sz 

        complex(dp), allocatable :: Sx_dense(:, :), Sy_dense(:, :), Sz_dense(:, :) 

        call getSpinop_dense(g_I, Sx_dense, Sy_dense, Sz_dense)

        call dense2coo(Sx_dense, Sx)
        call dense2coo(Sy_dense, Sy)
        call dense2coo(Sz_dense, Sz)

    end subroutine getSpinop_sp



end module spinop