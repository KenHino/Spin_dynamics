module spinop
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use utils
    use stdlib_sparse

    implicit none
    private 
    public :: get_spinop, spin_magnitude, spin_half_mag

    real(dp), parameter :: spin_half_mag = sqrt(0.75_dp) ! Magnitude of total spin vector of a spin-1/2 particle


    interface get_spinop
    ! Construct Sx,Sy,Sz,S+ and S- operators for particle with spin multiplicity g_I 
        procedure get_spinop_dense
        procedure get_spinop_sp
    end interface get_spinop

    contains

    subroutine get_spinop_dense(g_I, Sx, Sy, Sz)
        integer, intent(in)                     :: g_I
        complex(dp), allocatable, intent(inout) :: Sx(:, :), Sy(:, :), Sz(:, :) 
        
        real(dp)              :: S 
        real(dp), allocatable :: Ms(:)
        integer               :: i

        S = (real(g_I, kind=dp)-1.0_dp)/2.0_dp

        if (allocated(Sz)) then
            deallocate(Sx, Sy, Sz)
        end if 

        allocate(Sx(g_I, g_I), source=(0.0_dp, 0.0_dp))
        allocate(Sy(g_I, g_I), source=(0.0_dp, 0.0_dp))
        allocate(Sz(g_I, g_I), source=(0.0_dp, 0.0_dp))
        allocate(Ms(g_I))

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

    end subroutine get_spinop_dense

    subroutine get_spinop_sp(mult, Sz, Sp, Sm)
        integer, intent(in)               :: mult
        type(COO_cdp_type), intent(inout) :: Sz, Sp, Sm 

        real(dp)              :: S 
        real(dp), allocatable :: Ms(:)
        integer               :: i

        if (allocated(Sz%data)) then
            deallocate(Sz%data, Sp%data, Sm%data)
            deallocate(Sz%index, Sp%index, Sm%index)
        end if 

        S = (real(mult, kind=dp)-1.0_dp)/2.0_dp

        call Sz%malloc(num_rows=mult, num_cols=mult, nnz=mult)
        call Sp%malloc(num_rows=mult, num_cols=mult, nnz=mult-1)
        call Sm%malloc(num_rows=mult, num_cols=mult, nnz=mult-1)
        allocate(Ms(mult))

        Sz%data(1) = cmplx(S, 0.0_dp, kind=dp)
        Sz%index(:, 1) = 1 
        
        do i=2,mult
            Sz%data(i) = Sz%data(i-1) - (1.0_dp, 0.0_dp)
            Sz%index(:, i) = i

            Sp%data(i-1) = sqrt(S*(S+1) - Sz%data(i)*Sz%data(i-1))
            Sp%index(1,i-1) = i-1  
            Sp%index(2,i-1) = i

            Sm%data(i-1) = sqrt(S*(S+1) - Sz%data(i)*Sz%data(i-1))
            Sm%index(1,i-1) = i  
            Sm%index(2,i-1) = i-1
        end do


    end subroutine get_spinop_sp

    function spin_magnitude(mult) result(S)
        integer, intent(in) :: mult
        real(dp)            :: S

        ! S = sqrt(real((mult-1)*S, kind=dp))
        
    end function spin_magnitude

end module spinop