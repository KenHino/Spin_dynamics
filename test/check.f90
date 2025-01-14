program  check
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spinop
    use hamiltonian
    use utils
    implicit none

    complex(dp) :: w1(3), w2(3),D(3,3)
    integer :: g_I(2), i
    real(dp) :: kur = 0.5_dp
    complex(dp) :: a(2,2)
    complex(dp) :: Ahyp(2,3,3)
    ! w2 = cmplx([0.0_dp, 0.0_dp, 1.0_dp], 0.0_dp, dp)
    ! D = cmplx(0.0_dp, 0.0_dp, dp)

    real(dp) :: Iz(10)
    complex(dp) :: H(8,8)

    complex(dp), allocatable    :: S_x(:,:), S_y(:,:), S_z(:,:)

    call getSpinop(2, S_x, S_y, S_z)

    print*, S_x

    w1 = cmplx([0.0_dp, 0.0_dp, 1.0_dp], 0.0_dp)
    w2 = cmplx([0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp)
    g_I = [2,2]
    

    Ahyp = cmplx(1.0_dp, 0.0_dp)

    Iz = 0.0_dp
    
    ! a = getSz(kur)

    ! print*, 'Sz operator'
    ! do i = 1, 2
    !     print *, real(a(i,:), kind=dp)
    ! end do 

    H = getHi(w2, g_I, Ahyp)
    print*, 'Hyperfine Hamiltonian:', size(H)

    do i = 1, 8
        print *, real(H(i,:), kind=dp)
    end do

    ! print*, 'Electron Hamniltonian:', electronHamiltonian(w1,w2,cmplx(0.0_dp, 0.0_dp, dp),D)



end program  check

