program main
    use, intrinsic :: iso_fortran_env, only: dp => real64, i8 => int64
    use spinop
    use hamiltonian
    use variables
    use class_observables
    use utils
    use dynamics
    use stdlib_stats_distribution_normal
    ! use stdlib_random, only: random_seed
    use m_random
    use moments
    use symmetry
    implicit none

    ! integer :: N = 3
    ! real(dp) :: a(3), w(3)
    ! real(dp) :: a1(2), w2(2)

    ! a = 1.0_dp
    ! w = 1.0_dp
    ! call qrule(N, w, a, 2, a1, w2)


    character(:), allocatable :: input_file
    character(:), allocatable :: output_folder
    character(:), allocatable :: folder
    character(len=16) :: guz
    type(sys_param)           :: sys
    type(sim_param)           :: sim
    type(observables)         :: res
    type(RNG_t)               :: rng
    type(COO_cdp_type)        :: H_coo
    type(CSR_cdp_type)        :: H

    integer(8) :: seed(2)
    integer :: i

    real(dp), allocatable :: a(:)
    real(dp), allocatable :: w(:)
    integer :: N
    integer :: M
    real(dp), allocatable :: a_bar(:)
    real(dp), allocatable :: a2_bar(:)
    real(dp), allocatable :: w_bar(:)
    integer, allocatable :: n_bar(:)

    integer ::  n_trial(4)
    integer, allocatable ::  k_trial(:,:)

    n_trial = [4, 3, 3, 4]

    ! call cartesian_product(n_trial, k_trial)

    folder = '/home/damianko/fpm/spinchem/data/trial_sym'
    ! folder = '/home/damianko/fpm/spinchem/cpf_ini'
    input_file = '/home/damianko/fpm/spinchem/data/many.ini'
    call read_inp(input_file, sys, sim)

    call get_command_argument(1, input_file)
    call get_command_argument(2, folder)

    call rng%set_random_seed()
    seed = rng%s
    
    do i=1,size(sim%B)
!         print*, sim%B(i)

        call rng%set_seed(seed)

        ! Calculate Larmor frequencies of both electrons in mT units
        sys%e1%w = (sys%e1%g/g_e)*sim%B(i)
        sys%e2%w = (sys%e2%g/g_e)*sim%B(i)

        write(guz,'(F0.3)') sim%B(i) ! converting integer to string using a 'internal file'

        output_folder = folder // '/' // trim(guz)  

        call system(' mkdir ' // output_folder)

        call symmetrised_dynamics(sys, sim, rng, res, output_folder)
        ! call trace_sampling(sys, sim, rng, res, output_folder)
    end do

end program main
