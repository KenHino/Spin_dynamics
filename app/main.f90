program main
    use, intrinsic :: iso_fortran_env, only: dp => real64, i8 => int64
    use spinop
    use hamiltonian
    use variables
    use class_observables
    use utils
    use dynamics
    use class_radical_pair
    use sc_dynamics
    use m_random
    use moments
    use symmetry
    implicit none

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
    logical                   :: fileExists

    integer(8) :: seed(2)
    integer :: i


    input_file = '/home/damianko/fpm/spinchem/data/cpf.ini'

    inquire(file=input_file, exist=fileExists)
    if (.not. fileExists) stop 'Error: Input file does not exist.'

    ! call get_command_argument(1, guz)
    ! input_file = trim(guz)
    ! call get_command_argument(2, folder)
    
    call read_inp(input_file, sys, sim)

    call rng%set_random_seed()
    seed = rng%s
    

    
    do i=1,size(sim%B)
!         print*, sim%B(i)

        call rng%set_seed(seed)

        ! Calculate Larmor frequencies of both electrons in mT units
        sys%e1%w = (sys%e1%g/g_e)*sim%B(i)
        sys%e2%w = (sys%e2%g/g_e)*sim%B(i)

        write(guz,'(F0.3)') sim%B(i) ! converting integer to string using a 'internal file'

        folder = sim%output_folder // '/' // trim(guz)  

        call system(' mkdir ' // folder)

        call symmetrised_dynamics(sys, sim, rng, res)
        ! call trace_sampling(sys, sim, rng, res, output_folder)
    end do

end program main
