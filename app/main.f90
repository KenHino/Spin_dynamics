program main
    use, intrinsic :: iso_fortran_env, only: dp => real64, i8 => int64
    use spinop
    use hamiltonian
    use variables
    use utils
    use dynamics
    use stdlib_stats_distribution_normal
    use m_random
    use moments
    ! use class_radical_pair
    ! use sc_dynamics
    use class_observables
    use symmetry
    implicit none

    character(:), allocatable :: input_file
    character(:), allocatable :: output_folder
    character(len=600)        :: folder
    character(len=16)         :: guz
    type(sys_param)           :: sys
    type(sim_param)           :: sim
    type(RNG_t)               :: rng
    type(COO_cdp_type)        :: H_coo
    type(CSR_cdp_type)        :: H
    type(observables)         :: res 
    logical                   :: inpExists

    integer(i8) :: seed(2)
    integer :: i

    integer, dimension(3) :: start, finish, time
    
    
    call itime(start)
    print*, 'Program initialised'

    call get_command_argument(1, folder)
    input_file = trim(folder)

    inquire(file=input_file, exist=inpExists)
    if (.not. inpExists) stop 'Error: Input file does not exist.'

    ! input_file = '/home/sjoh5247/Spin_dynamics/input.ini'

    call read_inp(input_file, sys, sim)

    do i=1,size(sim%B)

        call rng%set_seed(sim%seed)

        ! Calculate Larmor frequencies of both electrons in mT units
        sys%e1%w = (sys%e1%g/g_e)*sim%B(i)
        sys%e2%w = (sys%e2%g/g_e)*sim%B(i)

        write(guz,'(F0.3)') sim%B(i) ! converting integer to string using a 'internal file'
        folder = sim%output_folder // '/' // trim(guz)  
        call system(' mkdir ' // folder)

        select case (sim%type)
            case('trace_sampling')
                call trace_sampling(sys, sim, rng, res, folder)
                ! call trace_sampling_para(sys, sim, rng, res, folder)
            case('symmetry_dynamics')
                print*, sim%B(i)
                call symmetrised_dynamics(sys, sim, rng, folder)
            case('semi_classical')
                print*, 'method not yet implemented'
            case default
                stop 'Error: Invalid simulation type entered.'
        end select
    end do

    call itime(finish)
    call get_time(start, finish, time)

    if (sim%type /= 'symmetry_dynamics') then
        print'(A15,I2,A1,I2,A1,I2)', 'Time elapsed:  ', time(1), ':',time(2), ':', time(3) 
    end if

end program main
