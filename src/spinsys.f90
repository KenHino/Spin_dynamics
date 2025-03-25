module variables
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use utils
    use stdlib_linalg, only: trace, eye
    use finer

    implicit none  
    private ::  check_inp, is_isotropic
    public  :: g_e, gamma_e, &
               sys_param,   &
               electron,    &
               sim_param,   & 
               read_inp,    & 
               electron_init
                
    real(dp), parameter :: g_e = 2.002319_dp
    real(dp), parameter :: gamma_e = 176.0_dp 

    type  electron
        ! Parameters of one electron spin Hamiltonian 
        real(dp)              :: g         ! g-factor
        real(dp)              :: w         ! Larmor frequency 
        integer, allocatable  :: g_I(:)    ! coupled nuclei spin multiplicities
        real(dp), allocatable :: A(:,:,:)  ! hyperfine copuling tensors    
        real(dp), allocatable :: a_iso(:)
        logical               :: isotropic ! true if all hyperfine couplings of electron are isotropic
    end type electron

    type sys_param
        ! Parameters of total spin Hamiltonian and recombination operator
        type(electron) :: e1, e2 ! electronic paramters
        real(dp)       :: J      ! exchange coupling
        real(dp)       :: D(3,3) ! dipolar tensor
        real(dp)       :: kS, kT ! recombination rate constants
        integer        :: Z1, Z2 ! sizes of individual nuclear Hilbert spaces
    end type sys_param
    
    type sim_param
        real(dp), allocatable     :: B(:)
        character(:), allocatable :: init_state       ! spin state in which radical pair is formed
        integer                   :: N_samples        ! # of Monte Carlo samples
        complex(dp), allocatable  :: SUZ_samples(:,:) ! S(UZ) Monte Carlo samples
        real(dp)                  :: t_end            ! duration of simulation
        real(dp)                  :: dt               ! integrator timestep
        integer                   :: N_krylov         ! size of Krylov space
        real(dp)                  :: tol              ! tolerance for recalculating the Krylov subspace
        integer                   :: M1               ! Number of symmetry blocks
        integer                   :: M2               ! Number of symmetry blocks
        real(dp)                  :: block_tol        ! Tolerance for discarding symmetry blocks
        character(:), allocatable :: output_folder    ! spin state in which radical pair is formed
    end type sim_param

    contains
 
    subroutine read_inp(filename, sys, sim) 
        ! Reads spin Hamiltonian paramters from input file
        character(:), allocatable, intent(in) :: filename ! Input file
        type(sys_param), intent(out)   :: sys  
        type(sim_param), intent(out)   :: sim  

        real(dp)            :: D_tmp(9)
        character(len=1000) :: tmp
        character(len=500)  :: tmp2
        type(file_ini)      :: fini
        logical             :: dirExists
        integer             :: err
        
        integer :: i

        call fini%load(filename=trim(filename))

        call check_inp(fini)

        call fini%get(section_name='system variables', option_name='J', val=sys%J, error=err)
        if (err /= 0) stop 'Error: Exchange coupling is missing'
        call fini%get(section_name='system variables', option_name='D', val=D_tmp)
        if (err /= 0) stop 'Error: Dipolar coupling is missing'
        sys%D = reshape(D_tmp, [3,3])
        call fini%get(section_name='system variables', option_name='kS', val=sys%kS)
        if (err /= 0) stop 'Error: Singlet recombination rate is missing'
        call fini%get(section_name='system variables', option_name='kT', val=sys%kT)
        if (err /= 0) stop 'Error: Triplet recombination rate is missing'
        sys%kS = sys%kS/gamma_e
        sys%kT = sys%kT/gamma_e

        allocate(sim%B(fini%count_values(section_name='simulation parameters', option_name='B')))
        call fini%get(section_name='simulation parameters', option_name='B', val=sim%B, error=err)
        if (err /= 0) stop 'Error: List of experiment magnetic fields is missing'
        call fini%get(section_name='simulation parameters', option_name='initial_state', val=tmp, error=err)
        sim%init_state = trim(tmp)
        if (err /= 0) stop 'Error: Initial state is missing'

        call fini%get(section_name='simulation parameters', option_name='dt', val=sim%dt, error=err)
        if (err /= 0) stop 'Error: Propagation timestep is missing'


        call fini%get(section_name='simulation parameters', option_name='N_samples', val=sim%N_samples)
        call fini%get(section_name='simulation parameters', option_name='simulation_time', val=sim%t_end)
        sim%t_end = (sim%t_end*gamma_e)/1000.0_dp
        sim%dt = (sim%dt*gamma_e)/1000.0_dp
        call fini%get(section_name='simulation parameters', option_name='block_tolerance', val=sim%block_tol)
        call fini%get(section_name='simulation parameters', option_name='N_krylov', val=sim%N_krylov)
        call fini%get(section_name='simulation parameters', option_name='integrator_tolerance', val=sim%tol)
        call fini%get(section_name='simulation parameters', option_name='M1', val=sim%M1)
        call fini%get(section_name='simulation parameters', option_name='M2', val=sim%M2)

        call fini%get(section_name='simulation parameters', option_name='output_folder', val=tmp2)
        sim%output_folder = trim(tmp2)
        inquire(file=sim%output_folder // '/.', exist=dirExists)
        if (.not. dirExists) stop 'Error: Output folder does not exist.'


        call electron_init(fini, 1, sys%e1)
        call electron_init(fini, 2, sys%e2)
        sys%Z1 = product(sys%e1%g_I) 
        sys%Z2 = product(sys%e2%g_I) 
        
        sys%e1%isotropic = is_isotropic(sys%e1) 
        sys%e2%isotropic = is_isotropic(sys%e2) 
   
    end subroutine read_inp

    subroutine check_inp(fini)
    ! Checks if all necessary parameters are specified in the input file
        type(file_ini), intent(in) :: fini

        character(len=24) :: sections(4)
        integer           :: i

        sections = [character(len=24) :: 'system variables', &
                                         'simulation parameters',&
                                         'electron 1', 'electron 2']

        if (fini%Ns < size(sections)) then
            print*, 'Error: Number of sections in input file is less than expected.'
            stop
        else if (fini%Ns > size(sections)) then
            print*, 'Error: Number of sections in input file is more than expected.'
            stop
        end if

        do i=1,size(sections)
            if (.not. fini%has_section(section_name=trim(sections(i)))) then
                print*, 'Error: Section "', trim(sections(i)), '" is missing from input file.'
                stop
            end if
        end do

    end subroutine check_inp

    subroutine electron_init(fini, el_number, e) 
        ! Reads electronic parameters from input file
        type(file_ini), intent(in)  :: fini      ! File handler object
        integer, intent(in)         :: el_number ! must 
        type(electron), intent(out) :: e         ! Electron must be 1 or 2
 
        integer, allocatable :: mult_I(:)
        integer, allocatable :: N_I(:)
        real(dp)             :: A_tmp(9)   ! Variable used for parsing individual hyperfines
        character(len=20)    :: el_section 
        character(len=10)    :: a_option 
        integer              :: err
        integer              :: start, finish
        integer              :: i

        write(el_section, '(I0)')  el_number
        el_section = 'electron '//el_section

        call fini%get(section_name=trim(el_section), option_name='g', val=e%g, error=err)
        if (err /= 0) then
            print'(A,A,I0, A)', 'Error: g-factor A'' for electron ', el_number, ' is missing.'
            stop
        end if
       
        allocate(mult_I(fini%count_values(section_name=trim(el_section), option_name='I')))
        call fini%get(section_name=trim(el_section), option_name='I',val=mult_I, error=err)

        allocate(N_I(fini%count_values(section_name=trim(el_section), option_name='N_I')))
        call fini%get(section_name=trim(el_section), option_name='N_I',val=N_I, error=err)

        ! allocate(e%g_I(fini%count_values(section_name=trim(el_section), option_name='g_I')))
        ! call fini%get(section_name=trim(el_section), option_name='g_I',val=e%g_I, error=err)
       
        ! If no electron couplings are specified, the g_I and A arrays have random sizes, so we need to reallocate them to 0 
        if (err /= 0) then
            deallocate(mult_I)
            deallocate(N_I)
            allocate(N_I(0))
            allocate(mult_I(0))
        end if

        if (size(N_I) /= size(mult_I)) stop 'Sections N_I and I must have the sam length.'

        allocate(e%g_I(sum(N_I)))
        allocate(e%A(size(e%g_I), 3, 3))
        allocate(e%a_iso(size(e%g_I)))

        if (size(e%g_I) == 0) return

        print*, sum(N_I)

        start = 1
        finish = N_I(1)
        do i=1,size(N_I)-1
                e%g_I(start:finish) = mult_I(i)
                start = start + N_I(i)
                finish = finish + N_I(i+1)
        end do
        e%g_I(start:finish) = mult_I(size(N_I))

        do i=1,size(e%g_I)
            write(a_option, '(I0)') i
            a_option = 'A'//a_option
            call fini%get(section_name=trim(el_section),option_name=trim(a_option),val=A_tmp, error=err)
            
            if (err /= 0) then
                print'(A,I0,A,I0, A)', 'Error: Hyperfine copuling A', i, ' for electron ', el_number, ' is missing.'
                stop
            end if
            
            e%A(i,:,:) = reshape(A_tmp, [3,3])
            e%a_iso(i) = (trace(e%A(i,:,:)))/3.0_dp
        end do

    end subroutine electron_init

    function is_isotropic(e) 
    ! Checks if spin system contains only isotropic hyperfine couplings
    type(electron), intent(in) :: e
    logical                    :: is_isotropic 

    real(dp)    :: identity(3,3)
    real(dp)    :: trA
    
    integer :: i

    is_isotropic = .true.
    identity = eye(3)

    do i=1,size(e%g_I)
        if(e%A(i,1,1) /= e%A(i,2,2)) is_isotropic = .false.
    end do

    end function is_isotropic

end module variables