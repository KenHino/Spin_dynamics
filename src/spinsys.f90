module system
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use utils
    use stdlib_linalg, only: trace
    use spinop
    use finer

    implicit none  
    private :: g_e, check_inp, is_isotropic
    public  :: spinsys,       &
               electron,      &
               sim_param,     & 
               sim_res,       & 
               spinsys_init,  & 
               electron_init, &
                

    complex(dp), parameter :: g_e = cmplx(2.002319_dp, 0.0_dp, kind=dp)

    type  electron
        complex(dp)              :: g         ! g-factor
        complex(dp)              :: w         ! larmor frequencies of electron, w =-g*B
        integer, allocatable     :: g_I(:)    ! spin multiplicities of nuclei coupled to electron
        complex(dp), allocatable :: A(:,:,:)  ! hyperfine tensors of all nuclei coupled to electron 1    
        logical                  :: isotropic ! true if all hyperfine couplings of system are isotropic
    end type electron
    
    type spinsys
        ! Contains all variables needed to specify a system's spin Hamiltonian
        complex(dp)    :: B            ! magnetic variables
        type(electron) :: e1, e2       ! elctronic paramters
        complex(dp)    :: J            ! exchange coupling
        complex(dp)    :: D(3,3)       ! dipolar tensor
    end type spinsys
    
    type sim_param
        integer                  :: N_sample       ! Total # of MC samples
        complex(dp), allocatable :: samples(:,:,:) ! Samples themselves (# sample, # nucleus, state parameters) 
        complex(dp)              :: dt             ! Integrator timestep
        integer                  :: N_steps        ! Total # of time steps
    end type sim_param

    contains
 
    subroutine initialise(filename, sys, sim) 
        ! Reads spin Hamiltonian paramters from input file
        character(len=50), intent(in) :: filename ! Input file
        type(spinsys), intent(out)    :: sys  
        type(sim_param), intent(out)  :: sim  

        type(file_ini) :: fini
        real(dp)       :: D_tmp(9)

        call fini%load(filename=trim(filename))

        call check_inp(fini)

        call fini%get(section_name='system variables', option_name='B', val=sys%B%re)
        sys%B%im = 0.0_dp
        call fini%get(section_name='system variables', option_name='J', val=sys%J%re)
        sys%J%im = 0.0_dp
        call fini%get(section_name='system variables', option_name='D', val=D_tmp)
        ! File parser cannot 2D arrays so we have to reshape the array we read
        sys%D = cmplx(reshape(D_tmp, [3,3]), 0.0_dp, kind=dp)

        call fini%get(section_name='simulation parameters', option_name='N_sample', val=sim%N_sample)
        call fini%get(section_name='simulation parameters', option_name='dt', val=sys%J%re)
        sym%dt%im = 0.0_dp
        call fini%get(section_name='simulation parameters', option_name='N_steps', val=sim%N_steps)

        sys%e1 = electron_init(fini, 1)
        sys%e2 = electron_init(fini, 2)

        sys%e1%isotropic = is_isotropic(sys%e1) 
        sys%e2%isotropic = is_isotropic(sys%e2) 

        ! Calculate Larmor frequencies of both electrons in mT units
        sys%e1%w = (sys%e1%g/g_e)*sys%B
        sys%e2%w = (sys%e2%g/g_e)*sys%B

    end subroutine spinsys_init

    function electron_init(fini, el_number) result(e)
        ! Reads electronic parameters from input file
        type(file_ini), intent(in) :: fini      ! File handler object
        integer, intent(in)        :: el_number ! must 
        type(electron)             :: e         ! Electron must be 1 or 2
 
        real(dp)          :: a_tmp(9)   ! Variable used for parsing individual hyperfines
        character(len=20) :: el_section 
        character(len=10) :: a_option 
        integer           :: err
        integer           :: i


        write(el_section, '(I0)') el_number
        el_section = 'electron '//el_section

        call fini%get(section_name=trim(el_section), option_name='g', val=e%g%re)
        e%g%im = 0.0_dp

        allocate(e%g_I(fini%count_values(section_name=trim(el_section), option_name='g_I')))
        call fini%get(section_name=trim(el_section), option_name='g_I',val=e%g_I, error=err)
        allocate(e%A(size(e%g_I), 3, 3))
        
        ! If no electron couplings are specified, the g_I and A arrays have random sizes, so we need to reallocate them to 0 
        if (err /= 0) then
            deallocate(e%g_I)
            deallocate(e%A)
            allocate(e%g_I(0))
            allocate(e%A(0,0,0))
        end if

        do i=1,size(e%g_I)
            ! Ai is the option name for each hyperfine
            write(a_option, '(I0)') i
            a_option = 'A'//a_option
            call fini%get(section_name=trim(el_section),option_name=trim(a_option),val=a_tmp, error=err)
            
            if (err /= 0) then
                print'(A,I0,A,I0, A)', 'Error: Hyperfine copuling A', i, ' for electron ', el_number, ' is missing.'
                stop
            end if
            
            ! File parser cannot 2D arrays so we have to reshape the array we read
            e%A(i,:,:) = cmplx(reshape(a_tmp, [3,3]), 0.0_dp, kind=dp)
        end do

    end function electron_init

    subroutine check_inp(fini)
    ! Checks if all necessary parameters are specified in the input file
        type(file_ini), intent(in) :: fini

        character(len=16) :: sections(3)
        integer           :: i

        sections = [character(len=16) :: 'system variables', 'electron 1', 'electron 2']

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

    function is_isotropic(e) 
    ! Checks if spin system contains only isotropic hyperfine couplings
    type(electron), intent(in) :: e
    logical                    :: is_isotropic 

    complex(dp) :: identity(3,3)
    complex(dp) :: trA
    integer :: i

    is_isotropic = .true.
    identity = eye(3)

    do i=1,size(e%g_I)
        trA = trace(e%A(i,:,:))

        if (any(e%A(i,:,:) /= (trA/3.0_dp)*identity)) then
            is_isotropic = .false.
        end if 
    end do

    end function is_isotropic

end module system