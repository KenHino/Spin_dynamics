module class_observables
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    public :: observables
    private 

    type observables
        real(dp), allocatable :: P_S(:)   ! Time-resolved |S> RP state population
        real(dp), allocatable :: P_T0(:)  ! Time-resolved |T0> RP state population
        real(dp), allocatable :: P_Tp(:)  ! Time-resolved |T+> RP state population
        real(dp), allocatable :: P_Tm(:)  ! Time-resolved |T-> RP state population
        real(dp), allocatable :: P_T(:)   ! Time-resolved |T-> RP state population
        real(dp), allocatable :: iden(:)  ! Trace of density operator - proportional to RP concentration
        real(dp), allocatable :: Phi_S(:) ! Singlet product yield
        real(dp), allocatable :: Phi_T(:) ! Triplet product yield
        
        contains
        procedure :: malloc
        procedure :: set
        procedure :: scale
        procedure :: update
        procedure :: get_kinetics
        procedure :: output

    end type observables


    contains

    subroutine malloc(self, N_steps)
        class(observables), intent(inout) :: self
        integer, intent(in)          :: N_steps

        allocate(self%P_S(N_steps))
        allocate(self%P_T0(N_steps))
        allocate(self%P_Tp(N_steps))
        allocate(self%P_Tm(N_steps))
        allocate(self%P_T(N_steps))
        allocate(self%iden(N_steps))
        allocate(self%Phi_S(N_steps/2))
        allocate(self%Phi_T(N_steps/2))

    end subroutine malloc

    subroutine set(self, N)
        class(observables), intent(inout) :: self
        real(dp), intent(in)              :: N

        self%P_S = N 
        self%P_Tp = N 
        self%P_T0 = N 
        self%P_Tm = N
        self%P_T = N

    end subroutine set

    subroutine scale(self, N) 
        class(observables), intent(inout) :: self
        real(dp), intent(in)              :: N

        self%P_S = self%P_S*N 
        self%P_Tp = self%P_Tp*N 
        self%P_T0 = self%P_T0*N 
        self%P_Tm = self%P_Tm*N 

    end subroutine scale

    subroutine update(self, obs) 
        class(observables), intent(inout) :: self
        class(observables), intent(in)    :: obs

        self%P_S = self%P_S + obs%P_S 
        self%P_Tp = self%P_Tp + obs%P_Tp 
        self%P_T0 = self%P_T0 + obs%P_T0 
        self%P_Tm = self%P_Tm + obs%P_Tm 
    
    end subroutine update

    subroutine get_kinetics(self, dt, kS, kT)

        class(observables), intent(inout) :: self
        real(dp), intent(in)              :: dt
        real(dp), intent(in)              :: kS, kT

        integer :: i

        self%Phi_S(1) = 0.0_dp 
        self%Phi_T(1) = 0.0_dp

        ! Evaluate yields with Simpson's rule
        do i = 1,size(self%Phi_S)-1
            self%Phi_S(i+1) = self%Phi_S(i) + (1.0_dp/3.0_dp)*dt*kS*(self%P_S(2*i-1) + 4*self%P_S(2*i) + self%P_S(2*i+1))
            self%Phi_T(i+1) = self%Phi_T(i) + (1.0_dp/3.0_dp)*dt*kT*(self%P_T(2*i-1) + 4*self%P_T(2*i) + self%P_T(2*i+1))
        end do

    end subroutine get_kinetics

    subroutine output(self, path)
        ! Subroutine saves all data from the time propagation of a tripletpair a_at a certain  field
        class(observables), intent(inout)        :: self
        character(:), allocatable, intent(inout) :: path ! Path to data Pbar_T_a

        character(:), allocatable :: filename

        filename = 'S_prob.data'
        call write_data(self%P_S, path, filename)
        filename = 'Tp_prob.data'
        call write_data(self%P_Tp, path, filename)
        filename = 'T0_prob.data'
        call write_data(self%P_T0, path, filename)
        filename = 'Tm_prob.data'
        call write_data(self%P_Tm, path, filename)
        filename = 'total_prob.data'
        call write_data(self%iden, path, filename)
        filename = 'S_yield.data'
        call write_data(self%Phi_S, path, filename)
        filename = 'T_yield.data'
        call write_data(self%Phi_T, path, filename)


    end subroutine output

    subroutine write_data(data, path, filename)
        ! Subroutine that saves an array to a file path/filename
        real(dp), dimension(:), intent(in) :: data
        character(:), allocatable, intent(inout) :: path, filename

        character(:), allocatable :: kur

        integer :: i

        ! We remove any trailing blanks from the file path
        path = path
        filename = filename

        open(1, file= path // '/' // filename)
        do i=1,size(data)
            write (1,*) data(i)
        end do

    end subroutine write_data

end module class_observables