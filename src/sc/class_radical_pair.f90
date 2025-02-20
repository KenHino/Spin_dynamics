module class_radical_pair
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use variables
    use stdlib_linalg, only: trace, eye, cross_product
    use m_random
    use spinop, only: spin_magnitude, spin_half_mag
    implicit none
    
    private
    public :: sc_system

    type sc_system
    ! Class handling dynamics of a semi-classical spin system
        real(dp), allocatable :: S12(:)           ! Vector of classical electronic variables 
        real(dp), allocatable :: I1(:,:), I2(:,:) ! Nuclear spins with electrons 1 and 2 
        real(dp)              :: P0_S, P0_T       ! Initial singlet/triplet probability

        contains

        procedure :: initialise
        procedure :: copy
        procedure :: evole_electron_spins
        procedure :: evole_nuclear_spins
        procedure :: propagate
        procedure :: sum_i1
        procedure :: sum_i2
        procedure :: evaluate_probs
        procedure :: deallocate_object
    end type sc_system
        
    contains
        
        subroutine initialise(self, sys, sim)
            class(sc_system), intent(inout) :: self
            type(sys_param), intent(in)     :: sys 
            type(sim_param), intent(inout)  :: sim 

            type(RNG_t) :: rng
            integer     :: i            
            
            call rng%set_random_seed()
            
            allocate(self%S12(16))
            allocate(self%I1(size(sys%e1%g_I), 3))
            allocate(self%I2(size(sys%e2%g_I), 3))
            
            ! The s_operators vector holds all operators of the electrons. The indices are as follows: S1: 1-3; S2: 4-6; T12: 7-15; Identity: 16 
            
            ! Sample electronic spins from a sphere
            self%S12(1:3) = rng%sphere(spin_half_mag) 
            self%S12(4:6) = rng%sphere(spin_half_mag)
            
            ! Now we sample the nuclear spins
            do i = 1, size(self%I1, dim=1)
                self%i1(i, :) = rng%sphere(spin_magnitude(sys%e1%g_I(i)))
            end do
            do i = 1, size(self%I2, dim=1)
                self%i2(i, :) = rng%sphere(spin_magnitude(sys%e2%g_I(i)))
            end do
            
            ! Initially the T12 tensor is the outer product of the electron spin vectors
            self%S12(7:15) = reshape(outer_product(self%S12(1:3), self%S12(4:6)), [9]) 
            
            ! Initially the identity is 1 
            self%S12(16) = 1.0_dp
            
            ! Save the initial values of the singlet/triplet probabilities for the integral evaluation
            call self%evaluate_probs(self%P0_S, self%P0_T)

        end subroutine initialise
        
        subroutine copy(self, RP)
            ! Function to create copy of a radical pair
            class(sc_system), intent(inout) :: self
            class(sc_system), intent(in)    :: RP

            integer :: i

            self%S12 = RP%S12
            self%I1 = RP%I1
            self%I2 = RP%I2
            self%P0_S = RP%P0_S
            self%P0_T = RP%P0_T

        end subroutine copy
        
        subroutine evole_electron_spins(self, w1, w2, a1, a2, K_exp, dt)
            ! Subroutine which precesses electon spin operators around corresponding total magnetic field vector
            class(sc_system), intent(inout) :: self
            real(dp), intent(in)            :: w1(3)
            real(dp), intent(in)            :: w2(3)
            real(dp), intent(in)            :: a1(:)
            real(dp), intent(in)            :: a2(:)
            real(dp), intent(in)            :: K_exp(16, 16)
            real(dp), intent(in)            :: dt

            real(dp):: w1_bar(3), w2_bar(3)

            ! Apply recombination operator for half a timestep
            self%S12 = matmul(K_exp, self%S12)

            ! First calculate the total magnetic (static + hyperfine + oscillating)
            w1_bar = w1 + self%sum_i1(a1) 
            w2_bar = w2 + self%sum_i2(a2)

            ! Rotate S1 and S2
            self%S12(1:3) = rotate(self%S12(1:3), w1_bar, dt)
            self%S12(4:6) = rotate(self%S12(1:3), w2_bar, dt)
            
            !7,8,9 ; 10,11,12; 13,14,15 are the colums of the T12 tensor 
            self%S12(7:9) = rotate(self%S12(7:9), w1_bar, dt)
            self%S12(10:12) = rotate(self%S12(10:12), w1_bar, dt)
            self%S12(13:15) = rotate(self%S12(13:15), w1_bar, dt)

            ! 7,10,13 ; 8,11,14; 9,12,15 are the rows of the T12 tensor 
            self%S12(7:13:3) = rotate(self%S12(7:13:3), w2_bar, dt)
            self%S12(8:14:3) = rotate(self%S12(8:14:3), w2_bar, dt)
            self%S12(9:15:3) = rotate(self%S12(9:15:3), w2_bar, dt)

        end subroutine evole_electron_spins


        subroutine evole_nuclear_spins(self, a1, a2, dt)
            ! Subroutine which precesses nuclear spin vectors around corresponding electronic spin vector
            class(sc_system)     :: self
            real(dp), intent(in) :: a1(3), a2(3)
            real(dp), intent(in) :: dt

            real(dp) :: s1_norm
            real(dp) :: s2_norm
            integer  :: i

            ! We need the factors to normalise the electronic spin vectors for the nuclear precession 
            s1_norm = spin_half_mag/norm2(self%S12(1:3))
            s2_norm = spin_half_mag/norm2(self%S12(4:6))

            ! Then rotate nuclear spins around electron spin for a whole timestep
            do i = 1, size(self%I1, dim=1)
                self%i1(i, :) = rotate(self%i1(i, :), s1_norm*a1(i)*self%S12(1:3), dt)
            end do
            do i = 1, size(self%I2, dim=1)
                self%i2(i, :) = rotate(self%i2(i, :), s2_norm*a2(i)*self%S12(4:6), dt)
            end do

        end subroutine

        subroutine propagate(self, w1, w2, a1, a2, K_exp, dt)
            ! The spins are propagated for a timestep dt
            class(sc_system), intent(inout) :: self
            real(dp), intent(in)            :: w1(3)
            real(dp), intent(in)            :: w2(3)
            real(dp), intent(in)            :: a1(:)
            real(dp), intent(in)            :: a2(:)
            real(dp), intent(in)            :: K_exp(16, 16)
            real(dp), intent(in)            :: dt

            ! First precess electron operators for half a timestep
            call self%evole_electron_spins(w1, w2, a1, a2, K_exp, dt*0.5_dp)
            ! Then rotate nuclear spins around electron spin for a whole timestep
            call self%evole_nuclear_spins(a1, a2, dt)
            ! Repeat electron operator evolution
            call self%evole_electron_spins(w1, w2, a1, a2, K_exp, dt*0.5_dp)

        end subroutine propagate

        subroutine evaluate_probs(self, P_S, P_T) 
            ! Subroutine to evaluate the singlet/triplet probabilities at time t
            class(sc_system), intent(in) :: self
            real(dp), intent(out)        :: P_S
            real(dp), intent(out)        :: P_T 
    
            ! S12(16) = identity
            ! self%S12(7) + self%S12(10) + self%S12(13) = tr(T12)
            P_S = 0.25_dp*self%S12(16) - (self%S12(7) + self%S12(10) + self%S12(13))
            P_T = self%S12(16) - P_S

        end subroutine evaluate_probs

        function sum_i1(self, a1) result(I_tot)
            ! Function to calculate total hyperfine field vector for electron 1
            class(sc_system)     :: self
            real(dp), intent(in) :: a1(:)
            real(dp)             :: I_tot(3)
            
            integer :: i

            I_tot = [0,0,0]
            do i = 1, size(self%I1, dim=1)
                I_tot = I_tot + self%i1(i, :) * a1(i)
            end do

        end function sum_i1

        function sum_i2(self, a2) result(I_tot)
            ! Function to calculate total hyperfine field vector for electron 2
            class(sc_system)     :: self
            real(dp), intent(in) :: a2(:)
            real(dp)             :: I_tot(3)
            
            integer :: i

            I_tot = [0,0,0]
            do i = 1, size(self%I2, dim=1)
                I_tot = I_tot + self%i2(i, :) * a2(i)
            end do

        end function sum_i2

        subroutine deallocate_object(self)
            class(sc_system), intent(inout) :: self

            deallocate(self%S12)
            deallocate(self%I1)
            deallocate(self%I2)

        end subroutine deallocate_object

        function rotate(spin, field, dt) result(new_spin)
        ! Function that precesses spin around field for time dt
            implicit none
            real(dp), intent(in) :: spin(3)
            real(dp), intent(in) :: field(3)
            real(dp), intent(in) :: dt
            real(dp)             :: new_spin(3)
            
            real(dp) :: freq
            real(dp) :: field_hat(3)
            real(dp) :: spin_parallel(3)
            real(dp) :: spin_perpendicular(3)
            real(dp) :: spin_cross(3)

            ! Find field magnitude
            freq = norm2(field)
            
            ! In zero field a radical with no nuclear spins may experinece a zero total magnetic field which will cause the normalisation to explode
            if (freq /= 0.0_dp) then
                new_spin = spin ! When there is no magnetic field the electronic spin remains unchanged
            else
                ! Normalise field vector
                 field_hat = field/freq
            
                ! Decompose spin vector into components 
                spin_parallel = dot_product(spin, field_hat) * field_hat
                spin_perpendicular = spin - spin_parallel
                spin_cross = cross_product(field_hat, spin)

                ! Calculate  spin vector
                new_spin = spin_parallel + cos(freq*dt)*spin_perpendicular + sin(freq*dt)*spin_cross 
            
            end if 
    end function rotate

end module class_radical_pair
