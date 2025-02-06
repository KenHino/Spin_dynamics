! module class_radical_pair
!     use, intrinsic :: iso_fortran_env, only: dp => real64
!     ! use system
!     ! use utils
!     ! use m_random
!     implicit none
    
!     private
!     public :: sc_var

!     type sc_var
!     ! Class handling dynamics of a semi-classical spin system
!         real(dp), dimension(:)                 :: S12(16) ! all of the electronic spin operators are kept in one vector so that we can define a recombination operator matrix
!         real(dp), dimension(:, :), allocatable :: I1, I2 ! Nuclear spins with electrons 1 and 2 
!         real(dp)                               :: P0_S, P0_T  ! Initial singlet/triplet probability
!         ! real(dp) :: P0_ST ! Mixed initial triplet and singlet populations

!         contains 
!         procedure :: initialise
!         ! procedure :: copy
!         ! procedure :: s1
!         ! procedure :: s2
!         ! procedure :: T12
!         ! procedure :: identity
!         ! procedure :: recombine
!         ! procedure :: evole_electron_spins
!         ! procedure :: evole_nuclear_spins
!         ! procedure :: propagate
!         ! procedure :: sum_i1
!         ! procedure :: sum_i2
!         ! procedure :: evaluate_probs
!         procedure :: deallocate_object
!     end type sc_var
        
!     contains
        
!         subroutine initialise(self, sys, theta, phi)
!             class(sc_var), intent(inout) :: self
!             real(dp)                     :: theta                     
!             real(dp)                     :: phi                     
!             ! integer, intent(in)          :: N_nuc1, N_nuc2
!             ! real(dp), intent(in)         :: lambda

!             integer :: i
            
!             ! Create a random number generator - code taken from https://github.com/jannisteunissen/rng_fortran
!             ! type(RNG_t) :: rng
            
!             ! ! Set the initial seed for the generator
!             ! call rng%set_random_seed()
            
!             ! ! Initialise all components of radical pair
!             ! self%N_nuc1 = N_nuc1
!             ! self%N_nuc2 = N_nuc2
            
!             ! Allocate all of the components of the system
!             ! allocate(self%s_operators(16))
!             ! allocate(self%i1(self%N_nuc1,3))
!             ! allocate(self%i2(self%N_nuc2,3))
            
!             ! ! The s_operators vector holds all operators of the electrons. The indices are as follows: S1: 1-3; S2: 4-6; T12: 7-15; Identity: 16 
            
!             ! ! Sample electronic spins from a sphere and save them in the first 6 elements of the component vector
!             ! self%s_operators(1:3) = rng%sphere(spin_magnitude) 
!             ! self%s_operators(4:6) = rng%sphere(spin_magnitude)
            
!             ! ! Now we sample the nuclear spins
!             ! do i = 1, self%N_nuc1
!             !     self%i1(i, :) = rng%sphere(spin_magnitude)
!             ! end do
!             ! do i = 1, self%N_nuc2
!             !     self%i2(i, :) = rng%sphere(spin_magnitude)
!             ! end do
            
!             ! ! Initially the T12 tensor is just the outer product of the electron spin vectors
!             ! self%s_operators(7:15) = reshape(outer_product(self%s_operators(1:3), self%s_operators(4:6)), [9]) 
            
!             ! ! Initially the identity is 1 and it will be the last element of the operator vector 
!             ! self%s_operators(16) = 1.0_dp
            
!             ! ! Save the initial values of the singlet/triplet probabilities for the integral evaluation
!             ! call self%evaluate_probs(self%P0_S, self%P0_T)

!             ! ! Initially we have can have a mixture of S and T so we mix the singlet/triplet populations
!             ! self%P0_ST = (1.0_dp-lambda)*self%P0_S + (lambda/3.0_dp)*self%P0_T
        
!         end subroutine initialise
        
!         ! subroutine copy(self, RP, lambda)
!         !     ! Function to create copy of a radical pair
!         !     class(radical_pair), intent(inout) :: self
!         !     class(radical_pair), intent(in) :: RP
!         !     real(dp), intent(in) :: lambda

!         !     integer :: i

!         !     self%N_nuc1 = RP%N_nuc1
!         !     self%N_nuc2 = RP%N_nuc2

!         !     allocate(self%s_operators(16))
!         !     allocate(self%i1(self%N_nuc1,3))
!         !     allocate(self%i2(self%N_nuc2,3))


!         !     ! Sample electronic spins from a sphere and save them in the first 6 elements of the component vector
!         !     self%s_operators(1:3) = RP%s_operators(1:3) 
!         !     self%s_operators(4:6) = RP%s_operators(4:6)
            
!         !     ! Now we sample the nuclear spins
!         !     do i = 1, self%N_nuc1
!         !         self%i1(i, :) = RP%i1(i, :)
!         !     end do
!         !     do i = 1, self%N_nuc2
!         !         self%i2(i, :) = RP%i2(i, :)
!         !     end do
            
!         !     ! Initially the T12 tensor is just the outer product of the electron spin vectors
!         !     self%s_operators(7:15) = reshape(outer_product(self%s_operators(1:3), self%s_operators(4:6)), [9]) 
            
!         !     ! Initially the identity is 1 and it will be the last element of the operator vector 
!         !     self%s_operators(16) = 1.0_dp
            
!         !     ! Save the initial values of the singlet/triplet probabilities for the integral evaluation
!         !     call self%evaluate_probs(self%P0_S, self%P0_T)

!         !     ! Initially we have can have a mixture of S and T so we mix the singlet/triplet populations
!         !     self%P0_ST = (1.0_dp-lambda)*self%P0_S + (lambda/3.0_dp)*self%P0_T
!         ! end subroutine copy
        

!         ! subroutine recombine(self, K_exp)
!         !     ! Applies recombination operator to electron spin operators of system
!         !     class(radical_pair), intent(inout) :: self
!         !     real(dp), dimension(16,16), intent(in) :: K_exp

!         !     self%s_operators = matmul(K_exp, self%s_operators)
        
!         ! end subroutine recombine
        
!         ! subroutine evole_electron_spins(self, w0, w1, w_rf, rf_dir, t, a1, a2, dt)
!         !     ! Subroutine which precesses electon spin operators around corresponding total magnetic field vector
!         !     class(radical_pair), intent(inout) :: self
!         !     real(dp), dimension(3), intent(in) :: w0
!         !     real(dp), intent(in) :: w1, w_rf, t
!         !     real(dp), dimension(:), intent(in) :: a1, a2
!         !     real(dp), intent(in) :: dt
!         !     character(1), intent(in) :: rf_dir


!         !     real(dp), dimension(3) :: w1_bar, w2_bar, w1_avrg


!         !     ! First calculate the total magnetic (static + hyperfine + oscillating)
!         !     w1_bar = w0 + self%sum_i1(a1) 
!         !     w2_bar = w0 + self%sum_i2(a2)

!         !     ! If there is an oscillating field we also add it
!         !     if (w1 /= 0.0_dp) then
!         !         if (rf_dir == 'x') then
!         !             w1_avrg = w1*(2.0_dp/(w_rf*dt))*[cos(w_rf*(t+dt*0.5_dp))*sin(w_rf*dt*0.5_dp), 0.0_dp, 0.0_dp]
!         !         end if
!         !         if (rf_dir == 'z') then
!         !             w1_avrg = w1*(2.0_dp/(w_rf*dt))*[0.0_dp, 0.0_dp, cos(w_rf*(t+dt*0.5_dp))*sin(w_rf*dt*0.5_dp)]
!         !         end if
!         !         w1_bar = w1_bar + w1_avrg
!         !         w2_bar = w2_bar + w1_avrg
!         !     end if

!         !     ! First we rorate the spin vectors
!         !     self%s_operators(1:3) = rotate(self%s1(), w1_bar, dt)
!         !     self%s_operators(4:6) = rotate(self%s2(), w2_bar, dt)
            
!         !     !7,8,9 ; 10,11,12; 13,14,15 are the colums of the T12 tensor - the numbers are indices in the s_operators vector
!         !     self%s_operators(7:9) = rotate(self%s_operators(7:9), w1_bar, dt)
!         !     self%s_operators(10:12) = rotate(self%s_operators(10:12), w1_bar, dt)
!         !     self%s_operators(13:15) = rotate(self%s_operators(13:15), w1_bar, dt)

!         !     ! 7,10,13 ; 8,11,14; 9,12,15 are the rows of the T12 tensor - the numbers are indices in the s_operators vector
!         !     self%s_operators(7:13:3) = rotate(self%s_operators(7:13:3), w2_bar, dt)
!         !     self%s_operators(8:14:3) = rotate(self%s_operators(8:14:3), w2_bar, dt)
!         !     self%s_operators(9:15:3) = rotate(self%s_operators(9:15:3), w2_bar, dt) 
        
!         ! end subroutine evole_electron_spins


!         ! subroutine evole_nuclear_spins(self, a1, a2, dt)
!         !     ! Subroutine which precesses nuclear spin vectors around corresponding electronic spin vector
!         !     class(radical_pair) :: self
!         !     real(dp), dimension(:), intent(in) :: a1, a2
!         !     real(dp), intent(in) :: dt

!         !     real(dp):: s1_norm, s2_norm ! Normalisation factors for electronic spin vectors
!         !     integer :: i

!         !     ! We need the factors to normalise the electronic spin vectors for the nuclear precession 
!         !     s1_norm = spin_magnitude/norm2(self%s1())
!         !     s2_norm = spin_magnitude/norm2(self%s2())

!         !     ! Then rotate nuclear spins around electron spin for a whole timestep
!         !     do i = 1, self%N_nuc1
!         !         self%i1(i, :) = rotate(self%i1(i, :), s1_norm*a1(i)*self%s1(), dt)
!         !     end do
!         !     do i = 1, self%N_nuc2
!         !         self%i2(i, :) = rotate(self%i2(i, :), s2_norm*a2(i)*self%s2(), dt)
!         !     end do

!         ! end subroutine
        

!         ! subroutine propagate(self, w0, w1, w_rf, rf_dir, t, a1, a2, dt, K_exp)
!         !     ! The spins are propagated for a timestep dt
!         !     class(radical_pair) :: self
!         !     real(dp), dimension(3), intent(in) :: w0
!         !     real(dp), intent(in) :: w1, w_rf, t
!         !     real(dp), dimension(:), intent(in) :: a1, a2
!         !     real(dp), intent(in) :: dt
!         !     real(dp), dimension(16,16), intent(in) :: K_exp
!         !     character(1), intent(in) :: rf_dir

            
!         !     ! We are gonna employ a split operator approach - recombination for dt/2 , precession for dt and recombination for dt/2 again

!         !     ! Apply recombination operator for half a timestep
!         !     call self%recombine(K_exp)
!         !     ! Now we will let the electronic and nuclear spins precess around each other
!         !     ! First precess electron operators for half a timestep
!         !     call self%evole_electron_spins(w0, w1, w_rf, rf_dir, t, a1, a2, dt*0.5_dp)
!         !     ! call self%evole_electron_spins(w0, a1, a2, dt*0.5_dp)
!         !     ! Then rotate nuclear spins around electron spin for a whole timestep
!         !     call self%evole_nuclear_spins(a1, a2, dt)
!         !     ! Repeat same evolution for electronic spins and T12
!         !     ! call self%evole_electron_spins(w0, a1, a2, dt*0.5_dp) 
!         !     call self%evole_electron_spins(w0, w1, w_rf, rf_dir, t+0.5_dp*dt, a1, a2, dt*0.5_dp) 
!         !     ! Apply recombination operator again
!         !     call self%recombine(K_exp)

!         ! end subroutine propagate

!         ! function s1(self) result(res)
!         !     ! Function that returns S1 vector
!         !     class(radical_pair) :: self
!         !     real(dp), dimension(3) :: res
!         !     res = self%s_operators(1:3)
!         ! end function s1

!         ! function s2(self) result(res)
!         !     ! Function that returns S2 vector
!         !     class(radical_pair) :: self
!         !     real(dp), dimension(3) :: res
!         !     res = self%s_operators(4:6)
!         ! end function s2

!         ! function T12(self) result(res)
!         !     ! Function that returns T12 tensor
!         !     class(radical_pair) :: self
!         !     real(dp), dimension(3,3) :: res
!         !     res = reshape(self%s_operators(7:15), [3,3])
!         ! end function T12
            
!         ! function identity(self) result(res)
!         !     ! Function that returns identity
!         !     class(radical_pair) :: self
!         !     real(dp) :: res
!         !     res = self%s_operators(16)
!         ! end function identity

!         ! subroutine evaluate_probs(self, P_S, P_T) 
!         !     ! Subroutine to evaluate the singlet/triplet probabilities at time t
!         !     class(radical_pair) :: self
!         !     real(dp), intent(out) :: P_S
!         !     real(dp), intent(out):: P_T 
    
!         !     P_S = 0.25_dp*self%identity() - tr(self%T12())
!         !     P_T = self%identity() - P_S

!         ! end subroutine evaluate_probs

!         ! function sum_i1(self,a1) result(I_tot)
!         !     ! Function to calculate total hyperfine field vector for electron 1
!         !     class(radical_pair) :: self
!         !     real(dp), dimension(:), intent(in) :: a1
!         !     real(dp), dimension(3) :: I_tot
!         !     integer :: i

!         !     I_tot = [0,0,0]
!         !     do i = 1, self%N_nuc1
!         !         I_tot = I_tot + self%i1(i, :) * a1(i)
!         !     end do

!         ! end function sum_i1

!         ! function sum_i2(self,a2) result(I_tot)
!         !     ! Function to calculate total hyperfine field vector for electron 2
!         !     class(radical_pair) :: self
!         !     real(dp), dimension(:), intent(in) :: a2
!         !     real(dp), dimension(3) :: I_tot
!         !     integer :: i

!         !     I_tot = [0,0,0]
!         !     do i = 1, self%N_nuc2
!         !         I_tot = I_tot + self%i2(i, :) * a2(i)
!         !     end do

!         ! end function sum_i2

!         subroutine deallocate_object(self)
!             class(sc_var), intent(inout) :: self

!             deallocate(self%s_operators)
!             deallocate(self%i1)
!             deallocate(self%i2)

!         end subroutine deallocate_object
    
! end module class_radical_pair