module sc_dynamics
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use dynamics
    use class_observables
    use variables
    use class_radical_pair
    use matrix_exponential
    use utils
    use m_random

    implicit none
    private 
    public :: semi_classical

    contains

    subroutine semi_classical(sim, sys, init_samples, res)
    ! Run semi-classical spin dynamics
        type(sys_param), intent(in)    :: sys 
        type(sim_param), intent(inout) :: sim 
        type(sc_system), intent(in)    :: init_samples(:) 
        type(observables), intent(out) :: res

        type(sc_system)       :: RP_samples(sim%N_samples) ! local copy of MC samples
        integer               :: N_steps
        real(dp)              :: w1(3)
        real(dp)              :: w2(3)
        real(dp)              :: k_bar                     ! Average recombination rate
        real(dp)              :: delta_k                   ! Difference in kS and kT
        real(dp)              :: K(16,16)                  ! Recombination operator
        real(dp)              :: J_exch(16,16)             ! Exchange operator
        real(dp), allocatable :: KJ_exp(:,:)               ! Exponential recombination operator
        real(dp)              :: int_S(sim%N_samples)
        real(dp)              :: int_T(sim%N_samples)
        real(dp)              :: P_S                       ! Singlet probabilities at time t
        real(dp)              :: P_T                       ! Triplet probability at time t

        integer :: i, j

        ! Generate local copy of initial RP samples
        do i=1, size(init_samples)
            call RP_samples(i)%copy(init_samples(i))
        end do

        N_steps = ceiling((sim%t_end+sim%dt)/sim%dt)
        call res%malloc(N_steps+1)

        ! Field vectors point in the z direction
        w1 = [0.0_dp, 0.0_dp, sys%e1%w]
        w2 = [0.0_dp, 0.0_dp, sys%e2%w]

        k_bar = (sys%kS + 3.0_dp*sys%kT) / 4.0_dp
        delta_k = (sys%kS - sys%kT) / 4.0_dp

        K = get_recomb_operator(k_bar, delta_k)
        J_exch = get_exchange_operator(sys%J)
        ! We get a combined decay operator, which is the exponent of K+R
        call expm(0.5_dp*sim%dt*(K+J_exch), KJ_exp)

        int_S = 0.0_dp
        int_T = 0.0_dp

        !$OMP PARALLEL DO SHARED(RP_samples, int_S, int_T_a, int_T_b, int_total_a,int_total_b) 
        do i=1, sim%N_samples
            int_S(i) = int_S(i) + RP_samples(i)%P0_S*RP_samples(i)%P0_S
            int_T(i) = int_T(i) + RP_samples(i)%P0_S*RP_samples(i)%P0_T
        end do
        !$OMP END PARALLEL DO

        ! We evaluate the probabilities at t=0
        res%P_S(1) = 4.0_dp*sum(int_S)/sim%N_samples 
        res%P_T(1) = 4.0_dp*sum(int_T)/sim%N_samples 

        ! No we let the system evolve, evaluating the probability correlations at each step
        do i = 1, N_steps

            print *, "Starting timestep ", i
            int_S = 0.0_dp
            int_T = 0.0_dp
            
            !$OMP PARALLEL DO SHARED (RP_samples, int_S, int_T_a, int_T_b, int_total_a,int_total_b)  PRIVATE(P_S, P_T)
            do j = 1, sim%N_samples
                ! Now we let the spin precess around each other
                call RP_samples(j)%propagate(w1, w2, sys%e1%a_iso, sys%e2%a_iso, KJ_exp, sim%dt)
                call RP_samples(j)%evaluate_probs(P_S, P_T)
                ! We collect the probability corelation funcitons of each sample seperately and then we sum it 
                int_S(j) = int_S(j) + (RP_samples(j)%P0_S)*P_S
                int_T(j) = int_T(j) + RP_samples(j)%P0_S*P_T
            end do
            !$OMP END PARALLEL DO

            res%P_S(i+1) = 4.0_dp*sum(int_S)/sim%N_samples 
            res%P_T(i+1) = 4.0_dp*sum(int_T)/sim%N_samples 

        end do    

        call res%get_kinetics(sim%dt, sys%kS, sys%kT)
        call res%output(sim%output_folder)

    end subroutine semi_classical

    function get_recomb_operator(k_avrg, delta_k) result(K)
        ! Function that constructs the recombination operator, K 
        real(dp) :: k_avrg, delta_k  
        real(dp) :: K(16, 16)

        integer :: i

        K = 0.0_dp

        ! All of the diagonal elements are k_avrg - corresponds to the underlying decay of each operator
        do i = 1,16
            K(i,i) = -k_avrg
        end do

        ! Then we add the elements that couple S1 and S2
        do i = 1,3
            K(i, i+3) = delta_k
            K(i+3, i) = delta_k
        end do

        ! Nowe we set the elements copuling pairs of symmetric off-diagonal elements of T12 (8-10, 9-13, 12-14)
        K(8,10) = delta_k
        K(10,8) = delta_k
        K(9,13) = delta_k
        K(9,13) = delta_k
        K(12,14) = delta_k
        K(14,12) = delta_k

        ! Then we add the coupling between the diagonal elements of the T12 are 7, 11 and 15 in the s_operators vector
        K(7,11) = -delta_k
        K(11,7) = -delta_k
        K(7,15) = -delta_k
        K(15,7) = -delta_k
        K(11,15) = -delta_k
        K(15,11) = -delta_k

        ! Finally we add the copulings between the diagonal elements and the identity
        K(7:15:4, 16) = 0.25_dp*delta_k
        K(16, 7:15:4) = 4*delta_k

    end function get_recomb_operator

    function get_exchange_operator(J_const) result(J)
        ! Function that constructs the exchange coupling operator, J
        real(dp), intent(in) :: J_const
        real(dp) :: twoJ, halfJ, J(16,16)

        ! We multiply J by 2 here for efficiency reasons
        twoJ = 2.0_dp*J_const 
        halfJ = 0.5_dp*J_const

        J = 0.0_dp

        ! Couplings of operator S1
        J(1, 12) = twoJ
        J(1, 14) = -twoJ
        J(2, 13) = twoJ
        J(2, 9) = -twoJ
        J(3, 8) = twoJ
        J(3, 10) = -twoJ

        ! Couplings of operator S2
        J(4, 12) = -twoJ
        J(4, 14) = twoJ
        J(5, 13) = -twoJ
        J(5, 9) = twoJ
        J(6, 8) = -twoJ
        J(6, 10) = twoJ

        ! Couplings of off-diagonal T12 elements (the diagonal ones commute with the exchange Hamiltonian)
        J(8, 3) = -halfJ
        J(8, 6) = halfJ
        J(9, 2) = halfJ
        J(9, 5) = -halfJ
        J(10, 3) = halfJ
        J(10, 6) = -halfJ
        J(12, 1) = -halfJ
        J(12, 4) = halfJ
        J(13, 2) = -halfJ
        J(13, 5) = halfJ
        J(14, 1) = halfJ
        J(14, 4) = -halfJ
        
    end function get_exchange_operator

    subroutine generate_init_samples(sys, sim, RP_samples)
        type(sys_param), intent(in)    :: sys 
        type(sim_param), intent(inout) :: sim  
        type(sc_system), intent(out)   :: RP_samples(sim%N_samples) 

        integer :: i

        !$OMP PARALLEL DO SHARED(RP_samples) 
        do i=1, sim%N_samples
            ! We randomly initialise each on sample
            call RP_samples(i)%initialise(sys, sim)
        end do
        !$OMP END PARALLEL DO
        print*, "Initial samples generated"

    end subroutine generate_init_samples

end module sc_dynamics