! module sc_dynamics
!     use, intrinsic :: iso_fortran_env, only: dp => real64
!     use dynamics
!     use system
!     use class_radical_pair
!     ! use matrix_exponential
!     ! use class_radical_pair
!     ! use utils
!     ! use m_random

!     implicit none
!     private 
!     public :: semi_classical

!     contains

!     subroutine semi_classical(sim, sys, res)
    
!         type(sys_param), intent(in) :: sys 
!         type(sim_param), intent(in) :: sim 
!         type(sim_res), intent(out)  :: res

!         type(sc_var), dimension(sim%N_samples) :: RP_samples ! local copy of MC samples
!         real(dp) :: k_bar, delta_k ! Recombination parameters
!         real(dp) :: K(16,16), R(16,16), J_exch(16,16), KRJ_exp(16,16) ! Exponential recombination operator


!         ! type(radical_pair), dimension(N_init)             :: RP_samples ! local copy of MC samples

!         k_bar = (sys%kS + 3.0_dp*sys%kT) / 4.0_dp
!         delta_k = (sys%kS - sys%kT) / 4.0_dp

!            ! Then we construct the recombination and relax operators
!         K = get_recomb_operator(k_bar, delta_k)
!         ! R = get_relax_operator(R1_1, R2_1, R1_2, R2_2)
!         J_exch = get_exchange_operator(sys%J)


!         ! Then we construct the recombination and relax operators
!         ! K = get_recomb_operator(k_avrg, delta_k)
!         ! R = get_relax_operator(R1_1, R2_1, R1_2, R2_2)
!         ! J_exch = get_exchange_operator(J_const)


!     end subroutine semi_classical

!     function get_recomb_operator(k_avrg, delta_k) result(K)
!         ! Function that constructs the recombination operator, K 
!         real(dp) :: k_avrg, delta_k  
!         real(dp) :: K(16, 16)

!         integer :: i

!         K = 0.0_dp

!         ! All of the diagonal elements are k_avrg - corresponds to the underlying decay of each operator
!         do i = 1,16
!             K(i,i) = -k_avrg
!         end do

!         ! Then we add the elements that couple S1 and S2
!         do i = 1,3
!             K(i, i+3) = delta_k
!             K(i+3, i) = delta_k
!         end do

!         ! Nowe we set the elements copuling pairs of symmetric off-diagonal elements of T12 (8-10, 9-13, 12-14)
!         K(8,10) = delta_k
!         K(10,8) = delta_k
!         K(9,13) = delta_k
!         K(9,13) = delta_k
!         K(12,14) = delta_k
!         K(14,12) = delta_k

!         ! Then we add the coupling between the diagonal elements of the T12 are 7, 11 and 15 in the s_operators vector
!         K(7,11) = -delta_k
!         K(11,7) = -delta_k
!         K(7,15) = -delta_k
!         K(15,7) = -delta_k
!         K(11,15) = -delta_k
!         K(15,11) = -delta_k

!         ! Finally we add the copulings between the diagonal elements and the identity
!         K(7:15:4, 16) = 0.25_dp*delta_k
!         K(16, 7:15:4) = 4*delta_k

!     end function get_recomb_operator

!         function get_relax_operator(R1_1, R2_1, R1_2, R2_2) result(R)
!         ! Function that constructs the phenomenological relaxation operator, R
!         real(dp), intent(in) :: R1_1, R2_1, R1_2, R2_2
!         real(dp) :: R(16,16)

!         R = 0.0_dp

!         ! The relaxation operator is just a diagonal operator that induces relaxation decay in the spin operators

!         ! First we set the relaxation times for the single electron operators (R1 - 1/T1 -> Sz decay, R2 - 1/T2 - Sx/Sy decay)
!         R(1,1) = -R2_1         
!         R(2,2) = -R2_1         
!         R(3,3) = -R1_1         
!         R(4,4) = -R2_2         
!         R(5,5) = -R2_2         
!         R(6,6) = -R1_2         
        
!         ! Now we need to set the relaxtion times for the T12 elements. For example S1xS2z decays with a rate of R2_1 + R1_2
!         R(7,7) =  -(R2_1 + R2_2)
!         R(8,8) =  -(R2_1 + R2_2)
!         R(9,9) = -(R1_1 + R2_2)
!         R(10,10) = -(R2_1 + R2_2)
!         R(11,11) =  -(R2_1 + R2_2)
!         R(12,12) =  -(R1_1 + R2_2)
!         R(13,13) =  -(R1_1 + R1_2)
!         R(14,14) =  -(R1_1 + R1_2)
!         R(15,15) =  -(R1_1 + R1_2)

!         ! Note that the identity opertor is not directly affected by spin relaxation

!     end function get_relax_operator 

!     function get_exchange_operator(J_const) result(J)
!         ! Function that constructs the exchange coupling operator, J
!         real(dp), intent(in) :: J_const
!         real(dp) :: twoJ, halfJ, J(16,16)

!         ! We multiply J by 2 here for efficiency reasons
!         twoJ = 2.0_dp*J_const 
!         halfJ = 0.5_dp*J_const

!         J = 0.0_dp

!         ! Couplings of operator S1
!         J(1, 12) = twoJ
!         J(1, 14) = -twoJ
!         J(2, 13) = twoJ
!         J(2, 9) = -twoJ
!         J(3, 8) = twoJ
!         J(3, 10) = -twoJ

!         ! Couplings of operator S2
!         J(4, 12) = -twoJ
!         J(4, 14) = twoJ
!         J(5, 13) = -twoJ
!         J(5, 9) = twoJ
!         J(6, 8) = -twoJ
!         J(6, 10) = twoJ

!         ! Couplings of off-diagonal T12 elements (the diagonal ones commute with the exchange Hamiltonian)
!         J(8, 3) = -halfJ
!         J(8, 6) = halfJ
!         J(9, 2) = halfJ
!         J(9, 5) = -halfJ
!         J(10, 3) = halfJ
!         J(10, 6) = -halfJ
!         J(12, 1) = -halfJ
!         J(12, 4) = halfJ
!         J(13, 2) = -halfJ
!         J(13, 5) = halfJ
!         J(14, 1) = halfJ
!         J(14, 4) = -halfJ

        
!     end function get_exchange_operator


!     ! subroutine generate_init_samples(N_init, N_nuc1, N_nuc2, lambda, RP_samples) 
!     !     integer, intent(in) :: N_nuc1, N_nuc2, N_init
!     !     real(dp), intent(in) :: lambda
!     !     type(radical_pair), dimension(N_init) :: RP_samples !  Monte Carlo samples

!     !     integer :: i

!     !     !$OMP PARALLEL DO SHARED(RP_samples) 
!     !     do i=1, N_init
!     !         ! We randomly initialise each on sample
!     !         call RP_samples(i)%initialise(N_nuc1, N_nuc2, lambda)
!     !         ! Then we create an off copy of it so that the on and off are run witrh exactly the same samples (eliminates statistical errors)
!     !     end do
!     !     !$OMP END PARALLEL DO
!     !     print*, "Initial samples generated"

!     ! end subroutine generate_init_samples


!     ! function get_recomb_operator(k_avrg, delta_k) result(K)
!     !     ! Function that constructs the recombination operator, K 
!     !     real(dp) :: k_avrg, delta_k  
!     !     real(dp) :: K(16, 16)

!     !     integer :: i

!     !     K = 0.0_dp

!     !     ! All of the diagonal elements are k_avrg - corresponds to the underlying decay of each operator
!     !     do i = 1,16
!     !         K(i,i) = -k_avrg
!     !     end do

!     !     ! Then we add the elements that couple S1 and S2
!     !     do i = 1,3
!     !         K(i, i+3) = delta_k
!     !         K(i+3, i) = delta_k
!     !     end do

!     !     ! Nowe we set the elements copuling pairs of symmetric off-diagonal elements of T12 (8-10, 9-13, 12-14)
!     !     K(8,10) = delta_k
!     !     K(10,8) = delta_k
!     !     K(9,13) = delta_k
!     !     K(9,13) = delta_k
!     !     K(12,14) = delta_k
!     !     K(14,12) = delta_k

!     !     ! Then we add the coupling between the diagonal elements of the T12 are 7, 11 and 15 in the s_operators vector
!     !     K(7,11) = -delta_k
!     !     K(11,7) = -delta_k
!     !     K(7,15) = -delta_k
!     !     K(15,7) = -delta_k
!     !     K(11,15) = -delta_k
!     !     K(15,11) = -delta_k

!     !     ! Finally we add the copulings between the diagonal elements and the identity
!     !     K(7:15:4, 16) = 0.25_dp*delta_k
!     !     K(16, 7:15:4) = 4*delta_k

!     ! end function get_recomb_operator

!     ! function get_relax_operator(R1_1, R2_1, R1_2, R2_2) result(R)
!     !     ! Function that constructs the phenomenological relaxation operator, R
!     !     real(dp), intent(in) :: R1_1, R2_1, R1_2, R2_2
!     !     real(dp) :: R(16,16)

!     !     R = 0.0_dp

!     !     ! The relaxation operator is just a diagonal operator that induces relaxation decay in the spin operators

!     !     ! First we set the relaxation times for the single electron operators (R1 - 1/T1 -> Sz decay, R2 - 1/T2 - Sx/Sy decay)
!     !     R(1,1) = -R2_1         
!     !     R(2,2) = -R2_1         
!     !     R(3,3) = -R1_1         
!     !     R(4,4) = -R2_2         
!     !     R(5,5) = -R2_2         
!     !     R(6,6) = -R1_2         
        
!     !     ! Now we need to set the relaxtion times for the T12 elements. For example S1xS2z decays with a rate of R2_1 + R1_2
!     !     R(7,7) =  -(R2_1 + R2_2)
!     !     R(8,8) =  -(R2_1 + R2_2)
!     !     R(9,9) = -(R1_1 + R2_2)
!     !     R(10,10) = -(R2_1 + R2_2)
!     !     R(11,11) =  -(R2_1 + R2_2)
!     !     R(12,12) =  -(R1_1 + R2_2)
!     !     R(13,13) =  -(R1_1 + R1_2)
!     !     R(14,14) =  -(R1_1 + R1_2)
!     !     R(15,15) =  -(R1_1 + R1_2)

!     !     ! Note that the identity opertor is not directly affected by spin relaxation

!     ! end function get_relax_operator 

!     ! function get_exchange_operator(J_const) result(J)
!     !     ! Function that constructs the exchange coupling operator, J
!     !     real(dp), intent(in) :: J_const
!     !     real(dp) :: twoJ, halfJ, J(16,16)

!     !     ! We multiply J by 2 here for efficiency reasons
!     !     twoJ = 2.0_dp*J_const 
!     !     halfJ = 0.5_dp*J_const

!     !     J = 0.0_dp

!     !     ! Couplings of operator S1
!     !     J(1, 12) = twoJ
!     !     J(1, 14) = -twoJ
!     !     J(2, 13) = twoJ
!     !     J(2, 9) = -twoJ
!     !     J(3, 8) = twoJ
!     !     J(3, 10) = -twoJ

!     !     ! Couplings of operator S2
!     !     J(4, 12) = -twoJ
!     !     J(4, 14) = twoJ
!     !     J(5, 13) = -twoJ
!     !     J(5, 9) = twoJ
!     !     J(6, 8) = -twoJ
!     !     J(6, 10) = twoJ

!     !     ! Couplings of off-diagonal T12 elements (the diagonal ones commute with the exchange Hamiltonian)
!     !     J(8, 3) = -halfJ
!     !     J(8, 6) = halfJ
!     !     J(9, 2) = halfJ
!     !     J(9, 5) = -halfJ
!     !     J(10, 3) = halfJ
!     !     J(10, 6) = -halfJ
!     !     J(12, 1) = -halfJ
!     !     J(12, 4) = halfJ
!     !     J(13, 2) = -halfJ
!     !     J(13, 5) = halfJ
!     !     J(14, 1) = halfJ
!     !     J(14, 4) = -halfJ

        
!     ! end function get_exchange_operator

! end module sc_dynamics