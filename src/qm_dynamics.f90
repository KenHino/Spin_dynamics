module dynamics
    use, intrinsic :: iso_fortran_env, only: dp => real64, i8 => int64
    use variables
    use kroenecker
    use class_observables
    use utils
    use stdlib_sparse
    use stdlib_stats_distribution_normal
    use hamiltonian
    use matrix_exponential
    use m_random

    implicit none

    public :: exact_dynamics_dense, &
              exact_dynamics, &
              trace_sampling, &
              trace_sampling_para, &
              sample_SUZ
    private

    contains

    subroutine exact_dynamics(sys, sim, res, out)
    ! Run quantum mechanical dynamics with full nuclear Hilbert space
        type(sys_param), intent(in)              :: sys
        type(sim_param), intent(inout)           :: sim
        type(observables), intent(out)           :: res
        character(len=600), optional, intent(inout) :: out


        ! type(observables), allocatable :: sample_obs(:)
        type(COO_cdp_type)             :: H_coo
        type(CSR_cdp_type)             :: H
        type(observables)              :: sample_obs
        integer                        :: N_steps
        complex(dp)                    :: theta(4) ! init electronic wavefunction
        complex(dp), allocatable       :: Z_ket(:) ! init nuclear wavefunction
        complex(dp), allocatable       :: psi(:)   ! theta X Z_ket
        integer                        :: Z        ! size of nuclear Hilbert space
        integer                        :: N        ! size of combined Hilbert space
        real(dp)                       :: t
        real(dp)                       :: init_pop(4)
        complex(dp)                    :: c(sim%N_krylov)
        real(dp)                       :: c2
        real(dp)                       :: c_m2
        complex(dp)                    :: A(sim%N_krylov, sim%N_krylov)
        complex(dp), allocatable       :: A_exp(:, :)
        complex(dp)                    :: P_S(sim%N_krylov, sim%N_krylov)
        complex(dp)                    :: P_Tp(sim%N_krylov, sim%N_krylov)
        complex(dp)                    :: P_T0(sim%N_krylov, sim%N_krylov)
        complex(dp)                    :: P_Tm(sim%N_krylov, sim%N_krylov)
        complex(dp), allocatable       :: Q(:, :)

        integer  :: i, j, k

        Z = sys%Z1*sys%Z2
        N = 4*Z

        if (sys%e1%isotropic .and. sys%e2%isotropic) then
            call getHiso_sp(sys, H_coo)
            call coo_to_csr(H_coo, H)
            ! print*, 'Hamiltonian generated.'
        else
            stop 'Anisotropic Hamiltoninans have not been implemented yet'
        end if

        call init_el_state(sim, theta, init_pop)

        N_steps = ceiling((sim%t_end+sim%dt)/sim%dt)
        call res%malloc(N_steps+1)
        call res%set(0.0_dp)
        ! allocate(sample_obs(Z))
        call sample_obs%malloc(N_steps+1)


        allocate(Z_ket(Z))
        allocate(psi(4*Z))
        allocate(Q(N, sim%N_krylov))


        do i=1,Z
            ! print*, 'start state', i
            Z_ket = (0.0_dp, 0.0_dp)
            Z_ket(i) = 1.0_dp
            psi = kron_vector(theta, Z_ket)
            ! call sample_obs(i)%malloc(N_steps+1)
            call sample_obs%set(0.0_dp)
            sample_obs%P_S(1) = init_pop(1)
            sample_obs%P_T0(1) = init_pop(2)
            sample_obs%P_Tp(1) = init_pop(3)
            sample_obs%P_Tm(1) = init_pop(4)

            t = 0.0_dp
            k = 2

            do while (t < sim%t_end)

                call get_krylov(N, psi, H, sim%N_krylov, Q, c, A)
                call expm(((0.0_dp, -1.0_dp)*sim%dt*A), A_exp)
                call convert_to_krylov(sim%N_krylov, Q, Z, P_S, P_T0, P_Tp, P_Tm)

                c2 = 1.0_dp
                c_m2 = 0.0_dp

                do while (c_m2 <= sim%tol*c2 .and. t < sim%t_end)
                    c = matmul(A_exp, c)

                    ! Collect observables
                    sample_obs%P_S(k) = real(dot_product(c, matmul(P_S, c)), kind=dp)
                    sample_obs%P_Tp(k) = real(dot_product(c, matmul(P_Tp, c)), kind=dp)
                    sample_obs%P_T0(k) = real(dot_product(c, matmul(P_T0, c)), kind=dp)
                    sample_obs%P_Tm(k) = real(dot_product(c, matmul(P_Tm, c)), kind=dp)

                    c2 = real(dot_product(c,c), kind=dp)
                    c_m2 = real(c(sim%N_krylov)*conjg(c(sim%N_krylov)), kind=dp)

                    t = t+sim%dt
                    k = k+1
                end do

                psi = matmul(Q, c)

            end do

            call res%update(sample_obs)

        end do


        call res%scale(1.0_dp/real(Z, kind=dp))
        if (present(out)) then
            ! call res%get_kinetics(sim%dt, sys%kS, sys%kT)
            call res%output(out)
        end if

    end subroutine exact_dynamics

    subroutine trace_sampling(sys, sim, rng, res, out)
    ! Run quantum mechanical dynamics with trace sampling
        type(sys_param), intent(in)                        :: sys
        type(sim_param), intent(inout)                     :: sim
        type(RNG_t), intent(inout)                         :: rng
        type(observables), intent(out)                     :: res
        character(len=600), optional, intent(inout)        :: out


        type(COO_cdp_type)       :: H_coo
        type(CSR_cdp_type)       :: H
        integer                  :: N_steps
        complex(dp)              :: theta(4) ! init electronic wavefunction
        complex(dp), allocatable :: Z_ket(:) ! init nuclear wavefunction
        complex(dp), allocatable :: psi(:)   ! theta X Z_ket
        integer                  :: Z        ! size of nuclear Hilbert space
        integer                  :: N        ! size of combined Hilbert space
        real(dp)                 :: t
        real(dp)                 :: init_pop(4)
        complex(dp)              :: c(sim%N_krylov)
        real(dp)                 :: c2
        real(dp)                 :: c_m2
        complex(dp)              :: A(sim%N_krylov, sim%N_krylov)
        complex(dp), allocatable :: A_exp(:, :)
        complex(dp)              :: P_S(sim%N_krylov, sim%N_krylov)
        complex(dp)              :: P_Tp(sim%N_krylov, sim%N_krylov)
        complex(dp)              :: P_T0(sim%N_krylov, sim%N_krylov)
        complex(dp)              :: P_Tm(sim%N_krylov, sim%N_krylov)
        complex(dp), allocatable :: Q(:, :)
        type(observables)        :: sample_obs
        ! type(observables), allocatable :: sample_obs(:)

        integer  :: i, j, k

        Z = sys%Z1*sys%Z2
        N = 4*Z

        if (sys%e1%isotropic .and. sys%e2%isotropic) then
            call getHiso_sp(sys, H_coo)
            call coo_to_csr(H_coo, H)
            ! print*, 'Hamiltonian generated.'
        else
            stop 'Anisotropic Hamiltoninans have not been implemented yet'
        end if

        call init_el_state(sim, theta, init_pop)

        N_steps = ceiling((sim%t_end+sim%dt)/sim%dt)
        call res%malloc(N_steps+1)
        call res%set(0.0_dp)
        call sample_obs%malloc(N_steps+1)


        allocate(Z_ket(Z))
        allocate(psi(4*Z))
        allocate(Q(N, sim%N_krylov))

        ! allocate(sample_obs(sim%N_samples))

        do i=1, sim%N_samples
            call sample_SUZ(rng, Z_ket)
            psi = kron_vector(theta, Z_ket)
            ! call sample_obs(i)%malloc(N_steps+1)
            call sample_obs%set(0.0_dp)
            sample_obs%P_S(1) = init_pop(1)
            sample_obs%P_T0(1) = init_pop(2)
            sample_obs%P_Tp(1) = init_pop(3)
            sample_obs%P_Tm(1) = init_pop(4)

            t = 0.0_dp
            k = 2

            ! print*, 'starting sample ', i
            do while (t < sim%t_end)

                ! print*, 'Generating Krylov subspace at time t =  ', t

                call get_krylov(N, psi, H, sim%N_krylov, Q, c, A)
                call expm(((0.0_dp, -1.0_dp)*sim%dt*A), A_exp)
                call convert_to_krylov(sim%N_krylov, Q, Z, P_S, P_T0, P_Tp, P_Tm)

                c2 = 1.0_dp
                c_m2 = 0.0_dp

                do while (c_m2 <= sim%tol*c2 .and. t < sim%t_end)
                    c = matmul(A_exp, c)

                    ! Collect observables
                    sample_obs%P_S(k) = real(dot_product(c, matmul(P_S, c)), kind=dp)
                    sample_obs%P_Tp(k) = real(dot_product(c, matmul(P_Tp, c)), kind=dp)
                    sample_obs%P_T0(k) = real(dot_product(c, matmul(P_T0, c)), kind=dp)
                    sample_obs%P_Tm(k) = real(dot_product(c, matmul(P_Tm, c)), kind=dp)

                    c2 = real(dot_product(c,c), kind=dp)
                    c_m2 = real(c(sim%N_krylov)*conjg(c(sim%N_krylov)), kind=dp)

                    t = t+sim%dt
                    k = k+1
                end do
                psi = matmul(Q, c)
            end do
            call res%update(sample_obs)

        end do

        ! do i=1,sim%N_samples
            ! call res%update(sample_obs)
        ! end do

        ! Average observables
        call res%scale(1.0_dp/real(sim%N_samples, kind=dp))
        if (present(out)) then
            ! call res%get_kinetics(sim%dt, sys%kS, sys%kT)
            call res%output(out)
        end if

    end subroutine trace_sampling

    subroutine trace_sampling_para(sys, sim, rng, res, out)
        ! Run quantum mechanical dynamics with trace sampling
            type(sys_param), intent(in)                 :: sys
            type(sim_param), intent(inout)              :: sim
            type(RNG_t), intent(inout)                  :: rng
            type(observables), intent(out)              :: res
            character(len=600), optional, intent(inout) :: out


            type(COO_cdp_type)       :: H_coo
            type(CSR_cdp_type)       :: H
            integer                  :: N_steps
            complex(dp)              :: theta(4) ! init electronic wavefunction
            complex(dp), allocatable :: Z_ket(:) ! init nuclear wavefunction
            complex(dp), allocatable :: psi(:)   ! theta X Z_ket
            integer                  :: Z        ! size of nuclear Hilbert space
            integer                  :: N        ! size of combined Hilbert space
            real(dp)                 :: t
            real(dp)                 :: init_pop(4)
            complex(dp)              :: c(sim%N_krylov)
            real(dp)                 :: c2
            real(dp)                 :: c_m2
            complex(dp)              :: A(sim%N_krylov, sim%N_krylov)
            complex(dp), allocatable :: A_exp(:, :)
            complex(dp)              :: P_S(sim%N_krylov, sim%N_krylov)
            complex(dp)              :: P_Tp(sim%N_krylov, sim%N_krylov)
            complex(dp)              :: P_T0(sim%N_krylov, sim%N_krylov)
            complex(dp)              :: P_Tm(sim%N_krylov, sim%N_krylov)
            complex(dp), allocatable :: Q(:, :)
            ! type(observables)        :: sample_obs
            type(observables), allocatable :: sample_obs(:)
            integer(i8), allocatable       :: seeds(:,:)

            integer  :: i, j, k

            Z = sys%Z1*sys%Z2
            N = 4*Z

            if (sys%e1%isotropic .and. sys%e2%isotropic) then
                call getHiso_sp(sys, H_coo)
                call coo_to_csr(H_coo, H)
                ! print*, 'Hamiltonian generated.'
            else
                stop 'Anisotropic Hamiltoninans have not been implemented yet'
            end if

            call init_el_state(sim, theta, init_pop)

            N_steps = ceiling((sim%t_end+sim%dt)/sim%dt)
            call res%malloc(N_steps+1)
            call res%set(0.0_dp)

            allocate(Z_ket(Z))
            allocate(psi(4*Z))
            allocate(Q(N, sim%N_krylov))

            allocate(sample_obs(sim%N_samples))

            allocate(seeds(sim%N_samples, 2))

            do i=1,size(seeds, dim=1)
                seeds(i,1) = rng%next()
                seeds(i,2) = rng%next()
            end do

            !$OMP PARALLEL DO SHARED(sys,sim,init_pop,N, sample_obs)&
            !$OMP& PRIVATE(t,k,c2,c_m2,c,psi,Q, P_S, P_T0, P_Tp, P_Tm,A,A_exp,Z_ket,rng)
            do i=1, sim%N_samples

                call rng%set_seed(seeds(i,:))
                call sample_SUZ(rng, Z_ket)
                psi = kron_vector(theta, Z_ket)
                call sample_obs(i)%malloc(N_steps+1)
                call sample_obs(i)%set(0.0_dp)
                sample_obs(i)%P_S(1) = init_pop(1)
                sample_obs(i)%P_T0(1) = init_pop(2)
                sample_obs(i)%P_Tp(1) = init_pop(3)
                sample_obs(i)%P_Tm(1) = init_pop(4)

                t = 0.0_dp
                k = 2

                ! print*, 'starting sample ', i
                do while (t < sim%t_end)

                    ! print*, 'Generating Krylov subspace at time t =  ', t

                    call get_krylov(N, psi, H, sim%N_krylov, Q, c, A)
                    call expm(((0.0_dp, -1.0_dp)*sim%dt*A), A_exp)
                    call convert_to_krylov(sim%N_krylov, Q, Z, P_S, P_T0, P_Tp, P_Tm)

                    c2 = 1.0_dp
                    c_m2 = 0.0_dp

                    do while (c_m2 <= sim%tol*c2 .and. t < sim%t_end)
                        c = matmul(A_exp, c)

                        ! Collect observables
                        sample_obs(i)%P_S(k) = real(dot_product(c, matmul(P_S, c)), kind=dp)
                        sample_obs(i)%P_Tp(k) = real(dot_product(c, matmul(P_Tp, c)), kind=dp)
                        sample_obs(i)%P_T0(k) = real(dot_product(c, matmul(P_T0, c)), kind=dp)
                        sample_obs(i)%P_Tm(k) = real(dot_product(c, matmul(P_Tm, c)), kind=dp)

                        c2 = real(dot_product(c,c), kind=dp)
                        c_m2 = real(c(sim%N_krylov)*conjg(c(sim%N_krylov)), kind=dp)
                        t = t+sim%dt
                        k = k+1
                    end do
                    psi = matmul(Q, c)
                end do
                if (mod(i, 100) == 0) then
                    print*, 'Sample ', i, '/', sim%N_samples, ' completed.'
                end if
            end do
            !$OMP END PARALLEL DO

            do i=1,sim%N_samples
                call res%update(sample_obs(i))
            end do

            ! Average observables
            call res%scale(1.0_dp/real(sim%N_samples, kind=dp))
            if (present(out)) then
                ! call res%get_kinetics(sim%dt, sys%kS, sys%kT)
                call res%output(out)
            end if

        end subroutine trace_sampling_para

        subroutine init_el_state(sim, theta, init_pop)
            type(sim_param), intent(in) :: sim
            complex(dp), intent(out)    :: theta(4) ! init electronic wavefunction
            real(dp), intent(out)       :: init_pop(4)

            init_pop=0.0_dp
            select case (sim%init_state)
                case('singlet')
                    theta%re = [0.0_dp, 1.0_dp/sqrt(2.0_dp), -1.0_dp/sqrt(2.0_dp), 0.0_dp]
                    theta%im = 0.0_dp
                    init_pop(1) = 1.0_dp
                case('T0')
                    theta%re = [0.0_dp, 1.0_dp/sqrt(2.0_dp), 1.0_dp/sqrt(2.0_dp), 0.0_dp]
                    theta%im = 0.0_dp
                    init_pop(2) = 1.0_dp
                case('Tp')
                    theta%re = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
                    theta%im = 0.0_dp
                    init_pop(3) = 1.0_dp
                case('Tm')
                    theta%re = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
                    theta%im = 0.0_dp
                    init_pop(4) = 1.0_dp
                case default
                    print*, 'Error: Invalid initial radical pair state entered.'
                    stop
            end select

    end subroutine init_el_state

    subroutine sample_SUZ(rng, Z_ket)
    ! Sample SU(Z) coherent state of dimension Z
        complex(dp), intent(inout) :: Z_ket(:)
        type(RNG_t), intent(inout) :: rng

        integer :: i

        do i=1,size(Z_ket)
            Z_ket(i)%re = rng%normal()
            Z_ket(i)%im = rng%normal()
        end do

        Z_ket= Z_ket/sqrt(dot_product(Z_ket, Z_ket))

    end subroutine sample_SUZ

    subroutine get_krylov(N, psi, H, m, Q, c, A)
    ! Constructs orthonormal basis of Krylov subspace (algorith adapted from Lewis )
        integer                        :: N
        complex(dp), intent(in)        :: psi(:)
        type(CSR_cdp_type), intent(in) :: H
        integer, intent(in)            :: m      ! Size of Krylov subspace
        complex(dp), intent(out)       :: Q(:,:) ! Krylov subspace basis vectors
        complex(dp), intent(out)       :: c(m)   ! Psi expressed in Krylov subspace basis
        complex(dp), intent(out)       :: A(m,m) ! Effective Hamiltonian in orthonormal Krylov subspace basis (upper Hessengerg form)

        complex(dp) :: p(N)
        integer     :: i, j

        c = (0.0, 0.0_dp)
        c(1) = sqrt(dot_product(psi, psi))
        A = (0.0, 0.0_dp)

        Q(:, 1) = psi/c(1)

        do i=1,m-1
            call spmv(H, Q(:, i), p)
            do j=1,i
                A(j, i) = dot_product(Q(:, j), p)
                p = p - A(j, i)*Q(:, j)
            end do
            A(i+1, i) = sqrt(dot_product(p, p))
            Q(:, i+1) = p/A(i+1, i)
        end do

        call spmv(H, Q(:, m), p)
        do j=1,m
            A(j, m) = dot_product(Q(:, j), p)
            p = p - A(j, m)*Q(:, j)
        end do

    end subroutine get_krylov

    subroutine convert_to_krylov(m, Q, Z, P_S, P_T0, P_Tp, p_Tm)
    ! Converts relevant projection operators to Krylov subspace
        integer, intent(in)     :: m
        complex(dp), intent(in) :: Q(:,:)
        integer, intent(in)     :: Z
        complex(dp), intent(out) :: P_S(m,m)
        complex(dp), intent(out) :: P_T0(m,m)
        complex(dp), intent(out) :: P_Tp(m,m)
        complex(dp), intent(out) :: P_Tm(m,m)

        complex(dp) :: Pab_ij, Pba_ij, ZQC_ij

        integer :: i, j

        do i=1,m
            do j=1,m
                P_Tp(i, j) = dot_product(Q(1:Z, i), Q(1:Z, j))

                Pab_ij = dot_product(Q(Z+1:2*Z, i), Q(Z+1:2*Z, j))
                Pba_ij = dot_product(Q(2*Z+1:3*Z, i), Q(2*Z+1:3*Z, j))
                ZQC_ij = dot_product(Q(2*Z+1:3*Z, i), Q(Z+1:2*Z, j))

                P_T0(i, j) = 0.5_dp*(Pab_ij + Pba_ij) + ZQC_ij
                P_S(i, j) = 0.5_dp*(Pab_ij + Pba_ij) - ZQC_ij

                p_Tm(i, J) = dot_product(Q(3*Z+1:4*Z, i), Q(3*Z+1:4*Z, j))
            end do
        end do

    end subroutine convert_to_krylov

    subroutine exact_dynamics_dense(sys, sim, res, out)
    ! Run quantum mechanical dynamics with dense matrices (only present for testing purposes)
        type(sys_param), intent(in)              :: sys
        type(sim_param), intent(in)              :: sim
        type(observables), intent(out)           :: res
        ! character(:), allocatable, intent(inout) :: output_folder ! Folder where all experimental data will be saved
        character(len=600), optional, intent(inout) :: out

        type(observables)          :: sample_obs
        integer                    :: N_steps
        complex(dp)                :: theta(4) ! init electronic wavefunction
        real(dp)                   :: init_pop(4)
        complex(dp), allocatable   :: Z_ket(:) ! init nuclear wavefunction
        complex(dp), allocatable   :: psi(:)   ! theta X Z_ket
        integer                    :: Z        ! size of nuclear Hilbert space
        integer                    :: N        ! size of combined Hilbert space
        complex(dp), allocatable   :: H(:,:)
        complex(dp), allocatable   :: propagator(:,:)
        complex(dp)                :: Pab, Pba, ZQC

        integer  :: i, j, k

        ! Singlet-born radical pair for now

        Z = sys%Z1*sys%Z2
        N = 4*Z

        call getH_dense(sys, H)
        call expm((0.0_dp, -1.0_dp)*sim%dt*H, propagator)

        call init_el_state(sim, theta, init_pop)

        allocate(Z_ket(Z))
        allocate(psi(4*Z))

        N_steps = floor(sim%t_end/sim%dt) + 1
        call sample_obs%malloc(N_steps+1)
        call res%malloc(N_steps+1)
        call res%set(0.0_dp)


        do i=1, Z
            Z_ket = (0.0_dp, 0.0_dp)
            Z_ket(i) = 1.0_dp
            psi = kron_vector(theta, Z_ket)
            call sample_obs%set(0.0_dp)
            sample_obs%P_S(1) = 1.0_dp

            do j=1,N_steps
                psi=matmul(propagator, psi)

                sample_obs%P_Tp(j+1) = real(dot_product(psi(1:Z), psi(1:Z)), kind=dp)

                Pab = dot_product(psi(Z+1:2*Z), psi(Z+1:2*Z))
                Pba = dot_product(psi(2*Z+1:3*Z), psi(2*Z+1:3*Z))
                ZQC = dot_product(psi(2*Z+1:3*Z), psi(Z+1:2*Z))

                sample_obs%P_S(j+1) = real(0.5_dp*(Pab + Pba) - ZQC, kind=dp)
                sample_obs%P_T0(j+1) = real(0.5_dp*(Pab + Pba) + ZQC, kind=dp)

                sample_obs%P_Tm(j+1) = real(dot_product(psi(3*Z+1:4*Z), psi(3*Z+1:4*Z)), kind=dp)

            end do

            call res%update(sample_obs)

        end do

        call res%scale(1.0_dp/real(Z, kind=dp))
        call res%get_kinetics(sim%dt, sys%kS, sys%kT)
        call res%output(out)

    end subroutine exact_dynamics_dense

end module dynamics
