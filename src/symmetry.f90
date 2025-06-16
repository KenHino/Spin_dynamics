module symmetry
    use, intrinsic :: iso_fortran_env, only: dp => real64, i8 => int64
    USE omp_lib
    use variables
    use hamiltonian
    use class_observables
    use dynamics
    use m_random
    use moments
    implicit none
    private
    public :: symmetrised_dynamics

    contains
    
    subroutine symmetrised_dynamics(sys, sim, rng, out)
    ! Run quantum mechanical dynamics with trace sampling
        type(sys_param), intent(in)        :: sys 
        type(sim_param), intent(inout)     :: sim 
        type(RNG_t), intent(inout)         :: rng
        character(len=600),  intent(inout) :: out

        integer(i8), allocatable       :: seeds(:,:) 
        type(observables)              :: res 
        integer(i8)                    :: Z
        integer                        :: N_steps
        real(dp), allocatable          :: a1_bar(:), a2_bar(:)
        integer, allocatable           :: n1_bar(:), n2_bar(:)
        integer, allocatable           :: N_bar(:)
        integer, allocatable           :: K_init(:,:), K(:,:) 
        integer, allocatable           :: gI1(:), gI2(:) 
        integer(i8), allocatable       :: Z_current(:)
        integer(i8)                    :: Z1, Z2
        integer(i8)                    :: w_ij
        real(dp)                       :: w_k
        type(sys_param)                :: sys_new
        type(observables), allocatable :: res_current(:)
        character(len=600)             :: folder
        character(len=16)              :: tmp
        integer                        :: start
        integer, dimension(3)          :: begin, finish, time, total

        integer :: maxi

        integer :: i, j

        ! Need to use 64-bit integer to calculate Z for large number of spins
        Z = product(int(sys%e1%g_I,kind=i8))*product(int(sys%e2%g_I,kind=i8))

        N_steps = ceiling((sim%t_end+sim%dt)/sim%dt)
        call res%malloc(N_steps+1)
        call res%set(0.0_dp)

        if (any(sys%e1%g_I /= 2) .or. any(sys%e2%g_I /= 2)) stop 'Method has not been implemented for nuclei with I>1/2.'

        call shrink(sys%e1%a_iso, sim%M1, a1_bar, n1_bar)
        call shrink(sys%e2%a_iso, sim%M2, a2_bar, n2_bar)

        allocate(n_bar(sim%M1 + sim%M2))

        n_bar(1:sim%M1) = n1_bar 
        n_bar(sim%M1+1:1:sim%M1+sim%M2) = n2_bar 

        call cartesian_product(n_bar, K_init)
        call sort_col_prod(K_init, Z_current, K)
        
        maxi = 0
        total = 0

        allocate(res_current(size(K, dim=2)))
        allocate(seeds(size(K, dim=2), 2))

        allocate(gI1(sim%M1))
        allocate(gI2(sim%M2))

        do i=1,size(seeds, dim=1)
            seeds(i,1) = rng%next()
            seeds(i,2) = rng%next()
        end do
        
        if (sim%isRestart) then
            start = sim%restart
        else
            start = 1
        end if    

        do i=start,size(K, dim=2)

            call itime(begin)
            call rng%set_seed(seeds(i,:))

            ! Make these two into pointers
            gI1 = K(1:sim%M1,i)
            gI2 = K(sim%M1+1:sim%M1+sim%M2,i)
            w_ij = weight(n1_bar, gI1) * weight(n2_bar, gI2)
            w_k = real(w_ij*Z_current(i), kind=dp)/real(Z, kind=dp) 

            if (w_k >= sim%block_tol) then
                write(tmp,'(I0)') i 
                folder = sim%output_folder // '/hamiltonian_' // trim(tmp)
                call system(' mkdir ' // folder // ' > /dev/null 2>&1')
                call reduce_system(sys, gI1, a1_bar, gI2, a2_bar, sys_new)

                if (Z_current(i) > maxi) then 
                    maxi = Z_current(i)
                end if
                print'(A, I0, A, I0)', 'starting trace sampling of hamiltonian ' , i,  '/' , size(K,dim=2)
                print'(A, I0)', 'Z = ', Z_current(i)   
                print'(A, I0, A, F0.16)', 'weight_', i ,' = ', w_k   
                call trace_sampling_para(sys_new, sim, rng, res_current(i), folder)

                ! Weight result of simulation 
                call res_current(i)%scale(real(w_k, kind=dp))
            else
                print'(A, I0, A, I0)', 'skipped hamiltonian ' , i
                print'(A, I0)', 'Z = ', Z_current(i)   
                print'(A, I0, A, F0.16)', 'weight_', i ,' = ', w_k   
                call res_current(i)%malloc(N_steps+1)
                call res_current(i)%set(0.0_dp)
            end if
            call itime(finish)

            call get_time(begin, finish, time)
            print'(A29,I0,A2,I2,A1,I2,A1,I2)', 'Time elapsed for hamiltonian ', i ,': ', time(1), ':',time(2), ':', time(3) 
            call add_time(total, time)

        end do

        do i=1,size(K, dim=2)
            call res%update(res_current(i))        
        end do

        print'(A20,I2,A1,I2,A1,I2)', 'Total time elapsed: ', total(1), ':',total(2), ':', total(3) 
        print'(A, I0)', 'Nuclear subspace of largest Hamiltonian: ' ,maxi
        
        call res%get_kinetics(sim%dt, sys%kS, sys%kT)
        call res%output(out)

    end subroutine symmetrised_dynamics

    subroutine reduce_system(sys, g_I1, a1_iso, g_I2, a2_iso, sys_new)
        type(sys_param), intent(in)  :: sys 
        integer, intent(in)          :: g_I1(:), g_I2(:)
        real(dp), intent(in)         :: a1_iso(:), a2_iso(:)
        type(sys_param), intent(out) :: sys_new 

        sys_new%J = sys%J
        sys_new%D = sys%D
        sys_new%kS = sys%kS
        sys_new%kT = sys%kT
        

        sys_new%e1%w = sys%e1%w 
        allocate(sys_new%e1%g_I(size(g_I1)))
        allocate(sys_new%e1%a_iso(size(g_I1)))
        sys_new%e1%g_I = g_I1
        sys_new%e1%a_iso = a1_iso
        sys_new%e1%isotropic = .true.
        sys_new%Z1 = product(sys_new%e1%g_I)
        
        sys_new%e2%w = sys%e2%w 
        allocate(sys_new%e2%g_I(size(g_I2)))
        allocate(sys_new%e2%a_iso(size(g_I2)))
        sys_new%e2%g_I = g_I2
        sys_new%e2%a_iso = a2_iso
        sys_new%e2%isotropic = .true.
        sys_new%Z2 = product(sys_new%e2%g_I)

    end subroutine reduce_system

    subroutine cartesian_product(n, k)
    ! Generate cartesian product of size(n) sets, where the i-th set has n(i) elements (not exactly because elements of na )
        integer, intent(in) :: n(:)
        integer, allocatable, intent(out) :: k(:,:)

        integer :: combinations 
        integer :: repeat
        integer :: length
        integer :: start, finish
        integer :: i, j

        combinations = product(n/2+1)
        allocate(k(size(n), combinations))

        repeat = 1
        do i=1,size(n)
            start = 1
            finish = repeat
            do j=1,n(i)/2+1
                k(i, start:finish) = 2*(j-1)+1+modulo(n(i), 2) 
                start = start + repeat
                finish = finish + repeat
            end do

            length = finish - repeat

            start = length+1
            finish = 2*length
            do j=2,combinations/length
                k(i, start:finish) = k(i, 1:length) 
                start = start + length
                finish = finish + length
            end do

            repeat = repeat*(n(i)/2+1)
        end do

    end subroutine cartesian_product

    function weight(n, k) result(w)
        integer, intent(in) :: n(:)
        integer, intent(in) :: k(:)
        integer(i8) :: w

        integer(i8) :: w_i
        integer  :: i

        w = 1
        do i=1,size(n)
            w_i = nCr(n(i), (n(i) + k(i) - 1)/2) * int(k(i)*2, kind=i8)
            w_i = w_i/int(n(i)+k(i)+1, kind=i8)                        
            w = w*w_i
        end do

    end function weight

    function nCr(n, r) result(C)
        integer, intent(in) :: n
        integer, intent(in) :: r
        integer(i8)         :: C

        integer(i8) :: numer
        integer(i8) :: denom
        integer(i8) :: i

        numer = 1_i8
        denom = 1_i8
        
        do i=1,int(min(n-r,r), kind=i8)
            numer = numer*(int(n, kind=i8)-i+1_i8)
            denom = denom*i
        end do

        C = numer/denom

    end function nCr

end module symmetry