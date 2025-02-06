module moments
    ! Updated code from ref. https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.120.220604
    use, intrinsic :: iso_fortran_env, only: dp => real64, i8 => int64
    ! use stdlib_linalg
    use stdlib_linalg_lapack, only: steqr
    use stdlib_linalg, only: solve_lu

    
    implicit none
    public :: shrink, quad_rule, count_set_bits
    private

    contains

    subroutine shrink (a_iso, M, N, a_bar, n_bar)
    !  ------------------------------------------------------------------
    !  Optimally approximates N inequivalent nuclei with hyperfine
    !  coupling constants a(j) by M < N sets of equivalent nuclei
    !  with hyperfine constants abar(j), with nbar(j) nuclei in set j.
    !  ------------------------------------------------------------------
    !  We use integer*8 so that we can deal with all M < 64:
 
        real(dp), intent(in)  :: a_iso(N)  ! Isotropic hyperfine couplings 
        integer, intent(in)   :: N         ! Number of nuclear spins
        integer, intent(in)   :: M         ! Number of symmetry blocks
        real(dp), intent(out) :: a_bar(M) ! Isotropic hyperfine couplings 
        integer, intent(out)  :: n_bar(M) ! Isotropic hyperfine couplings 
        
        real(dp) :: W(N)
        real(dp) :: w_bar(M)
        real(dp) :: aw_bar(M)
        real(dp) :: atemp(M)
        integer :: nftot
        integer :: ntemp(M)
        integer :: nfloor(M)
        integer :: ndiff
        real(dp) :: dmom
        real(dp) :: emom
        integer(i8) :: nset
        real(dp) :: exact
        real(dp) :: approx
        integer :: i, j, k
        integer(i8) :: ib
        integer :: icode
        integer(i8) :: nind

        ! It only makes any sense to call this subroutine with M < N:
        if (M > N) stop 'M > N in shrink?'
        if (M == N) stop 'M = N in shrink?'
    
        W = 1.0_dp

        ! Optimum solution with non-integer weights:
        call quad_rule (n,w,a_iso,m,w_bar,aw_bar)
       
        nfloor = int(w_bar)
        nftot = sum(nfloor)
        ndiff = n-nftot

        ! Optimum solution with integer weights:
        dmom = huge(0.0_dp)
        do ib=0,LSHIFT(1,m)-1
            call count_set_bits(ib,nset)
            if (nset == ndiff) then
                do i=1,M
                    nind = RSHIFT(IAND(ib,LSHIFT(1_i8,i-1_i8)),i-1_i8)
                    ntemp(i) = nfloor(i)+nind
                end do
                atemp = aw_bar
                call newton (aw_bar,w_bar,atemp,ntemp,m,icode)
                if (icode == 0) then
                    emom = 0.0_dp
                    do k = 1,m+1
                        exact = sum(a_iso**k)
                        approx = sum(atemp**k)
                        emom = emom+abs(approx-exact)
                    end do
                    if (emom < dmom) then
                        dmom = emom  
                        n_bar = ntemp
                        a_bar = atemp
                    endif
                endif
            endif
        end do
    end subroutine shrink

    subroutine quad_rule(N_init, W_init, X_init, N, W, X)
        ! -----------------------------------------------------------------
        ! Uses a discrete Stieltjes procedure to construct an n-point
        ! contracted quadrature rule from an np-point primitive quadrature
        ! rule with the same (non-negative) weight function.
        ! -----------------------------------------------------------------
        integer, intent(in)   :: N_init    ! # of points in initial quadrature
        real(dp), intent(in)  :: W_init(N_init) ! weights in initial quadrature
        real(dp), intent(in)  :: X_init(N_init) ! nodes in initial quadrature
        integer, intent(in)  :: N         ! # of points in converted quadrature
        real(dp), intent(out) :: W(N)      ! weights in converted quadrature
        real(dp), intent(out) :: X(N)      ! nodes in converted quadrature

        real(dp) :: P(N_init, N)
        real(dp) :: q(N_init)
        real(dp), allocatable :: Z(:,:)
        real(dp) :: q2
        real(dp) :: pq
        real(dp) :: weight
        integer :: ldz
        real(dp) :: work(2*n-2)
        real(dp) :: E(n)
        integer :: info

        integer :: j, k

        ldz = N
        allocate(Z(ldz, N))

        q = sqrt(W_init)
        Z = 0.0_dp

        do k=1,N
            q2 = norm2(q)
            if (q2 == 0.0_dp) stop 'quad_rule 1'
            p(:, k) = q/q2
            q = X_init*P(:, k)

            do j = k,1,-1
                pq = dot_product(p(:, j), q)
                if(j == k) then
                    w(k) = q2
                    x(k) = pq
                end if
                q = q - pq*p(:, j)
            end do
        end do

        weight = W(1)
        Z(1,1) = 1.0_dp
        E = cshift(W, 1)
        
        
        call steqr('I', N, X, E(1:size(E)-1), Z, ldz, work, info)

        if (info /= 0) stop 'Eigenvalue algorithm in Stieltjes procedure did not work properly.'

        do j = 1,N
            w(j) = (weight*Z(1,j))**2
        end do

    end subroutine quad_rule

    subroutine count_set_bits(i,nset)
    ! ------------------------------------------------------------------
    ! Determines the total number of non-zero bits in a 64 bit integer
    ! using a 64 bit implementation of the Hamming Weight algorithm.
    ! ------------------------------------------------------------------
        integer(i8), intent(in)  :: i
        integer(i8), intent(out) :: nset
        
        integer(i8) :: v
        integer(i8) :: m1,m2,m3,m4
        
        DATA m1 /Z'5555555555555555'/, m2/Z'3333333333333333'/
        DATA m3 /Z'f0f0f0f0f0f0f0f'/, m4/Z'101010101010101'/

        v = i
        v = v-IAND(RSHIFT(v,1),m1)
        v = IAND(V,m2)+IAND(RSHIFT(V,2),m2)
        v = IAND(v+RSHIFT(v,4),m3)
        nset = RSHIFT(v*m4,(SIZEOF(nset)-1)*8)
    
    end subroutine count_set_bits

    subroutine newton(aw_bar,w_bar,a_bar,n_bar,m,icode)
    ! ------------------------------------------------------------------
    ! Uses Newton’s method to solve M non-linear moment equations in M
    ! unknowns (the hyperfine constants of M sets of equivalent nuclei),
    ! starting from an appropriate initial guess.
    ! ------------------------------------------------------------------
        real(dp), intent(in) :: aw_bar(M)
        real(dp), intent(out) :: a_bar(M)
        real(dp), intent(in) :: w_bar(M)
        integer, intent(in) :: n_bar(M)
        integer, intent(in) :: M
        integer, intent(out) :: icode
        
        real(dp) :: f(M)
        real(dp) :: q(M)
        real(dp) :: g(M,M)
        integer :: indx(M)
        integer :: ierr

        real(dp) :: a_abs, d_abs
        real(dp) :: da, dap
        integer :: iter
        integer, parameter :: maxit = 30
        integer :: j, k

        do iter = 1,maxit
            do k = 1,m
                f(k) = dot_product(n_bar, a_bar**k) - dot_product(w_bar, aw_bar**k)
                g(k, :) = k*n_bar*a_bar**(k-1)
            end do

            call solve_lu(G, f, q)

            a_abs = sum(abs(a_bar))
            d_abs = sum(abs(q))
            g(:,1) = a_bar
            a_bar = a_bar-q

            da = d_abs/a_abs
            
            if (da > 1.0_dp) then
                icode = 1
                return
            endif
            if (1.0_dp+da == 1.0_dp) then
                icode = 0
                return
            endif
            if (iter > 1 .and. da >= dap) then
                a_bar = g(:,1)
                icode = 0
                return
            endif
            
            dap = da
        end do
        icode = 2

    end subroutine newton

end module moments