module hamiltonian
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use utils
    use spinop
    use variables
    use kroenecker
    use stdlib_sparse
    use stdlib_linalg, only: eye


    implicit none
    private 
    public :: getHiso_sp, getH_dense

    contains

    subroutine getHiso_sp(sys, H_coo) 
        type(sys_param), intent(in)     :: sys ! spin system paraWmeters
        type(COO_cdp_type), intent(out) :: H_coo
        ! type(CSR_cdp_type), intent(out) :: H_csr
        
        integer            :: N      ! Size of Hilbert space
        integer            :: nnz    ! Number of non-zero elements in spin Hamiltonian
        type(COO_cdp_type) :: H_tmp
        real(dp)           :: k_bar, delta_k 

        ! complex(dp), allocatable :: H(:,:) ! keep for testing

        type(COO_cdp_type) :: Sz, Sp, Sm    
        type(COO_cdp_type) :: S1p, S1m, S1z 
        type(COO_cdp_type) :: S2p, S2m, S2z 
        type(COO_cdp_type) :: S1pS2m, S1mS2p, S1zS2z 
        type(COO_cdp_type) :: Hel_diag, Hel_off 
        integer            :: Z_left, Z_right 
        type(COO_cdp_type) :: Ix, Iy, Iz, Ip, Im

        integer :: i, j
        integer :: start, end

        call get_spinop(2, Sz, Sp, Sm)

        S1p = kron_iden(Sp, 2)
        S1m = kron_iden(Sm, 2)
        S1z = kron_iden(Sz, 2)
        S2p = kron_iden(2, Sp)
        S2m = kron_iden(2, Sm)
        S2z = kron_iden(2, Sz)

        S1pS2m = kron(Sp, Sm) 
        S1mS2p = kron(Sm, Sp) 
        S1zS2z = kron(Sz, Sz) 
    
        N = 4*sys%Z1*sys%Z2

        k_bar = (sys%kS + 3.0_dp*sys%kT)/4.0_dp
        delta_k = (sys%kS - sys%kT)/4.0_dp

        nnz = 0
        ! Diagonal elements
        nnz = nnz + N
        ! Off-diagonal terms from S1xS2x and S2xS2y from exchange and recombination
        nnz = nnz + N/2
        ! Off-diagonal terms from SxIx and SyIy operators
        do i=1,size(sys%e1%g_I)  
            nnz = nnz + (sys%e1%g_I(i)-1)*(N/sys%e1%g_I(i))            
        end do
        do i=1,size(sys%e2%g_I)  
            nnz = nnz + (sys%e2%g_I(i)-1)*(N/sys%e2%g_I(i))            
        end do

        call H_coo%malloc(num_rows=N, num_cols=N, nnz=nnz)

        ! Initialise indices which keep track of off-diagonal elements
        start = N + 1
        
        ! Add diagonal terms of electronic Hamiltonian (S1z, S2z, S1zS2z, 1)
        Hel_diag = sys%e1%w*S1z + sys%e2%w*S2z & ! Zeeman
                    + (2*sys%J + (0.0_dp, 2.0_dp)*delta_k)*S1zS2z ! exchange and recombination

        ! There is also an identity term arising from K
        Hel_diag%data = Hel_diag%data - (0.0_dp, 0.5_dp)*k_bar
        
        H_tmp = kron_iden(Hel_diag, sys%Z1*sys%Z2)
        H_coo%data(1:N) = H_tmp%data
        H_coo%index(:, 1:N) = H_tmp%index

        ! Add off-diagonal terms of electronic Hamiltonian (S1+S2-, S1-S2+)
        H_tmp = kron_iden((sys%J + (0.0_dp, 1.0_dp)*delta_k)*S1pS2m, sys%Z1*sys%Z2)
        end = start+H_tmp%nnz-1
        H_coo%data(start:end) = H_tmp%data
        H_coo%index(:, start:end) = H_tmp%index
        start = end + 1

        H_tmp = kron_iden((sys%J + (0.0_dp, 1.0_dp)*delta_k)*S1mS2p, sys%Z1*sys%Z2)
        end = start+H_tmp%nnz-1
        H_coo%data(start:end) = H_tmp%data
        H_coo%index(:, start:end) = H_tmp%index
        start = end + 1

        ! Add isotropic hyperfine couplings for electron 1
        if (size(sys%e1%g_I)  > 0) then
            call get_spinop(sys%e1%g_I(1), Iz, Ip, Im)
        end if

        do i=1,size(sys%e1%g_I)

            ! The g_I array should be sorted so we need to recompute spin operators only if the gurrent g_I is different from the preivous one   
            if (i > 1) then 
                if (sys%e1%g_I(i) /= sys%e1%g_I(i-1)) then
                    ! Maybe deallocate before getting spinoperators again  
                    call get_spinop(sys%e1%g_I(i), Iz, Ip, Im)
                end if
            end if

            ! Calculate identity matrices needed to construct S.A.I operators
            Z_left = product(sys%e1%g_I(:i-1)) 
            Z_right = product(sys%e1%g_I(i+1:))*sys%Z2

            ! Add SzIz (diagonal operator)
            H_tmp = kron(kron_iden(sys%e1%a_iso(i)*S1z, Z_left), kron_iden(Iz, Z_right)) 
            H_coo%data(1:N) = H_coo%data(1:N) + H_tmp%data 

            ! Add SpIm and SmIp (off-diagonal elements)  
            H_tmp = kron(kron_iden(0.5*sys%e1%a_iso(i)*S1p, Z_left), kron_iden(Im, Z_right)) 
            end = start+H_tmp%nnz-1
            H_coo%data(start:end) = H_tmp%data
            H_coo%index(:, start:end) = H_tmp%index
            start = end + 1

            H_tmp = kron(kron_iden(0.5*sys%e1%a_iso(i)*S1m, Z_left), kron_iden(Ip, Z_right))
            end = start+H_tmp%nnz-1
            H_coo%data(start:end) = H_tmp%data
            H_coo%index(:, start:end) = H_tmp%index
            start = end + 1

        end do

        ! Repeat same procedure for electron 2
        if (size(sys%e2%g_I) > 0) then
            call get_spinop(sys%e2%g_I(1), Iz, Ip, Im)
        end if

        do i=1,size(sys%e2%g_I)

            ! The g_I array should be sorted so we need to recompute spin operators only if the gurrent g_I is different from the preivous one   
            if (i > 1) then 
                if (sys%e2%g_I(i) /= sys%e2%g_I(i-1)) then
                    ! Maybe deallocate before getting spinoperators again  
                    call get_spinop(sys%e2%g_I(i), Iz, Ip, Im)
                end if
            end if

            ! Calculate dimensions of identity matrices needed to construct S.A.I operators
            Z_left = sys%Z1*product(sys%e2%g_I(:i-1)) 
            Z_right = product(sys%e2%g_I(i+1:))

            ! Add SzIz (diagonal operator)
            H_tmp = kron(kron_iden(sys%e2%a_iso(i)*S1z, Z_left), kron_iden(Iz, Z_right)) 
            H_coo%data(1:N) = H_coo%data(1:N) + H_tmp%data 

            ! Remeber to change indixes back

            ! Add SyIy (off-diagonal, we need to add new set of indices) 
            H_tmp = kron(kron_iden(0.5*sys%e2%a_iso(i)*S1p, Z_left), kron_iden(Im, Z_right)) 
            end = start+H_tmp%nnz-1
            H_coo%data(start:end) = H_tmp%data
            H_coo%index(:, start:end) = H_tmp%index
            start = end + 1

            ! Add SxIx 
            H_tmp = kron(kron_iden(0.5*sys%e2%a_iso(i)*S1m, Z_left), kron_iden(Ip, Z_right)) 
            end = start+H_tmp%nnz-1
            H_coo%data(start:end) = H_tmp%data
            H_coo%index(:, start:end) = H_tmp%index
            start = end + 1

        end do

    end subroutine getHiso_sp

    subroutine getH_dense(sys, H) 
        type(sys_param), intent(in) :: sys ! spin system paraWmeters
        complex(dp), allocatable, intent(out) :: H(:,:)
        
        integer            :: N      ! Size of Hilbert space
        real(dp)           :: k_bar, delta_k 
        complex(dp)        :: H_el(4,4)

        complex(dp), allocatable :: Sx(:,:), Sy(:,:), Sz(:,:)    
        complex(dp)              :: S1x(4,4), S1y(4,4), S1z(4,4) 
        complex(dp)              :: S2x(4,4), S2y(4,4), S2z(4,4) 
        complex(dp)              :: D_full(3,3)         
        integer                  :: Z_left, Z_right 
        complex(dp), allocatable :: Ix(:,:), Iy(:,:), Iz(:,:), Ip(:,:), Im(:,:)    

        integer :: i
         
        N = 4*sys%Z1*sys%Z2

        k_bar = (sys%kS + 3.0_dp*sys%kT)/4.0_dp
        delta_k = (sys%kS - sys%kT)/4.0_dp

        allocate(H(N,N), source=(0.0_dp, 0.0_dp))

        call get_spinop(2, Sx, Sy, Sz)

        S1x = kron_iden(Sx, 2)
        S1y = kron_iden(Sy, 2)
        S1z = kron_iden(Sz, 2)
        S2x = kron_iden(2, Sx)
        S2y = kron_iden(2, Sy)
        S2z = kron_iden(2, Sz)

        ! Construct 4x4 electronic spin Hamiltonian 
        H_el = (0.0_dp, 0.0_dp)
        ! Add Zeeman terms
        H_el = H_el + sys%e1%w*S1z + sys%e2%w*S2z
        ! Combine isotropic and anisotropic electron coupling
        D_full = sys%D + 2.0_dp*sys%J*eye(3)
        ! Add coupling 
        H_el = H_el + kron(Sx, D_full(1,1)*Sx + D_full(1,2)*Sy + D_full(1,3)*Sz) &
                    + kron(Sy, D_full(2,1)*Sx + D_full(2,2)*Sy + D_full(2,3)*Sz) &
                    + kron(Sz, D_full(3,1)*Sx + D_full(3,2)*Sy + D_full(3,3)*Sz)
        ! Add Haberkorn recombination operator
        H_el = H_el - (0.0_dp, 1.0_dp)*k_bar*eye_cp(4)
        H_el = H_el + (0.0_dp, 1.0_dp)*delta_k*(kron(Sx, Sx) + kron(Sy, Sy) + kron(Sz, Sz))

        H = H + kron_iden(H_el, sys%Z1*sys%Z2)

        ! Add isotropic hyperfine couplings for electron 1
        if (size(sys%e1%g_I)  > 0) then
            call get_spinop(sys%e1%g_I(1), Ix, Iy, Iz)
        end if

        do i=1,size(sys%e1%g_I)

            ! The g_I array should be sorted so we need to recompute spin operators only if the gurrent g_I is different from the preivous one   
            if (i > 1) then 
                if (sys%e1%g_I(i) /= sys%e1%g_I(i-1)) then
                    call get_spinop(sys%e1%g_I(i), Ix, Iy, Iz)
                end if
            end if

            ! Calculate identity matrices needed to construct S.A.I operators
            Z_left = product(sys%e1%g_I(:i-1)) 
            Z_right = product(sys%e1%g_I(i+1:))*sys%Z2

            
            H = H + kron(sys%e1%A(i,1,1)*S1x + sys%e1%A(i,2,1)*S1y + sys%e1%A(i,3,1)*S1z, &
                         kron_iden(Z_left, kron_iden(Ix, Z_right)))                       &
                  + kron(sys%e1%A(i,1,2)*S1x + sys%e1%A(i,2,2)*S1y + sys%e1%A(i,3,2)*S1z, &
                         kron_iden(Z_left, kron_iden(Iy, Z_right)))                       &  
                  + kron(sys%e1%A(i,1,3)*S1x + sys%e1%A(i,2,3)*S1y + sys%e1%A(i,3,3)*S1z, &
                         kron_iden(Z_left, kron_iden(Iz, Z_right)))

        end do

        ! Repeat same procedure for electron 2
        if (size(sys%e2%g_I) > 0) then
            call get_spinop(sys%e2%g_I(1), Ix, Iy, Iz)
        end if

        do i=1,size(sys%e2%g_I)

            ! The g_I array should be sorted so we need to recompute spin operators only if the gurrent g_I is different from the preivous one   
            if (i > 1) then 
                if (sys%e2%g_I(i) /= sys%e2%g_I(i-1)) then
                    call get_spinop(sys%e2%g_I(i), Ix, Iy, Iz)
                end if
            end if

            ! Calculate dimensions of identity matrices needed to construct S.A.I operators
            Z_left = sys%Z1*product(sys%e2%g_I(:i-1)) 
            Z_right = product(sys%e2%g_I(i+1:))

            H = H + kron(sys%e2%A(i,1,1)*S2x + sys%e2%A(i,2,1)*S2y + sys%e2%A(i,3,1)*S2z, &
                          kron_iden(Z_left, kron_iden(Ix, Z_right)))                      &
                  + kron(sys%e2%A(i,1,2)*S2x + sys%e2%A(i,2,2)*S2y + sys%e2%A(i,3,2)*S2z, &
                          kron_iden(Z_left, kron_iden(Iy, Z_right)))                      &  
                  + kron(sys%e2%A(i,1,3)*S2x + sys%e2%A(i,2,3)*S2y + sys%e2%A(i,3,3)*S2z, &
                         kron_iden(Z_left, kron_iden(Iz, Z_right)))
        end do

    end subroutine getH_dense


end module hamiltonian