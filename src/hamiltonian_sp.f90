module sph
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use system
    use stdlib_sparse
    use spinop
    use kroenecker
    use utils

    implicit none
    private 
    public :: getHsp 

    contains

    function getHsp(sys) result(H)
        type(spinsys), intent(in) :: sys ! spin system parameters
        type(COO_cdp_type)        :: H

        integer            :: Z1, Z2 ! Size of individual nuclear Hilbert spaces
        integer            :: N      ! Size of Hilberrt space
        integer            :: nnz    ! Number of non-zero elements in spin Hamiltonian
        type(COO_cdp_type) :: H_tmp

        type(COO_cdp_type) :: Sx, Sy, Sz    
        type(COO_cdp_type) :: S1x, S1y, S1z 
        type(COO_cdp_type) :: S2x, S2y, S2z 
        integer            :: Z_left, Z_right 
        type(COO_cdp_type) :: Ix, Iy, Iz
    
        integer :: i
        
        call getSpinop(2, Sx, Sy, Sz)

        S1x = kron_iden(Sx, 2)
        S1y = kron_iden(Sy, 2)
        S1z = kron_iden(Sz, 2)
        S2x = kron_iden(2, Sx)
        S2y = kron_iden(2, Sy)
        S2z = kron_iden(2, Sz)

        Z1 = product(sys%e1%g_I) 
        Z2 = product(sys%e2%g_I) 
        N = 4*Z1*Z2

        nnz = 0

        if (sys%e1%isotropic .and. sys%e2%isotropic) then

            ! Diagonal elements (=  to total size of Hilbert space)
            nnz = nnz + N
            ! Off-diagonal terms from SxIx and SyIy operators  
            nnz = nnz + N*(size(sys%e1%g_I) + size(sys%e2%g_I))
            ! Off-diagonal terms from S1xS2x and S2xS2y (only present for non-zero J)
            if (sys%J%re /= 0.0_dp) then
                nnz = nnz + N
            end if

            call H%malloc(num_rows=N, num_cols=N, nnz=nnz)

            ! Add Zeeman interaction
            H_tmp = kron_iden(sys%e1%w*S1z + sys%e2%w*S2z, Z1*Z2)
            H%data(1:N) = H_tmp%data
            H%index(:, 1:N) = H_tmp%index

            ! Add isotropic hyperfine couplings 

            call getSpinop(sys%e1%g_I(1), Ix, Iy, Iz)

            do i=1,size(sys%e1%g_I)

                ! The g_I array should be sorted so we need to recompute spin operators only if the gurrent g_I is different from the preivous one   
                if (i > 1) then 
                    if (sys%e1%g_I(i) /= sys%e1%g_I(i-1)) then
                        call getSpinop(sys%e1%g_I(i), Ix, Iy, Iz)
                    end if
                end if

                ! Calculate identity matrices needed to construct S.A.I operators
                Z_left = product(sys%e1%g_I(:i-1)) 
                Z_right = product(sys%e1%g_I(i+1:))*Z2

                ! Add SzIz (diagonal operator)
                H_tmp = kron(kron_iden(sys%e1%A(i,1,1)*S1z, Z_left), kron_iden(Iz, Z_right)) 
                H%data(1:N) = H%data(1:N) + H_tmp%data 

                ! Add SyIy (off-diagonal, we need to add new set of indices) 
                H_tmp = kron(kron_iden(sys%e1%A(i,1,1)*S1y, Z_left), kron_iden(Iy, Z_right)) 
                H%data((i-1)*N+1:i*N) = H_tmp%data
                H%index(:, (i-1)*N+1:i*N) = H_tmp%index

                ! Add SxIx - remeber to try the faster algorithm 
                H_tmp = kron(kron_iden(sys%e1%A(i,1,1)*S1x, Z_left), kron_iden(Ix, Z_right)) 
                H%data((i-1)*N+1:i*N) = H%data((i-1)*N+1:i*N) + H_tmp%data

            end do

            call getSpinop(sys%e2%g_I(1), Ix, Iy, Iz)

            do i=1,size(sys%e2%g_I)

                ! The g_I array should be sorted so we need to recompute spin operators only if the gurrent g_I is different from the preivous one   
                if (i > 1) then 
                    if (sys%e2%g_I(i) /= sys%e2%g_I(i-1)) then
                        call getSpinop(sys%e2%g_I(i), Ix, Iy, Iz)
                    end if
                end if

                ! Calculate identity matrices needed to construct S.A.I operators
                Z_left = product(sys%e2%g_I(:i-1)) 
                Z_right = product(sys%e2%g_I(i+1:))*Z2

                ! Add SzIz (diagonal operator)
                H_tmp = kron(kron_iden(sys%e2%A(i,1,1)*S2z, Z_left), kron_iden(Iz, Z_right)) 
                H%data(1:N) = H%data(1:N) + H_tmp%data 

                ! Add SyIy (off-diagonal, we need to add new set of indices) 
                H_tmp = kron(kron_iden(sys%e2%A(i,1,1)*S2y, Z_left), kron_iden(Iy, Z_right)) 
                H%data((i-1)*N+1:i*N) = H_tmp%data
                H%index(:, (i-1)*N+1:i*N) = H_tmp%index

                ! Add SxIx - remeber to try the faster algorithm 
                H_tmp = kron(kron_iden(sys%e2%A(i,1,1)*S2x, Z_left), kron_iden(Ix, Z_right)) 
                H%data((i-1)*N+1:i*N) = H%data((i-1)*N+1:i*N) + H_tmp%data

            end do

            ! We may not need to recompute SxIx

            ! Sort H in the end

        end if

    end function getHsp

end module sph