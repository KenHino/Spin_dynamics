module hamiltonian
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use utils
    use spinop
    use system
    use kroenecker

    implicit none
    private 
    public :: getHi, getHfull

    contains

    function getHi(e) result(H)
    ! Create one electron Hamiltonian

        type(electron), intent(in) :: e      ! Electron input  
        complex(dp), allocatable   :: H(:,:) ! Output hamiltonian

        integer                     :: Z 
        complex(dp), allocatable    :: S_x(:,:), S_y(:,:), S_z(:,:)
        integer                     :: Z_left, Z_right 
        complex(dp), allocatable    :: I_x(:,:), I_y(:,:), I_z(:,:)

        integer :: i

        ! We will need the electon spin operators to construct the Hamiltonian
        call getSpinop(2, S_x, S_y, S_z)

        ! Calculate dimension of nuclear subspace
        Z = product(e%g_I) 
        ! Total Hamiltonian subspace will be 2Z because we include the electron spin as well 
        allocate(H(2*Z, 2*Z))
        H = cmplx(0.0_dp, 0.0_dp, kind=dp)

        ! We first include the Zeeman terms in the interaction (the field is always in the z direction)
        H = H + kron_iden(e%w*S_z, Z)

        ! Add hyperfine coupling terms S.A.I for all nuclei
        call getSpinop(e%g_I(1), I_x, I_y, I_z)
        
        if (e%isotropic) then
            do i=1,size(e%g_I)

                ! The g_I array is sorted so we need to recompute spin operators only if the gurrent g_I is different from the preivous one   
                if (i > 1) then 
                    if (e%g_I(i) /= e%g_I(i-1)) then
                        deallocate(I_x)
                        deallocate(I_y)
                        deallocate(I_z)
                        call getSpinop(e%g_I(i), I_x, I_y, I_z)
                    end if
                end if

                Z_left = product(e%g_I(:i-1)) 
                Z_right = product(e%g_I(i+1:))

                ! Evaluation of the Hamiltonian is simplified in the case where all elecron couplings are isotropic
                H = H + (kron(kron_iden(e%A(i,1,1)*S_x, Z_left), kron_iden(I_x, Z_right)) &
                      + kron(kron_iden(e%A(i,1,1)*S_y, Z_left), kron_iden(I_y, Z_right)) &
                      + kron(kron_iden(e%A(i,1,1)*S_z, Z_left), kron_iden(I_z, Z_right)))

            end do
        else 
            do i=1,size(e%g_I)

                ! The g_I array is sorted so we need to recompute spin operators only if the gurrent g_I is different from the preivous one   
                if (i > 1) then 
                    if (e%g_I(i) /= e%g_I(i-1)) then
                        deallocate(I_x)
                        deallocate(I_y)
                        deallocate(I_z)
                        call getSpinop(e%g_I(i), I_x, I_y, I_z)
                    end if
                end if

                ! Calculate identity matrices needed to construct S.A.I operators
                Z_left = product(e%g_I(:i-1)) 
                Z_right = product(e%g_I(i+1:))

                ! This factorisation reduces the number of Kroenencker products that need to be evaluated
                H = H + kron(kron_iden(S_x, Z_left), kron_iden(e%A(i,1,1)*I_x + e%A(i,1,2)*I_y +e%A(i,1,3)*I_z, Z_right)) &
                    + kron(kron_iden(S_y, Z_left), kron_iden(e%A(i,2,1)*I_x + e%A(i,2,2)*I_y +e%A(i,2,3)*I_z, Z_right)) &
                    + kron(kron_iden(S_z, Z_left), kron_iden(e%A(i,3,1)*I_x + e%A(i,3,2)*I_y +e%A(i,3,3)*I_z, Z_right))

            end do
        end if

    end function getHi

    function getHfull(sys) result(H)
    ! Calculate fulll hamiltonian of radical pair 
        type(spinsys), intent(in) :: sys ! data strucuture containing all parameters needed to specify hamiltonian
        complex(dp), allocatable  :: H(:,:)

        integer                  :: Z1, Z2
        complex(dp), allocatable :: S_x(:,:), S_y(:,:), S_z(:,:) 
        complex(dp)              :: D_full(3,3)         

        ! We will need the electon spin operators to construct the Hamiltonian
        call getSpinop(2, S_x, S_y, S_z)

        ! Calculate dimension of nuclear subspaces 1 and 2
        Z1 = product(sys%e1%g_I) 
        Z2 = product(sys%e2%g_I) 

        allocate(H(4*Z1*Z2, 4*Z1*Z2))
        H = cmplx(0.0_dp, 0.0_dp, kind=dp)

        ! Combine isotropic and anisotropic 
        D_full = sys%D + 2.0_dp*sys%J*eye(3)

        H = H + kron(kron_iden(S_x, Z1), kron_iden(D_full(1,1)*S_x + D_full(1,2)*S_y + D_full(1,3)*S_z, Z2)) &
              + kron(kron_iden(S_y, Z1), kron_iden(D_full(2,1)*S_x + D_full(2,2)*S_y + D_full(2,3)*S_z, Z2)) &
              + kron(kron_iden(S_z, Z1), kron_iden(D_full(3,1)*S_x + D_full(3,2)*S_y + D_full(3,3)*S_z, Z2))

        ! Finally we add the electron couplings
        
        if (size(sys%e1%g_I) > 0) then
            H = H + kron_iden(getHi(sys%e1), 2*Z2)
        end if

        if (size(sys%e2%g_I) > 0) then
            H = H + kron_iden(2*Z1, getHi(sys%e2))
        end if

    end function 

end module hamiltonian