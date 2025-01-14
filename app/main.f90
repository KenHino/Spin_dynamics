program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spinop
    use hamiltonian
    use utils
    use finer
    use system
    use stdlib_sparse
    use sph
    use kroenecker
    implicit none

    character(len=50)     :: input_file
    type(spinsys) :: sys
    type(spinsys) :: sim
    real(dp), allocatable :: a(:)
    complex(dp) :: j, H(8,8)
    integer :: i 
    complex(dp), dimension(:,:), allocatable :: s_x, s_y, s_z, id
    complex(dp), allocatable :: sprod(:,:)
    type(COO_cdp_type)  :: Hsp
    type(COO_cdp_type)  :: Sxsp, Szsp, prod

    ! j%re=1.0_dp
    ! print*, sys%B

    ! call getSpinop(2, S_x, S_y, S_z)

    ! call dense2coo(S_x, Sxsp)
    ! call dense2coo(S_z, Szsp)

    ! prod = kron_iden(Szsp, 2*)

    ! call coo2dense(prod, sprod)    

    ! call print_matrix(sprod)
    ! print('(F0.0,SP,F0.0,"i")'), 



    ! call print_matrix(kron(s_z, s_x+s_y+s_z))

    ! H = kron_iden(S_x, 2)
    ! H = kron_iden(S_x, 1)

    input_file = '/home/damianko/fpm/spinchem/input.ini'
    sys = spinsys_init(input_file)

    Hsp = getHsp(sys)

    ! print*, Hsp%data(1:8)
    ! print*, 
    ! call print_matrix(Hsp%index(:,1:8))


    ! Hsp = getHsp(sys)

    ! call print_matrix(sys%e1%A(1,:,:)) 

    ! print*, size(sys%e2%g_I)

    ! H = getHfull(sys)

    ! H = getHi(sys%e1)

    ! call print_matrix(H)
    ! call print_matrix(getHi(sys%e1))

end program main
