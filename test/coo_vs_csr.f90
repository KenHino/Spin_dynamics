program coo_vs_csr      
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use hamiltonian
    use variables
    use utils
    use stdlib_sparse
    use stdlib_linalg_blas, only : dotc
    implicit none

    character(:), allocatable :: input_file
    type(sys_param)           :: sys
    type(sim_param)           :: sim
    type(COO_cdp_type)        :: H_coo
    type(CSR_cdp_type)        :: H_csr
    complex(dp), allocatable  :: psi(:), psi0(:)
    complex(dp), allocatable  :: H(:, :)
    real                      :: start, finish
    integer                   :: i
    integer :: N
    complex(dp) :: prods

    input_file = '/home/damianko/fpm/spinchem/test/coo_vs_csr.ini'
    call read_inp(input_file, sys, sim)

    call getHiso_sp(sys, H_coo)
    call getH_dense(sys, H)

    ! allocate(psi0(H_coo%ncols))
    ! allocate(psi(H_coo%ncols))

   allocate(psi0(H_coo%ncols), source=(1.0_dp, 0.0_dp))
   allocate(psi(H_coo%ncols), source=(1.0_dp, 0.0_dp))


    call cpu_time(start)
    do i=1,100000
        call spmv(H_coo, psi0, psi)
    end do
    call cpu_time(finish)
    print*, "COO Time = ", finish-start, ' s'


    ! call coo_to_csr(H_coo, H_csr)

    print*, 'Size of Hilbert space', H_csr%ncols
    
    call cpu_time(start)
        do i=1,100000
            psi = matmul(H, psi0)
            ! call spmv(H_csr, psi0, psi)
            ! call spmv(H_csr, psi0, psi)
        end do
    call cpu_time(finish)
    ! print*, "CSR Time = ", finish-start, ' s'
    print*, "Dense Time = ", finish-start, ' s'


end program  coo_vs_csr

