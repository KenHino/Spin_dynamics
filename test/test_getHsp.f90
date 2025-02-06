! program test_getHsp      
!     use, intrinsic :: iso_fortran_env, only: dp => real64
!     use hamiltonian
!     use variables
!     use utils
!     use stdlib_sparse
!     implicit none

!     character(:), allocatable :: input_file
!     type(sys_param)           :: sys
!     type(sim_param)           :: sim
!     type(COO_cdp_type)        :: H_coo
!     complex(dp), allocatable  :: H_dense(:, :)
!     complex(dp), allocatable  :: H_conv(:, :)
!     complex(dp), allocatable  :: psi_sp(:), psi0_sp(:)
!     complex(dp), allocatable  :: psi_dense(:), psi0_dense(:)
!     complex(dp), allocatable  :: delta_psi(:)
!     real                      :: start, finish
!     integer                   :: i


!     input_file = '/home/damianko/fpm/spinchem/test/test_getHsp.ini'
!     call read_inp(input_file, sys, sim)



!     ! call getHiso_sp(sys, H_coo)
        
!     ! call getH_dense(sys, H_dense)

!     ! allocate(psi0_sp(H_coo%ncols), source=(0.0_dp, 0.0_dp))
!     ! psi0_sp(1) =(1.0_dp, 0.0_dp)
!     ! allocate(psi_sp(H_coo%ncols))
!     ! allocate(psi0_dense(H_coo%ncols), source=(0.0_dp, 0.0_dp))
!     ! psi0_dense(1) =(1.0_dp, 0.0_dp)
!     ! allocate(psi_dense(H_coo%ncols))
!     ! allocate(delta_psi(H_coo%ncols))

!     ! psi_dense = matmul(H_dense, psi0_dense)
!     ! call spmv(H_coo, psi0_sp, psi_sp)

!     ! delta_psi = psi_dense - psi_sp
    
!     ! if (any(delta_psi /= (0.0_dp, 0.0_dp))) then
!     !     print*, 'Sparse Hamiltonian is different from dense one.'
!     ! else
!     !     print*, 'Sparse Hamiltonian is same as dense one.'
!     ! end if 

! end program  test_getHsp

