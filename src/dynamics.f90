module dynamics
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use system
    use stdlib_sparse

    implicit none
    public  
    private 

    type sim_res
        complex(dp), allocatable :: P_S(:)  ! Time-resolved |S> RP state population
        complex(dp), allocatable :: P_T0(:) ! Time-resolved |T0> RP state population
        complex(dp), allocatable :: P_Tp(:) ! Time-resolved |T+> RP state population
        complex(dp), allocatable :: P_Tm(:) ! Time-resolved |T-> RP state population
    end type sim_res

    contains

    subroutine monte_carlo(sys)
        type(spinsys), intent(in) :: sys ! data strucuture containing all parameters needed to specify hamiltonian

        type(COO_cdp_type)  :: Hsp
    

    end subroutine monte_carlo

end module dynamics