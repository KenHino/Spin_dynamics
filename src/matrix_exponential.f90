module matrix_exponential
    ! Code adapted from Expokit software package - https://www.maths.uq.edu.au/expokit/
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use utils
    use stdlib_linalg, only: inv

    implicit none

    public :: expm
    private :: log_2, mat_norm_linf

    interface expm
      procedure expm_dp
      procedure expm_cdp
    end interface expm

    interface mat_norm_linf
      procedure mat_norm_linf_dp
      procedure mat_norm_linf_cdp
    end interface mat_norm_linf

    contains

    subroutine expm_dp(A, E)
    ! Computes exponential of matrix
      real(dp), intent(in)               :: A(:,:)
      real(dp), allocatable, intent(out) :: E(:,:)
      
      
      integer                  :: n
      complex(dp), allocatable :: A2(:,:)
      real(dp), allocatable    :: a_norm
      real(dp)                 :: c
      real(dp)                 :: c8mat_norm_li
      complex(dp), allocatable :: D(:,:)
      integer                  :: ee   
      integer                  :: k   
      logical                  :: p
      integer , parameter      :: q = 6
      real(dp)                 :: r8_log_2
      integer                  :: s
      complex(dp), allocatable :: X(:,:)
    
      n = size(A,dim=1)

      allocate(E(n,n))
      allocate(A2(n,n))
      allocate(D(n,n))
      allocate(X(n,n))

      ! Make a copy of the matrix.
      A2 = A

      ! Compute the L-infinity norm.
      a_norm = mat_norm_linf(A2)

      ! Determine a scaling factor for the matrix.
      ee = int(log_2(a_norm)) + 1

      s = max(0, ee + 1)
      A2 = A2/2.0_dp**s
      x = A2

      c = 0.5_dp
      E = eye_cp(n) 
      E = E + c*A2

      D = eye_cp(n)
      D = D - c*A2
      p = .true.

      do k = 2, q

        c = c * real(q-k+1, kind=dp) &
          / real(k*(2*q-k+1), kind=dp)
        
        X = matmul(A2, x)
        E = E + c*X
        
        if (p) then
          D = D + c*X
        else
          D = D - c*X
        end if
        p = .not. p

      end do

      ! E -> inverse(D) * E
      E = matmul(inv(D), E)

      ! E -> E^(2*S)
      do k = 1, s
        E = matmul(E,E)
      end do

    end subroutine expm_dp

    subroutine expm_cdp(A, E)
    ! Computes exponential of matrix
      complex(dp), intent(in)               :: A(:,:)
      complex(dp), allocatable, intent(out) :: E(:,:)
      
      
      integer                  :: n
      complex(dp), allocatable :: A2(:,:)
      real(dp), allocatable    :: a_norm
      real(dp)                 :: c
      real(dp)                 :: c8mat_norm_li
      complex(dp), allocatable :: D(:,:)
      integer                  :: ee   
      integer                  :: k   
      logical                  :: p
      integer , parameter      :: q = 6
      real(dp)                 :: r8_log_2
      integer                  :: s
      complex(dp), allocatable :: X(:,:)
    
      n = size(A,dim=1)

      allocate(E(n,n))
      allocate(A2(n,n))
      allocate(D(n,n))
      allocate(X(n,n))

      ! Make a copy of the matrix.
      A2 = A

      ! Compute the L-infinity norm.
      a_norm = mat_norm_linf(A2)

      ! Determine a scaling factor for the matrix.
      ee = int(log_2(a_norm)) + 1

      s = max(0, ee + 1)
      A2 = A2/2.0_dp**s
      x = A2

      c = 0.5_dp
      E = eye_cp(n) 
      E = E + c*A2

      D = eye_cp(n)
      D = D - c*A2
      p = .true.

      do k = 2, q

        c = c * real(q-k+1, kind=dp) &
          / real(k*(2*q-k+1), kind=dp)
        
        X = matmul(A2, x)
        E = E + c*X
        
        if (p) then
          D = D + c*X
        else
          D = D - c*X
        end if
        p = .not. p

      end do

      ! E -> inverse(D) * E
      E = matmul(inv(D), E)

      ! E -> E^(2*S)
      do k = 1, s
        E = matmul(E,E)
      end do

    end subroutine expm_cdp

    function mat_norm_linf_dp(a) result(norm_linf)
        ! returns the matrix L_inf norm of a matrix
          real(dp)    :: A(:,:) ! A(M,N), the matrix whose Loo norm is desired.
          real(dp)    :: norm_linf ! the Loo norm of A.
          
          real(dp) :: row_sum
          integer  :: i
        
          norm_linf = 0.0_dp
        
          do i = 1, size(A, dim=1)
            row_sum = sum(abs(A(i,:)))
            norm_linf = max(norm_linf, row_sum)
          end do
        
    end function mat_norm_linf_dp

    function mat_norm_linf_cdp(a) result(norm_linf)
        ! returns the matrix L_inf norm of a matrix
          complex(dp) :: A(:,:) ! A(M,N), the matrix whose Loo norm is desired.
          real(dp)    :: norm_linf ! the Loo norm of A.
          
          real(dp) :: row_sum
          integer  :: i
        
          norm_linf = 0.0_dp
        
          do i = 1, size(A, dim=1)
            row_sum = sum(abs(A(i,:)))
            norm_linf = max(norm_linf, row_sum)
          end do
        
    end function mat_norm_linf_cdp

    function log_2(x) result(log_2x)
  
          real(dp), intent(in) :: x
          real(dp)             :: log_2x
        
          if ( x == 0.0_dp ) then
            log_2x = -huge(x)
          else
            log_2x = log(abs(x))/log(2.0_dp)
          end if
        
    end function log_2


end module matrix_exponential