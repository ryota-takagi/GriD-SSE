module m_sort

  implicit none
  private
  public :: sort__ascend
  public :: sort__descend

contains

  !--------------------------------------------------------------------------------------
  subroutine sort__ascend(n,x,ix)
  !--------------------------------------------------------------------------------------
    integer,  intent(in)    :: n
    real(8),  intent(inout) :: x(:)
    integer,  intent(inout) :: ix(:)
    integer                 :: iminloc(1),iswap,i,ixswap
    real(8)                 :: xswap
    do i=1,n-1
      ! search minloc witnin x(i:n)
      iminloc = minloc(x(i:n))
      iswap   = iminloc(1) + i - 1
      ! swap
      xswap     = x(iswap)
      x(iswap)  = x(i)
      x(i)      = xswap
      ixswap    = ix(iswap)
      ix(iswap) = ix(i)
      ix(i)     = ixswap
    enddo
  end subroutine sort__ascend

  !--------------------------------------------------------------------------------------
  subroutine sort__descend(n,x,ix)
  !--------------------------------------------------------------------------------------
    integer,  intent(in)    :: n
    real(8),  intent(inout) :: x(:)
    integer,  intent(inout) :: ix(:)
    integer                 :: imaxloc(1),iswap,i,ixswap
    real(8)                 :: xswap
    do i=1,n-1
      ! search maxloc witnin x(i:n)
      imaxloc = maxloc(x(i:n))
      iswap   = imaxloc(1) + i - 1
      ! swap
      xswap     = x(iswap)
      x(iswap)  = x(i)
      x(i)      = xswap
      ixswap    = ix(iswap)
      ix(iswap) = ix(i)
      ix(i)     = ixswap
    enddo
  end subroutine sort__descend

end module m_sort
