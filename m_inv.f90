module m_inv

  implicit none
  private
  public :: inv

contains

  subroutine inv(n,A)
    integer, intent(in)    :: n
    real(8), intent(inout) :: A(:,:)
    integer                :: IPIV(n),INFO,LWORK
    real(8), allocatable   :: WORK(:)
    real(8)                :: LWORK0
    call DGETRF(n,n,A,n,IPIV,INFO)
    call DGETRI(n,A,n,IPIV,LWORK0,-1,INFO)
    LWORK = int(LWORK0)
    allocate(WORK(LWORK))
    call DGETRI(n,A,n,IPIV,WORK,LWORK,INFO)
    deallocate(WORK)
    return
  end subroutine inv

end module m_inv
