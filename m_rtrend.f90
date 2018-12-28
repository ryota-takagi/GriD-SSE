module m_rtrend

  implicit none
  private

  public :: rtrend__rtrend
  public :: rtrend__rtrendp

contains

  !----------------------------------------------------------------------------------------------------------------
  subroutine rtrend__rtrend(nd,d,t)
  !----------------------------------------------------------------------------------------------------------------
    integer, intent(in)             :: nd
    real(8), intent(inout)          :: d(nd)
    real(8), intent(in),   optional :: t(nd)
    real(8) :: x(nd),tscale
    real(8) :: sum_x,sum_xx,sum_y,sum_xy,delta,a,b
    integer :: id
    if(nd<2) return
    if(present(t))then
      tscale = 1.0d0/(t(nd)-t(1))
      do id=1,nd
        x(id) = (t(id)-t(1))*tscale
      enddo
    else
      do id=1,nd
        x(id) = dble(id)/dble(nd)
      enddo
    endif
    sum_x  = sum(x(1:nd))
    sum_xx = sum(x(1:nd)*x(1:nd))
    sum_y  = sum(d(1:nd))
    sum_xy = sum(x(1:nd)*d(1:nd))
    delta = nd * sum_xx - sum_x * sum_x
    a = ( nd * sum_xy - sum_y * sum_x ) / delta
    b = ( sum_y * sum_xx - sum_x * sum_xy ) / delta
    do id=1,nd
      d(id) = d(id) - (a*x(id)+b)
    enddo
  end subroutine rtrend__rtrend

  !----------------------------------------------------------------------------------------------------------------
  subroutine rtrend__rtrendp(nd,d,t,t1,t2)
  !----------------------------------------------------------------------------------------------------------------
    integer, intent(in)    :: nd
    real(8), intent(inout) :: d(nd)
    real(8), intent(in)    :: t(nd)
    real(8), intent(in)    :: t1,t2
    real(8) :: x(nd),tscale
    real(8) :: sum_x,sum_xx,sum_y,sum_xy,delta,a,b
    integer :: id,nsum
    nsum = 0
    sum_x  = 0.0d0
    sum_xx = 0.0d0
    sum_y  = 0.0d0
    sum_xy = 0.0d0
    do id=1,nd
      if(t(id)>=t1.and.t(id)<=t2)then
        sum_x  = sum_x  + t(id)
        sum_xx = sum_xx + t(id)*t(id)
        sum_y  = sum_y  + d(id)
        sum_xy = sum_xy + t(id)*d(id)
        nsum = nsum + 1
      endif
    enddo
    if(nsum==0)return
    delta = nsum * sum_xx - sum_x * sum_x
    a = ( nsum * sum_xy - sum_y * sum_x ) / delta
    b = ( sum_y * sum_xx - sum_x * sum_xy ) / delta
    do id=1,nd
      d(id) = d(id) - (a*t(id)+b)
    enddo
  end subroutine rtrend__rtrendp

end module m_rtrend
