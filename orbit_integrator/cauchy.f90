module cauchy
  implicit none
  interface
    function f_rn_rn(u, t) result(f)
      real (kind = 8), intent(in) :: u(0:),  t
      real (kind = 8) :: f(0:size(u))
    end function
  end interface
contains
  subroutine cauchy_problem(u, tf, f)
    real (kind = 8), intent(inout) :: u(0:,0:) !0:M (time steps), 0:N(variables)
    real (kind = 8), intent(in) :: tf !0:M
    procedure(f_rn_rn) :: f

    integer :: m
    real (kind = 8) :: delta_t
    integer :: n
    real (kind = 8) :: tn

    write(*,*) u(0,:)

    m = size(u(:,0))
    delta_t = tf / m

    do n = 0, m - 1
      tn = n * delta_t
      u(n + 1, :) = u(n, :) + delta_t * f(u(n, :), tn)
    end do
  end subroutine
end module cauchy
