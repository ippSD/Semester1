module cauchy
  use interfaces
  implicit none

contains
  subroutine cauchy_problem(time_domain, differential_operator, temporal_scheme, solution)
    real, intent(in) :: time_domain(0:)  !0:M
    procedure(odes) :: differential_operator
    procedure(scheme) :: temporal_scheme
    real, intent(inout) :: solution(0:,:) !0:M (time steps), N(variables)
    integer :: m, i

    m = size(time_domain) - 1
    do i = 0, m - 1
        call temporal_scheme(differential_operator, time_domain(i), time_domain(i+1), solution(i,:), solution(i+1,:))
    end do
  end subroutine
end module cauchy
