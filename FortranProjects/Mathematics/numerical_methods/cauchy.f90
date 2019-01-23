module cauchy
  use ode_interfaces, only: ode_function, propagator => temporal_scheme
  implicit none
    
    contains
    
    !-----------------------------------------------------------------------!
    !   ( SUBROUTINE ) cauchy_problem                                       !
    !-----------------------------------------------------------------------!
    !   Solves the Cauchy Problem defined by a given time domain,           !
    !   differential operator, temporal scheme and initial conditions.      !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(0:M)) time_domain                               !
    !               ! Ordered vector containing the times at which the      !
    !               ! solution of the cauchy problem is recorded.           !
    !---------------!-------------------------------------------------------!
    !   IN          ! (ode_function) differential_operator ( u, t )         !
    !               ! Derivative function of the Cauchy Problem where       !
    !               ! u is the state vector ant t the evaluation time.      !
    !---------------!-------------------------------------------------------!
    !   IN          ! (propagator) temporal scheme ( f, t1, t2, u1, u2 )    !
    !               ! Temporal scheme used for solving the Cauchy problem.  !
    !---------------!-------------------------------------------------------!
    !   INOUT       ! (real(0:M,N)) solution                                !
    !               ! Solution matrix of the Cauchy Problem. Rows contains  !
    !               ! the state vector at the time given by time_domain,    !
    !               ! while columns store the corresponding component of    !
    !               ! the state vector. On row 0, the initial condition of  !
    !               ! the Cauchy Problem must be provided.                  !
    !---------------!-------------------------------------------------------!
    subroutine cauchy_problem(     &
            time_domain          , &
            differential_operator, &
            temporal_scheme      , &
            solution               &
            )
        real, intent(in) :: time_domain(0:)  !0:M
        procedure(ode_function) :: differential_operator
        procedure(propagator) :: temporal_scheme
        real, intent(inout) :: solution(0:,:) !0:M (time steps), N(variables)
        integer :: m, i

        m = size(time_domain) - 1
        do i = 0, m - 1
            call temporal_scheme(          &
                f = differential_operator, &
                t1 = time_domain(i)      , &
                t2 = time_domain(i+1)    , &
                u1 = solution(i,:)       , &
                u2 = solution(i+1,:)       &
            )
        end do
      end subroutine
end module cauchy
