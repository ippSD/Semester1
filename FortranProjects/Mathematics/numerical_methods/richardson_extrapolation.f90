!---------------------------------------------------------------------------!
!   ( MODULE ) richardson_extrapolation                                     !
!---------------------------------------------------------------------------!
!   Implementation of Richardson extrapolation.                             !
!---------------------------------------------------------------------------!
    
module richardson_extrapolation
    use ode_interfaces, propagator => temporal_scheme
    implicit none
    
    contains
    
    !-----------------------------------------------------------------------!
    !   ( SUBROUTINE ) richardson_extrapolator                              !
    !-----------------------------------------------------------------------!
    !   Solves the Cauchy Problem for a time domain, a refined time domain  !
    !   (DT_refined = DT / 2) and returns an extrapolation of the solution  !
    !   based on the included temporal scheme's order as well as the        !
    !   extrapolation of the commited error.                                !
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
    !   IN          ! (integer) order                                       !
    !               ! Order of the temporal scheme.                         !
    !---------------!-------------------------------------------------------!
    !   INOUT       ! (real(0:M,N)) solution                                !
    !               ! Solution matrix of the Cauchy Problem with            !
    !               ! Richardson extrapolation:                             !
    !               !                                                       !
    !               !            order                                      !
    !               !          2       * Refined_Solution  +  1 * Solution  !
    !               !    W = ---------------------------------------------- !
    !               !            order                                      !
    !               !          2                           +  1             !
    !---------------!-------------------------------------------------------!
    !   INOUT       ! (real(1:M)) error                                     !
    !               ! Error vector of the Cauchy Problem with               !
    !               ! Richardson extrapolation:                             !
    !               !                                                       !
    !               !      |    order                (2n              (n |  !
    !               !   (n |  2       * Ref._Solution  -  1 * Solution   |  !
    !               !  W = ------------------------------------------------ !
    !               !           order                                       !
    !               !         2                        -  1                 !
    !---------------!-------------------------------------------------------!
    subroutine richardson_extrapolator( &
        time_domain                   , &
        differential_operator         , &
        temporal_scheme               , &
        order                         , &
        solution                      , &
        error)
    
        real, intent(in) :: time_domain(0:)
        procedure(ode_function) :: differential_operator
        procedure(propagator) :: temporal_scheme
        integer, intent(in) :: order
        real, intent(inout) :: solution(0:,:)
        real, intent(out) :: error(:)
        
        real, allocatable :: time_domain_refined(:)
        procedure(ode_function), pointer :: f
        procedure(propagator), pointer :: integrate
        
        integer :: m, n, i
        real :: time_current, time_half, time_next
        real, allocatable :: y(:), y_refined(:)
        
        f   => differential_operator
        integrate => temporal_scheme
        m = size(time_domain) - 1
        n = size(solution(0,:))
        allocate(y(n), y_refined(n))
        y = solution(0,:)
        y_refined = solution(0,:)

        do i = 0, m - 1
            time_current = time_domain(i)
            time_next = time_domain(i+1)
            time_half = (time_current + time_next) / 2d0
            call integrate(f, time_current, time_next, y        , y)
            call integrate(f, time_current, time_half, y_refined, y_refined)
            call integrate(f, time_half   , time_next, y_refined, y_refined)
            
            solution(i+1,:) = (2d0**order*y_refined-y)/(2d0**order - 1d0)
            error(i+1) = norm2((y_refined - y) / (1d0 - 2d0 ** -order))
        end do
            
    end subroutine richardson_extrapolator
    
end module richardson_extrapolation