program richardson
    use cauchy
    use numerical_methods
    use dislin
    implicit none
    
    logical, parameter :: DO_PLOT = .false.
    integer, parameter :: M = 1000, M0 = 1
    
    integer :: i, l
    real, parameter :: TF = 1d2
    real :: time(0:M), u_c(0:M, 1), u_r(0:M, 1), error(M), real_error_c(M), real_error_r(M)
    
    time = [(i*TF/M, i = 0, M)]
    u_c(0,1) = 1d0
    u_r(0,1) = 1d0
    
    call cauchy_problem(time, cauchy_neg_exponential, runge_kutta, u_c)
    call richardson_extrapolation(time, cauchy_neg_exponential, runge_kutta, 4, u_r, error)
    
    real_error_c = [(abs(neg_exponential(time(i)) - u_c(i,1)), i = 1, M)]
    real_error_r = [(abs(neg_exponential(time(i)) - u_r(i,1)), i = 1, M)]
    
    call qplot(time(M0:M), log(error(M0:M)), M-M0+1)
    call qplot(time(M0:M), log(real_error_c(M0:M)), M-M0+1)
    call qplot(time(M0:M), log(real_error_r(M0:M)), M-M0+1)
    call qplot(time(M0:M), log(error(M0:M)/real_error_c(M0:M)), M-M0+1)
    call qplot(time(M0:M), log(error(M0:M))/time(M0:M), M-M0+1)
    call qplot(time(M0:M), log(real_error_c(M0:M))/time(M0:M), M-M0+1)
    call qplot(time(M0:M), log(real_error_r(M0:M))/time(M0:M), M-M0+1)
    contains
    
    function cauchy_neg_exponential(u,t) result(f)
        real, intent(in) :: u(:), t
        real :: f(size(u))
        
        f(1) = - u(1)
    end function
    
    function neg_exponential(t) result(f)
        real, intent(in) :: t
        real :: f(1)
        
        f(1) = exp(-t)
    end function
    
    subroutine richardson_extrapolation(time_domain, differential_operator, temporal_scheme, order, solution, error, component)
        real, intent(in) :: time_domain(0:)
        procedure(odes) :: differential_operator
        procedure(scheme) :: temporal_scheme
        integer, intent(in) :: order
        real, intent(inout) :: solution(0:, :)
        real, intent(out) :: error(:)
        integer, optional, intent(in) :: component
        
        real, allocatable :: time_domain_refined(:)
        
        integer :: m, n, i, c
        real :: t_i, t_half, t_ipp
        real, allocatable :: y_i(:), y_i_refined(:)
        
        c = 1
        m = size(time_domain) - 1
        n = size(solution(0,:))
        allocate(y_i(n), y_i_refined(n))
        
        if(present(component)) then
            if(abs(component) > 0 .and. abs(component) <= n) then
                c = component
            end if
        end if

        y_i = solution(0,:)
        y_i_refined = solution(0,:)
        do i = 0, m - 1
            t_i = time_domain(i)
            t_ipp = time_domain(i+1)
            t_half = (t_i + t_ipp) / 2d0
            call temporal_scheme(differential_operator, t_i   , t_ipp , y_i        , y_i)
            call temporal_scheme(differential_operator, t_i   , t_half, y_i_refined, y_i_refined)
            call temporal_scheme(differential_operator, t_half, t_ipp , y_i_refined, y_i_refined)
            
            solution(i+1,:) = (2d0**order * y_i_refined - y_i) / (2d0**order - 1d0)
            error(i+1) = norm2((y_i_refined - y_i) / (1d0 - 2d0 ** -order))
        end do
            
    end subroutine
    
end program