!---------------------------------------------------------------------------!
!   richardson                                                              !
!---------------------------------------------------------------------------!
!   Use Richardson extrapolation on a sample function.                      !
!---------------------------------------------------------------------------!
    
program richardson
    use richardson_extrapolation
    use cauchy_problem_solver
    use temporal_schemes
    use dislin_mod
    implicit none
    
    integer, parameter :: M = 1000, M0 = 1
    
    integer :: i, l
    real, parameter :: TF = 1d2
    real :: time(0:M)
    real :: u_cauchy(0:M, 1), u_richard(0:M, 1)
    real :: richard_error(M), real_error_cauchy(M), real_error_richard(M)
    real :: richard_lyapunov(M), real_lyapunov(M)
    
    time = [(i*TF/M, i = 0, M)]
    u_cauchy(0,1) = 1d0
    u_richard(0,1) = 1d0
    
    call cauchy_problem(        &
        time                  , &
        cauchy_neg_exponential, &
        runge_kutta_4         , &
        u_cauchy                &
    )
    call richardson_extrapolator (                      &
        time_domain           = time                  , &
        differential_operator = cauchy_neg_exponential, &
        temporal_scheme       = runge_kutta_4         , &
        order                 = 4                     , &
        solution              = u_richard             , &
        error                 = richard_error           &
    )
    
    real_error_cauchy  = abs( exp(-time(1:M)) -  u_cauchy(1:M,1) )
    real_error_richard = abs( exp(-time(1:M)) - u_richard(1:M,1) )
    
    richard_lyapunov = log(richard_error/richard_error(1))
    real_lyapunov = log(real_error_richard/real_error_richard(1))
    
    !Plot Richardson's extrapolation error.
    call plot(                          &
        x = time(M0:M)                , &
        y = log10(richard_error(M0:M)), &
        color = "RED"                 , &
        ylabel = "$Log10(ERR)$"       , &
        hold_on = .true.                &
    )
    call plot_title("Richardson Error.")
    call plot_end()
    
    !Plot Richardson VS Cauchy analytic error.
    call plot(                               &
        x = time(M0:M)                     , &
        y = log10(real_error_cauchy(M0:M)) , &
        color = "RED"                      , &
        ylabel = "$Log10(ERR)$"            , &
        hold_on = .true.                     &
    )
    call plot(                               &
        x = time(M0:M)                     , &
        y = log10(real_error_richard(M0:M)), &
        color = "BLUE",                      &
        hold_on = .true.                     &
    )
    call plot_legend(["Cauchy", "Richardson"])
    call plot_title("Committed real error:", 1)
    call plot_subtitle("$Cauchy \quad VS \quad Richardson$")
    call plot_end()
    
    !Plot Richardson's expected error VS real error.
    call plot(                               &
        x = time(M0:M)                     , &
        y = log10(richard_error(M0:M))     , &
        color = "RED"                      , &
        ylabel = "$Log10(ERR)$"            , &
        hold_on = .true.                     &
    )
    call plot(                               &
        x = time(M0:M),                      &
        y = log10(real_error_richard(M0:M)), &
        color = "BLUE"                     , &
        hold_on = .true.                     &
    )
    call plot_legend(["Richardson", "Real"])
    call plot_title("Expected Error:", 1)
    call plot_subtitle("$Richardson \quad VS \quad Real$")
    call plot_end()
    
    ! Plot Lyapunov Exponents form real and Richardson's errors.
    call plot(                  &
        x = time(1:M)         , &
        y = real_lyapunov     , &
        color = "RED"         , &
        xlabel = "$t_N$"      , &
        ylabel = "$\lambda_N$", &
        hold_on = .true.        &
    )
    call plot(                  &
        x = time(1:M)         , &
        y = richard_lyapunov  , &
        color = "BLUE"        , &
        xlabel = "$t_N$"      , &
        ylabel = "$\lambda_N$", &
        hold_on = .true.        &
    )
    call plot_title("Evolution of Lyapunov Exponent")
    call plot_title("$\epsilon_N = \epsilon_0^{\lambda_N\cdot t_N}$", 3)
    call plot_legend(["Real Lyapunov Exp.", "Richardson's Lyapunov Exp."])
    call plot_end();
    
    ! Plot Deviation between Lyapunov Exp.s of real and Richardson errors.
    call plot(                                        &
        x = time(1:M)                               , &
        y = (richard_lyapunov - real_lyapunov) * 1d4, &
        xlabel = "$t_N$"                            , &
        ylabel = "$\lambda_N\cdot 10^{-4}$"         , &
        hold_on = .true.                              &
    )
    call plot_title("Lyapunov Exponent Deviation")
    call plot_end();
    
    
    contains
    
    !-----------------------------------------------------------------------!
    !   ( FUNCTION ) cauchy_neg_exponential                                 !
    !-----------------------------------------------------------------------!
    !   Derivative of the negative exponential for Cauchy Problem.          !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(1)) u                                           !
    !               ! State vector of the cauchy problem at time t.         !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real) t                                              !
    !               ! Time at which the state vector are given.             !
    !---------------!-------------------------------------------------------!
    !   OUT         ! (real(1)) f                                           !
    !               ! Derivative of the state vector at time t.             !
    !---------------!-------------------------------------------------------!
    function cauchy_neg_exponential(u,t) result(f)
        real, intent(in) :: u(:), t
        real :: f(size(u))
        
        f(1) = - u(1)
    end function
    
    !-----------------------------------------------------------------------!
    !   ( FUNCTION ) neg_exponential                                        !
    !-----------------------------------------------------------------------!
    !   Derivative of the negative exponential for Cauchy Problem.          !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real) t                                              !
    !               ! Time argument.                                        !
    !---------------!-------------------------------------------------------!
    !   OUT         ! (real(1)) f                                           !
    !               ! Negative exponential of time argument (exp(-t)).      !
    !---------------!-------------------------------------------------------!
    function neg_exponential(t) result(f)
        real, intent(in) :: t
        real :: f(1)
        
        f(1) = exp(-t)
    end function
    
end program richardson