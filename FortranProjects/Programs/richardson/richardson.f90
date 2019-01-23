!---------------------------------------------------------------------------!
!   richardson                                                              !
!---------------------------------------------------------------------------!
!   Use Richardson extrapolation on a sample function.                      !
!---------------------------------------------------------------------------!
    
program richardson
    use richardson_extrapolation
    use cauchy
    use numerical_methods
    use dislin_mod
    implicit none
    
    logical, parameter :: DO_PLOT = .false.
    integer, parameter :: M = 1000, M0 = 1
    
    integer :: i, l
    real, parameter :: TF = 1d2
    real :: time(0:M)
    real :: u_cauchy(0:M, 1), u_richard(0:M, 1)
    real :: richard_error(M), real_error_cauchy(M), real_error_richard(M)
    
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
    
    call plot(time(M0:M), log10(richard_error(M0:M)), "RED", ylabel = "$Log10(ERR)$", hold_on = .true.)
    call plot_title("Richardson Error.", 1)
    call plot_end()
    
    call plot(time(M0:M), log10(real_error_cauchy(M0:M)), "RED", ylabel = "$Log10(ERR)$", hold_on = .true.)
    call plot(time(M0:M), log10(real_error_richard(M0:M)), "BLUE", hold_on = .true.)
    call plot_legend(["Cauchy", "Richardson"])
    call plot_title("Committed real error:", 1)
    call plot_subtitle("$Cauchy \quad VS \quad Richardson$")
    call plot_end()
    
    call plot(time(M0:M), log10(richard_error(M0:M)), "RED", ylabel = "$Log10(ERR)$", hold_on = .true.)
    call plot(time(M0:M), log10(real_error_richard(M0:M)), "BLUE", hold_on = .true.)
    call plot_legend(["Richardson", "Real"])
    call plot_title("Expected Error:", 1)
    call plot_subtitle("$Richardson \quad VS \quad Real$")
    call plot_end()
    
    
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