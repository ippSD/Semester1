module temporal_schemes
    use rkf45
    use Linear_Systems
    use Non_Linear_Systems
    use ode_interfaces, only: ode_function
    use dop853_module
    use module_ode, only: ode => ode_int
    use module_odex, only: odex => odex_int
    implicit none
    
    private
    public :: euler_explicit, euler_implicit, runge_kutta_4, dopri853, ode, odex

    contains
    
    !-----------------------------------------------------------------------!
    !   ( SUBROUTINE ) euler_explicit                                       !
    !-----------------------------------------------------------------------!
    !   Extrapolates the solution of the Cauchy Problem at a given time     !
    !   and initial conditions by means of the Explicit Euler               !
    !   numerical scheme.                                                   !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (ode_function) f(u,t)                                 !
    !               ! f: RN x R => RN                                       !
    !               !     u , t => du_dt                                    !
    !               ! Derivative function of the Cauchy Problem where       !
    !               ! u is the state vector ant t the evaluation time.      !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real) t1                                             !
    !               ! Time at which the initial conditions are given.       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real) t2                                             !
    !               ! Time at which the solution is evaluated.              !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(N)) u1                                          !
    !               ! Initial conditions' state vector at time 't1'.        !
    !---------------!-------------------------------------------------------!
    !   OUT         ! (real(N)) u2                                          !
    !               ! Solution's state vector at time 't2'.                 !
    !---------------!-------------------------------------------------------!
    recursive subroutine euler_explicit(f, t1, t2, u1, u2, tol)
        procedure(ode_function) :: f
        real, intent(in) :: t1, t2, u1(:)
        real, intent(out) :: u2(size(u1))
        real, intent(in), optional :: tol
        real :: dt = 1e-1, u_half(size(u1))
                
        u2 = u1 + (t2 - t1) * f( u1, t1)
    end subroutine euler_explicit
    
    !-----------------------------------------------------------------------!
    !   ( SUBROUTINE ) euler_implicit                                       !
    !-----------------------------------------------------------------------!
    !   Extrapolates the solution of the Cauchy Problem at a given time     !
    !   and initial conditions by means of the Implicit/Inverse Euler       !
    !   numerical scheme.                                                   !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (ode_function) f(u,t)                                 !
    !               ! f: RN x R => RN                                       !
    !               !     u , t => du_dt                                    !
    !               ! Derivative function of the Cauchy Problem where       !
    !               ! u is the state vector ant t the evaluation time.      !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real) t1                                             !
    !               ! Time at which the initial conditions are given.       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real) t2                                             !
    !               ! Time at which the solution is evaluated.              !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(N)) u1                                          !
    !               ! Initial conditions' state vector at time 't1'.        !
    !---------------!-------------------------------------------------------!
    !   OUT         ! (real(N)) u2                                          !
    !               ! Solution's state vector at time 't2'.                 !
    !---------------!-------------------------------------------------------!
    recursive subroutine euler_implicit(f, t1, t2, u1, u2, tol)
        procedure(ode_function) :: f
        real, intent(in) :: t1, t2, u1(:)
        real, intent(out) :: u2(size(u1))
        real, intent(in), optional :: tol
        real :: dt = 1e-1, u_half(size(u1))
        
        call euler_explicit( f, t1, t2, u1, u2 )
        call newton( F = euler_implicit_equation, x0 = u2 )
        
        contains
    
        function euler_implicit_equation(x) result(g)
            real, intent(in) :: x(:) 
            real :: g( size(x) )
        
            g = x - u1 - (t2 - t1) * f( x, t2 )
        end function
        
        end subroutine euler_implicit
    
    !-----------------------------------------------------------------------!
    !   ( SUBROUTINE ) runge_kutta_4                                        !
    !-----------------------------------------------------------------------!
    !   Extrapolates the solution of the Cauchy Problem at a given time     !
    !   and initial conditions by means of the 4th order Runge Kutta        !
    !   numerical scheme.                                                   !
    !   This is a wrapper of the code written by John Burkardt's rk4vec.    !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (ode_function) f(u,t)                                 !
    !               ! f: RN x R => RN                                       !
    !               !     u , t => du_dt                                    !
    !               ! Derivative function of the Cauchy Problem where       !
    !               ! u is the state vector ant t the evaluation time.      !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real) t1                                             !
    !               ! Time at which the initial conditions are given.       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real) t2                                             !
    !               ! Time at which the solution is evaluated.              !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(N)) u1                                          !
    !               ! Initial conditions' state vector at time 't1'.        !
    !---------------!-------------------------------------------------------!
    !   OUT         ! (real(N)) u2                                          !
    !               ! Solution's state vector at time 't2'.                 !
    !---------------!-------------------------------------------------------!
    recursive subroutine runge_kutta_4(f, t1, t2, u1, u2, tol)
        procedure(ode_function) :: f
        real, intent(in) :: t1, t2, u1(:)
        real, intent(out) :: u2(size(u1))
        real, intent(in), optional :: tol
        real :: dt = 1e-1, u_half(size(u1))
        integer :: m
        
        m = size(u1)
        call rk4vec( t1, m, u1, t2 - t1, f, u2 )
    end subroutine runge_kutta_4
    
    !-----------------------------------------------------------------------!
    !   ( SUBROUTINE ) dopri853                                       !
    !-----------------------------------------------------------------------!
    !   Extrapolates the solution of the Cauchy Problem at a given time     !
    !   and initial conditions by means of the Explicit Euler               !
    !   numerical scheme.                                                   !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (ode_function) f(u,t)                                 !
    !               ! f: RN x R => RN                                       !
    !               !     u , t => du_dt                                    !
    !               ! Derivative function of the Cauchy Problem where       !
    !               ! u is the state vector ant t the evaluation time.      !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real) t1                                             !
    !               ! Time at which the initial conditions are given.       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real) t2                                             !
    !               ! Time at which the solution is evaluated.              !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(N)) u1                                          !
    !               ! Initial conditions' state vector at time 't1'.        !
    !---------------!-------------------------------------------------------!
    !   OUT         ! (real(N)) u2                                          !
    !               ! Solution's state vector at time 't2'.                 !
    !---------------!-------------------------------------------------------!
    recursive subroutine dopri853(f, t1, t2, u1, u2, tol)
        procedure(ode_function) :: f
        real, intent(in) :: t1, t2, u1(:)
        real, intent(out) :: u2(size(u1))
        real :: dt = 1e-1, u_half(size(u1)), t_end
        real, intent(in), optional :: tol
        type(dop853_class) :: dop_class
        
        integer :: n, iout, idi
        logical :: output
        
        real :: tols(size(u1));
        
        n = size(u1)
        t_end = t1
        u2 = u1
        !if ( present(tol) ) then
        !    tols = tol
        !else
        !    tols = 1d-3
        !endif
        tols = 1d-3
        iout = 0
        idi = 1
        
        call dop_class%initialize(n = n, fcn = sub_f, status_ok = output)
        call dop_class%integrate(t_end, u2, t2, tols, tols, iout, idi)
        
        contains
        
        subroutine sub_f(me,x,y,yp)
            class(dop853_class),intent(inout) :: me
            real, intent(in) :: x
            real, intent(in) :: y(:)
            real, intent(out) :: yp(:)
            yp = f(y,x)
        end subroutine
    end subroutine
    
end module temporal_schemes
 