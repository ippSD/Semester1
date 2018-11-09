module euler_explicit
    !Module with Explicit Euler numerical Scheme
    use interfaces
    implicit none
    
    contains
    
    subroutine euler_explicit_s( t0, u0, dt, f, u )
        !  Parameters:
        !
        !    Input, real ( kind = 8 ) T0, the current time.
        !
        !    Input, real ( kind = 8 ) U0, the solution estimate at the current time.
        !
        !    Input, real ( kind = 8 ) DT, the time step.
        !
        !    Input, external F, a subroutine of the form 
        !      subroutine f ( t, u, uprime ) 
        !    which evaluates the derivative UPRIME given the time T and
        !    solution scalar U.
        !
        !    Output, real ( kind = 8 ) U, the fourth-order Runge-Kutta solution 
        !    estimate at time T0+DT.
        double precision, intent(in) :: t0
        double precision, intent(in) :: u0
        double precision, intent(in) :: dt
        procedure(f_s_schemes) :: f
        double precision, intent(out) :: u
        
        double precision :: u_inter
        
        call f(t0, u0, u_inter)
        u = u_inter  + u0
    end subroutine euler_explicit_s
    
    subroutine euler_explicit_v( t0, m, u0, dt, f, u )
        !  Parameters:
        !
        !    Input, real ( kind = 8 ) T0, the current time.
        !
        !    Input, integer ( kind = 4 ) M, the dimension of the system.
        !
        !    Input, real ( kind = 8 ) U0(M), the solution estimate at the current time.
        !
        !    Input, real ( kind = 8 ) DT, the time step.
        !
        !    Input, external F, a subroutine of the form 
        !      subroutine f ( t, m, u, uprime ) 
        !    which evaluates the derivative UPRIME(1:M) given the time T and
        !    solution vector U(1:M).
        !
        !    Output, real ( kind = 8 ) U(M), the fourth-order Runge-Kutta solution 
        !    estimate at time T0+DT.
        double precision, intent(in) :: t0
        integer, intent(in) :: m
        double precision, intent(in) :: u0(m)
        double precision, intent(in) :: dt
        procedure(f_v_schemes) :: f
        double precision, intent(out) :: u(m)
        
        double precision :: u_inter(m)
        
        call f(t0, m, u0, u_inter)
        u = u_inter  + u0
    end subroutine euler_explicit_v
    
end module euler_explicit