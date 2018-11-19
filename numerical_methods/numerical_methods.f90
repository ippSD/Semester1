module numerical_methods
    use rkf45
    use interfaces
    implicit none

    contains
    
    recursive subroutine euler_explicit(f, t1, t2, u1, u2)
        !  Parameters:
        !    F, a function of the form 
        !       f ( u, t ) result(uprime)
        !    which evaluates the derivative UPRIME given the time T and
        !    solution scalar U.
        !    Input, real ( kind = 8 ) T1, the current time.
        !    Input, real ( kind = 8 ) T2, the end time.
        !    Input, real ( kind = 8 ) U1, the solution estimate at the current time.
        !    Output, real ( kind = 8 ) U2, the solution estimate at time T2.
    
        procedure(odes) :: f
        double precision, intent(in) :: t1, t2, u1(:)
        double precision, intent(out) :: u2(size(u1))
        double precision :: dt = 1e-1, u_half(size(u1))
        
        !allocate(u2(lbound(u1):ubound(u1))
        
        u2 = u1 + (t2 - t1) * f( u1, t1)
    end subroutine euler_explicit
    
    recursive subroutine runge_kutta(f, t1, t2, u1, u2)
        procedure(odes) :: f
        double precision, intent(in) :: t1, t2, u1(:)
        double precision, intent(out) :: u2(size(u1))
        double precision :: dt = 1e-1, u_half(size(u1))
        integer :: m
        
        m = size(u1)
        call rk4vec( t1, m, u1, t2 - t1, f, u2 )
    end subroutine runge_kutta
    
end module
 