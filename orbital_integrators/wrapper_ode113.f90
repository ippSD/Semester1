module wrapper_ode113
    use ode_interfaces
    implicit none

    contains
    
    subroutine ode113(f, t1, t2, u1, u2)
        procedure(ode_function) :: f
        real, intent(in) :: t1, t2, u1(:)
        real, intent(out) :: u2(size(u1))
        real :: dt = 1e-1, u_half(size(u1))
          
        u2 = u1
        u3 = u2(1:3*N)
        call ode(f_sub,n,u3,t1,t2)
        u2(1:3*N) = u3(1:3*N)
        
        contains
    
        subroutine f_sub(t, y, yp)
        real, intent(in) :: t, y
        real, intent(out) :: yp
        
        yp = f(y,t)
        end subroutine f_sub
        
    end subroutine ode113
    
end module wrapper_ode113