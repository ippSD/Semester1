module interfaces
    implicit none
    
    abstract interface
    
        function odes(u, t) result(f)
          real, intent(in) :: u(:),  t
          real :: f(size(u))
        end function
        
        subroutine f_s_schemes ( t, u, uprime )
            real, intent(in) :: t
            real, intent(in) :: u
            real, intent(out) :: uprime
        end subroutine f_s_schemes
        
        subroutine f_v_schemes ( t, m, u, uprime )
            real, intent(in) :: t
            integer, intent(in) :: m
            real, intent(in) :: u(m)
            real, intent(out) :: uprime(m)
        end subroutine f_v_schemes
    end interface
    
    abstract interface

        subroutine scheme(f, t1, t2, u1, u2)
            real, intent(in) :: t1, t2, u1(:)
            real, intent(out) :: u2(size(u1))
            interface
                function f(u, t)
                  real, intent(in) :: u(:),  t
                  real :: f(size(u))
                end function
            end interface
        end subroutine scheme
        
    end interface
    
end module interfaces