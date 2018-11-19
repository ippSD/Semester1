module interfaces
    implicit none
    
    abstract interface
    
        function odes(u, t) result(f)
          double precision, intent(in) :: u(:),  t
          double precision :: f(size(u))
        end function
        
        subroutine f_s_schemes ( t, u, uprime )
            double precision,intent(in) :: t
            double precision, intent(in) :: u
            double precision, intent(out) :: uprime
        end subroutine f_s_schemes
        
        subroutine f_v_schemes ( t, m, u, uprime )
            double precision,intent(in) :: t
            integer, intent(in) :: m
            double precision, intent(in) :: u(m)
            double precision, intent(out) :: uprime(m)
        end subroutine f_v_schemes
    end interface
    
    abstract interface

        subroutine scheme(f, t1, t2, u1, u2)
            double precision, intent(in) :: t1, t2, u1(:)
            double precision, intent(out) :: u2(size(u1))
            interface
                function f(u, t)
                  double precision, intent(in) :: u(:),  t
                  double precision :: f(size(u))
                end function
            end interface
        end subroutine scheme
        
    end interface
    
end module interfaces