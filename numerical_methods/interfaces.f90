module interfaces
    implicit none
    
    interface
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
    
end module