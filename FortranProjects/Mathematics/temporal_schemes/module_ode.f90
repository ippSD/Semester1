module module_ode
    use ode_interfaces, only: ode_function
    private
    public :: ode_int
    contains
    include "External/ode.f"
    
    subroutine ode_int(f, t1, t2, u1, u2)
        procedure(ode_function) :: f
        real, intent(in) :: t1, t2, u1(:)
        real, intent(out) :: u2(size(u1))
        real :: dt = 1e-1, u_half(size(u1)), work(30)
        real :: atol, rtol
        integer :: n, ipar(size(u1)), iflag = 1, liwork = 30, iwork(30), lwork = 30
        real :: rpar(size(u1))
        
        atol = 1d-5
        rtol = 1d-5
        n = size(u1)
        u2 = u1
        
        call ode(f_sub,n,u2,t1,t2,rtol,atol,iflag,work,iwork)
        write(*,*) iflag
        write(*,*)u2(2) - u1(2)
        
        
        contains
        
        subroutine f_sub(t,y,yp)
            real, intent(in) :: t, y(:)
            real, intent(out) :: yp(:)
            
            yp = f(y,t);
        end subroutine f_sub
    end subroutine
end module module_ode