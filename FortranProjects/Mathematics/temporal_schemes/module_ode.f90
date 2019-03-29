module module_ode
    use ode_interfaces, only: ode_function
    private
    public :: ode_int
    contains
    include "External/ode.f"
    
    subroutine ode_int(f, t1, t2, u1, u2, tol)
        procedure(ode_function) :: f
        real, intent(in) :: t1, t2, u1(:)
        real, intent(out) :: u2(size(u1))
        real ::work(100+21*size(u1)+1)
        real, intent(in), optional :: tol
        real :: t, atol, rtol
        integer :: n, iflag, iwork(5)
        real :: rpar(size(u1)), u_temp(size(u1))
        
        if ( present(tol) ) then
            atol = tol
            rtol = tol
        else
            atol = 1d-3
            rtol = 1d-3
        endif
        !atol = 1d-3
        !rtol = 1d-3
        t = t1
        iflag = 1
        n = size(u1)
        u2 = u1
        
        call ode(sub_f,n,u2,t,t2,rtol,atol,iflag,work,iwork)
        !if ( abs(iflag) /= 2 ) write(*,*) "Error: bad flag: ", iflag
        if ( abs(iflag) == 6 ) then
            write(*,*) "Error: iflag = 6"
            write(*,*) "T - TEND: ", t1 - t2
        endif
        
        contains
        
        subroutine sub_f(t,u,du_dt)
            real, intent(in) :: t, u(n)
            real, intent(out) :: du_dt(size(u))
            du_dt = f(u,t)
        end subroutine sub_f
        
        
    end subroutine
end module module_ode