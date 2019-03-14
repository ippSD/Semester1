module module_odex
    use ode_interfaces, only: ode_function
    
    private
    public :: odex_int
    contains
    
    include "External/odex.f"
    
    subroutine odex_int(f, t1, t2, u1, u2)
        procedure(ode_function) :: f
        real, intent(in) :: t1, t2, u1(:)
        real, intent(out) :: u2(size(u1))
        real :: dt = 1e-1, u_half(size(u1)), work(30)
        real :: atol(size(u1)), rtol(size(u1))
        integer :: n, ipar(size(u1)), idid, liwork = 30, iwork(30), lwork = 30
        real :: rpar(size(u1))
        
        atol = 1d-5
        rtol = 1d-5
        n = size(u1)
        u2 = u1
        
        call ODEX(n,f_sub,t1,u2,t2,1d-3,&
                rtol,atol,0,&
                SOLOUT,0,&
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
        
        contains
    
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,NCON,ICOMP,ND, RPAR,IPAR,IRTRN)
            integer :: nr, n, ncon, nd, ipar, irtrn
            integer :: con(N), icomp(nd)
            real :: xold, x, y(N), rpar
        END SUBROUTINE
        
        subroutine f_sub(n,x,y,yp,rpar,ipar)
            integer, intent(in) :: n, ipar
            real, intent(in) :: x, y(:), rpar
            real, intent(out) :: yp(:)
            
            yp = f(y,x);
        end subroutine f_sub
    end subroutine
end module module_odex