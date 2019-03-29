module module_odex
    use ode_interfaces, only: ode_function
    implicit none
    
    private
    public :: odex_int
    contains
    
    include "External/odex.f"
    
    subroutine odex_int(f, t1, t2, u1, u2, ttol)
        procedure(ode_function) :: f
        real, intent(in) :: t1, t2, u1(:)
        real, intent(out) :: u2(size(u1))
        real, intent(in), optional :: ttol
        integer :: liwork, lwork
        real, allocatable :: work(:)
        integer, allocatable :: iwork(:)
        real :: dt = 1e-1, u_half(size(u1))
        real :: atol(size(u1)), rtol(size(u1))
        integer :: n, ipar(1), idid, km, nrdens, iout, itol
        real :: rpar(1)
        
        iout = 0
        itol = 0
        n = size(u1)
        km = 9
        nrdens = 0
        liwork = 2*km+21+nrdens +1
        lwork = n*(km+5)+5*km+20+(2*km*(km+2)+5)*nrdens +1
        allocate(work(lwork), iwork(liwork))
        iwork = 0
        work = 0d0
        !write(*,*) present(ttol)
        !if ( present(ttol) ) then
        !    write(*,*) ttol
        !    atol = ttol
        !    rtol = ttol
        !else
        !    atol = 1d-3
        !    rtol = 1d-3
        !endif
        atol = 1d-3
        rtol = 1d-3
        
        n = size(u1)
        u2 = u1
        
        call odex(&
            n,&
            f_sub,&
            t1,&
            u2,&
            t2,&
            dt,&
            rtol,&
            atol,&
            itol,&
            SOLOUT,&
            iout,&
            work,&
            lwork,&
            iwork,&
            liwork,&
            rpar,&
            ipar,&
            idid)
        
        contains
    
        subroutine solout (nr,xold,x,y,n,con,ncon,ICOMP,ND, RPAR,IPAR,IRTRN)
            integer :: nr, n, ncon, nd, ipar, irtrn
            integer :: con(N), icomp(nd)
            real :: xold, x, y(N), rpar
        END SUBROUTINE
        
        subroutine f_sub(n,x,y,yp,rpar,ipar)
            integer, intent(in) :: n, ipar
            real, intent(in) :: x, y(n), rpar
            real, intent(out) :: yp(n)
            
            yp = f(y,x);
        end subroutine f_sub
    end subroutine odex_int
end module module_odex