module modulos
    implicit none
    contains
    
    subroutine Newton_Solution()
    
        real :: x0(3) = [1., 1., 1. ]
        
        call Newton( F, x0 )
        
        write(*,*) 'Zeroes of F(x) are x = ', x0
        
    end subroutine
end module