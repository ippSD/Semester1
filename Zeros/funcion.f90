module funcion
    implicit none
    contains
    function F(xv)
    
        real , intent(in) :: xv(:)
        real:: F(size(xv))
        
        real :: x, y, z
        
        x = xv(1)
        y = xv(2)
        z = xv(3)
        
        F(1) = x**2 - y**3 - 2
        F(2) = 3 * x * y - z
        F(3) = z**2 - x
        
    end function
end module