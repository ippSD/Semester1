module hamiltonian_functions
  implicit none
  contains
    
    function example(u,t) result(f)
        real, intent(in) :: u(:), t
        real :: f(size(u))
    
        real :: x , p
        
        x = u(1)
        p = u(2)

        f = [p, 3d0*x-x**3d0+2d0]
    end function
    
    function kepler(u,t) result(f)
        real, intent(in) :: u(:), t
        real :: f(size(u))
    
        real :: x , p
        
        x = u(1)
        p = u(2)

        f = [p, -x]
    end function

end module hamiltonian_functions
