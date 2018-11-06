module numerical_methods
    use rkf45
    implicit none
    
    interface runge_kutta
        module procedure rk4, rk4vec
    end interface

    contains
end module
 