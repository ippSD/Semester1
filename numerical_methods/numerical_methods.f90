module numerical_methods
    use rkf45
    use euler_explicit
    implicit none
    
    interface runge_kutta
        module procedure rk4, rk4vec
    end interface
    
    interface e_exp
        module procedure euler_explicit_s, euler_explicit_v
    end interface e_exp

    contains
end module
 