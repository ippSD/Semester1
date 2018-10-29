module functions
  implicit none
contains

  function oscilador(u,t) result(f)
    real, intent(in) :: u(0:), t
    real :: f(0:size(u))

    f = [u(1), -u(0)]
  end function

  function kepler(u,t) result(f)
    real, intent(in) :: u(0:), t
    real :: f(0:size(u))

    f = [u(2), u(3), -u(0)/norm2(u(0:1)), -u(1)/norm2(u(0:1))]
  end function
end module
