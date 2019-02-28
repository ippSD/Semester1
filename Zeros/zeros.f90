program zeros
    use funcion
    use modulos
    use Non_Linear_Systems
    implicit none
    
    real :: xv(3) = [1d0, 2d0, 1d0];
    
    call Newton_Solution()
    
end program