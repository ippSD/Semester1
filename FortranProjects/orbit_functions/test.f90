program test
    use numerical_methods
    use cauchy
    use hamiltonian_functions 
    
    integer, parameter :: N = 10, M = 1000
    real, parameter :: DT = 1d-2
    
    integer :: i, j
    real :: u0(2,N), u(0:M,2,N), time(0:M)
    
    u0(1,:) = [(3d0*(i-1-N/2)/(N/2-1d0),i=1,N)]
    u0(2,:) = 0d0
    
    u(0,:,:) = u0
    time = [(i*DT,i=0,M)]
    
    do i=1,N
        call cauchy_problem(time_domain = time, differential_operator = example, temporal_scheme = runge_kutta, solution = u(0:M,:,i))
    end do
    
    open(unit=13,file="B1.dat")
    
    do i=1,N
        do j=0,M
            write(13,"(2(F6.3,A1),F6.3)") time(j), ',', u(j,1,i), ',', u(j,2,i)
        end do
        write(13,*) ""
    end do
    close(13)
        
end program