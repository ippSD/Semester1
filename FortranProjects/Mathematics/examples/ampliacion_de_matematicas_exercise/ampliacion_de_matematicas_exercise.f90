!---------------------------------------------------------------------------!
!   ampliacion_de_matematicas_exercise                                      !
!---------------------------------------------------------------------------!
!   Teacher told us to sketch the solution of a hamiltonian function        !
!   based on its stability.  Well, why don't we solve it numerically?       !
!---------------------------------------------------------------------------!

program ampliacion_de_matematicas_exercise
    use temporal_schemes, only: runge_kutta_4
    use cauchy_problem_solver, only: cauchy_problem
    use hamiltonian_functions 
    use dislin_mod
    implicit none
    
    integer, parameter :: N = 10, M = 1000
    real, parameter :: DT = 1d-2
    
    integer :: i, j
    real :: u0(2,N), u(0:M,2,N), time(0:M)
    
    u0(1,:) = [(3d0*(i-1-N/2)/(N/2-1d0),i=1,N)]
    u0(2,:) = 0d0
    
    u(0,:,:) = u0
    time = [(i*DT,i=0,M)]
    
    do i=1,N
        call cauchy_problem( &
            time_domain = time, &
            differential_operator = example, &
            temporal_scheme = runge_kutta_4, &
            solution = u(0:M,:,i) &
        )
    end do
    
    
    
    open(unit=13,file="B1.dat")
    
    do i=1,N
        call plot(u(:,1,i), u(:,2,i), hold_on = .true.)
        do j=0,M
            write(13,"(2(F6.3,A1),F6.3)") time(j), ',', u(j,1,i), ',', u(j,2,i)
        end do
        write(13,*) ""
    end do
    close(13)
    
    call plot_title("Bombilla")
    call plot_end()
        
end program ampliacion_de_matematicas_exercise