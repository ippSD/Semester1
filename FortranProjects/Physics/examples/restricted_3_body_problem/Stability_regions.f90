module Stability_regions 

    use problema_cauchy
    use tools

implicit none

    contains

!***************************************************************************
!   Stability region                                        
!                                                                       
!It determines the Stability of absolute stability
! of an ODE numerical scheme    
!*************************************************************************

subroutine Region_Estabilidad_Absoluta(Esquema_temporal, Dominio)

interface    !  especificacion de parametros que son procedimentos

subroutine Esquema_temporal(Sistema_EDOS, t1, t2,  U,ierr)
    interface 
    subroutine Sistema_EDOS(t, U, F)  ! Calcula F(U, t) 
        real, intent(in)  :: t
        real, intent(inout)  ::  U(:)
        real, intent(out) :: F(:)
    end subroutine Sistema_EDOS   
    end interface 
    real, intent(inout)    :: t1, t2, U(:)
    integer, intent(out)   :: ierr                
end subroutine Esquema_temporal 

subroutine Dominio(x,y)
    real, pointer   :: x(:,:), y(:,:)
endsubroutine

end interface

integer :: i, j, k, l, ierr
integer :: N(2), pasos
complex :: unidad_img=(0.d0,1.d0), w
real    :: t0 = 0.d0, t1 = 1.d0, pi,escala
real, allocatable   :: Est_Val(:, :)
complex, allocatable    :: AutoVal(:), A(:)
real, pointer   :: U(:), x(:,:), y(:,:)

open(10, file = "Region_estabilidad.img")
write(10,*) "x,   y,   z"
write(10,*) "tiempo=", 0.d0

call Dominio(x,y)

N(1) = size( x, dim=1 )-1
N(2) = size( x, dim=2 )-1

call determina_pasos(pasos,Esquema_temporal)

allocate (U(2), A(0:pasos), Est_Val(0:N(1),0:N(2)),Autoval(pasos) )

pi = dacos(-1.d0)
A = -1.d0
do i = 0,N(1)
    do j = 0, N(2)

        w = cmplx ( x(i,j), y(i,j) )  
   
        do k = 1, pasos
            do l = 1, pasos

                if (k == l) then                        
                    U(1) = 1.0d0
                    U(2) = 0.0d0
                else
                    U = 0.d0
                endif

                call Esquema_temporal(EDO_RA, t0, t1,  U, ierr)

            enddo

            A(k-1) = cmplx ( U(1), U(2) ) 
        enddo
    
        ! Inicializamos el vector de autovalores
    
        
        do k=1, pasos     
            autoVal(k) = exp( 2*pi*(k-1)*unidad_img/pasos )
        end do     
                    
        call raices_polinomio(A, AutoVal)
            
        Est_Val(i, j) = maxval(abs(AutoVal))

        write(10,*) x(i,j), y(i,j),  Est_Val(i, j)
        
    enddo   
enddo

close(10)

contains

!********************************************************************************************
    subroutine EDO_RA(t,U,F)
!********************************************************************************************

        real, intent(in)  :: t
        real, intent(inout)  ::  U(:)
        real, intent(out) :: F(:)

        complex ::  Fc

            call EDO_RAc( cmplx(U(1), U(2)), Fc ) 

            F(1:2) = (/ real(Fc), imag(Fc) /) 

    endsubroutine

!*****************************************************************************************
    subroutine EDO_RAc(U, F)
!*****************************************************************************************
          complex, intent(in)  ::  U
          complex, intent(out) :: F 

        F = w * U 

    endsubroutine

endsubroutine









!*****************************************************************************************
subroutine determina_pasos(pasos,Esquema_temporal)
!*****************************************************************************************

interface

subroutine Esquema_temporal(Sistema_EDOS, t1, t2,  U,ierr)
    interface 
    subroutine Sistema_EDOS(t, U, F)  ! Calcula F(U, t) 
        real, intent(in)  :: t
        real, intent(inout)  ::  U(:)
        real, intent(out) :: F(:)
    end subroutine Sistema_EDOS   
    end interface 
    real, intent(inout)    :: t1, t2, U(:)
    integer, intent(out)   :: ierr                
end subroutine Esquema_temporal 

end interface

integer, intent(out)    :: pasos

real, pointer   :: U(:)
real    :: t0 = 0.d0, t1=1.d0
integer :: i, ierr

allocate(U(1))

i = 0

U = 1.d0


do      
    call Esquema_temporal(EDO_Lineal, t0, t1,  U, ierr)
    
    if ( abs(U(1)) <= epsilon(1d0) ) then

        pasos = i
        exit

    endif

    U = 0.d0
    i = i+1
enddo

contains

!********************************************************************************************
    subroutine EDO_Lineal(t,U,F)
!********************************************************************************************

        real, intent(in)  :: t
        real, intent(inout)  ::  U(:)
        real, intent(out) :: F(:)

        
            F = -U 

    end subroutine

end subroutine

end module
