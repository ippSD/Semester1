!***********************************************************************
module Problema_Cauchy
!***********************************************************************


implicit none 

public :: Problema_evolucion 

contains    

!***************************************************************************
! Integra  el siguiente problema de Cauchy 
!
!       dU/dt = F(U, t, parametros),      U(0) = U^0 (CI) 
!
!       U^{n+1} = G( U^n... U^{n-1+p}, F^n... F^{n-1+p}, dt, parametros ) 
!
!    donde  
!              Sistema            :   F(U,t, parametros) 
!              Condicion_Inicial  :   U(0) 
!              Paso_temporal      :   dt 
!              Esquema_temporal   :   G(....)  
!
!***************************************************************************
subroutine Problema_evolucion( Dominio,  Condicion_Inicial, Sistema_EDOS,   &
                                Esquema_temporal, Salidas) 
        
     interface    !  especificacion de parametros que son procedimentos

        subroutine  Dominio( Intervalo) ! dominio de integracion  
                       real, pointer :: Intervalo(:) 
        end subroutine Dominio  

        subroutine  Condicion_Inicial( t0, U) ! condicion inicial 
                       real, intent(in) :: t0 
                       real, pointer :: U(:)
        end subroutine Condicion_Inicial 

        subroutine  Sistema_EDOS(t, U,  F)            ! Calcula F(U)  
              real, intent(in)  :: t
              real, intent(inout)  ::  U(:)
              real, intent(out) :: F(:)
        end subroutine Sistema_EDOS  

        
        subroutine Esquema_temporal(Sistema_EDOS, t1, t2,  U, ierr )
                interface 
                      subroutine Sistema_EDOS(t, U, F)  ! Calcula F(U, t) 
                           real, intent(in)  :: t
                           real, intent(inout)  ::  U(:)
                           real, intent(out) :: F(:)
                      end subroutine Sistema_EDOS   
                end interface 
                real, intent(inout)    :: t1, t2, U(:)
                integer, intent(out) :: ierr  
              
       end subroutine Esquema_temporal 

       subroutine Salidas(Sistema_EDOS, t, U) ! Salidas   
                 interface 
                      subroutine Sistema_EDOS(t, U, F)  ! Calcula F(U, t) 
                           real, intent(in)  :: t
                           real, intent(inout)  ::  U(:) 
                           real, intent(out) :: F(:)
                      end subroutine Sistema_EDOS   
                 end interface 

                 real, intent(in) :: t, U(:)
           
       end subroutine Salidas  
     end interface 
 

!  *** especificacion de  variables
       real :: t0, t1, t2 
       integer ::  i, n_pasos, ierr  

       real, pointer :: U(:), Intervalo(:)    

!  *** condicion inicial y dominio de integracion 
       call Dominio( Intervalo ) 
       t0 = Intervalo(1) 
       n_pasos = size(Intervalo) - 1 
        
       call Condicion_inicial ( t0, U )

       call Salidas(  Sistema_EDOS, t0, U )

!  *** bucle para la integracion temporal
       do i=1, n_pasos 

          t1 = Intervalo(i) 
          t2 = Intervalo(i+1)
         
          call Esquema_Temporal( Sistema_EDOS, t1, t2, U,  ierr )  ! integra entre t1 y t2 
          call Salidas(  Sistema_EDOS, t2, U ) 

          if (ierr>0) exit 
             
             
       enddo 
        

       deallocate( U ) 

       deallocate( Intervalo )

                 
   end subroutine 

  
end module 
