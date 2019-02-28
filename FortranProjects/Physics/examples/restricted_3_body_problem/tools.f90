!*********************************************************************
module tools
!*********************************************************************

implicit none


contains


!*************************************************************************
!*   Subrutina: closeall                                                 *
!*                                                                       *
!*       Cierra todas las unidades                                       *
!*************************************************************************

 subroutine create_file(file_name) 
   character(len=*), intent(in) :: file_name


     open(103, file=file_name ) 
     close(103)


 endsubroutine

!*************************************************************************
!*   Subrutina: closeall                                                 *
!*                                                                       *
!*       Cierra todas las unidades                                       *
!*************************************************************************

 subroutine closeall

 integer :: i

 do i = 1, 100
     close(i)
 end do

 endsubroutine

!*************************************************************************
!*   Subrutina: get_name                                                 *
!*                                                                       *
!*   Ofrece el name del programa que se está ejecutando sin extension  *
!*************************************************************************

 subroutine get_name(name)

 character, intent(out)  :: name*(*)

 character   :: name_temporal*400
 integer     :: i, j

 call getarg(0,name_temporal)

 i = SCAN (name_temporal, '\' , .true.)
 
 j = len_trim(name_temporal)

 name = name_temporal(i+1:j-4)

 endsubroutine

!*************************************************************************
!*   Subrutina: intcha                                                   *
!*                                                                       *
!*       Convierte un entero en un caracter                              *
!*                                                                       *
!*   Entrada                                                             *
!*       -numero: numero a convertir                                     *
!*   Salida                                                              *
!*       -numero_caracteres: espacios que ocupa el numero                *
!*       -texto: el caracter resultante                                  *
!*************************************************************************

 subroutine intcha(numero, numero_caracteres, texto)
 
 integer, intent(in) :: numero

 integer, intent(out)    :: numero_caracteres
 character, intent(out)  :: texto*(*)

 integer     :: i, j, k, numero_alterado

 if (numero == 0) then

     numero_caracteres = 1
     texto(1:1) = "0"

 elseif (numero == 1) then

     numero_caracteres = 1
     texto(1:1) = "1"

 else

     numero_caracteres = int( log10(real(numero)) ) + 1

     j = 1

     numero_alterado = numero

     do i = numero_caracteres, 1, -1
         k = int( numero_alterado / (10.d0**real(i-1)) )
         selectcase (k)

             case(1)
             texto(j:j) = "1"

             case(2)
             texto(j:j) = "2"

             case(3)
             texto(j:j) = "3"

             case(4)
             texto(j:j) = "4"

             case(5)
             texto(j:j) = "5"

             case(6)
             texto(j:j) = "6"

             case(7)
             texto(j:j) = "7"

             case(8)
             texto(j:j) = "8"

             case(9)
             texto(j:j) = "9"

             case(0)
             texto(j:j) = "0"

         endselect
         j = j+1
         numero_alterado = numero_alterado - k*int(10.d0**real(i-1))
     enddo

 endif

!*************************************************************************
!*   Subrutina: crea_lista                                               *
!*                                                                       *
!*       Crea un caracter con una lista con separacion por comas         *
!*************************************************************************

 endsubroutine

 subroutine crea_lista(dimension_lista, lista_variables, texto1, texto2)
 
 integer, intent(in)     :: dimension_lista
 character, intent(in)   :: texto1*(*), texto2*(*)
 character, intent(out)  :: lista_variables*(*)

 integer     :: i, numero_caracteres, numero_texto1, numero_texto2
 character   :: indice_variable*2

 numero_texto1 = len(texto1)
 numero_texto2 = len(texto2)
 lista_variables = texto1

 do i = 1, dimension_lista
     lista_variables = adjustl(lista_variables)

     if (i < 10 ) then
         call intcha(i, numero_caracteres, indice_variable)
         indice_variable = adjustl(indice_variable)
         lista_variables = lista_variables(1:(numero_texto1+(i-1)*(6+numero_texto2)))//&
                           ",    "//texto2//indice_variable(1:numero_caracteres)
     else
         call intcha(i, numero_caracteres, indice_variable)
         indice_variable = adjustl(indice_variable)
         lista_variables = lista_variables(1:(numero_texto1+9*(6+numero_texto2)+(i-10)*(7+numero_texto2)))//&
                           ",    "//texto2//indice_variable(1:numero_caracteres)
     endif

 enddo

 endsubroutine

!*************************************************************************
!*   Subrutina: busca_name                                             *
!*                                                                       *
!*       Busca el primer fichero no existente de una lista que tiene la  *
!*   siguiente forma: AAA.ext, AAA1.ext, AAA2.ext,...                    *
!*************************************************************************

 subroutine busca_name(name_fichero_, name_busca, extension)

 character,intent(inout) :: name_fichero_*(*)
 character,intent(in)    :: extension*(*)
 character,intent(out)   :: name_busca*(*)

 logical             :: existe
 integer             :: indice_archivo = 0 

 do
     call crea_name(name_fichero_, indice_archivo, extension, name_busca)
     inquire (file = name_busca, Exist=existe)
     if (existe == .false. ) exit
     indice_archivo = indice_archivo + 1
 enddo

 end subroutine

!*************************************************************************
!*   Subrutina: crea_name                                              *
!*                                                                       *
!*       Crea un nobre a partir de un name general(AAA) un numero(n) y *
!*   una extension(ext) de la siguiente forma: AAAn.ext; para n = 0, el  *
!*   name sera AAA.ext                                                 *
!*************************************************************************

 subroutine crea_name(name_fichero_, numero, extension, name_creado)

     character,intent(in)    :: name_fichero_*(*)
     character,intent(in)    :: extension*(*)
     integer, intent(in)     :: numero
     character, intent(out)  :: name_creado*(*)

     integer     :: numero_caracteres
     character   :: indice_texto*2

     if (numero == 0) then
         name_creado = trim(name_fichero_)//extension
     else
         call intcha(numero, numero_caracteres, indice_texto)
         indice_texto = adjustl(indice_texto)
         name_creado = trim(name_fichero_)//&
                         indice_texto(1:numero_caracteres)//extension
     endif


 endsubroutine






!*************************************************************************
!*   Subrutina: lee_error                                                *
!*                                                                       *
!*       Lee el error de una solucion del problemas de evolucion         *
!*   devolviendo el vector error                                         *
!*************************************************************************

 subroutine lee_error_paso(name_fichero,extension,n_frames_inicial,i,error)

 implicit none

 character, intent(inout)    :: name_fichero*(*)
 character, intent(in)   :: extension*(*)
 integer, intent(in)     :: i, n_frames_inicial
 real, intent(out)       :: error(:)
 
 character               :: name_fichero_hst*26
 integer                 :: k, l, m, unidad_lee, &
                            numero_caracteres, numero
 real,allocatable        :: datos(:,:,:), vector_error(:)
 real                    :: tiempo

 allocate (datos(size(error),n_frames_inicial+1,2),&
           vector_error(n_frames_inicial+1) )

 unidad_lee = 45
 
 name_fichero = adjustr(name_fichero)
 
 do k = 1, 2
     call crea_name(name_fichero, (i-1+k), extension, name_fichero_hst)
     open(unidad_lee, file=name_fichero_hst)
     read(unidad_lee, *)
     read(unidad_lee, *) tiempo, datos(:,1,k)

     do l = 2, n_frames_inicial + 1
         numero = (2**(i-1+k))-1

         do m = 1, numero
             read(unidad_lee,*)
         enddo
         
         read(unidad_lee, *) tiempo, datos(:,l,k)
     enddo

     close(unidad_lee)
 enddo
 
 do k = 1, size(error)
     
     vector_error(:) = datos(k,:,1) - datos(k,:,2)
     error(k) = dsqrt( sum( vector_error(:) ** 2.d0) )
     
 enddo 

 endsubroutine

 subroutine reduce_polinomio(A, b)

     complex, intent(inout)  :: A(:), b
     
     complex :: C,D
     integer :: n,i

     n = size(A)

     C = A(n-1) + A(n)*b
     A(n-1) = A(n)

     do i = n-2,1,-1

         D = A(i) + C*b
         A(i) =  C
         C = D

     enddo

 endsubroutine

 subroutine calcula_polinomio(A,x,y)

     complex, intent(in) :: A(:), x
     complex, intent(out):: y

     integer :: n,i

     n = size(A)

     y = 0.d0

     do i = 1,n

         y = y + A(i) * (x**(i-1.d0))

     enddo

 endsubroutine

 subroutine deriva_polinimio(A,x,dydx)

     complex, intent(in) :: A(:), x
     complex, intent(out):: dydx

     integer :: n, i

     n = size(A)

     dydx = 0.d0

     do i = 1,n

         dydx = dydx + (i-1.d0) * A(i) * (x**(i-2.d0))

     enddo

 endsubroutine

 subroutine raices_polinomio(A, raices)

     complex, intent(inout) :: A(:)
     complex, intent(inout):: raices(:)

     integer :: n, i
     
     n = size(A)

     do i = 1, n-1

         call resuelve_raiz(A(1:n+1-i), raices(i))
         call reduce_polinomio(A(1:n+1-i),raices(i))

     enddo

 endsubroutine

 subroutine resuelve_raiz(A,raiz)

     complex, intent(in) :: A(:)
     complex, intent(inout):: raiz


     complex :: y, dydx
     integer :: n, i
     
     n = size(A)

     call calcula_polinomio(A,raiz,y)

     do i=1, 1000

         call deriva_polinimio(A,raiz,dydx)

         raiz = raiz - y/dydx

         call calcula_polinomio(A,raiz,y)

         if (abs(y) < 1.d-8) exit

         if (i==1000) then

             write(*,*) "Numero maximo de iteraciones alcanzado"
             read(*,*)
             stop

         endif

     enddo

 endsubroutine


 subroutine busca_en_lista_cmplx(A,b,flag)

     complex, intent(in) :: A(:), b
     logical, intent(out):: flag

     integer :: n,i

     n = size(A)
     flag = .FALSE.

     do i = 1, n

         if (b==A(i)) flag=.TRUE.

     enddo

 endsubroutine

end module
