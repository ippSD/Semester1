module objects
    implicit none
    
    interface Body
        module procedure body_constructor
    end interface Body
    
    interface Galaxy
        module procedure galaxy_constructor
    end interface Galaxy
    
    type Body
        double precision :: mass
        double precision, dimension(3) :: r, v, a
    contains
    end type Body
    
    type :: BodyPointer
        double precision, pointer :: mass
        double precision, pointer :: r(:), v(:)
        double precision, dimension(3) :: a
    end type BodyPointer
    
    type Galaxy
        type(Body), allocatable :: bodies(:)
    contains
        procedure, pass :: calc_accelerations
    end type Galaxy
    
    type ::  GalaxyPointer
        type(BodyPointer), allocatable :: bodies(:)
    contains
        procedure, pass :: set_pointer
        procedure, pass :: calc_accelerations_pointer
    end type GalaxyPointer
    
    contains
    
    function body_constructor(init_vector) result(this)
        double precision, dimension(7), intent(in) :: init_vector
        type(Body) :: this
        
        this%mass = init_vector(1)
        this%r = init_vector(2:4)
        this%v = init_vector(5:7)
        this%a = [0d0, 0d0, 0d0]
        
    end function body_constructor
    
    function galaxy_constructor(init_vector) result(this)
        double precision, intent(in) :: init_vector(:)
        type(Galaxy) :: this
        
        integer :: i, n, s
        double precision :: mass_i
        double precision, dimension(3) :: r_i, v_i
        
        s = size(init_vector)
        
        if (mod(s,7) /= 0) then
            write(*,*) "Error: init vector of size /= 7N."
            call exit(1)
        end if
        
        n = s / 7
        allocate(this%bodies(n))
        
        do i = 1, n
            mass_i = init_vector(i)
            r_i = init_vector(1*n+3*(i-1)+1:1*n+3*(i-1)+3)
            v_i = init_vector(4*n+3*(i-1)+1:4*n+3*(i-1)+3)
            this%bodies(i) = Body([mass_i, r_i, v_i])
        end do
        
    end function galaxy_constructor
    
    
    subroutine calc_accelerations(this)
        class(Galaxy) :: this
        
        integer :: i, j, n
        double precision, dimension(3) :: a_i
        type(Body) :: body_i, body_j
        
        n = size(this%bodies)/7
        do i = 1, n
            body_i = this%bodies(i)
            a_i = 0d0
            do j = 1, n
                if (j /= i) then
                    body_j = this%bodies(j)
                    a_i = a_i - body_j%mass * (body_j%r - body_i%r) / norm2(body_j%r - body_i%r) ** 3d0
                end if
            end do
            body_i%a = a_i
        end do
    end subroutine calc_accelerations
    
    subroutine set_pointer(this, t)
        class(GalaxyPointer) :: this
        double precision, target, intent(in) :: t(:)
        
        integer :: i, n
        
        n = size(t) / 7
        allocate(this%bodies(n))
        
        do i=1, n
            this%bodies(i)%mass => t(i)
            this%bodies(i)%r => t(1*n+3*(i-1)+1:1*n+3*(i-1)+3)
            this%bodies(i)%v => t(4*n+3*(i-1)+1:4*n+3*(i-1)+3)
        end do
    end subroutine set_pointer
    
    subroutine calc_accelerations_pointer(this)
        class(GalaxyPointer) :: this
        
        integer :: i, j, n
        double precision, dimension(3) :: a_i
        type(BodyPointer) :: body_i, body_j
        
        n = size(this%bodies)/7
        do i = 1, n
            body_i = this%bodies(i)
            a_i = 0d0
            do j = 1, n
                if (j /= i) then
                    body_j = this%bodies(j)
                    a_i = a_i - body_j%mass * (body_j%r - body_i%r) / norm2(body_j%r - body_i%r) ** 3d0
                end if
            end do
            body_i%a = a_i
        end do
    end subroutine calc_accelerations_pointer
    
end module objects
