module dislin_mod
    use dislin
    implicit none
    
    contains
    
    subroutine plot_end()
        call disfin()
    end subroutine
    
    subroutine plot_legend(curves, legend_title)
        character(len=*), intent(in) :: curves(:)
        character(len=*), optional :: legend_title
        character(len=:), allocatable :: buffer
        
        integer :: n_lines, string_length, i
        character(len=6) :: default_legend = "Legend"
        
        n_lines = size(curves)
        string_length = len(curves(1))
        write(*,*) n_lines, string_length
        allocate(character(n_lines * string_length) :: buffer)
        
        call legini(buffer, n_lines, string_length)
        
        do i = 1, n_lines
            write(*,*) buffer
            call leglin(buffer, curves(i), i)
        end do
        
        if(present(legend_title)) then
            call legtit(legend_title)
        else
            call legtit(default_legend)
        end if
        
        call color('FORE')
        call legend(buffer, 7)
        
    end subroutine
    
    subroutine plot_title(tit, plot_level)
        character(len=*), intent(in) :: tit
        integer, optional :: plot_level
        
        integer :: level, lvl = 1
        call getlev(level)
        if(present(plot_level)) lvl = plot_level
        if(level > 0) call titlin(tit, lvl)
        
        call color('FORE')
        call title()
    end subroutine plot_title
    
    subroutine plot_subtitle(subtit)
        character(len=*), intent(in) :: subtit
        
        integer :: level
        call getlev(level)
        if(level > 0) call titlin(subtit, 2)
        
        call color('FORE')
        call title()
    end subroutine plot_subtitle
        
    
    subroutine plot(x, y, plotcolor, xlabel, ylabel, subtitle, plottitle, file_name, hold_on)
        real :: x(:), y(:)
        character(len=*), optional :: plotcolor, xlabel, ylabel, plottitle, subtitle, file_name
        logical, optional :: hold_on
        
        integer :: i, ic, level, filename_len
        character(len=3) :: extension
        character(len=:), allocatable :: trimmed_filename
        real :: xmin, xmax, ymin, ymax, xstep, ystep
        
        call getlev(level)
        xmin = minval(x)
        ymin = minval(y)
        xmax = maxval(x)
        ymax = maxval(y)
        xstep = (xmax - xmin) / 10e0
        ystep = (ymax - ymin) / 10e0
        
        if(level == 0) then
            if(present(file_name)) then
                filename_len = len(trim(file_name))
                allocate(character(len=filename_len) :: trimmed_filename)
                trimmed_filename = trim(file_name)
                extension = trimmed_filename(filename_len - 2:filename_len)
                do i=1,3
                    ic = ichar(extension(i:i))
                    if (ic >= 65 .and. ic < 90) extension(i:i) = char(ic+32) 
                end do
                call metafl(extension)
                call setfil(file_name)
            else
                call metafl('CONS')
            endif
            !call metafl('CONS')
            call scrmod('REVERS')
            call disini()
            call pagera()
            call complx()
            call texmod('ON')
            call axspos(450,1800)
            call axslen(2200,1200)
            
            ic=intrgb(0.95,0.95,0.95)
            call axsbgd(ic)
            
            call graf(xmin,xmax,xmin,xstep,ymin,ymax,ymin,ystep)
            call setrgb(0.7,0.7,0.7)
            call grid(1,1)
            
            if(present(xlabel)) call name(xlabel, 'X')
            if(present(ylabel)) call name(ylabel, 'Y')
            
            !call labdig(-1,'X')
            !call tics(10,'XY')
        
            if(present(plottitle)) call plot_title(plottitle)
            if(present(subtitle)) call plot_subtitle(subtitle)

        end if
        
        if(present(plotcolor)) call color(plotcolor)
        call curve(x,y,size(x))
        
        if(present(hold_on)) then
            if(.not. hold_on) call disfin()
        else
            call disfin()
        end if
    
    end subroutine plot
end module dislin_mod