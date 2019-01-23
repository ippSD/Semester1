module dislin_mod
    use dislin
    implicit none
    
    contains
    
    !-----------------------------------------------------------------------!
    !   ( SUBROUTINE ) plot_end                                             !
    !-----------------------------------------------------------------------!
    !   Ends the currently running dislin plot.                             !
    !-----------------------------------------------------------------------!
    subroutine plot_end()
        call disfin()
    end subroutine
    
    !-----------------------------------------------------------------------!
    !   ( SUBROUTINE ) plot_legend                                          !
    !-----------------------------------------------------------------------!
    !   Plots the legend on the currently running dislin plot.              !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character(N)) legends                                !
    !               ! character vector containing the legend's entries.   !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character) legend_title                              !
    !   OPTIONAL    ! Legend title, "Legend" by default.                    !
    !---------------!-------------------------------------------------------!
    subroutine plot_legend(legends, legend_title)
        character(len=*), intent(in) :: legends(:)
        character(len=*), optional :: legend_title
        character(len=:), allocatable :: buffer
        
        integer :: n_lines, string_length, i
        character(len=6) :: default_legend = "Legend"
        
        n_lines = size(legends)
        string_length = len(legends(1))
        write(*,*) n_lines, string_length
        allocate(character(n_lines * string_length) :: buffer)
        
        call legini(buffer, n_lines, string_length)
        
        do i = 1, n_lines
            write(*,*) buffer
            call leglin(buffer, legends(i), i)
        end do
        
        if(present(legend_title)) then
            call legtit(legend_title)
        else
            call legtit(default_legend)
        end if
        
        call color('FORE')
        call legend(buffer, 7)
        
    end subroutine
    
    !-----------------------------------------------------------------------!
    !   ( SUBROUTINE ) plot_title                                           !
    !-----------------------------------------------------------------------!
    !   Plots the title on the currently running dislin plot.               !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character) tit                                       !
    !               ! Title of the plot.                                    !
    !---------------!-------------------------------------------------------!
    !   IN          ! (integer) plot_level                                  !
    !   OPTIONAL    ! Row on which the title is plotted  starting from      !
    !               ! the top: 1 -> top, 2 -> under the top, etc.           !
    !---------------!-------------------------------------------------------!
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
    
    !-----------------------------------------------------------------------!
    !   ( SUBROUTINE ) plot_subtitle                                        !
    !-----------------------------------------------------------------------!
    !   Plots the subtitle on the currently running dislin plot.            !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character) subtit                                    !
    !               ! Subtitle of the plot.                                 !
    !---------------!-------------------------------------------------------!
    subroutine plot_subtitle(subtit)
        character(len=*), intent(in) :: subtit
        
        call plot_title(subtit, 2)
    end subroutine plot_subtitle 
    
    !-----------------------------------------------------------------------!
    !   ( SUBROUTINE ) plot                                                 !
    !-----------------------------------------------------------------------!
    !   Plots (x,y) curves with Matlab style.                               !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(N)) x                                           !
    !               ! X-axis plot data.                                     !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(N)) y                                           !
    !               ! Y-axis plot data.                                     !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character) plotcolor                                 !
    !   OPTIONAL    ! Color of the curve. Available colors are:             !
    !               ! 'BLACK', 'RED', 'GREEN', 'BLUE', 'CYAN', 'YELLOW',    !
    !               ! 'ORANGE', 'MAGENTA', 'WHITE', 'FORE' (default),       !
    !               ! 'BACK' (background), 'GRAY' and                       !
    !               ! 'HALF' (half intensity of foreground)                 !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character) xlabel                                    !
    !   OPTIONAL    ! Title of the X-label.                                 !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character) ylabel                                    !
    !   OPTIONAL    ! Title of the Y-label.                                 !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character) ylabel                                    !
    !   OPTIONAL    ! Title of the Y-label.                                 !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character) file_name                                 !
    !   OPTIONAL    ! Filename in which the plot is saved.                  !
    !               ! File type goes according to filename extension.       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (logical) hold_on                                     !
    !   OPTIONAL    ! Does not close the plot on exit,                      !
    !               ! allowing multiple plots.                              !
    !---------------!-------------------------------------------------------!
    subroutine plot(x, y, plotcolor, xlabel, ylabel, file_name, hold_on)
        real :: x(:), y(:)
        character(len=*), optional :: plotcolor, xlabel, ylabel, file_name
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
                if ( extension(1:1) == '.' ) then
                    call metafl(trim(extension(2:3)))
                else
                    call metafl(extension)
                end if
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
            
            if(present(xlabel)) call name(xlabel, 'X')
            if(present(ylabel)) call name(ylabel, 'Y')
            
            !ic=intrgb(0.95,0.95,0.95)
            !call axsbgd(ic)
            
            call graf(xmin,xmax,xmin,xstep,ymin,ymax,ymin,ystep)
            call setrgb(0.7,0.7,0.7)
            call grid(1,1)
                        
            !call labdig(-1,'X')
            !call tics(10,'XY')
        
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