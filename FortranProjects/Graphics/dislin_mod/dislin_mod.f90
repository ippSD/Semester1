module dislin_mod
    use dislin, color_raw => color
    implicit none
    
    interface plot_title
       module procedure plot_title_default, plot_title_position
    end interface plot_title
    
    public
    private :: plot_title_default, plot_title_position
    
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
        allocate(character(n_lines * string_length) :: buffer)
        
        call legini(buffer, n_lines, string_length)
        
        do i = 1, n_lines
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
    !   ( SUBROUTINE ) plot_title_default                                   !
    !-----------------------------------------------------------------------!
    !   Plots the title on the currently running dislin plot.               !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character) tit                                       !
    !               ! Title of the plot.                                    !
    !---------------!-------------------------------------------------------!
    subroutine plot_title_default(tit)
        character(len=*), intent(in) :: tit
        
        integer :: level, lvl = 1
        call getlev(level)
        if(level > 0) call titlin(tit, lvl)
        
        call color_raw('FORE')
        call title()
    end subroutine plot_title_default
    
    !-----------------------------------------------------------------------!
    !   ( SUBROUTINE ) plot_title_position                                  !
    !-----------------------------------------------------------------------!
    !   Plots the title on an specific line.                                !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character) tit                                       !
    !               ! Title of the plot.                                    !
    !---------------!-------------------------------------------------------!
    !   IN          ! (integer) plot_level                                  !
    !               ! Row on which the title is plotted  starting from      !
    !               ! the top: 1 -> top, 2 -> under the top, etc.           !
    !---------------!-------------------------------------------------------!
    subroutine plot_title_position(tit, plot_level)
        character(len=*), intent(in) :: tit
        integer :: plot_level
        
        integer :: level
        call getlev(level)
        if(level > 0) call titlin(tit, plot_level)
        
        call color_raw('FORE')
        call title()
    end subroutine plot_title_position
    
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
    !   ( SUBROUTINE ) scatter                                              !
    !-----------------------------------------------------------------------!
    !   Plots (x,y) scatter points with Matlab style.                       !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(N)) x                                           !
    !               ! X-axis plot data.                                     !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(N)) y                                           !
    !               ! Y-axis plot data.                                     !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character) color                                     !
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
    !   IN          ! (character) filename                                  !
    !   OPTIONAL    ! Filename in which the plot is saved.                  !
    !               ! File type goes according to filename extension.       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (logical) hold_on                                     !
    !   OPTIONAL    ! Does not close the plot on exit,                      !
    !               ! allowing multiple plots.                              !
    !---------------!-------------------------------------------------------!
    subroutine scatter(x, y, color, xlabel, ylabel, filename, hold_on)
        real, intent(in) :: x(:), y(:)
        character(len=*), optional :: color, xlabel, ylabel, filename
        logical, optional :: hold_on
        
        logical :: do_color  = .false., &
                   do_xlabel = .false., &
                   do_ylabel = .false., &
                   do_save   = .false., &
                   do_hold  = .false.
        
        character( len = 127 ) :: opt_color    = "FORE", &
                                  opt_xlabel   = ""    , &
                                  opt_ylabel   = ""    , &
                                  opt_filename = ""
        
        if ( present( color ) ) then
            do_color = .true.
            opt_color = color
        end if
        
        if ( present( xlabel ) ) then
            do_xlabel = .true.
            opt_xlabel = xlabel
        end if
        
        if ( present( ylabel ) ) then
            do_ylabel = .true.
            opt_ylabel = ylabel
        end if
        
        if ( present( filename ) ) then
            do_save = .true.
            opt_filename = filename
        end if
        
        if ( present( hold_on ) ) do_hold = hold_on
        
        call plot_hub( &
            x = x, &
            y = y, &
            color = opt_color, &
            plot_type = 1, &
            hold_on = do_hold, &
            do_xlabel = do_xlabel, &
            do_ylabel = do_ylabel, &
            do_save = do_save, &
            xlabel = opt_xlabel, &
            ylabel = opt_ylabel, &
            filename = opt_filename, &
            logaxis = .false. &
        )
    
    end subroutine scatter
    
    !-----------------------------------------------------------------------!
    !   ( SUBROUTINE ) semilogy                                             !
    !-----------------------------------------------------------------------!
    !   Plots (x,y) curves on Log10 y axis with Matlab style.               !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(N)) x                                           !
    !               ! X-axis plot data.                                     !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(N)) y                                           !
    !               ! Y-axis plot data.                                     !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character) color                                     !
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
    !   IN          ! (character) filename                                  !
    !   OPTIONAL    ! Filename in which the plot is saved.                  !
    !               ! File type goes according to filename extension.       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (logical) hold_on                                     !
    !   OPTIONAL    ! Does not close the plot on exit,                      !
    !               ! allowing multiple plots.                              !
    !---------------!-------------------------------------------------------!
    subroutine semilogy(x, y, color, xlabel, ylabel, filename, hold_on)
        real, intent(in) :: x(:), y(:)
        character(len=*), optional :: color, xlabel, ylabel, filename
        logical, optional :: hold_on
        
        logical :: do_color  = .false., &
                   do_xlabel = .false., &
                   do_ylabel = .false., &
                   do_save   = .false., &
                   do_hold  = .false.
        
        character( len = 127 ) :: opt_color    = "FORE", &
                                  opt_xlabel   = ""    , &
                                  opt_ylabel   = ""    , &
                                  opt_filename = ""
        
        if ( present( color ) ) then
            do_color = .true.
            opt_color = color
        end if
        
        if ( present( xlabel ) ) then
            do_xlabel = .true.
            opt_xlabel = xlabel
        end if
        
        if ( present( ylabel ) ) then
            do_ylabel = .true.
            opt_ylabel = ylabel
        end if
        
        if ( present( filename ) ) then
            do_save = .true.
            opt_filename = filename
        end if
        
        if ( present( hold_on ) ) do_hold = hold_on
        
        call plot_hub( &
            x = x, &
            y = y, &
            color = opt_color, &
            plot_type = 0, &
            hold_on = do_hold, &
            do_xlabel = do_xlabel, &
            do_ylabel = do_ylabel, &
            do_save = do_save, &
            xlabel = opt_xlabel, &
            ylabel = opt_ylabel, &
            filename = opt_filename, &
            logaxis = .true. &
        )
    end subroutine semilogy
    
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
    !   IN          ! (character) color                                     !
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
    !   IN          ! (character) filename                                  !
    !   OPTIONAL    ! Filename in which the plot is saved.                  !
    !               ! File type goes according to filename extension.       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (logical) hold_on                                     !
    !   OPTIONAL    ! Does not close the plot on exit,                      !
    !               ! allowing multiple plots.                              !
    !---------------!-------------------------------------------------------!
    subroutine plot(x, y, color, xlabel, ylabel, filename, hold_on)
        real, intent(in) :: x(:), y(:)
        character(len=*), optional :: color, xlabel, ylabel, filename
        logical, optional :: hold_on
        
        logical :: do_color  = .false., &
                   do_xlabel = .false., &
                   do_ylabel = .false., &
                   do_save   = .false., &
                   do_hold  = .false.
        
        character( len = 127 ) :: opt_color    = "FORE", &
                                  opt_xlabel   = ""    , &
                                  opt_ylabel   = ""    , &
                                  opt_filename = ""
        
        if ( present( color ) ) then
            do_color = .true.
            opt_color = color
        end if
        
        if ( present( xlabel ) ) then
            do_xlabel = .true.
            opt_xlabel = xlabel
        end if
        
        if ( present( ylabel ) ) then
            do_ylabel = .true.
            opt_ylabel = ylabel
        end if
        
        if ( present( filename ) ) then
            do_save = .true.
            opt_filename = filename
        end if
        
        if ( present( hold_on ) ) do_hold = hold_on
        
        call plot_hub( &
            x = x, &
            y = y, &
            color = opt_color, &
            plot_type = 0, &
            hold_on = do_hold, &
            do_xlabel = do_xlabel, &
            do_ylabel = do_ylabel, &
            do_save = do_save, &
            xlabel = opt_xlabel, &
            ylabel = opt_ylabel, &
            filename = opt_filename, &
            logaxis = .false. &
        )
    
    end subroutine plot
    
    !-----------------------------------------------------------------------!
    !   ( SUBROUTINE ) plot_hub                                             !
    !-----------------------------------------------------------------------!
    !   Plots (x,y) curves, scatters, ... with Matlab style.                !
    !-----------------------------------------------------------------------!
    !   Parameters: !                                                       !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(N)) x                                           !
    !               ! X-axis plot data.                                     !
    !---------------!-------------------------------------------------------!
    !   IN          ! (real(N)) y                                           !
    !               ! Y-axis plot data.                                     !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character) color                                     !
    !               ! Colour of the curve. Available colors are:            !
    !               ! 'BLACK', 'RED', 'GREEN', 'BLUE', 'CYAN', 'YELLOW',    !
    !               ! 'ORANGE', 'MAGENTA', 'WHITE', 'FORE' (default),       !
    !               ! 'BACK' (background), 'GRAY' and                       !
    !               ! 'HALF' (half intensity of foreground)                 !
    !---------------!-------------------------------------------------------!
    !   IN          ! (integer) plot_type                                   !
    !               ! Integer of what to plot:                              !
    !               ! 0 => plot, 1 => scatter                               !
    !---------------!-------------------------------------------------------!
    !   IN          ! (logical) hold_on                                     !
    !               ! If true, the current plot is not closed on exit,      !
    !               ! allowing multiple plots.                              !
    !---------------!-------------------------------------------------------!
    !   IN          ! (logical) do_xlabel                                   !
    !               ! If true, X-Axis label is set to xlabel.               !
    !---------------!-------------------------------------------------------!
    !   IN          ! (logical) do_ylabel                                   !
    !               ! If true, Y-Axis label is set to ylabel.               !
    !---------------!-------------------------------------------------------!
    !   IN          ! (logical) do_save                                     !
    !               ! If true, plot outputs to filename.                    !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character) xlabel                                    !
    !               ! Title of the X-label if do_xlabel is true.            !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character) ylabel                                    !
    !               ! Title of the Y-label if do_ylabel is true.            !
    !---------------!-------------------------------------------------------!
    !   IN          ! (character) filename                                  !
    !               ! Filename in which the plot is saved if do_save is     !
    !               ! true. File type goes according to filename extension. !
    !---------------!-------------------------------------------------------!
    !   IN          ! (logical) logaxis                                     !
    !               ! If true, use Log Y-Axis.                              !
    !---------------!-------------------------------------------------------!
    subroutine plot_hub(x, y, color, plot_type, hold_on, do_xlabel, do_ylabel, do_save, xlabel, ylabel, filename, logaxis)
        real, intent(in) :: x(:), y(:)
        integer, intent(in) :: plot_type
        logical, intent(in) :: hold_on, do_xlabel, do_ylabel, do_save
        character(len=*), intent(in) :: color, xlabel, ylabel, filename
        logical, intent(in) :: logaxis
        
        integer :: i, ic, level, filename_len
        character(len=3) :: extension
        character(len=:), allocatable :: trimmed_filename
        real :: xmin, xmax, ymin, ymax, xstep, ystep, yminlog10, ymaxlog10, ysteplog10
        
        call getlev(level)
        xmin = minval(x)
        ymin = minval(y)
        xmax = maxval(x)
        ymax = maxval(y)
        xstep = (xmax - xmin) / 10e0
        ystep = (ymax - ymin) / 10e0
        !yminlog10 = 10e0 ** (int(log10(ymin)))
        !ymaxlog10 = 10e0 ** (int(log10(ymax)))
        !ysteplog10  = 10e0 ** (int(log10(ymax) - log10(ymin)))
        yminlog10 = int(log10(ymin)) - 1
        ymaxlog10 = int(log10(ymax)) + 1
        ysteplog10  = (ymaxlog10 - yminlog10) / 10e0
        
        if(level == 0) then
            if( do_save ) then
                filename_len = len(trim(filename))
                allocate(character(len=filename_len) :: trimmed_filename)
                trimmed_filename = trim(filename)
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
                call setfil(trim(filename))
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
            
            !ic=intrgb(0.95,0.95,0.95)
            !call axsbgd(ic)
            !call ticks(5,'XY')
            if( do_xlabel ) call name(xlabel, 'X')
            if( do_ylabel ) call name(ylabel, 'Y')
            
            if ( xmax - xmin < 1d-2 ) call labels('EXP', 'X')
            if ( xmax - xmin >= 1d5 ) call labels('EXP', 'X')
            if ( ymax - ymin < 1d-2 ) call labels('EXP', 'Y')
            if ( ymax - ymin >= 1d5 ) call labels('EXP', 'Y')
            if ( logaxis )  then
                call axsscl('LOG', 'Y')
                call labels('EXP', 'Y')
                call graf(xmin,xmax,xmin,xstep,yminlog10,ymaxlog10,yminlog10,ysteplog10)
            else
                call graf(xmin,xmax,xmin,xstep,ymin,ymax,ymin,ystep)
            end if
            
            !if ( plot_type == 1 ) call graf3(xmin,xmax,xmin,xstep,ymin,ymax,ymin,ystep, 0e0, 1e-3, 0e0, 1e-3)
            call setrgb(0.7,0.7,0.7)
            call grid(1,1)
                        
            !call labdig(-1,'X')
            
        
        end if
        call color_raw(color)
        if ( plot_type == 0 ) call curve( x, y, size(x) )
        if ( plot_type == 1 ) then
            plot_each_point: do i = 1, size(x)
                call rlsymb( 3, x(i), y(i) )
            end do plot_each_point
        end if
        
        if( .not. hold_on ) call disfin()
    
    end subroutine plot_hub
    
end module dislin_mod