!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file pyplot_module.f90
!>
!> @brief Interface for plots in python (x-y plots)
!>
!> @details All subroutines and functions of this file are related to the
!> frequency domain analysis.
!>
!> @author Guenter Glanzer, SAP SE
!>
!> @date Feb 2019

!===============================================================================
!> @brief Interface for plots in python (x-y plots).
module pyplot_module

  use kindModule, only : dp

  implicit none

  private

  integer, parameter :: max_int_len  = 10 !< max string length for integers
  integer, parameter :: max_real_len = 30 !< max string length for reals

  !> @brief Main python plot class.
  type, public :: pyplot

     private

     character(len=:), allocatable :: python   !< the python executable
     character(len=:), allocatable :: str      !< String buffer
     character(len=:), allocatable :: real_fmt !< Real number formatting

     logical :: show_legend  = .false. !< Show legend into plot
     logical :: use_numpy    = .true.  !< Use numpy python module
     logical :: use_oo_api   = .false. !< Use OO interface of matplotlib (incompatible with showfig subroutine)
     logical :: mplot3d      = .false. !< It is a 3d plot
     logical :: polar        = .false. !< It is a polar plot
     logical :: axis_equal   = .false. !< Equal scale on each axis
     logical :: tight_layout = .false. !< Tight layout option
     logical :: initialized  = .false. !< Set to .true. in first initialization

   contains

     procedure, public :: initialize  !< Initialize pyplot instance
     procedure, public :: add_plot    !< Add a 2d plot to pyplot instance
     procedure, public :: add_3d_plot !< Add a 3d plot to pyplot instance
     procedure, public :: add_sphere  !< Add a 3d sphere to pyplot instance
     procedure, public :: add_contour !< Add a contour plot to pyplot instance
     procedure, public :: add_bar     !< Add a barplot to pyplot instance
     procedure, public :: add_imshow  !< Add an image plot (using `imshow`)
     procedure, public :: add_hist    !< Add a histogram plot to pyplot instance
     procedure, public :: savefig     !< Save plots of pyplot instance
     procedure, public :: showfig     !< Show plots of pyplot instance
     procedure, public :: destroy     !< Destroy pyplot instance

     procedure         :: execute     !< Execute pyplot commands
     procedure         :: add_str     !< Add string to pytplot instance buffer
     procedure         :: finish_ops  !< Some final ops before saving

  end type pyplot


contains

  !=============================================================================
  !> @brief Destroys an instance of pyplot
  !> @param[in] me    class handler

  subroutine destroy(me)

    class(pyplot), intent(inout)       :: me

    if (allocated(me%str))       deallocate(me%str)
    if (allocated(me%real_fmt))  deallocate(me%real_fmt)

  end subroutine destroy


  !=============================================================================
  !> @brief Add a string to the string buffer
  !> @param[in] me    class handler
  !> @param[in] str   string buffer

  subroutine add_str(me, str)

    class(pyplot),    intent(inout) :: me
    character(len=*), intent(in)    :: str

    me%str = me%str//str//new_line(' ')

  end subroutine add_str


  !=============================================================================
  !> @brief Initializes curve plots
  !> @param[inout] me    class handler
  !> @param[in] grid               flag for grid drawing
  !> @param[in] xlabel             x axis label
  !> @param[in] ylabel             y axis label
  !> @param[in] zlabel             z axis label
  !> @param[in] title              drawing title
  !> @param[in] legend             drawing legend
  !> @param[in] use_numpy          flag for using numpy
  !> @param[in] figsize            figure dimesnion
  !> @param[in] font_size          font size
  !> @param[in] axes_labelsize     axis label size
  !> @param[in] xtick_labelsize    size of x axis tick lables
  !> @param[in] ytick_labelsize    size of y axis tick lables
  !> @param[in] ztick_labelsize    size of z axis tick lables
  !> @param[in] legend_fontsize    size of legend font
  !> @param[in] mplot3d            set true for 3d plots
  !> @param[in] axis_equal         set true for axis = 'equal'
  !> @param[in] polar              set true for polar plots
  !> @param[in] real_fmt           format string for real numbers
  !> @param[in] use_oo_api         avoid matplotlib's GUI by using the OO interface
  !> @param[in] axisbelow          to put the grid lines below the other chart elements
  !> @param[in] tight_layout       enable tight layout
  !> @param[out] istat             Non-zero if python cannot be invoked

  subroutine initialize (me, grid, xlabel, ylabel, zlabel, title, legend, &
       &                 use_numpy, figsize, font_size, axes_labelsize, &
       &                 xtick_labelsize, ytick_labelsize, ztick_labelsize, &
       &                 legend_fontsize, mplot3d, axis_equal, polar, real_fmt, &
       &                 use_oo_api, axisbelow, tight_layout, istat)

    use reportErrorModule, only : error_p, reportError

    class(pyplot),         intent(inout)        :: me
    logical,               intent(in), optional :: grid
    character(len=*),      intent(in), optional :: xlabel
    character(len=*),      intent(in), optional :: ylabel
    character(len=*),      intent(in), optional :: zlabel
    character(len=*),      intent(in), optional :: title
    logical,               intent(in), optional :: legend
    logical,               intent(in), optional :: use_numpy
    integer, dimension(2), intent(in), optional :: figsize
    integer,               intent(in), optional :: font_size
    integer,               intent(in), optional :: axes_labelsize
    integer,               intent(in), optional :: xtick_labelsize
    integer,               intent(in), optional :: ytick_labelsize
    integer,               intent(in), optional :: ztick_labelsize
    integer,               intent(in), optional :: legend_fontsize
    logical,               intent(in), optional :: mplot3d
    logical,               intent(in), optional :: axis_equal
    logical,               intent(in), optional :: polar
    character(len=*),      intent(in), optional :: real_fmt
    logical,               intent(in), optional :: use_oo_api
    logical,               intent(in), optional :: axisbelow
    logical,               intent(in), optional :: tight_layout
    integer,               intent(out)          :: istat

    character(len=max_int_len)                  :: width_str                    !< figure width dummy string
    character(len=max_int_len)                  :: height_str                   !< figure height dummy string
    character(len=max_int_len)                  :: font_size_str                !< font size dummy string
    character(len=max_int_len)                  :: axes_labelsize_str           !< size of axis labels dummy string
    character(len=max_int_len)                  :: xtick_labelsize_str          !< size of x axis tick labels dummy string
    character(len=max_int_len)                  :: ytick_labelsize_str          !< size of x axis tick labels dummy string
    character(len=max_int_len)                  :: ztick_labelsize_str          !< size of z axis tick labels dummy string
    character(len=max_int_len)                  :: legend_fontsize_str          !< size of legend font dummy string

    character(len=:),      allocatable          :: python_fig_func              !< Python's function for creating a new Figure instance

#if defined(win32) || defined(win64)
    character(len=*), parameter :: devnull = 'nul'       !< Null device on Windows
#else
    character(len=*), parameter :: devnull = '/dev/null' !< Null device on Linux
#endif

    call me%destroy()

    if (.not. allocated(me%python)) then
       if (me%initialized) then
          istat = -1
          return
       else
          me%initialized = .true.
       end if
       !! Find valid python executable
       call execute_command_line ('python3 --version 2>'//devnull, &
            &                     exitstat=istat)
       if (istat == 0) then
          me%python = 'python3 '
       else
          call execute_command_line ('python --version 2>'//devnull, &
               &                     exitstat=istat)
          if (istat == 0) then
             me%python = 'python '
          else
             call reportError (error_p,'python is not available, no plotting', &
                  &            addString='pyplot_module::initialize')
             return
          end if
       end if
    end if

    if (present(legend)) then
        me%show_legend = legend
    else
        me%show_legend = .false.
    end if
    if (present(use_numpy)) then
        me%use_numpy = use_numpy
    else
        me%use_numpy = .true.
    end if
    if (present(use_oo_api)) then
        me%use_oo_api = use_oo_api
    else
        me%use_oo_api = .false.
    end if
    if (present(figsize)) then
        call integer_to_string(figsize(1), width_str)
        call integer_to_string(figsize(2), height_str)
    end if
    if (present(mplot3d)) then
        me%mplot3d = mplot3d
    else
        me%mplot3d = .false.
    end if
    if (present(polar)) then
        me%polar = polar
    else
        me%polar = .false.
    end if
    if (present(axis_equal)) then
        me%axis_equal = axis_equal
    else
        me%axis_equal = .false.
    end if
    if (present(real_fmt)) then
        me%real_fmt = trim(adjustl(real_fmt))
    else
        me%real_fmt = '(E30.16)' ! default real number format
    end if
    if (present(tight_layout)) then
        me%tight_layout = tight_layout
    else
        me%tight_layout = .false.
    end if

    call optional_int_to_string(font_size, font_size_str, '15')
    call optional_int_to_string(axes_labelsize, axes_labelsize_str, '15')
    call optional_int_to_string(xtick_labelsize, xtick_labelsize_str, '15')
    call optional_int_to_string(ytick_labelsize, ytick_labelsize_str, '15')
    call optional_int_to_string(ztick_labelsize, ztick_labelsize_str, '15')
    call optional_int_to_string(legend_fontsize, legend_fontsize_str, '15')

    me%str = ''

    call me%add_str('#!/usr/bin/env python3')
    call me%add_str('')

    call me%add_str('import matplotlib')
    if (me%use_oo_api) then
        call me%add_str('from matplotlib.figure import Figure')
        call me%add_str('from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas')
    else
        call me%add_str('import matplotlib.pyplot as plt')
    end if
    if (me%mplot3d) call me%add_str('from mpl_toolkits.mplot3d import Axes3D')
    if (me%use_numpy) call me%add_str('import numpy as np')
    call me%add_str('')

    call me%add_str('matplotlib.rcParams["font.family"] = "Serif"')
    call me%add_str('matplotlib.rcParams["font.size"] = '//trim(font_size_str))
    call me%add_str('matplotlib.rcParams["axes.labelsize"] = '//trim(axes_labelsize_str))
    call me%add_str('matplotlib.rcParams["xtick.labelsize"] = '//trim(xtick_labelsize_str))
    call me%add_str('matplotlib.rcParams["ytick.labelsize"] = '//trim(ytick_labelsize_str))
    call me%add_str('matplotlib.rcParams["legend.fontsize"] = '//trim(legend_fontsize_str))

    call me%add_str('')

    if (me%use_oo_api) then
        allocate(python_fig_func, source='Figure')
    else
        allocate(python_fig_func, source='plt.figure')
    end if
    if (present(figsize)) then  !if specifying the figure size
        call me%add_str('fig = '//python_fig_func//'(figsize=('//trim(width_str)//','//trim(height_str)//'),facecolor="white")')
    else
        call me%add_str('fig = '//python_fig_func//'(facecolor="white")')
    end if

    if (me%mplot3d) then
        call me%add_str('ax = fig.gca(projection=''3d'')')
    elseif (me%polar) then
        call me%add_str('ax = fig.gca(projection=''polar'')')
    else
        call me%add_str('ax = fig.gca()')
    end if

    if (present(grid)) then
        if (grid) call me%add_str('ax.grid()')
    end if

    if (.not. present(axisbelow)) then
       call me%add_str('ax.set_axisbelow(True)') ! default
    else if (axisbelow) then
       call me%add_str('ax.set_axisbelow(True)')
    end if

    if (present(xlabel)) call me%add_str('ax.set_xlabel("'//trim(xlabel)//'")')
    if (present(ylabel)) call me%add_str('ax.set_ylabel("'//trim(ylabel)//'")')
    if (present(zlabel)) call me%add_str('ax.set_zlabel("'//trim(zlabel)//'")')
    if (present(title))  call me%add_str('ax.set_title("' //trim(title) //'")')

    call me%add_str('')

  end subroutine initialize


  !=============================================================================
  !> @cond NO_DOCUMENTAION
  function initError(prg) result(istat)
    use reportErrorModule, only : error_p, reportError
    character(len=*), intent(in) :: prg
    integer :: istat
    call reportError (error_p,'pyplot class not properly initialized', &
         &            addString='pyplot_module::'//prg)
    istat = -1
  end function initError
  !> @endcond


  !=============================================================================
  !> @brief Add an x,y plot
  !> @param[inout] me      class handler
  !> @param[in] x          x values
  !> @param[in] y          y values
  !> @param[in] label      plot label
  !> @param[in] linestyle  style of the plot line
  !> @param[in] markersize size of the plot markers
  !> @param[in] linewidth  width of the plot line
  !> @param[in] xlim       x-axis range
  !> @param[in] ylim       y-axis range
  !> @param[in] xscale     example: 'linear' (default), 'log'
  !> @param[in] yscale     example: 'linear' (default), 'log'
  !> @param[in] color      RGB color tuple
  !> @param[out] istat      status output

  subroutine add_plot(me, x, y, label, linestyle, markersize, linewidth, xlim, ylim, xscale, yscale, color, istat)

    class(pyplot),          intent (inout)        :: me
    real(dp), dimension(:), intent (in)           :: x
    real(dp), dimension(:), intent (in)           :: y
    character(len=*),       intent (in)           :: label
    character(len=*),       intent (in)           :: linestyle
    integer,                intent (in), optional :: markersize
    integer,                intent (in), optional :: linewidth
    real(dp),dimension(2),  intent (in), optional :: xlim
    real(dp),dimension(2),  intent (in), optional :: ylim
    character(len=*),       intent (in), optional :: xscale
    character(len=*),       intent (in), optional :: yscale
    real(dp),dimension(:),  intent (in), optional :: color
    integer,                intent (out)          :: istat

    character(len=:),                 allocatable :: arg_str      !< the arguments to pass to `plot`
    character(len=:),                 allocatable :: xstr         !< x values stringified
    character(len=:),                 allocatable :: ystr         !< y values stringified
    character(len=:),                 allocatable :: xlimstr      !< xlim values stringified
    character(len=:),                 allocatable :: ylimstr      !< ylim values stringified
    character(len=:),                 allocatable :: color_str    !< color values stringified
    character(len=max_int_len)                    :: imark        !< actual markers size
    character(len=max_int_len)                    :: iline        !< actual line width
    character(len=*), parameter                   :: xname = 'x'  !< x variable name for script
    character(len=*), parameter                   :: yname = 'y'  !< y variable name for script


    if (allocated(me%str)) then

        istat = 0

        !axis limits (optional):
        if (present(xlim)) call vec_to_string(xlim, me%real_fmt, xlimstr, me%use_numpy)
        if (present(ylim)) call vec_to_string(ylim, me%real_fmt, ylimstr, me%use_numpy)

        !convert the arrays to strings:
        call vec_to_string(x, me%real_fmt, xstr, me%use_numpy)
        call vec_to_string(y, me%real_fmt, ystr, me%use_numpy)

        !get optional inputs (if not present, set default value):
        call optional_int_to_string(markersize, imark, '3')
        call optional_int_to_string(linewidth, iline, '3')

        !write the arrays:
        call me%add_str(trim(xname)//' = '//xstr)
        call me%add_str(trim(yname)//' = '//ystr)
        call me%add_str('')

        !main arguments for plot:
        allocate(arg_str, source=trim(xname)//','//&
                  trim(yname)//','//&
                  '"'//trim(linestyle)//'",'//&
                  'linewidth='//trim(adjustl(iline))//','//&
                  'markersize='//trim(adjustl(imark))//','//&
                  'label="'//trim(label)//'"')

        ! optional arguments:
        if (present(color)) then
            if (size(color)<=3) then
                call vec_to_string(color(1:3), '*', color_str, use_numpy=.false., is_tuple=.true.)
                arg_str = arg_str//',color='//trim(color_str)
            end if
        end if

        !write the plot statement:
        call me%add_str('ax.plot('//arg_str//')')

        !axis limits:
        if (allocated(xlimstr)) call me%add_str('ax.set_xlim('//xlimstr//')')
        if (allocated(ylimstr)) call me%add_str('ax.set_ylim('//ylimstr//')')

        !axis scales:
        if (present(xscale)) call me%add_str('ax.set_xscale("'//xscale//'")')
        if (present(yscale)) call me%add_str('ax.set_yscale("'//yscale//'")')

        call me%add_str('')

    else
       istat = initError('add_plot')
    end if

  end subroutine add_plot


  !=============================================================================
  !> @brief Add an histogram plot
  !> @param[inout] me      class handler
  !> @param[in] x          x values
  !> @param[in] label      plot label
  !> @param[in] xlim       x-axis range
  !> @param[in] ylim       y-axis range
  !> @param[in] xscale     example: 'linear' (default), 'log'
  !> @param[in] yscale     example: 'linear' (default), 'log'
  !> @param[in] bins       number of bins
  !> @param[in] normed     boolean flag that determines whether bin counts are normalized
  !> @param[in] cumulative boolean flag that determines whether histogram represents the cumulative density of dataset
  !> @param[out] istat     status output

  subroutine add_hist(me, x, label, xlim, ylim, xscale, yscale, bins, normed, cumulative, istat)

    class(pyplot),          intent (inout)        :: me
    real(dp), dimension(:), intent (in)           :: x
    character(len=*),       intent (in)           :: label
    real(dp),dimension(2),  intent (in), optional :: xlim
    real(dp),dimension(2),  intent (in), optional :: ylim
    character(len=*),       intent (in), optional :: xscale
    character(len=*),       intent (in), optional :: yscale
    integer,                intent (in), optional :: bins
    logical,                intent (in), optional :: normed
    logical,                intent (in), optional :: cumulative
    integer,                intent (out)          :: istat

    character(len=*),                 parameter   :: xname = 'x'      !< x variable name for script
    character(len=:),                 allocatable :: xstr             !< x values stringified
    character(len=:),                 allocatable :: xlimstr          !< xlim values stringified
    character(len=:),                 allocatable :: ylimstr          !< ylim values stringified
    character(len=:),                 allocatable :: normedstr        !< optional stuff
    character(len=:),                 allocatable :: cumulativestr    !< string for cumulative
    character(len=max_int_len)                    :: binsstr          !< string for bins

    if (allocated(me%str)) then

        istat = 0

        !axis limits (optional):
        if (present(xlim)) call vec_to_string(xlim, me%real_fmt, xlimstr, me%use_numpy)
        if (present(ylim)) call vec_to_string(ylim, me%real_fmt, ylimstr, me%use_numpy)

        !convert the arrays to strings:
        call vec_to_string(x, me%real_fmt, xstr, me%use_numpy)

        !write the arrays:
        call me%add_str(trim(xname)//' = '//xstr)
        call me%add_str('')

        !get optional inputs (if not present, set default value):
        call optional_int_to_string(bins, binsstr, '10')
        call optional_logical_to_string(normed, normedstr, 'False')
        call optional_logical_to_string(cumulative, cumulativestr, 'False')

        !write the plot statement:
        call me%add_str('ax.hist('//&
                        trim(xname)//','//&
                        'label="'//trim(label)//'",'//&
                        'bins='//trim(binsstr)//','//&
                        'cumulative='//trim(cumulativestr)//','//&
                        'normed='//trim(normedstr)//')')

        !axis limits:
        if (allocated(xlimstr)) call me%add_str('ax.set_xlim('//xlimstr//')')
        if (allocated(ylimstr)) call me%add_str('ax.set_ylim('//ylimstr//')')

        !axis scales:
        if (present(xscale)) call me%add_str('ax.set_xscale("'//xscale//'")')
        if (present(yscale)) call me%add_str('ax.set_yscale("'//yscale//'")')

        call me%add_str('')

    else
       istat = initError('add_hist')
    end if

  end subroutine add_hist


  !=============================================================================
  !> @brief Add a contour plot
  !> @note This requires `use_numpy` to be True.
  !> @param[inout] me      class handler
  !> @param[in] x          x values
  !> @param[in] y          y values
  !> @param[in] z          z values
  !> @param[in] label      plot label
  !> @param[in] linestyle  style of the plot line
  !> @param[in] linewidth  width of the plot line
  !> @param[in] levels     contour levels to plot
  !> @param[in] color      color of the contour line
  !> @param[in] filled     use filled control (default=False)
  !> @param[in] cmap       colormap if filled=True (examples: 'jet', 'bone')
  !> @param[in] colorbar   add a colorbar (default=False)
  !> @param[out] istat     status output

  subroutine add_contour(me, x, y, z, label, linestyle, linewidth, levels, color, &
                           filled, cmap, colorbar, istat)

    class(pyplot),           intent (inout)        :: me
    real(dp),dimension(:),   intent (in)           :: x
    real(dp),dimension(:),   intent (in)           :: y
    real(dp),dimension(:,:), intent (in)           :: z
    character(len=*),        intent (in)           :: label
    character(len=*),        intent (in)           :: linestyle
    integer,                 intent (in), optional :: linewidth
    real(dp),dimension(:),   intent (in), optional :: levels
    character(len=*),        intent (in), optional :: color
    logical,                 intent (in), optional :: filled
    character(len=*),        intent (in), optional :: cmap
    logical,                 intent (in), optional :: colorbar
    integer,                 intent (out)          :: istat

    character(len=:),                  allocatable :: xstr          !< x values stringified
    character(len=:),                  allocatable :: ystr          !< y values stringified
    character(len=:),                  allocatable :: zstr          !< z values stringified
    character(len=:),                  allocatable :: levelstr      !< levels vector stringified
    character(len=max_int_len)                     :: iline         !< actual line width
    character(len=*), parameter                    :: xname = 'x'   !< x variable name for script
    character(len=*), parameter                    :: yname = 'y'   !< y variable name for script
    character(len=*), parameter                    :: zname = 'z'   !< z variable name for script
    character(len=*), parameter                    :: xname_ = 'X'  !< X variable name for contour
    character(len=*), parameter                    :: yname_ = 'Y'  !< Y variable name for contour
    character(len=*), parameter                    :: zname_ = 'Z'  !< Z variable name for contour
    character(len=:), allocatable                  :: extras        !< optional stuff
    character(len=:), allocatable                  :: contourfunc   !< 'contour' or 'contourf'
    logical                                        :: is_filled     !< if it is a filled contour plot


    if (allocated(me%str)) then

        istat = 0

        !convert the arrays to strings:
        call vec_to_string(x, me%real_fmt, xstr, me%use_numpy)
        call vec_to_string(y, me%real_fmt, ystr, me%use_numpy)
        call matrix_to_string(z, me%real_fmt, zstr, me%use_numpy)
        if (present(levels)) call vec_to_string(levels, me%real_fmt, levelstr, me%use_numpy)

        !get optional inputs (if not present, set default value):
        call optional_int_to_string(linewidth, iline, '3')

        !write the arrays:
        call me%add_str(trim(xname)//' = '//xstr)
        call me%add_str(trim(yname)//' = '//ystr)
        call me%add_str(trim(zname)//' = '//zstr)
        call me%add_str('')

        !convert inputs for contour plotting:
        call me%add_str(xname_//', '//yname_//' = np.meshgrid('//trim(xname)//', '//trim(yname)//')')
        call me%add_str(zname_//' = np.transpose('//zname//')')

        !optional arguments:
        allocate(extras, source='')
        if (present(levels))     extras = extras//','//'levels='//levelstr
        if (present(color))      extras = extras//','//'colors="'//color//'"'
        if (present(linewidth))  extras = extras//','//'linewidths='//trim(adjustl(iline))
        if (present(cmap))       extras = extras//','//'cmap="'//cmap//'"'

        !filled or regular:
        allocate(contourfunc, source='contour')  !default
        is_filled = .false.
        if (present(filled)) then
            is_filled = filled
            if (filled) contourfunc = 'contourf'  !filled contour
        end if

        !write the plot statement:
        call me%add_str('CS = ax.'//contourfunc//'('//xname_//','//yname_//','//zname_//','//&
                                        'label="'//trim(label)//'",'//&
                                        'linestyles="'//trim(adjustl(linestyle))//'"'//&
                                        extras//')')

        if (present(colorbar)) then
            if (colorbar) call me%add_str('fig.colorbar(CS)')
        end if

        if (.not. is_filled) call me%add_str('ax.clabel(CS, fontsize=9, inline=1)')
        call me%add_str('')

    else
       istat = initError('add_contour')
    end if

  end subroutine add_contour


  !=============================================================================
  !> @brief Add a 3D x, y, z plot
  !> @note Must initialize the class with ```mplot3d=.true.```
  !> @param[inout] me        class handler
  !> @param[in] x            x values
  !> @param[in] y            y values
  !> @param[in] z            z values
  !> @param[in] label        plot label
  !> @param[in] linestyle    style of the plot line
  !> @param[in] markersize   size of the plot markers
  !> @param[in] linewidth    width of the plot line
  !> @param[out] istat       status output

  subroutine add_3d_plot(me, x, y, z, label, linestyle, markersize, linewidth, istat)

    class(pyplot),          intent (inout)        :: me
    real(dp), dimension(:), intent (in)           :: x
    real(dp), dimension(:), intent (in)           :: y
    real(dp), dimension(:), intent (in)           :: z
    character(len=*),       intent (in)           :: label
    character(len=*),       intent (in)           :: linestyle
    integer,                intent (in), optional :: markersize
    integer,                intent (in), optional :: linewidth
    integer,                intent (out)          :: istat

    character(len=:),                 allocatable :: xstr         !< x values stringified
    character(len=:),                 allocatable :: ystr         !< y values stringified
    character(len=:),                 allocatable :: zstr         !< z values stringified
    character(len=max_int_len)                    :: imark        !< actual markers size
    character(len=max_int_len)                    :: iline        !< actual line width
    character(len=*), parameter                   :: xname = 'x'  !< x variable name for script
    character(len=*), parameter                   :: yname = 'y'  !< y variable name for script
    character(len=*), parameter                   :: zname = 'z'  !< z variable name for script

    if (allocated(me%str)) then

        istat = 0

        !convert the arrays to strings:
        call vec_to_string(x, me%real_fmt, xstr, me%use_numpy)
        call vec_to_string(y, me%real_fmt, ystr, me%use_numpy)
        call vec_to_string(z, me%real_fmt, zstr, me%use_numpy)

        !get optional inputs (if not present, set default value):
        call optional_int_to_string(markersize, imark, '3')
        call optional_int_to_string(linewidth, iline, '3')

        !write the arrays:
        call me%add_str(trim(xname)//' = '//xstr)
        call me%add_str(trim(yname)//' = '//ystr)
        call me%add_str(trim(zname)//' = '//zstr)
        call me%add_str('')

        !write the plot statement:
        call me%add_str('ax.plot('//&
                        trim(xname)//','//&
                        trim(yname)//','//&
                        trim(zname)//','//&
                        '"'//trim(linestyle)//'",'//&
                        'linewidth='//trim(adjustl(iline))//','//&
                        'markersize='//trim(adjustl(imark))//','//&
                        'label="'//trim(label)//'")')
        call me%add_str('')

    else
       istat = initError('add_3d_plot')
    end if

  end subroutine add_3d_plot


  !=============================================================================
  !> @brief Add a sphere to a 3D x,y,z plot
  !> @note Must initialize the class with `mplot3d=.true.` and `use_numpy=.true.`.
  !> @param[inout] me        class handler
  !> @param[in] r            radius of the sphere
  !> @param[in] xc           x value of sphere center
  !> @param[in] yc           y value of sphere center
  !> @param[in] zc           z value of sphere center
  !> @param[in] n_facets     [default is 100]
  !> @param[in] linewidth    line width
  !> @param[in] antialiased  enabled anti-aliasing
  !> @param[in] color        color of the contour line
  !> @param[out] istat       status output (0 means no problems)

  subroutine add_sphere(me, r, xc, yc, zc, n_facets, linewidth, antialiased, color, istat)

    implicit none

    class(pyplot),    intent (inout)        :: me
    real(dp),         intent (in)           :: r
    real(dp),         intent (in)           :: xc
    real(dp),         intent (in)           :: yc
    real(dp),         intent (in)           :: zc
    integer,          intent (in), optional :: n_facets
    integer,          intent (in), optional :: linewidth
    logical,          intent (in), optional :: antialiased
    character(len=*), intent (in), optional :: color
    integer,          intent (out)          :: istat

    character(len=:),           allocatable :: rstr               !< `r` value stringified
    character(len=:),           allocatable :: xcstr              !< `xc` value stringified
    character(len=:),           allocatable :: ycstr              !< `yc` value stringified
    character(len=:),           allocatable :: zcstr              !< `zc` value stringified
    character(len=*),           parameter   :: xname = 'x'        !< `x` variable name for script
    character(len=*),           parameter   :: yname = 'y'        !< `y` variable name for script
    character(len=*),           parameter   :: zname = 'z'        !< `z` variable name for script

    character(len=max_int_len)              :: linewidth_str      !< `linewidth` input stringified
    character(len=:), allocatable           :: antialiased_str    !< `antialised` input stringified
    character(len=max_int_len)              :: n_facets_str       !< `n_facets` input stringified
    character(len=:), allocatable           :: extras             !< optional stuff string


    if (allocated(me%str)) then

        !get optional inputs (if not present, set default value):
        call optional_int_to_string(n_facets, n_facets_str, '100')
        allocate(extras, source='')
        if (present(linewidth)) then
            call optional_int_to_string(linewidth, linewidth_str, '1')
            extras = extras//','//'linewidth='//linewidth_str
        end if
        if (present(antialiased)) then
            call optional_logical_to_string(antialiased, antialiased_str, 'False')
            extras = extras//','//'antialiased='//antialiased_str
        end if
        if (present(color)) then
            extras = extras//','//'color="'//trim(color)//'"'
        end if

        istat = 0

        !convert the arrays to strings:
        call real_to_string(r , me%real_fmt, rstr)
        call real_to_string(xc, me%real_fmt, xcstr)
        call real_to_string(yc, me%real_fmt, ycstr)
        call real_to_string(zc, me%real_fmt, zcstr)

        ! sphere code:
        call me%add_str('u = np.linspace(0, 2 * np.pi, '//n_facets_str//')')
        call me%add_str('v = np.linspace(0, np.pi, '//n_facets_str//')')
        call me%add_str(xname//' = '//xcstr//' + '//rstr//' * np.outer(np.cos(u), np.sin(v))')
        call me%add_str(yname//' = '//ycstr//' + '//rstr//' * np.outer(np.sin(u), np.sin(v))')
        call me%add_str(zname//' = '//zcstr//' + '//rstr//' * np.outer(np.ones(np.size(u)), np.cos(v))')
        call me%add_str('ax.plot_surface('//xname//', '//yname//', '//zname//extras//')')
        call me%add_str('')

    else
       istat = initError('add_sphere')
    end if

  end subroutine add_sphere


  !=============================================================================
  !> @brief Add a bar plot
  !> @param[inout] me      class handler
  !> @param[in] x          x bar values
  !> @param[in] height     height bar values
  !> @param[in] label      plot label
  !> @param[in] width      width values
  !> @param[in] bottom     bottom values
  !> @param[in] color      plot color
  !> @param[in] yerr       yerr values
  !> @param[in] align      default: 'center'
  !> @param[in] xlim       x-axis range
  !> @param[in] ylim       y-axis range
  !> @param[in] xscale     example: 'linear' (default), 'log'
  !> @param[in] yscale     example: 'linear' (default), 'log'
  !> @param[out] istat     status output (0 means no problems)

  subroutine add_bar(me, x, height, label, width, bottom, color, &
                        yerr, align, xlim, ylim, xscale, yscale, istat)

    class(pyplot),          intent(inout)         :: me
    real(dp), dimension(:), intent(in)            :: x
    real(dp), dimension(:), intent(in)            :: height
    character(len=*),       intent(in)            :: label
    real(dp), dimension(:), intent(in),  optional :: width
    real(dp), dimension(:), intent(in),  optional :: bottom
    character(len=*),       intent(in),  optional :: color
    real(dp), dimension(:), intent(in),  optional :: yerr
    character(len=*),       intent(in),  optional :: align
    real(dp),dimension(2),  intent (in), optional :: xlim
    real(dp),dimension(2),  intent (in), optional :: ylim
    character(len=*),       intent (in), optional :: xscale
    character(len=*),       intent (in), optional :: yscale
    integer,                intent (out)          :: istat

    character(len=:), allocatable                 :: xstr               !< x axis values stringified
    character(len=:), allocatable                 :: ystr               !< y axis values stringified
    character(len=:), allocatable                 :: xlimstr            !< xlim values stringified
    character(len=:), allocatable                 :: ylimstr            !< ylim values stringified
    character(len=:), allocatable                 :: wstr               !< width values stringified
    character(len=:), allocatable                 :: bstr               !< bottom values stringified
    character(len=:), allocatable                 :: plt_str            !< plot string
    character(len=:), allocatable                 :: yerr_str           !<  yerr values stringified
    character(len=*), parameter                   :: xname = 'x'        !< x axis name
    character(len=*), parameter                   :: yname = 'y'        !< y axis name
    character(len=*), parameter                   :: wname = 'w'        !< width name
    character(len=*), parameter                   :: bname = 'b'        !< bottom name
    character(len=*), parameter                   :: yerrname = 'yerr'  !< yerr name

    if (allocated(me%str)) then

        istat = 0

        !axis limits (optional):
        if (present(xlim)) call vec_to_string(xlim, me%real_fmt, xlimstr, me%use_numpy)
        if (present(ylim)) call vec_to_string(ylim, me%real_fmt, ylimstr, me%use_numpy)

        !convert the arrays to strings:
                             call vec_to_string(x,      me%real_fmt, xstr,     me%use_numpy)
                             call vec_to_string(height, me%real_fmt, ystr,     me%use_numpy)
        if (present(width))  call vec_to_string(width,  me%real_fmt, wstr,     me%use_numpy)
        if (present(bottom)) call vec_to_string(bottom, me%real_fmt, bstr,     me%use_numpy)
        if (present(yerr))   call vec_to_string(yerr,   me%real_fmt, yerr_str, me%use_numpy)

        !write the arrays:
                             call me%add_str(trim(xname)//' = '//xstr)
                             call me%add_str(trim(yname)//' = '//ystr)
        if (present(width))  call me%add_str(trim(wname)//' = '//wstr)
        if (present(bottom)) call me%add_str(trim(bname)//' = '//bstr)
        if (present(yerr))   call me%add_str(trim(yerrname)//' = '//yerr_str)
        call me%add_str('')

        !create the plot string:
        allocate(plt_str, source = 'ax.bar('//&
                  'x='//trim(xname)//','//&
                  'height='//trim(yname)//',')
        if (present(yerr))   plt_str=plt_str//'yerr='//trim(yerrname)//','
        if (present(width))  plt_str=plt_str//'width='//trim(wname)//','
        if (present(bottom)) plt_str=plt_str//'bottom='//trim(bstr)//','
        if (present(color))  plt_str=plt_str//'color="'//trim(color)//'",'
        if (present(align))  plt_str=plt_str//'align="'//trim(align)//'",'
        plt_str=plt_str//'label="'//trim(label)//'")'

        !write the plot statement:
        call me%add_str(plt_str)

        !axis limits:
        if (allocated(xlimstr)) call me%add_str('ax.set_xlim('//xlimstr//')')
        if (allocated(ylimstr)) call me%add_str('ax.set_ylim('//ylimstr//')')

        !axis scales:
        if (present(xscale)) call me%add_str('ax.set_xscale("'//xscale//'")')
        if (present(yscale)) call me%add_str('ax.set_yscale("'//yscale//'")')

        call me%add_str('')

    else
       istat = initError('add_bar')
    end if

  end subroutine add_bar


  !=============================================================================
  !> @brief Add an image plot using `imshow`
  !> @param[inout] me      class handler
  !> @param[in] x          x values
  !> @param[in] xlim       x-axis range
  !> @param[in] ylim       y-axis range
  !> @param[in] exValues   array values to extent
  !> @param[out] istat     status output (0 means no problems)

  subroutine add_imshow(me, x, xlim, ylim, exValues, istat)

    class(pyplot),          intent (inout)        :: me
    real(dp),dimension(:,:),intent (in)           :: x
    real(dp),dimension(2),  intent (in), optional :: xlim, ylim
    real(dp),dimension(4),  intent (in), optional :: exValues
    integer,                intent (out)          :: istat

    character(len=:), allocatable                 :: leftstr      !< string for bounding box definition, pos left
    character(len=:), allocatable                 :: rightstr     !< string for bounding box definition, pos right
    character(len=:), allocatable                 :: bottomstr    !< string for bounding box definition, pos bottom
    character(len=:), allocatable                 :: topstr       !< string for bounding box definition, pos top

    character(len=:), allocatable                 :: xstr         !< x values stringified
    character(len=*), parameter                   :: xname = 'x'  !< x variable name for script

    !axis limits (optional):
    character(len=:), allocatable                 :: xlimstr      !< xlim values stringified
    character(len=:), allocatable                 :: ylimstr      !< ylim values stringified

    if (allocated(me%str)) then

        istat = 0

        if (present(xlim)) call vec_to_string(xlim, me%real_fmt, xlimstr, me%use_numpy)
        if (present(ylim)) call vec_to_string(ylim, me%real_fmt, ylimstr, me%use_numpy)

        !convert the arrays to strings:
        call matrix_to_string(x, me%real_fmt, xstr, me%use_numpy)

        !write the arrays:
        call me%add_str(trim(xname)//' = '//xstr)
        call me%add_str('')

        if(present(exValues)) then
          call real_to_string(exValues(1), me%real_fmt, leftstr)
          call real_to_string(exValues(2), me%real_fmt, rightstr)
          call real_to_string(exValues(3), me%real_fmt, bottomstr)
          call real_to_string(exValues(4), me%real_fmt, topstr)
          call me%add_str('ax.imshow(' //trim(xname)//  ',' // 'origin="' // 'lower"' // &
                ',' // 'cmap="' // 'viridis"' // &
                ',' // 'extent=(' // leftstr // ',' // rightstr // ',' // bottomstr // ',' // topstr  // ')' // &
                ',' // 'aspect="' // 'auto"' // ')')
          call me%add_str('')
        else
          ! write the plot statement:
          call me%add_str('ax.imshow('//trim(xname)//')')
          call me%add_str('')
        end if

        !axis limits:
        if (allocated(xlimstr)) call me%add_str('ax.set_xlim('//xlimstr//')')
        if (allocated(ylimstr)) call me%add_str('ax.set_ylim('//ylimstr//')')

    else
       istat = initError('add_imshow')
    end if

  end subroutine add_imshow


  !=============================================================================
  !> @brief Integer to string, specifying the default value if the optional argument is not present
  !> @param[in] int_value        integer
  !> @param[out] string_value    integer value stringified
  !> @param[in] default_value    default integer value

  subroutine optional_int_to_string(int_value, string_value, default_value)

    integer,          intent(in), optional :: int_value
    character(len=*), intent(out)          :: string_value
    character(len=*), intent(in)           :: default_value

    if (present(int_value)) then
        call integer_to_string(int_value, string_value)
    else
        string_value = default_value
    end if

  end subroutine optional_int_to_string


  !=============================================================================
  !> @brief Logical to string, specifying the default value if the optional argument is not present
  !> @param[in] logical_value        logical value
  !> @param[out] string_value        integer value stringified
  !> @param[in] default_value        default integer value

  subroutine optional_logical_to_string(logical_value, string_value, default_value)

    logical,intent(in),optional              :: logical_value
    character(len=:),allocatable,intent(out) :: string_value
    character(len=*),intent(in)              :: default_value

    if (present(logical_value)) then
        if (logical_value) then
            string_value = 'True'
        else
            string_value = 'False'
        end if
    else
        string_value = default_value
    end if

  end subroutine optional_logical_to_string


  !=============================================================================
  !> @brief  Integer to string conversion.
  !> @param[in] i       integer
  !> @param[out] s      integer value stringified

  subroutine integer_to_string(i, s)

    integer,          intent(in), optional  :: i
    character(len=*), intent(out)           :: s

    integer                                 :: istat !< IO status

    write(s,'(I10)',iostat=istat) i
    if (istat /= 0) then
        s = '****'
    else
        s = adjustl(s)
    end if

  end subroutine integer_to_string


  !=============================================================================
  !> @brief  Integer to string conversion.
  !> @param[in]  v      real
  !> @param[in]  fmt    real format string
  !> @param[out] str    real values stringified

  subroutine real_to_string(v, fmt, str)

    real(dp),                      intent(in)  :: v
    character(len=*),              intent(in)  :: fmt
    character(len=:), allocatable, intent(out) :: str

    integer                                    :: istat     !< IO status
    character(len=max_real_len)                :: tmp       !< dummy string

    if (fmt=='*') then
        write(tmp, *, iostat=istat) v
    else
        write(tmp, fmt, iostat=istat) v
    end if
    if (istat /= 0) then
        str = '****'
    else
        str = trim(adjustl(tmp))
    end if

  end subroutine real_to_string


  !=============================================================================
  !> @brief  Real vector to string
  !> @param[in]  v           real
  !> @param[out] fmt         real format string
  !> @param[out] str         real values stringified
  !> @param[in]  use_numpy   activate numpy python module usage
  !> @param[in]  is_tuple    if true [default], use '()', if false use '[]'

  subroutine vec_to_string(v, fmt, str, use_numpy, is_tuple)

    real(dp), dimension(:),        intent(in)  :: v
    character(len=*),              intent(in)  :: fmt
    character(len=:), allocatable, intent(out) :: str
    logical,                       intent(in)  :: use_numpy
    logical,intent(in),optional                :: is_tuple

    integer                                    :: i         !< counter
    integer                                    :: istat     !< IO status
    character(len=max_real_len)                :: tmp       !< dummy string
    logical :: tuple

    if (present(is_tuple)) then
        tuple = is_tuple
    else
        tuple = .false.
    end if

    if (tuple) then
        str = '('
    else
        str = '['
    end if

    do i=1, size(v)
        if (fmt=='*') then
            write(tmp, *, iostat=istat) v(i)
        else
            write(tmp, fmt, iostat=istat) v(i)
        end if
        if (istat /= 0) then
            str = '****'
            return
        end if
        str = str//trim(adjustl(tmp))
        if (i<size(v)) str = str // ','
    end do

    if (tuple) then
        str = str // ')'
    else
        str = str // ']'
    end if

    !convert to numpy array if necessary:
    if (use_numpy) str = 'np.array('//str//')'

  end subroutine vec_to_string


  !=============================================================================
  !> @brief  Real matrix (rank 2) to string
  !> @param[in] v            real
  !> @param[in] fmt          real format string
  !> @param[out] str         real values stringified
  !> @param[in] use_numpy    activate numpy python module usage

  subroutine matrix_to_string(v, fmt, str, use_numpy)

    real(dp), dimension(:,:),      intent(in)  :: v
    character(len=*),              intent(in)  :: fmt
    character(len=:), allocatable, intent(out) :: str
    logical,                       intent(in)  :: use_numpy

    integer                                    :: i         !> counter
    character(len=:),              allocatable :: tmp       !> dummy string

    str = '['
    do i=1, size(v,1)  !rows
        call vec_to_string(v(i,:), fmt, tmp, use_numpy)  !one row at a time
        str = str//trim(adjustl(tmp))
        if (i<size(v,1)) str = str // ','
    end do
    str = str // ']'

    !convert to numpy array if necessary:
    if (use_numpy) str = 'np.array('//str//')'

  end subroutine matrix_to_string


  !=============================================================================
  !> @brief Write the buffer to a specified file, then execute it with Python.
  !> @param[in] me     python plot handler
  !> @param[in] pyfile name of the python script to generate
  !> @param[out] istat status output (0 means no problems)

  subroutine execute(me, pyfile, istat)

    use reportErrorModule, only : error_p, reportError

    class(pyplot),    intent(in)  :: me
    character(len=*), intent(in)  :: pyfile
    integer,          intent(out) :: istat

    integer :: iunit !< File unit number for the generated python script

    istat = 0
    if (.not.(allocated(me%python) .and. allocated(me%str))) return

    !open the file:
    open(newunit=iunit, file=trim(pyfile), status='REPLACE', iostat=istat)
    if (istat /= 0) then
       call reportError (error_p,'Failed to open file '//pyfile, &
            &            addString='pyplot_module::execute',ierr=istat)
       return
    end if

    !write to the file:
    write(iunit, '(A)') me%str

    !to ensure that the file is there for the next
    !command line call, we have to close it here.
    close(iunit, iostat=istat)
    if (istat /= 0) then
       call reportError (error_p,'Failed to close file '//pyfile, &
            &            addString='pyplot_module::execute',ierr=istat)
       return
    end if

    !run the file using python:
    call execute_command_line(me%python//trim(pyfile),exitstat=istat)
    if (istat /= 0) then
       call reportError (error_p,'Failed to execute python script '//pyfile, &
            &            addString='pyplot_module::execute',ierr=istat)
    end if

  end subroutine execute


  !=============================================================================
  !> @brief Some final things to add before saving or showing the figure
  !> @param[inout] me         pytplot handler

  subroutine finish_ops(me)

    class(pyplot),intent(inout) :: me  !! pyplot handler

    if (me%show_legend) then
        call me%add_str('ax.legend(loc="best")')
        call me%add_str('')
    end if
    if (me%axis_equal) then
        if (me%mplot3d) then
            call me%add_str('ax.set_aspect("equal")')
            call me%add_str('')

            call me%add_str('def set_axes_equal(ax):')
            call me%add_str('    x_limits = ax.get_xlim3d()')
            call me%add_str('    y_limits = ax.get_ylim3d()')
            call me%add_str('    z_limits = ax.get_zlim3d()')
            call me%add_str('    x_range = abs(x_limits[1] - x_limits[0])')
            call me%add_str('    x_middle = np.mean(x_limits)')
            call me%add_str('    y_range = abs(y_limits[1] - y_limits[0])')
            call me%add_str('    y_middle = np.mean(y_limits)')
            call me%add_str('    z_range = abs(z_limits[1] - z_limits[0])')
            call me%add_str('    z_middle = np.mean(z_limits)')
            call me%add_str('    plot_radius = 0.5*max([x_range, y_range, z_range])')
            call me%add_str('    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])')
            call me%add_str('    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])')
            call me%add_str('    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])')
            call me%add_str('set_axes_equal(ax)')

        else
            call me%add_str('ax.axis("equal")')
        end if
        call me%add_str('')
    end if
    if (me%tight_layout) then
        call me%add_str('fig.tight_layout()')
        call me%add_str('')
    end if

  end subroutine finish_ops


  !=============================================================================
  !> @brief Save the figure
  !> @param[inout] me         pytplot handler
  !> @param[in] figfile       file name for the figure
  !> @param[in] pyfile        name of the Python script to generate
  !> @param[in] dpi           resolution of the figure for png
  !> @param[in] transparent   transparent background (T/F)
  !> @param[in] facecolor     the colors of the figure rectangle
  !> @param[in] edgecolor     the colors of the figure rectangle
  !> @param[in] orientation   'landscape' or 'portrait'
  !> @param[out] istat        status output (0 means no problems)

  subroutine savefig(me, figfile, pyfile, dpi, transparent, facecolor, edgecolor, orientation, istat)

    use reportErrorModule, only : debugFileOnly_p, reportError

    class(pyplot),    intent(inout)        :: me
    character(len=*), intent(in)           :: figfile
    character(len=*), intent(in), optional :: pyfile
    character(len=*), intent(in), optional :: dpi

    logical, intent(in), optional          :: transparent
    character(len=*), intent(in), optional :: facecolor
    character(len=*), intent(in), optional :: edgecolor
    character(len=*), intent(in), optional :: orientation
    integer,          intent (out)         :: istat

    character(len=:),          allocatable :: tmp  !! for building the `savefig` arguments.

    if (allocated(me%str)) then

        !finish up the string:
        call me%finish_ops()
        !build the savefig arguments:
        allocate(tmp, source='"'//trim(figfile)//'"')
        if (present(dpi)) tmp = tmp//', dpi='//trim(dpi)
        if (present(transparent)) then
            if (transparent) then
                tmp = tmp//', transparent=True'
            else
                tmp = tmp//', transparent=False'
            end if
        end if
        if (present(facecolor)) tmp = tmp//', facecolor="'//trim(facecolor)//'"'
        if (present(edgecolor)) tmp = tmp//', edgecolor="'//trim(edgecolor)//'"'
        if (present(orientation)) tmp = tmp//', orientation="'//trim(orientation)//'"'
        if (me%use_oo_api) then
            call me%add_str('canvas = FigureCanvas(fig)')
            call me%add_str('canvas.print_figure('//tmp//')')
        else
            call me%add_str('plt.savefig('//tmp//')')
        end if
        deallocate(tmp)

        !run it:
        call me%execute(pyfile, istat=istat)
        if (istat /= 0) then
           call reportError (debugFileOnly_p,'pyplot_module::savefig')
        end if

    else
       istat = initError('savefig')
    end if

  end subroutine savefig


  !=============================================================================
  !> @brief Show the figure
  !> @param[inout] me         pytplot handler
  !> @param[in] pyfile        name of the Python script to generate
  !> @param[out] istat        status output (0 means no problems)

  subroutine showfig(me, pyfile, istat)

    use reportErrorModule, only : error_p, debugFileOnly_p, reportError

    class(pyplot),    intent(inout)        :: me
    character(len=*), intent(in), optional :: pyfile
    integer,          intent(out)          :: istat

    if (.not. allocated(me%str)) then
       istat = initError('showfig')
    else if (me%use_oo_api) then
       istat = -2
       call reportError (error_p,'Not compatible with "use_oo_api" option', &
            &            addString='pyplot_module::showfig')
    else

       !finish up the string:
       call me%finish_ops()
       !show figure:
       call me%add_str('plt.show()')

       !run it:
       call me%execute(pyfile, istat=istat)
       if (istat /= 0) then
          call reportError (debugFileOnly_p,'pyplot_module::showfig')
       end if

    end if

  end subroutine showfig

end module pyplot_module
