!> @brief
!> Module of elfe3d containing subroutines
!> to define initial model parameters (from input file)
!> 'in/elfe3D_input.txt'
! 
!> written by Paula Rulff, 20/06/2019
!> modified by Paula Rulff, 25/07/2024
!!
!> Copyright (C) Paula Rulff 2024
!>
!>  This file is part of elfe3D.
!> 
!>  Licensed under the Apache License, Version 2.0 (the "License"); 
!>  you may not use this file except in compliance with the License.  
!>  You may obtain a copy of the License at
!> 
!>      https://www.apache.org/licenses/LICENSE-2.0
!> 
!>  Unless required by applicable law or agreed to in writing, software
!>  distributed under the License is distributed on an "AS IS" BASIS, 
!>  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or 
!>  implied.  See the
!>  License for the specific language governing permissions and 
!>  limitations under the License. 

module define_model

  use mod_util

  implicit none

  ! LOCAL VARIABLES
  ! length of names of I/O files
  integer, parameter, public :: fnl = 200
  ! length of lines in I/O files
  integer, parameter, public :: lil = 200
  ! format specifier for input lines of length lil
  character(len=6), parameter, public :: lfm = "(a200)"
  ! length of method strings
  integer, parameter, public :: msl=20
  ! elfe3D start file
  character(len=fnl) :: FileName = 'in/elfe3D_input.txt'
  ! variables to read file
  integer :: in_unit = 20
  integer :: opening
  integer :: skipcol
  integer :: OpenCode, ReadCode, ctmpCode
  integer :: bwrd,ewrd
  character(len=lil) :: ctmp, ccom

contains

  !---------------------------------------------------------------------
  !> @brief'
  !> subroutine for defining solver
  !> Define the following in input file:
  !> solver
  !> 1: PARDISO (currently not available!)
  !> 2: MUMPS
  !---------------------------------------------------------------------
  subroutine define_solver (solver)
  
    ! OUTPUT
    integer, intent(out) :: solver

    !-------------------------------------------------------------------
    ! Read from elfe3D_input.txt
    ! open the file
    open (in_unit, file = trim(FileName), status='old', &
                   action = 'read', iostat = opening)

    ! was opening successful?
    if (opening /= 0) then
        call Write_Error_Message(log_unit, &
        'define_solver: file '//trim(FileName)//' could not be opened')
    else
       ! read parameters
       read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       line_read_loop: do while (ReadCode == 0)
          if (index(ctmp,'solver') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (unit=ctmp(bwrd:ewrd),fmt=*,iostat=ctmpCode) solver
             if ((ctmpCode /= 0) .or. (solver <= 0) .or. &
                                      (solver> 2)) then
                call Check_Input(log_unit, 'solver')
             end if
          end if
          ! read next line
          read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       end do line_read_loop
       ! close input file
       close (unit = in_unit)
    end if
    
   !--------------------------------------------------------------------
  end subroutine define_solver

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for defining refinement parameters
  !> Define the following in input file:
  !> maxRefSteps, maxUnknowns, betaRef, accuracyTol, vtk , 
  !> errorEst_method , refStrategy             0
  !>
  !> Maximum number of refinement steps:
  !> maxRefSteps = 0 means no refinement
  !> run several frequencies only without refinement
  !> Maximum number of unknowns:
  !> maxUnknowns  
  !> threshold for number of elements to be refined:
  !> betaRef
  !> (e.g. 0.9 = those with 90% highest error estimator)
  !> desired accuracy tolerance < 1:
  !> accuracyTol
  !> write .vtk files for paraview with error estimates (yes/no) = (1/0)
  !> vtk
  !> Method for error estimation:
  !> errorEst_method
  !> residuals (1)
  !> residuals and face jumps J (2)
  !> residuals and face jumps J and H (3)
  !> face jumps J (4)
  !> face jumps H (5)
  !> face jumps J & H  (6)
  !> Refinement strategy:
  !> refStrategy
  !> constant quality factor (0)
  !> maxRefSteps-1 on low quality mesh, last step high quality mesh (1)
  !> increasing quality factor (2)
  !> increasing quality factor on mesh with detailled subsurface anomaly
  !> (-T and -d  option added) (3)
  !---------------------------------------------------------------------
  subroutine define_refinement (maxRefSteps, maxUnknowns, accuracyTol, &
                                betaRef, vtk, errorEst_method, &
                                refStrategy)
  
    ! OUTPUT
    integer, intent(out) :: maxRefSteps, maxUnknowns, vtk, &
                            errorEst_method, refStrategy
    real(kind=dp), intent(out) :: accuracyTol, betaRef 
    !-------------------------------------------------------------------
    ! Read from elfe3D_input.txt
    ! open the file
    open (in_unit, file = trim(FileName), status='old', &
                   action = 'read', iostat = opening)

    ! was opening successful?
    if (opening /= 0) then
        call Write_Error_Message(log_unit, &
        'define_solver: file '//trim(FileName)//' could not be opened')
    else
       ! read parameters
       read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       line_read_loop: do while (ReadCode == 0)
          ! maxRefSteps
          if (index(ctmp,'maxRefSteps') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (unit=ctmp(bwrd:ewrd),fmt=*,iostat=ctmpCode) &
                                                             maxRefSteps
             if ((ctmpCode /= 0) .or. (maxRefSteps < 0)) then
                call Check_Input(log_unit, 'maxRefSteps')
             end if
          ! maxUnknowns
          else if (index(ctmp,'maxUnknowns') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (unit=ctmp(bwrd:ewrd),fmt=*,iostat=ctmpCode) &
                                                             maxUnknowns
             if ((ctmpCode /= 0) .or. (maxUnknowns < 0)) then
                call Check_Input(log_unit, 'maxUnknowns')
             end if
          ! betaRef
          else if (index(ctmp,'betaRef') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (unit=ctmp(bwrd:ewrd),fmt=*,iostat=ctmpCode) betaRef
             if ((ctmpCode /= 0) .or. (betaRef < 0) &
                                 .or. (betaRef > 1)) then
                call Check_Input(log_unit, 'betaRef')
             end if
          ! accuracyTol
          else if (index(ctmp,'accuracyTol') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (unit=ctmp(bwrd:ewrd),fmt=*,iostat=ctmpCode) &
                                                             accuracyTol
             if ((ctmpCode /= 0) .or. (accuracyTol <= 0) &
                                 .or. (accuracyTol >= 1)) then
                call Check_Input(log_unit, 'accuracyTol')
             end if
          ! vtk
          else if (index(ctmp,'vtk') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (unit=ctmp(bwrd:ewrd),fmt=*,iostat=ctmpCode) vtk
             if ((ctmpCode /= 0) .or. (vtk < 0) &
                                 .or. (vtk > 1)) then
                call Check_Input(log_unit, 'vtk')
             end if
          ! errorEst_method
          else if (index(ctmp,'errorEst_method') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (unit=ctmp(bwrd:ewrd),fmt=*,iostat=ctmpCode) &
                                                         errorEst_method
             if ((ctmpCode /= 0) .or. (errorEst_method < 0) &
                                 .or. (errorEst_method > 6)) then
                call Check_Input(log_unit, 'errorEst_method')
             end if
          ! refStrategy
          else if (index(ctmp,'refStrategy') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (unit=ctmp(bwrd:ewrd),fmt=*,iostat=ctmpCode) &
                                                         refStrategy
             if ((ctmpCode /= 0) .or. (refStrategy < 0) &
                                 .or. (refStrategy > 3)) then
                call Check_Input(log_unit, 'erefStrategy')
             end if
          end if
          ! read next line
          read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       end do line_read_loop
       ! close input file
       close (unit = in_unit)
    end if
    
   !--------------------------------------------------------------------
  end subroutine define_refinement

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for defining mesh
  !> Define the following in input file: model_file_name
  !---------------------------------------------------------------------
  subroutine define_mesh (NodeFile, EdgeFile, ElementFile, NeighFile, &
                         refStep, StringName)
  
    ! INPUT
    ! current refinement step
    integer, intent(in) :: refStep
    ! OUTPUT
    character(len = 255), intent(out) :: NodeFile, EdgeFile, &
                                         ElementFile, NeighFile
    character(len = 50), intent(out) :: StringName
    ! number of nodes
    !integer, intent(out) :: N, E, M
    
    ! LOCAL variables:
    character(len = 50) :: StringStep, StringEnding
    !-------------------------------------------------------------------
    ! Set the Filenames for the current refinement step:
    write(StringStep , *) refStep

    ! Read from elfe3D_input.txt
    ! open the file
    open (in_unit, file = trim(FileName), status='old', &
                   action = 'read', iostat = opening)

    ! was opening successful?
    if (opening /= 0) then
        call Write_Error_Message(log_unit, &
        'define_mesh: file '//trim(FileName)//' could not be opened')
    else
       ! read model file name
       read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       line_read_loop: do while (ReadCode == 0)
          ! model_file_name 
          if (index(ctmp,'model_file_name') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (ctmp(bwrd:ewrd),fmt='(a)', iostat=ctmpCode) &
                                                              StringName
             if ((ctmpCode /= 0)) then
                call Check_Input(log_unit, 'model_file_name')
             end if
          end if
          ! read next line
          read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       end do line_read_loop
       ! close input file
       close (unit = in_unit)
    end if
    
    ! How to build the filenames: // = string concatenation operator
    ! adjustl() -> align string left, 
    ! remove leading blanks, e.g. ("     text" -> "text     ")
    ! trim() -> remove trailing blanks, e.g. ("text     " -> "text")
    
    ! Filenames are always: "Name.Step.Ending" 
    ! e.g. 'in/nodes_refined.5.txt'
    
    ! File containing nodes with coordinates:
    StringEnding = ".node"
    NodeFile = trim(adjustl(StringName))//trim(adjustl(StringStep))// &
               trim(adjustl(StringEnding))
    
    ! File containing edges withy nodenumbers
    StringEnding = ".edge"
    EdgeFile = trim(adjustl(StringName))//trim(adjustl(StringStep))// &
               trim(adjustl(StringEnding))
    
    !File containing elements with node numbers
    StringEnding = ".ele"
    ElementFile = trim(adjustl(StringName))//trim(adjustl(StringStep)) &
                //trim(adjustl(StringEnding))

    ! File containing elements with their neighbour-elements
    StringEnding = ".neigh"
    NeighFile = trim(adjustl(StringName))//trim(adjustl(StringStep))// &
                trim(adjustl(StringEnding))
   !--------------------------------------------------------------------
  end subroutine define_mesh

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for defining size of the source
  !>
  !> Define the following in input file:
  !> source_type
  !> source start coordinates
  !> source end coordinates
  !> current_direction
  !> source_moment
  !>
  !> Options for source_type:
  !> HED_x (0), HED_y (1), loop_source (2), arbitrary HED_x (3), 
  !> arbitrary HED_y (4), straight_source_segment in any direction (5)
  !> segmented_line_source (6), segmented_loop_source (7)
  !> PR comment:  only 6 or 7 would be sufficient

  !> HED_x (0), HED_y (1) are straight HEDs in exactly either x- or 
  !> y- direction
  !> loop_source (2) is a squared loop source defined by northern
  !> line-source and midpoint, segments exactly on x- and y- direction
  !> if you choose arbitrary HED, you have to use node-boundary 
  !> markers = 3 for all source nodes in your .poly-inputfile
  !> this will be the fastest option, but it does probably not work 
  !> for refinement!
  !> straight_source_segment in any direction (5) every edge between 
  !> source start and source endpoint will be a source edge
  !> line/loop source split into straight source segments - provide 
  !> input file with coordinates of segment nodes! (6/7)
  !> Loop source coordinates always in clockwise direction!


  !> Define coordinates of the horizontal source:
  !> if z negative downwards:
  !> x positive to E
  !> y positive to N
  !> z positive upwards
  !> in case of a horizontal loop source, define north line-source 
  !> start & endpoints

  !> if z positive downwards:
  !> x positive to N
  !> y positive to E
  !> z positive downwards
  !> in case of a horizontal loop source, define east line-source 
  !> start & endpoint

  !> Start (starting-point/smaller x,y,z-coordinate) of
  !> electric dipole source
  !> End (end-point/larger x,y,z-coordinate) of electric dipole source


  !> Define direction:
  !> if z negative downwards:
  !> for CSTYPE = HED:
  !> current in positive direction: direction = 0
  !> current in negative direction: direction = 1
  !> for CSTYPE = loop_source:
  !> clockwise current: direction  = 0
  !> anticlocwise current: direction = 1

  !> if z positive downwards:
  !> for CSTYPE = HED:
  !> current in positive direction: direction = 1
  !> current in negative direction: direction = 0
  !> for CSTYPE = loop_source:
  !> clockwise current: direction  = 1
  !> anticlocwise current: direction = 0
  !---------------------------------------------------------------------
  subroutine define_source (sx_start, sy_start, sz_start, sx_end, &
                            sy_end, sz_end, CSTYPE, direction, p, &
                            midp_source)

    ! OUTPUT
    real(kind=dp), intent(out) :: sx_start, sy_start, sz_start, &
                                  sx_end, sy_end, sz_end
    integer, intent(out) :: CSTYPE, direction
    real(kind=dp), intent(out) :: p ! source moment
    ! source midpoint (x,y,z)
    real(kind=dp), dimension(3), intent(out) :: midp_source 
    !-------------------------------------------------------------------
    ! Read from elfe3D_input.txt
    ! open the file
    open (in_unit, file = trim(FileName), status='old', &
                   action = 'read', iostat = opening)

    ! was opening successful?
    if (opening /= 0) then
        call Write_Error_Message(log_unit, &
        'define_solver: file '//trim(FileName)//' could not be opened')
    else
       ! read parameters
       read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       line_read_loop: do while (ReadCode == 0)
          ! source_type
          if (index(ctmp,'source_type') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (unit=ctmp(bwrd:ewrd),fmt=*,iostat=ctmpCode) CSTYPE
             if ((ctmpCode /= 0) .or. (CSTYPE < 0)) then
                call Check_Input(log_unit, 'source_type')
             end if
             ! read source start and end coordinates
             read (in_unit,*) sx_start, sy_start, sz_start 
             read (in_unit,*) sx_end, sy_end, sz_end
          ! current_direction 
          else if (index(ctmp,'current_direction') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (unit=ctmp(bwrd:ewrd),fmt=*,iostat=ctmpCode) direction
             if ((ctmpCode /= 0) .or. (direction < 0) &
                                 .or. (direction > 1)) then
                call Check_Input(log_unit, 'current_direction')
             end if
          ! source_moment 
          else if (index(ctmp,'source_moment ') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (unit=ctmp(bwrd:ewrd),fmt=*,iostat=ctmpCode) p
             if ((ctmpCode /= 0) .or. (p < 0.0_dp)) then
                call Check_Input(log_unit, 'source_moment ')
             end if
          end if  
          ! read next line
          read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       end do line_read_loop
       ! close input file
       close (unit = in_unit)
    end if

    !-------------------------------------------------------------------
    select case (CSTYPE)

      case (0,1) ! HED
        ! Calculate source midpoint (for line-source)
        midp_source(1) = 0.5 * (sx_start + sx_end)
        midp_source(2) = 0.5 * (sy_start + sy_end)
        midp_source(3) = 0.5 * (sz_start + sz_end)

      case (2) ! loop source
        ! Define coordinates of loop midpoint here:
        midp_source(1) = -0.625_dp ! x-coordinate
        midp_source(2) = -0.625_dp ! y-coordinate
        midp_source(3) = 0.0_dp ! z-coordinate

    end select
   !--------------------------------------------------------------------
  end subroutine define_source



  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for defining existence and size of a PEC
  !> to model metallic borehole casing
  !> Define the following in input file:
  !> no PEC present: 
  !> PEC_present = 0
  !> PEC present: 
  !> PEC_present = 1
  !> define number of PEC:
  !> num_PEC
  !> followed by coordinates of start points
  !> PEC1 start x,y,z
  !> PEC2 start x,y,z
  !> PEC.. start x,y,z
  !> PEC1 end x,y,z
  !> PEC2 end x,y,z
  !> PEC.. end x,y,z
  !---------------------------------------------------------------------
  subroutine define_PEC (PEC, num_PEC, PEC_start, PEC_end)

    ! OUTPUT
    integer, intent(out) :: PEC, num_PEC
    real(kind=dp), allocatable, intent(out), dimension(:,:) :: &
                                                      PEC_start, PEC_end

    ! LOCAL variables
    integer :: allo_stat
    integer :: i
    !-------------------------------------------------------------------
    ! Define PEC presence

    ! 1 (yes)
    ! 0 (no) - default
    ! default
    PEC = 0

    ! Read from elfe3D_input.txt
    ! open the file
    open (in_unit, file = trim(FileName), status='old', &
                   action = 'read', iostat = opening)

    ! was opening successful?
    if (opening /= 0) then
        call Write_Error_Message(log_unit, &
        'define_solver: file '//trim(FileName)//' could not be opened')
    else
       ! read parameters
       read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       line_read_loop: do while (ReadCode == 0)
          ! PEC 
          if (index(ctmp,'PEC_present') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (unit=ctmp(bwrd:ewrd),fmt=*,iostat=ctmpCode) PEC
             if ((ctmpCode /= 0) .or. (PEC < 0) .or. (PEC > 1)) then
                call Check_Input(log_unit, 'PEC_present')
             end if
          ! num_PEC 
          else if (index(ctmp,'num_PEC') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (unit=ctmp(bwrd:ewrd),fmt=*,iostat=ctmpCode) num_PEC
             if ((ctmpCode /= 0) .or. (num_PEC < 0)) then
                call Check_Input(log_unit, 'num_PEC')
             end if
             ! allocate space fuer PEC coordinates
             if (PEC == 1) then
               ! allocate
               allocate (PEC_start(num_PEC,3), PEC_end(num_PEC,3), &
                         stat = allo_stat)
               call allocheck(log_unit, allo_stat, &
                         "define_PEC: error allocating PEC variables")
               ! read PEC start and end coordinates
               do i = 1,num_PEC
                 read (in_unit,*) PEC_start(i,1), &
                                  PEC_start(i,2), &
                                  PEC_start(i,3)
               end do 
               do i = 1,num_PEC
                 read (in_unit,*) PEC_end(i,1), &
                                  PEC_end(i,2), &
                                  PEC_end(i,3)
               end do 
             end if
          end if  
          ! read next line
          read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       end do line_read_loop
       ! close input file
       close (unit = in_unit)
    end if
   !--------------------------------------------------------------------
  end subroutine define_PEC


  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for defining output files
  !> Define the following in input file: output_E_file, output_E_file
  !---------------------------------------------------------------------
  subroutine define_output (EFile, HFile, num_rec)
  
    ! OUTPUT
    character(len = 500), intent(out) :: EFile, HFile
    integer, intent(out) :: num_rec
    
    !-------------------------------------------------------------------
    ! Read from elfe3D_input.txt
    ! open the file
    open (in_unit, file = trim(FileName), status='old', &
                   action = 'read', iostat = opening)

    ! was opening successful?
    if (opening /= 0) then
        call Write_Error_Message(log_unit, &
        'define_output: file '//trim(FileName)//' could not be opened')
    else
       ! read in elfe3D_input file
       read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       line_read_loop: do while (ReadCode == 0)
          ! num_rec
          if (index(ctmp,'num_rec') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (unit=ctmp(bwrd:ewrd),fmt=*,iostat=ctmpCode) num_rec
             if ((ctmpCode /= 0) .or. (num_rec < 0)) then
                call Check_Input(log_unit, ' num_rec')
             end if
          ! read output_E_file name
          else if (index(ctmp,'output_E_file') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (ctmp(bwrd:ewrd),fmt='(a)', iostat=ctmpCode) EFile
             if ((ctmpCode /= 0)) then
                call Check_Input(log_unit, 'output_E_file')
             end if
          ! read output_H_file name
          else if (index(ctmp,'output_H_file') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (ctmp(bwrd:ewrd),fmt='(a)', iostat=ctmpCode) HFile
             if ((ctmpCode /= 0)) then
                call Check_Input(log_unit, 'output_H_file')
             end if
          end if
          ! read next line
          read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       end do line_read_loop
       ! close input file
       close (unit = in_unit)
    end if

   !--------------------------------------------------------------------
  end subroutine define_output

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for defining frequencies
  !> Define the following in input file: 
  !> num_freq                1
  !> followed by a list of frequencies
  !---------------------------------------------------------------------

  subroutine define_freq (Nfreq, freq)
  
    ! OUTPUT
    integer, intent(out) :: Nfreq
    real (kind=dp), allocatable, dimension(:) :: freq

    ! LOCAL variables
    integer :: allo_stat
    integer :: i
    !-------------------------------------------------------------------
    ! Read from elfe3D_input.txt
    ! open the file
    open (in_unit, file = trim(FileName), status='old', &
                   action = 'read', iostat = opening)

    ! was opening successful?
    if (opening /= 0) then
        call Write_Error_Message(log_unit, &
        'define_output: file '//trim(FileName)//' could not be opened')
    else
       ! read in elfe3D_input file
       read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       line_read_loop: do while (ReadCode == 0)
          ! num_freq
          if (index(ctmp,'num_freq') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (unit=ctmp(bwrd:ewrd),fmt=*,iostat=ctmpCode) Nfreq
             if ((ctmpCode /= 0) .or. (Nfreq < 0)) then
                call Check_Input(log_unit, 'num_freq')
             end if
             allocate (freq(Nfreq),stat = allo_stat)
               call allocheck(log_unit, allo_stat, &
                        "read_model_param: error allocating array freq")
              ! initialise
              freq = 0.0_dp
              ! read frequencies
              do i = 1,Nfreq
                read (in_unit,*) freq(i)
              end do    
          end if
          ! read next line
          read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       end do line_read_loop
       ! close input file
       close (unit = in_unit)
    end if

   !--------------------------------------------------------------------
  end subroutine define_freq

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for defining model_size
  !> Define the following in input file below keyword model_size
  !> minimum x,y,z
  !> maximum x,y,z
  !---------------------------------------------------------------------
  subroutine define_model_size (x_min, x_max, &
                                y_min, y_max, &
                                z_min, z_max)

  
    ! OUTPUT
    real(kind=dp), intent(out) :: x_min, x_max, y_min, &
                                  y_max, z_min, z_max
    !------------------------------------------------------------------- 
    ! Read from elfe3D_input.txt
    ! open the file
    open (in_unit, file = trim(FileName), status='old', &
                   action = 'read', iostat = opening)

    ! was opening successful?
    if (opening /= 0) then
        call Write_Error_Message(log_unit, &
        'define_output: file '//trim(FileName)//' could not be opened')
    else
       ! read in elfe3D_input file
       read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       line_read_loop: do while (ReadCode == 0)
          ! find model_size keyword
          if (index(ctmp,'model_size') > 0) then
              ! read size
              read (in_unit,*) x_min, y_min, z_min  
              read (in_unit,*) x_max, y_max, z_max
          end if
          ! read next line
          read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       end do line_read_loop
       ! close input file
       close (unit = in_unit)
    end if

   !--------------------------------------------------------------------
  end subroutine define_model_size

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for defining receiver locations
  !> Define the following in input file: 
  !> num_rec
  !> followed by a list of receiver coordinates x,y,z
  !---------------------------------------------------------------------
  subroutine define_rec (M, num_rec, a, b, c, d, Ve, u1, v1, w1, &
                         rec1_el)
  
    ! Input variables
    integer, intent(in) :: M, num_rec
    !real(kind=dp), dimension(:,:), intent(in) :: nd
    !integer, dimension(:,:), intent(in) :: el2nd
    !real(kind=dp), dimension(:), intent(in) :: rho
    real(kind=dp), dimension(:,:), intent(in) :: a,b,c,d
    real(kind=dp), dimension(:), intent(in) :: Ve

    ! OUTPUT
    real(kind=dp), allocatable, dimension(:), intent(out) :: u1,v1,w1
    integer, allocatable, dimension(:), intent(out) :: rec1_el

    ! LOCAL variables
    integer :: allo_stat, l, i
    real(kind=dp) :: bc_nd1, bc_nd2, bc_nd3, bc_nd4, factor
    integer :: num_rec_dummy
    !------------------------------------------------------------------- 
    ! Read from elfe3D_input.txt
    ! open the file
    open (in_unit, file = trim(FileName), status='old', &
                   action = 'read', iostat = opening)

    ! was opening successful?
    if (opening /= 0) then
        call Write_Error_Message(log_unit, &
        'define_output: file '//trim(FileName)//' could not be opened')
    else
       ! read in elfe3D_input file
       read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       line_read_loop: do while (ReadCode == 0)
          ! find num_rec
          if (index(ctmp,'num_rec') > 0) then
             bwrd = BegWrd(ctmp,2)
             ewrd = EndWrd(ctmp,2)
             read (unit=ctmp(bwrd:ewrd),fmt=*,iostat=ctmpCode) &
                                                           num_rec_dummy
             if ((ctmpCode /= 0) .or. (num_rec_dummy < 0)) then
                call Check_Input(log_unit, 'num_rec')
             end if
             allocate (u1(num_rec), v1(num_rec), w1(num_rec), &
                       rec1_el(num_rec), stat = allo_stat)
                 call allocheck(log_unit, allo_stat, &
                     "define_rec: error allocating receiver arrays")
              ! initialise
              u1 = 0.0_dp
              v1 = 0.0_dp
              w1 = 0.0_dp
              rec1_el = 0.0_dp
              ! read locations
              do i = 1,num_rec
                read (in_unit,*) u1(i),v1(i), w1(i)
              end do    
          end if
          ! read next line
          read (unit=in_unit, fmt=lfm, iostat=ReadCode) ctmp
       end do line_read_loop
       ! close input file
       close (unit = in_unit)
    end if


    !------------------------------------------------------------------- 
    ! Find receiver-earth elements using linear interpolation 
    ! functions a,b,c,d

    receiver_loop: do l = 1, num_rec

      ! Element loop
      do i = 1,M

        ! If you insert the point coordinates x, y, z in the linear
        ! interp functions a + bx + cy+ dz (at nodes) and
        ! their values are all in between 0 and 1, 
        ! then this point is inside the tetrahedron.
        ! If one of the barycentric coordinates is larger then 1 or 
        ! smaller then 0, the point is outside.
        factor = (1.0_dp/ (6.0_dp * Ve(i)))

        bc_nd1 = factor * &
             (a(i,1) + b(i,1) * u1(l) + c(i,1) * v1(l) + d(i,1) * w1(l))
        bc_nd2 = factor * &
             (a(i,2) + b(i,2) * u1(l) + c(i,2) * v1(l) + d(i,2) * w1(l))
        bc_nd3 = factor * &
             (a(i,3) + b(i,3) * u1(l) + c(i,3) * v1(l) + d(i,3) * w1(l))
        bc_nd4 = factor * &
            (a(i,4) + b(i,4) * u1(l) + c(i,4) * v1(l) + d(i,4) * w1(l))

        if (bc_nd1 .gt. 0.0_dp .and. bc_nd1 .lt. 1.0_dp .and. & ! node1

            bc_nd2 .gt. 0.0_dp .and. bc_nd2 .lt. 1.0_dp .and. & ! node2
 
            bc_nd3 .gt. 0.0_dp .and. bc_nd3 .lt. 1.0_dp .and. & ! node3

            bc_nd4 .gt. 0.0_dp .and. bc_nd4 .lt. 1.0_dp) then  ! node4

            rec1_el(l) = i

        end if

      end do 

      write (*,'(a)') "Receiver number and its Earth element:"
      print *, l, rec1_el(l)
    end do receiver_loop
   !--------------------------------------------------------------------
  end subroutine define_rec

end module define_model

