!> @brief
!> Module of elfe3D containing subroutines to obtain source arrays
!> for dipole and loop sources
!!
!> written by Paula Rulff, 23/07/2019, extended 04/2020
!!
!> Copyright (C) Paula Rulff 2020
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

module calculate_global_source

  use mod_util
  use vector_products

  implicit none

contains

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for detecting source edges and their directions
  !>
  !> for x-directed electric dipole source 
  !>
  !> store it in global vector "source" with +/- 1
  !---------------------------------------------------------------------
  subroutine HED_x (E, direction, nd, ed2nd, &
                    sx_start, sy_start, sz_start, sx_end, source)

    ! INPUT
    integer, intent(in) :: E, direction
    real(kind=dp), dimension(:,:), intent(in) :: nd 
    integer, dimension(:,:), intent(in) :: ed2nd
    real(kind=dp), intent(in) :: sx_start, sy_start, sz_start, sx_end


    ! OUTPUT
    real(kind=dp), dimension(:), intent(inout) :: source

    ! LOCAL variables
    integer :: i

    !-------------------------------------------------------------------
    global_edge_loop: do i = 1,E  

      if(abs(nd(ed2nd(i,1),2) - sy_start) <= eps_dp .and. & ! y
         abs(nd(ed2nd(i,1),3) - sz_start) <= eps_dp .and. & ! z
         ! node1 x on source edge
         nd(ed2nd(i,1),1) .ge. sx_start .and. &
         nd(ed2nd(i,1),1) .le. sx_end .and. &! x

         abs(nd(ed2nd(i,2),2) - sy_start) <= eps_dp .and. & ! y
         abs(nd(ed2nd(i,2),3) - sz_start) <= eps_dp .and. & ! z
         ! node2 x on source edge
         nd(ed2nd(i,2),1) .ge. sx_start .and. &
         nd(ed2nd(i,2),1) .le. sx_end) then


            ! "source sign" similar to local edge sign 
            ! with increasing node-number start & end node of 
            ! source edge i
            source(i) = ed2nd(i,1)-ed2nd(i,2)
            source(i) = source(i)/abs(source(i))
            ! ensure that current will point in correct direction
            ! on the source edge
            if (nd(ed2nd(i,1),1) < nd(ed2nd(i,2),1)) then
              select case (direction)
                case (0) ! positive curent
                  source(i) = source(i)*(-1.0_dp)
                case (1) ! negative curent
                  source(i) = source(i)*(1.0_dp)
               end select
            else if (nd(ed2nd(i,1),1) > nd(ed2nd(i,2),1)) then
              select case (direction)
                case (0) ! positive curent
                  source(i) = source(i)*(1.0_dp)
                case (1) ! negative curent
                  source(i) = source(i)*(-1.0_dp)
                end select
            end if
      end if
    end do global_edge_loop 

  end subroutine HED_x

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for detecting source edges and their directions
  !>
  !> for y-directed electric dipole source
  !>
  !> stores it in global vector "source" with +/- 1
  !---------------------------------------------------------------------
  subroutine HED_y (E, direction, nd, ed2nd, &
                    sx_start, sy_start, sz_start, sy_end, source)

    ! INPUT
    integer, intent(in) :: E, direction
    real(kind=dp), dimension(:,:), intent(in) :: nd 
    integer, dimension(:,:), intent(in) :: ed2nd
    real(kind=dp), intent(in) :: sx_start, sy_start, sz_start, sy_end

    ! OUTPUT
    real(kind=dp), dimension(:), intent(inout) :: source

    ! LOCAL variables
    integer :: i

    !-------------------------------------------------------------------
    global_edge_loop: do i = 1,E  

      if(abs(nd(ed2nd(i,1),1) - sx_start) <= eps_dp .and. & ! x
         abs(nd(ed2nd(i,1),3) - sz_start) <= eps_dp .and. & ! z
         ! node1 y on source edge
         nd(ed2nd(i,1),2) .ge. sy_start .and. &
         nd(ed2nd(i,1),2) .le. sy_end .and. &! x

         abs(nd(ed2nd(i,2),1) - sx_start) <= eps_dp .and. & ! x
         abs(nd(ed2nd(i,2),3) - sz_start) <= eps_dp .and. & ! z
         ! node2 y on source edge
         nd(ed2nd(i,2),2) .ge. sy_start .and. &
         nd(ed2nd(i,2),2) .le. sy_end) then

            ! "source sign" similar to local edge sign 
            ! with increasing node-number start & end node of 
            ! source edge i
            source(i) = ed2nd(i,1)-ed2nd(i,2)
            source(i) = source(i)/abs(source(i))
            ! ensure that current will point in correct direction on
            ! the source edge
            if (nd(ed2nd(i,1),2) < nd(ed2nd(i,2),2)) then
              select case (direction)
                case (0) ! positive current
                  source(i) = source(i)*(-1.0_dp)
                case (1) ! negative curent
                  source(i) = source(i)*(1.0_dp)
              end select
            else if (nd(ed2nd(i,1),2) > nd(ed2nd(i,2),2)) then
              select case (direction)
                case (0) ! positive current
                  source(i) = source(i)*(1.0_dp)
                case (1) ! negative curent
                  source(i) = source(i)*(-1.0_dp)
                end select
            end if

      end if

    end do global_edge_loop 

  end subroutine HED_y

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for detecting source edges and their directions
  !>
  !> for a horizontal square-loop source 
  !>
  !> stores it in global vector "source" with +/- 1
  !---------------------------------------------------------------------
  subroutine loop_source (E, direction, midp_source, nd, ed2nd, &
                          sx_start, sy_start, sz_start, sx_end, source)

    ! INPUT
    integer, intent(in) :: E, direction
    ! source midpoint (x,y,z)
    real(kind=dp), dimension(3), intent(in) :: midp_source 
    real(kind=dp), dimension(:,:), intent(in) :: nd 
    integer, dimension(:,:), intent(in) :: ed2nd
    real(kind=dp), intent(in) :: sx_start, sy_start, sz_start, sx_end


    ! OUTPUT
    real(kind=dp), dimension(:), intent(inout) :: source

    ! LOCAL variables
    integer :: current_direction
    real(kind=dp) :: source_diameter, sx_start_dummy, &
                     sy_start_dummy, sy_end_dummy

    !-------------------------------------------------------------------
    ! Initialise dummy source-corner-coordinates to zero
    sx_start_dummy = 0.0_dp
    sy_start_dummy = 0.0_dp
    sy_end_dummy = 0.0_dp

    source_diameter = 2*abs(sy_start - midp_source(2))
    !print *, source_diameter,'loop source diameter'

    select case (direction)
      case (0) ! clockwise current

          ! source edge 1 (N)
          current_direction = 0
          call HED_x (E, current_direction, nd, ed2nd, &
                      sx_start, sy_start, sz_start, sx_end, source)

          ! source edge 4 (W)
          current_direction = 0
          sy_end_dummy = sy_start
          sy_start_dummy = sy_start-source_diameter
          call HED_y (E, current_direction, nd, ed2nd, &
                      sx_start, sy_start_dummy, sz_start, &
                      sy_end_dummy, source)

          ! source edge 3 (S)
          current_direction = 1
          call HED_x (E, current_direction, nd, ed2nd, &
                      sx_start, sy_start_dummy, sz_start, &
                      sx_end, source)

          ! source edge 2 (E)
          current_direction = 1
          sx_start_dummy = sx_end
          call HED_y (E, current_direction, nd, ed2nd, &
                      sx_start_dummy, sy_start_dummy, sz_start, &
                      sy_end_dummy, source)


      case (1) ! anticlockwise current

          ! source edge 1 (N)
          current_direction = 1
          call HED_x (E, current_direction, nd, ed2nd, &
                      sx_start, sy_start, sz_start, sx_end, source)

          ! source edge 4 (W)
          current_direction = 1
          sy_end_dummy = sy_start
          sy_start_dummy = sy_start-source_diameter
          call HED_y (E, current_direction, nd, ed2nd, &
                      sx_start, sy_start_dummy, sz_start, &
                      sy_end_dummy, source)

          ! source edge 3 (S)
          current_direction = 0
          call HED_x (E, current_direction, nd, ed2nd, &
                      sx_start, sy_start_dummy, sz_start, &
                      sx_end, source)

          ! source edge 2 (E)
          current_direction = 0
          sx_start_dummy = sx_end
          call HED_y (E, current_direction, nd, ed2nd, &
                      sx_start_dummy, sy_start_dummy, &
                      sz_start, sy_end_dummy, source)

    end select

  end subroutine loop_source

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for detecting source edges and their directions
  !>
  !> for x-directed arbitrary electric dipole source, store it in
  !> global vector "source" with +/- 1
  !>
  !> ! might not work for refinement as node markers are used !
  !---------------------------------------------------------------------
  subroutine arbitrary_HED_x (E, direction, nodemarker, nd, ed2nd, &
                              source)

    ! INPUT
    integer, intent(in) :: E, direction
    integer, dimension(:), intent(in) :: nodemarker
    real(kind=dp), dimension(:,:), intent(in) :: nd 
    integer, dimension(:,:), intent(in) :: ed2nd


    ! OUTPUT
    real(kind=dp), dimension(:), intent(inout) :: source

    ! LOCAL variables
    integer :: i

    !-------------------------------------------------------------------
    global_edge_loop: do i = 1,E  
      ! source-nodes have the nodemarker 3
      if (nodemarker(ed2nd(i,1)) .eq. 3 .and. &
          nodemarker(ed2nd(i,2)) .eq. 3 ) then
            ! "source sign" similar to local edge sign 
            ! with increasing node-number start & end node of source 
            ! edge i
            source(i) = ed2nd(i,1)-ed2nd(i,2)
            source(i) = source(i)/abs(source(i))
            ! ensure that current will point in correct direction
            ! on the source edge
            if (nd(ed2nd(i,1),1) < nd(ed2nd(i,2),1)) then
              select case (direction)
                case (0) ! positive curent
                  source(i) = source(i)*(-1.0_dp)
                case (1) ! negative curent
                  source(i) = source(i)*(1.0_dp)
               end select
            else if (nd(ed2nd(i,1),1) > nd(ed2nd(i,2),1)) then
              select case (direction)
                case (0) ! positive curent
                  source(i) = source(i)*(1.0_dp)
                case (1) ! negative curent
                  source(i) = source(i)*(-1.0_dp)
                end select
            end if

      end if

    end do global_edge_loop 

  end subroutine arbitrary_HED_x

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for detecting source edges and their directions
  !>
  !> for y-directed electric dipole source, store it in
  !> global vector "source" with +/- 1
  !>
  !> ! might not work for refinement as node markers are used !
  !---------------------------------------------------------------------
  subroutine arbitrary_HED_y (E, direction, nodemarker, nd, ed2nd, &
                              source)

    ! INPUT
    integer, intent(in) :: E, direction
    integer, dimension(:), intent(in) :: nodemarker
    real(kind=dp), dimension(:,:), intent(in) :: nd 
    integer, dimension(:,:), intent(in) :: ed2nd

    ! OUTPUT
    real(kind=dp), dimension(:), intent(inout) :: source

    ! LOCAL variables
    integer :: i

    !-------------------------------------------------------------------
    global_edge_loop: do i = 1,E  

      ! source-nodes have the nodemarker 3
      if (nodemarker(ed2nd(i,1)) .eq. 3 .and. &
          nodemarker(ed2nd(i,2)) .eq. 3 ) then
            ! "source sign" similar to local edge sign 
            ! with increasing node-number start & end node of source
            ! edge i
            source(i) = ed2nd(i,1)-ed2nd(i,2)
            source(i) = source(i)/abs(source(i))
            ! ensure that current will point in correct direction on
            ! the source edge
            if (nd(ed2nd(i,1),2) < nd(ed2nd(i,2),2)) then
              select case (direction)
                case (0) ! positive current
                  source(i) = source(i)*(-1.0_dp)
                case (1) ! negative curent
                  source(i) = source(i)*(1.0_dp)
              end select
            else if (nd(ed2nd(i,1),2) > nd(ed2nd(i,2),2)) then
              select case (direction)
                case (0) ! positive current
                  source(i) = source(i)*(1.0_dp)
                case (1) ! negative curent
                  source(i) = source(i)*(-1.0_dp)
                end select
            end if

      end if

    end do global_edge_loop 

  end subroutine arbitrary_HED_y

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for detecting source edges and their directions
  !> for one straight source segment
  !>
  !> global vector "source" with +/- 1
  !---------------------------------------------------------------------
  subroutine straight_source_segment (E, direction, nd, ed2nd, &
                                      sx_start, sy_start, sz_start, &
                                      sx_end, sy_end, sz_end, source)

    ! INPUT
    integer, intent(in) :: E, direction
    real(kind=dp), dimension(:,:), intent(in) :: nd 
    integer, dimension(:,:), intent(in) :: ed2nd
    real(kind=dp), intent(in) :: sx_start, sy_start, sz_start, &
                                 sx_end, sy_end, sz_end


    ! OUTPUT
    real(kind=dp), dimension(:), intent(inout) :: source

    ! LOCAL variables
    integer :: i, ORIENTATION
    real(kind=dp), dimension(3) :: p1, p2, cp1, cp2
    real(kind=dp) :: length_x, length_y, length_z

    !-------------------------------------------------------------------
    ! Detect if source segment is more x, y, or z x-oriented
    length_x = abs(sx_end - sx_start) 
    length_y = abs(sy_end - sy_start) 
    length_z = abs(sz_end - sz_start) 

    if (length_x .ge. length_y .and. length_x .ge. length_z) then
      ORIENTATION = 0 ! x-dipole
    else if (length_y .gt. length_x .and. length_y .ge. length_z) then
      ORIENTATION = 1 ! y-dipole
    else if (length_z .gt. length_x .and. length_z .gt. length_y) then
      ORIENTATION = 2 ! z-dipole
    else
      ORIENTATION = 999 ! undefined
    end if

    ! Define start and endpoint of line segment
    ! (starting-point: smaller x,y,z-coordinate) - might be unnecesary?
    select case (ORIENTATION)
      case(0)
        if (sx_start .lt. sx_end) then
          p1 = (/sx_start, sy_start, sz_start/)
          p2 = (/sx_end, sy_end, sz_end/)
        else
          p2 = (/sx_start, sy_start, sz_start/)
          p1 = (/sx_end, sy_end, sz_end/)
        end if
      case(1)
        if (sy_start .lt. sy_end) then
          p1 = (/sx_start, sy_start, sz_start/)
          p2 = (/sx_end, sy_end, sz_end/)
        else
          p2 = (/sx_start, sy_start, sz_start/)
          p1 = (/sx_end, sy_end, sz_end/)
        end if
      case(2)
        if (sz_start .lt. sz_end) then
          p1 = (/sx_start, sy_start, sz_start/)
          p2 = (/sx_end, sy_end, sz_end/)
        else
          p2 = (/sx_start, sy_start, sz_start/)
          p1 = (/sx_end, sy_end, sz_end/)
        end if
      case default
          call Write_Message (log_unit, &
            'Straight source sgment: Wrong source orientation')
    end select


    global_edge_loop: do i = 1,E  

      ! egde start and endpoint
      cp1 = nd(ed2nd(i,1),:)
      cp2 = nd(ed2nd(i,2),:)

         ! node1 on source segment
      if(dist2segment(p1,p2,cp1) .le. eps_dp .and. &  
         ! node2 on source segment
         dist2segment(p1,p2,cp2) .le. eps_dp) then 

            ! "source sign" similar to local edge sign 
            ! with increasing node-number start & end node of source
            ! edge i
            source(i) = ed2nd(i,1)-ed2nd(i,2)
            source(i) = source(i)/abs(source(i))

            select case (ORIENTATION)
              ! case HED more x-oriented:
              case(0)
                ! ensure that current will point in correct direction
                ! on the source edge
                if (nd(ed2nd(i,1),1) < nd(ed2nd(i,2),1)) then
                  select case (direction)
                    case (0) ! positive curent
                      source(i) = source(i)*(-1.0_dp)
                    case (1) ! negative curent
                      source(i) = source(i)*(1.0_dp)
                   end select
                else if (nd(ed2nd(i,1),1) > nd(ed2nd(i,2),1)) then
                  select case (direction)
                    case (0) ! positive curent
                      source(i) = source(i)*(1.0_dp)
                    case (1) ! negative curent
                      source(i) = source(i)*(-1.0_dp)
                    end select
                end if

              ! case HED more y-oriented:
              case(1)
                ! ensure that current will point in correct direction
                ! on the source edge
                if (nd(ed2nd(i,1),2) < nd(ed2nd(i,2),2)) then
                  select case (direction)
                    case (0) ! positive current
                      source(i) = source(i)*(-1.0_dp)
                    case (1) ! negative curent
                      source(i) = source(i)*(1.0_dp)
                  end select
                else if (nd(ed2nd(i,1),2) > nd(ed2nd(i,2),2)) then
                  select case (direction)
                    case (0) ! positive current
                      source(i) = source(i)*(1.0_dp)
                    case (1) ! negative curent
                      source(i) = source(i)*(-1.0_dp)
                    end select
                end if

              ! case HED more z-oriented:
              case(2)
                ! ensure that current will point in correct direction
                ! on the source edge
                if (nd(ed2nd(i,1),3) < nd(ed2nd(i,2),3)) then
                  select case (direction)
                    case (0) ! positive current
                      source(i) = source(i)*(-1.0_dp)
                    case (1) ! negative curent
                      source(i) = source(i)*(1.0_dp)
                  end select
                else if (nd(ed2nd(i,1),3) > nd(ed2nd(i,2),3)) then
                  select case (direction)
                    case (0) ! positive current
                      source(i) = source(i)*(1.0_dp)
                    case (1) ! negative curent
                      source(i) = source(i)*(-1.0_dp)
                    end select
                end if
            end select


      end if

    end do global_edge_loop 

  end subroutine straight_source_segment

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for detecting  the minimum distance of a point (cp)
  !> to a line-segment defined by two points (p1,p2)
  !---------------------------------------------------------------------
  function dist2segment (p1, p2, cp)

    ! INPUT
     real(kind=dp), dimension(3), intent(in) :: p1, p2, cp

    ! OUTPUT
     real(kind=dp) :: dist2segment

    ! LOCAL variables
    real(kind=dp) :: p1p2, p1cp, p2cp ! distances

    !-------------------------------------------------------------------
    ! calculate the distance of point cp to the straight line
    ! through p1,p2
    dist2segment = sqrt(sum((cross_product_real(p2-p1,p1-cp))**2.0_dp))& 
                   / sqrt(sum((p2-p1)**2.0_dp))

    ! if the point lies on the line...
    if (dist2segment .le. eps_dp) then

      ! Check if the point lies within the segment:
      ! Find the distance of point cp from both the line end points
      ! p1, p2. 
      ! If p1p2 = p1cp + p2cp, then cp lies on the line segment p1p2
      p1p2 = sqrt((p2(1)-p1(1))*(p2(1)-p1(1))+ &
                  (p2(2)-p1(2))*(p2(2)-p1(2))+ &
                  (p2(3)-p1(3))*(p2(3)-p1(3)))
      p1cp = sqrt((cp(1)-p1(1))*(cp(1)-p1(1))+ &
                  (cp(2)-p1(2))*(cp(2)-p1(2))+ &
                  (cp(3)-p1(3))*(cp(3)-p1(3)))
      p2cp = sqrt((p2(1)-cp(1))*(p2(1)-cp(1))+ &
                  (p2(2)-cp(2))*(p2(2)-cp(2))+ &
                  (p2(3)-cp(3))*(p2(3)-cp(3)))

      ! if point is outside segment (p1p2 /= p1cp + p2cp), 
      ! set dist2segment == 999
      if(p1p2 .ne. p1cp + p2cp) then
         dist2segment = 999.0_dp
      end if

    end if

  end function dist2segment


  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for detecting source edges and their directions
  !> for a segmented line or loop source
  !>
  !> global vector source with +/- 1
  !>
  !> The file: 'in/source.txt' has to be provided
  !>
  !> best imeplemented option 
  !---------------------------------------------------------------------
  subroutine segmented_source (E, direction, nd, ed2nd, CSTYPE, source)

    ! INPUT
    integer, intent(in) :: E, direction
    real(kind=dp), dimension(:,:), intent(in) :: nd 
    integer, dimension(:,:), intent(in) :: ed2nd
    integer, intent(in) :: CSTYPE


    ! OUTPUT
    real(kind=dp), dimension(:), intent(inout) :: source

    ! LOCAL variables
    character(len=20):: FileName
    integer :: i, j, num_segments, num_corners
    integer :: in_unit = 1
    integer :: opening
    integer :: allo_stat
    integer :: current_direction, comp_direction
    real(kind=dp) :: sx_start, sy_start, sz_start, &
                     sx_end, sy_end, sz_end
    real(kind=dp), allocatable, dimension(:,:) :: segments
    real(kind=dp) :: length_x, length_y, length_z
    !-------------------------------------------------------------------
    ! initialise
    num_corners = 0
    comp_direction = 9

    ! Read source segment coordinates from in/source.txt and write them 
    ! into array segments
    ! First entry is the NW-most source coordinate
    FileName = 'in/source.txt'

    ! open the file
    open (in_unit, file = trim(FileName), status='old', &
          action = 'read', iostat = opening)

    ! was opening successful?
    if (opening /= 0) then
        call Write_Error_Message(log_unit, &
                'segmented_source: file ' // trim(FileName) // &
                ' could not be opened')
    else


    
       ! read number of nodes from first line
       ! (only first column, ignore second column)
       read (in_unit, *) num_corners

        
       ! allocate array
       allocate (segments(num_corners,3), stat = allo_stat)
       call allocheck(log_unit, allo_stat, &
                   "segmented_source: error allocating array segments" )
       ! initialise
       segments = D0
            
    
       do i = 1,num_corners
            read (in_unit,*) (segments(i,j),j=1,3)
       end do
            
       ! test output
       ! write (*,'(a)') "segment coordinates are ="
       ! write (*,*) segments
        
       close (unit = in_unit)

    end if


    select case(CSTYPE)

      case(6) ! segmented line source
        num_segments = num_corners-1

        do i = 1,num_segments
          sx_start = segments(i,1)
          sy_start = segments(i,2)
          sz_start = segments(i,3)
          sx_end = segments(i+1,1)
          sy_end = segments(i+1,2)
          sz_end = segments(i+1,3)

          call straight_source_segment (E, direction, nd, ed2nd, &
                                        sx_start, sy_start, sz_start, &
                                        sx_end, sy_end, sz_end, source)
        end do

      case(7) ! segmented horizontal loop source
        num_segments = num_corners

        ! loop segments
        do i = 1,num_segments-1
          sx_start = segments(i,1)
          sy_start = segments(i,2)
          sz_start = segments(i,3)
          sx_end = segments(i+1,1)
          sy_end = segments(i+1,2)
          sz_end = segments(i+1,3)

          length_x = abs(sx_end - sx_start) 
          length_y = abs(sy_end - sy_start) 


          ! do which compass direction does the segment belong?
          ! N: comp_direction = 0
          ! E: comp_direction = 1
          ! S: comp_direction = 2
          ! W: comp_direction = 3

          if (sx_start .lt. sx_end .and. &
              length_x .ge. length_y) then
            comp_direction = 0 ! N 
          else if (sy_start .gt. sy_end .and. &
                   length_y .gt. length_x) then
            comp_direction = 1 ! E
          else if (sx_start .gt. sx_end .and. &
                   length_x .ge. length_y) then
            comp_direction = 2 ! S
          else if (sy_start .lt. sy_end .and. &
                   length_y .gt. length_x) then
            comp_direction = 3 ! W
          end if

          ! Depending on the source-current direction, 
          ! the compass direction, 
          ! decide which current_direction is needed for the segment i
          select case (direction)
            case (0) ! clockwise current

              select case (comp_direction)
                case(0) ! N
                  current_direction = 0
                case(1) ! E
                  current_direction = 1
                case(2) ! S
                  current_direction = 1
                case(3) ! W
                  current_direction = 0
              end select

            case (1) ! anticlockwise current  

              select case (comp_direction)
                case(0) ! N
                  current_direction = 1
                case(1) ! E
                  current_direction = 0
                case(2) ! S
                  current_direction = 0
                case(3) ! W
                  current_direction = 1
              end select

          end select

          call straight_source_segment (E, current_direction, nd, &
                                        ed2nd, &
                                        sx_start, sy_start, sz_start, &
                                        sx_end, sy_end, sz_end, source)

         end do

        ! segmented vertical electric loop source (resulting in HMD)
        case(9) 
            num_segments = num_corners

            ! loop segments
            do i = 1,num_segments-1
              sx_start = segments(i,1)
              sy_start = segments(i,2)
              sz_start = segments(i,3)
              sx_end = segments(i+1,1)
              sy_end = segments(i+1,2)
              sz_end = segments(i+1,3)

              length_x = abs(sx_end - sx_start) 
              length_y = abs(sy_end - sy_start)
              length_z = abs(sz_end - sz_start) 


              ! do which compass direction does the segment belong?
              ! N: comp_direction = 0
              ! E: comp_direction = 1
              ! S: comp_direction = 2
              ! W: comp_direction = 3

              if (sx_start .lt. sx_end .and. &
                  length_x .ge. length_y) then
                comp_direction = 0 ! N 
              else if (sy_start .gt. sy_end .and. &
                       length_y .gt. length_x) then
                comp_direction = 1 ! E
              else if (sx_start .gt. sx_end .and. &
                       length_x .ge. length_y) then
                comp_direction = 2 ! S
              else if (sy_start .lt. sy_end .and. &
                       length_y .gt. length_x) then
                comp_direction = 3 ! W
              end if

              ! Depending on the source-current direction, 
              ! the compass direction,
              ! decide which current_direction is needed for the
              ! segment i
              select case (direction)
                case (0) ! clockwise current

                  select case (comp_direction)
                    case(0) ! N
                      current_direction = 0
                    case(1) ! E
                      current_direction = 1
                    case(2) ! S
                      current_direction = 1
                    case(3) ! W
                      current_direction = 0
                  end select

                case (1) ! anticlockwise current  

                  select case (comp_direction)
                    case(0) ! N
                      current_direction = 1
                    case(1) ! E
                      current_direction = 0
                    case(2) ! S
                      current_direction = 0
                    case(3) ! W
                      current_direction = 1
                  end select

              end select

              call straight_source_segment (E, current_direction, nd, &
                                            ed2nd, &
                                            sx_start, sy_start, &
                                            sz_start, &
                                            sx_end, sy_end, sz_end, &
                                            source)

           end do

        

        ! last segment to close the loop
        i = num_segments
        sx_start = segments(i,1)
        sy_start = segments(i,2)
        sz_start = segments(i,3)
        sx_end = segments(1,1)
        sy_end = segments(1,2)
        sz_end = segments(1,3)

        length_x = abs(sx_end - sx_start) 
        length_y = abs(sy_end - sy_start) 
        length_z = abs(sz_end - sz_start) 


        ! do which compass direction does the segment belong?
        ! N: comp_direction = 0
        ! E: comp_direction = 1
        ! S: comp_direction = 2
        ! W: comp_direction = 3

        if (sx_start .lt. sx_end .and. length_x .ge. length_y) then
          comp_direction = 0 ! N 
        else if (sy_start .gt. sy_end .and. length_y .gt. length_x) then
          comp_direction = 1 ! E
        else if (sx_start .gt. sx_end .and. length_x .ge. length_y) then
          comp_direction = 2 ! S
        else if (sy_start .lt. sy_end .and. length_y .gt. length_x) then
          comp_direction = 3 ! W
        end if

        ! Depending on the source-current direction, 
        ! the compass direction, 
        ! decide which current_direction is needed for the segment i
        select case (direction)
          case (0) ! clockwise current

            select case (comp_direction)
              case(0) ! N
                current_direction = 0
              case(1) ! E
                current_direction = 1
              case(2) ! S
                current_direction = 1
              case(3) ! W
                current_direction = 0
            end select

          case (1) ! anticlockwise current  

            select case (comp_direction)
              case(0) ! N
                current_direction = 1
              case(1) ! E
                current_direction = 0
              case(2) ! S
                current_direction = 0
              case(3) ! W
                current_direction = 1
            end select

        end select

        
        call straight_source_segment (E, current_direction, nd, &
                                      ed2nd, &
                                      sx_start, sy_start, sz_start, &
                                      sx_end, sy_end, sz_end, source)

    end select

    if (allocated(segments)) deallocate(segments)


  end subroutine segmented_source

!-----------------------------------------------------------------------
end module calculate_global_source
