!> @brief
!> Module of elfe3d containing subroutines to read mesh from ASCII files
!> created with TetGen
!!
!> written by Paula Rulff, 28/08/2018
!> Last change: March 2024
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
module read_mesh

  use mod_util

  implicit none

contains

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for reading *.nodes file
  !---------------------------------------------------------------------
  subroutine read_nodes (FileName,N,nd,nodemarker)

    ! INPUT
    character(len=*), intent(in) :: FileName

    ! OUTPUT
    real(kind=dp), allocatable, dimension(:,:), intent(out) :: nd
    ! number of nodes
    integer, intent(out) :: N
    !node boundary marker
    integer, allocatable, dimension(:), intent(out) :: nodemarker

    ! LOCAL variables
    integer :: in_unit = 1
    integer :: opening
    integer :: i,j
    integer :: allo_stat
    integer :: skipcol

    !-----------------------------------------------------------------
    ! number of nodes (number of lines to be read in from in node-file)
    N = 0 
    
    ! open the file
    open (in_unit, file = trim(FileName), status='old', &
                   action = 'read', iostat = opening)

    ! was opening successful?
    if (opening /= 0) then
        call Write_Error_Message(log_unit, &
        'read_nodes: file ' // trim(FileName) // ' could not be opened')
    else
        ! read number of nodes from first line 
        ! (only first column, ignore second column)
        read (in_unit, *) N
        ! allocate space for nd and nodemarker
        allocate (nd(N,3), stat = allo_stat)
            call allocheck(log_unit, allo_stat, &
                           "read_nodes: error allocating array nd" )
        allocate (nodemarker(N), stat = allo_stat)
            call allocheck(log_unit, allo_stat, &
                       "read_nodes: error allocating array nodemarker" )

        do i = 1,N
           read (in_unit,*) skipcol,(nd(i,j),j=1,3)!,nodemarker(i)
        end do

       close (unit = in_unit)

    end if

  end subroutine read_nodes
  
  !---------------------------------------------------------------------	
  !> @brief
  !> subroutine for reading *.edges file
  !---------------------------------------------------------------------
  subroutine read_edges (FileName,E,ed2nd, edgemarker)

    ! INPUT
    character(len=*), intent(in) :: FileName

    ! OUTPUT
    integer, allocatable, dimension(:,:), intent(out) :: ed2nd
    ! number of edges
    integer, intent(out) :: E
    ! edge boundary marker
    integer, allocatable, dimension(:), intent(out) :: edgemarker

    ! LOCAL variables
    integer :: in_unit = 2
    integer :: opening
    integer :: i,j
    integer :: allo_stat
    integer :: skipcol
    
    !-------------------------------------------------------------------
    ! determine number of edges (= number of lines in edge-file)
    E = 0 
    
    ! open the file
    open (in_unit, file = trim(FileName), status='old', &
                   action = 'read', iostat = opening)

    ! was opening successful?
    if (opening /= 0) then
        call Write_Error_Message(log_unit, &
        'read_edges: file ' // trim(FileName) // ' could not be opened')
    else
        ! read number of edges from first line
        read (in_unit, *) E
        ! allocate space for ed2nd and edgemarker
        allocate (ed2nd(E,2), stat = allo_stat)
        call allocheck(log_unit, allo_stat, &
                             "read_edges: error allocating array ed2nd")
        allocate (edgemarker(E), stat = allo_stat)
            call allocheck(log_unit, allo_stat, &
                       "read_edges: error allocating array edgemarker" )

       do i = 1,E
            read (in_unit,*) skipcol,(ed2nd(i,j),j=1,2),edgemarker(i)
       end do

       close (unit = in_unit)

    end if

  end subroutine read_edges

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for reading *.elements file
  !---------------------------------------------------------------------
  subroutine read_elements(FileName,M,el2nd,eleattr)

    ! INPUT
    character(len=*), intent(in) :: FileName

    ! OUTPUT
    integer, allocatable, dimension(:,:), intent(out) :: el2nd
    ! number of elements
    integer, intent(out) :: M
    ! element attributes
    integer, allocatable, dimension(:), intent(out) :: eleattr
    
    ! LOCAL variables
    integer :: in_unit = 3
    integer :: opening
    integer :: i,j
    integer :: allo_stat
    integer :: skipcol

    !-------------------------------------------------------------------
    ! determine number of elements (= number of lines in elements-file)
    M = 0 
    
    ! open the file
    open (in_unit, file = trim(FileName), status='old', &
                   action = 'read', iostat = opening)

    ! was opening successful?
    if (opening /= 0) then
        call Write_Error_Message(log_unit, &
     'read_elements: file ' // trim(FileName) // ' could not be opened')
    else
        ! read number of elements from first line
        read (in_unit, *) M
        ! allocate space for el2nd and eleattr
        allocate (el2nd(M,4), stat = allo_stat)
        call allocheck(log_unit, allo_stat, &
                       "read_elements: error allocating array el2nd")  
        allocate (eleattr(M), stat = allo_stat)
        call allocheck(log_unit, allo_stat, &
                       "read_elements: error allocating array eleattr" )


        do i = 1,M
            read (in_unit,*) skipcol,(el2nd(i,j),j=1,4),eleattr(i)
        end do

       close (unit = in_unit)

    end if

  end subroutine read_elements


  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for reading *.neigh file
  !---------------------------------------------------------------------
  subroutine read_neigh(FileName,M,el2neigh)

    ! INPUT
    character(len=*), intent(in) :: FileName
    ! number of elements
    integer, intent(in) :: M

    ! OUTPUT
    integer, allocatable, dimension(:,:), intent(out) :: el2neigh

    ! LOCAL variables
    integer :: in_unit = 3
    integer :: opening
    integer :: i,j, M_check
    integer :: allo_stat
    integer :: skipcol

    
    !-------------------------------------------------------------------
    ! determine number of elements (= number of lines in elements-file)
    M_check = 0 
    
    ! open the file
    open (in_unit, file = trim(FileName), status='old', &
                   action = 'read', iostat = opening)

    ! was opening successful?
    if (opening /= 0) then
        call Write_Error_Message(log_unit, &
        'read_neigh: file ' // trim(FileName) // ' could not be opened')
    else
        ! read number of elements from first line
        read (in_unit, *) M_check
        
        ! check if number of elements in .ele file is consistent with
        ! number of elements in .neigh file
        if (M /= M_check) then
            call Write_Message (log_unit, &
               'Element numbers in .ele and .neigh files are different')
        end if
        
        ! allocate el2neigh map (Mx4)
        allocate (el2neigh(M,4), stat = allo_stat)
        call allocheck(log_unit, allo_stat, &
                       "read_elements: error allocating array el2neigh")
            
       do i = 1,M
          read (in_unit,*) skipcol,(el2neigh(i,j),j=1,4)
       end do

       close (unit = in_unit)

    end if
    
  end subroutine read_neigh
  !---------------------------------------------------------------------
end module read_mesh
