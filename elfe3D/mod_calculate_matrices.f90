!> @brief
!> Module of elfe3D containing subroutines to calculate 
!> connectivity matrix elements2edges, coordinate matrices, edge signs
!!
!> written by Paula Rulff, 24/08/2018
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

module calculate_matrices

  use mod_util

  implicit none

contains

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for calculating connectivity matrix el2ed 
  !---------------------------------------------------------------------
  subroutine calc_connect_matrix (el2nd, ed2nd, E, M, el2ed)

    ! INPUT
    ! total numbers of edges and elements
    integer, intent(in) :: E,M
    ! declaration of array for element nodes
    integer, dimension(:,:), intent(in) :: el2nd ! (1:M,1:4)
    ! declaration of array for edge nodes
    integer, dimension(:,:), intent(in) :: ed2nd ! (1:E,1:2)

    ! OUTPUT
    ! declaration of array for element edges and element edge lengths
    integer, allocatable, dimension(:,:), intent(out) :: el2ed

    ! LOCAL variables
    integer :: i,j
    integer :: allo_stat
    !-------------------------------------------------------------------
    allocate (el2ed(M,6), stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
                    "calc_connect_matrix: error allocating array el2ed")

    ! searching for the combination of two nodes belonging to one edge 
    ! in ed2nd and assigning the edge to the element matrix as in table

    ! initialise
    el2ed = 0

    !$OMP PARALLEL DO
      element_loop: do i = 1,M
        
         edge_loop: do j = 1,E

            if (((el2nd(i,1) .eq. ed2nd(j,1)) .and. &
                 (el2nd(i,2) .eq. ed2nd(j,2))) .or. &
                 ((el2nd(i,1) .eq. ed2nd(j,2)) .and. &
                  (el2nd(i,2) .eq. ed2nd(j,1)))) then

                  el2ed(i,1) = j


            else if (((el2nd(i,1) .eq. ed2nd(j,1)) .and. &
                      (el2nd(i,3) .eq. ed2nd(j,2))) .or. &
                 ((el2nd(i,1) .eq. ed2nd(j,2)) .and. &
                  (el2nd(i,3) .eq. ed2nd(j,1)))) then
                  el2ed(i,2) = j


            else if (((el2nd(i,1) .eq. ed2nd(j,1)) .and. &
                      (el2nd(i,4) .eq. ed2nd(j,2))) .or. &
                 ((el2nd(i,1) .eq. ed2nd(j,2)) .and. &
                  (el2nd(i,4) .eq. ed2nd(j,1)))) then
                  
                  el2ed(i,3) = j 


            else if (((el2nd(i,2) .eq. ed2nd(j,1)) .and. &
                      (el2nd(i,3) .eq. ed2nd(j,2))) .or. &
                 ((el2nd(i,2) .eq. ed2nd(j,2)) .and. &
                  (el2nd(i,3) .eq. ed2nd(j,1)))) then
                  
                  el2ed(i,4) = j


            else if (((el2nd(i,4) .eq. ed2nd(j,1)) .and. &
                      (el2nd(i,2) .eq. ed2nd(j,2))) .or. &
                 ((el2nd(i,4) .eq. ed2nd(j,2)) .and. &
                  (el2nd(i,2) .eq. ed2nd(j,1)))) then

                  el2ed(i,5) = j


            else if (((el2nd(i,3) .eq. ed2nd(j,1)) .and. &
                      (el2nd(i,4) .eq. ed2nd(j,2))) .or. &
                 ((el2nd(i,3) .eq. ed2nd(j,2)) .and. &
                  (el2nd(i,4) .eq. ed2nd(j,1)))) then

                  el2ed(i,6) = j
            end if

         end do edge_loop

      end do element_loop
    !$OMP END PARALLEL DO

  end subroutine calc_connect_matrix

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for calculating coordinate arrays
  !---------------------------------------------------------------------
  subroutine calc_coord_matrices (M, nd, el2nd, x, y, z, &
       x_start, x_end, y_start, y_end, z_start, z_end)

    ! INPUT
    ! total number of elements
    integer, intent(in) :: M
    real(kind=dp), dimension(:,:), intent(in) :: nd ! (1:N,1:3)
    integer, dimension(:,:), intent(in) :: el2nd ! (1:M,1:4)

    ! OUTPUT
    real(kind=dp), allocatable, dimension (:,:), intent(out) :: x,y,z
    real(kind=dp), allocatable, dimension (:,:), intent(out) :: &
                                                        x_start,x_end, &
                                                        y_start,y_end, &
                                                        z_start,z_end

    ! LOCAL variables
    integer :: i
    integer :: allo_stat

    !-------------------------------------------------------------------
    allocate (x(M,4),y(M,4),z(M,4), stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
                    "calc_coord_matrices: error allocating array x,y,z")

    allocate (x_start(M,6),x_end(M,6), &
              y_start(M,6),y_end(M,6), &
              z_start(M,6),z_end(M,6), &
              stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
             "calc_coord_matrices: error allocating array x_start etc.")
       
    ! loop over nodes of all elements
    do i = 1,4 
       x(:,i) = nd(el2nd(:,i),1) 
       y(:,i) = nd(el2nd(:,i),2) 
       z(:,i) = nd(el2nd(:,i),3) 
    end do


    x_start(:,1) = x(:,1)
    x_start(:,2) = x(:,1)
    x_start(:,3) = x(:,1)
    x_start(:,4) = x(:,2)
    x_start(:,5) = x(:,4)
    x_start(:,6) = x(:,3)

    x_end(:,1) = x(:,2)
    x_end(:,2) = x(:,3)
    x_end(:,3) = x(:,4)
    x_end(:,4) = x(:,3)
    x_end(:,5) = x(:,2)
    x_end(:,6) = x(:,4)

    y_start(:,1) = y(:,1)
    y_start(:,2) = y(:,1)
    y_start(:,3) = y(:,1)
    y_start(:,4) = y(:,2)
    y_start(:,5) = y(:,4)
    y_start(:,6) = y(:,3)

    y_end(:,1) = y(:,2)
    y_end(:,2) = y(:,3)
    y_end(:,3) = y(:,4)
    y_end(:,4) = y(:,3)
    y_end(:,5) = y(:,2)
    y_end(:,6) = y(:,4)

    z_start(:,1) = z(:,1)
    z_start(:,2) = z(:,1)
    z_start(:,3) = z(:,1)
    z_start(:,4) = z(:,2)
    z_start(:,5) = z(:,4)
    z_start(:,6) = z(:,3)

    z_end(:,1) = z(:,2)
    z_end(:,2) = z(:,3)
    z_end(:,3) = z(:,4)
    z_end(:,4) = z(:,3)
    z_end(:,5) = z(:,2)
    z_end(:,6) = z(:,4)

  end subroutine calc_coord_matrices

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine to calculate matrices for signs that are assigned to 
  !> the vector shape
  !> functions of the edges of all elements to compensate for local and
  !> global directions
  !---------------------------------------------------------------------
  subroutine calc_edge_signs(M, el2nd, ed_sign)

    ! INPUT
    ! total number of elements
    integer, intent(in) :: M
    !integer, dimension(:,:), intent(in) :: ed2nd ! (1:E,1:2)
    !integer, dimension(:,:), intent(in) :: el2ed ! (1:M,1:6)
    integer, dimension(:,:), intent(in) :: el2nd ! (1:M,1:4)

    ! OUTPUT
    ! array of local edge signs
    real(kind=dp), allocatable, dimension(:,:), intent(out) :: ed_sign

    ! LOCAL variables
    integer :: allo_stat

    allocate (ed_sign(M,6), stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
                      "calc_edge_signd: error allocating array ed_sign")

    ! PR calculate edge signs as in matlab routine signs_edges.m of 
    ! Anjam et al (2015).

    ed_sign(:,1) = el2nd(:,1)-el2nd(:,2)
    ed_sign(:,1) = ed_sign(:,1)/abs(ed_sign(:,1))

    ed_sign(:,2) = el2nd(:,1)-el2nd(:,3)
    ed_sign(:,2) = ed_sign(:,2)/abs(ed_sign(:,2))

    ed_sign(:,3) = el2nd(:,1)-el2nd(:,4)
    ed_sign(:,3) = ed_sign(:,3)/abs(ed_sign(:,3))

    ed_sign(:,4) = el2nd(:,2)-el2nd(:,3)
    ed_sign(:,4) = ed_sign(:,4)/abs(ed_sign(:,4))

    ed_sign(:,5) = el2nd(:,4)-el2nd(:,2)
    ed_sign(:,5) = ed_sign(:,5)/abs(ed_sign(:,5))

    ed_sign(:,6) = el2nd(:,3)-el2nd(:,4)
    ed_sign(:,6) = ed_sign(:,6)/abs(ed_sign(:,6))

  end subroutine calc_edge_signs

  !---------------------------------------------------------------------
end module calculate_matrices
