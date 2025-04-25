!> @brief
!> Module of elfe3d containing subroutines to to calculate coefficients 
!> for interpolation functions and element characteristics
!!
!> written by Paula Rulff, 27/08/2018
!!
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

module interp_functions

  use mod_util

  implicit none

contains

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for calculating coefficients for linear interpolation 
  !> functions beween nodes
  !> and abcd_start & abcd_end matrices M times 6 each
  !---------------------------------------------------------------------
  subroutine calc_abcd (M, x, y, z, a, b, c, d, &
       a_start, a_end, b_start, b_end, c_start, c_end, d_start, d_end)

    ! INPUT
    ! total number of elements
    integer, intent(in) :: M
    ! coordinates of 4 nodes for all M elements
    real(kind=dp), dimension(:,:), intent(in) :: x,y,z ! (1:M,1:4)

    ! OUTPUT
    ! coefficients for linear shape functions determined by 4 nodes of 
    ! all M elements
    real(kind=dp), allocatable, dimension(:,:), intent(out) :: a,b,c,d
    real(kind=dp), allocatable, dimension(:,:), intent(out) :: &
          a_start, a_end, b_start, b_end, c_start, c_end, d_start, d_end

    ! LOCAL variables
    real(kind=dp), allocatable, dimension(:) :: a1, a2, a3, a4, &
                                                b1, b2, b3, b4, &
                                                c1, c2, c3, c4, &
                                                d1, d2, d3, d4
    integer :: allo_stat
    !-------------------------------------------------------------------
    allocate (a(M,4), b(M,4), c(M,4), d(M,4), stat = allo_stat)
    call allocheck &
    (log_unit, allo_stat, "calc_abcd: error allocating arrays a,b,c,d" )

    allocate (a_start(M,6), a_end(M,6), b_start(M,6), b_end(M,6), &
              c_start(M,6), c_end(M,6), d_start(M,6), d_end(M,6), &
              stat = allo_stat)
    call allocheck &
      (log_unit, allo_stat, "calc_abcd: error allocating array a_start")

    ! allocate local arrays
    allocate (a1(M), a2(M), a3(M), a4(M), &
              b1(M), b2(M), b3(M), b4(M), &
              c1(M), c2(M), c3(M), c4(M), &
              d1(M), d2(M), d3(M), d4(M), stat = allo_stat)
    call allocheck &
    (log_unit, allo_stat, "calc_abcd: error allocating local arrays")

    
    a1 = x(:,2)*y(:,3)*z(:,4)-x(:,2)*y(:,4)*z(:,3)+x(:,3)*y(:,4)*z(:,2) &
        -x(:,3)*y(:,2)*z(:,4)+x(:,4)*y(:,2)*z(:,3)-x(:,4)*y(:,3)*z(:,2)
    a2 = -x(:,1)*y(:,3)*z(:,4)+x(:,1)*y(:,4)*z(:,3)-x(:,3)*y(:,4)*z(:,1) &
         +x(:,3)*y(:,1)*z(:,4)-x(:,4)*y(:,1)*z(:,3)+x(:,4)*y(:,3)*z(:,1)
    a3 = x(:,1)*y(:,2)*z(:,4)-x(:,1)*y(:,4)*z(:,2)+x(:,2)*y(:,4)*z(:,1) &
        -x(:,2)*y(:,1)*z(:,4)+x(:,4)*y(:,1)*z(:,2)-x(:,4)*y(:,2)*z(:,1)
    a4 = -x(:,1)*y(:,2)*z(:,3)+x(:,1)*y(:,3)*z(:,2)-x(:,2)*y(:,3)*z(:,1) &
         +x(:,2)*y(:,1)*z(:,3)-x(:,3)*y(:,1)*z(:,2)+x(:,3)*y(:,2)*z(:,1)


    b1 = -y(:,3)*z(:,4) + y(:,4)*z(:,3) &
        - y(:,4)*z(:,2) + y(:,2)*z(:,4) &
        - y(:,2)*z(:,3) + y(:,3)*z(:,2)
    b2 = y(:,3)*z(:,4) - y(:,4)*z(:,3) &
       + y(:,4)*z(:,1) - y(:,1)*z(:,4) &
       + y(:,1)*z(:,3) - y(:,3)*z(:,1)
    b3 = -y(:,2)*z(:,4) + y(:,4)*z(:,2) &
        - y(:,4)*z(:,1) + y(:,1)*z(:,4) &
        - y(:,1)*z(:,2) + y(:,2)*z(:,1)
    b4 = y(:,2)*z(:,3) - y(:,3)*z(:,2) &
       + y(:,3)*z(:,1) - y(:,1)*z(:,3) &
       + y(:,1)*z(:,2) - y(:,2)*z(:,1)

    c1 = -(-x(:,3)*z(:,4) + x(:,4)*z(:,3) - x(:,4)*z(:,2) &
          + x(:,2)*z(:,4) - x(:,2)*z(:,3) + x(:,3)*z(:,2))
    c2 = -(x(:,3)*z(:,4) - x(:,4)*z(:,3) + x(:,4)*z(:,1) &
         - x(:,1)*z(:,4) + x(:,1)*z(:,3) - x(:,3)*z(:,1))
    c3 = -(-x(:,2)*z(:,4) + x(:,4)*z(:,2) - x(:,4)*z(:,1) &
          + x(:,1)*z(:,4) - x(:,1)*z(:,2) + x(:,2)*z(:,1))
    c4 = -(x(:,2)*z(:,3) - x(:,3)*z(:,2) + x(:,3)*z(:,1) &
         - x(:,1)*z(:,3) + x(:,1)*z(:,2) - x(:,2)*z(:,1))

    d1 = -x(:,3)*y(:,4) + x(:,4)*y(:,3) - x(:,4)*y(:,2) &
        + x(:,2)*y(:,4) - x(:,2)*y(:,3) + x(:,3)*y(:,2)
    d2 = x(:,3)*y(:,4) - x(:,4)*y(:,3) + x(:,4)*y(:,1) &
       - x(:,1)*y(:,4) + x(:,1)*y(:,3) - x(:,3)*y(:,1)
    d3 = -x(:,2)*y(:,4) + x(:,4)*y(:,2) - x(:,4)*y(:,1) &
        + x(:,1)*y(:,4) - x(:,1)*y(:,2) + x(:,2)*y(:,1)
    d4 = x(:,2)*y(:,3) - x(:,3)*y(:,2) + x(:,3)*y(:,1) &
       - x(:,1)*y(:,3) + x(:,1)*y(:,2) - x(:,2)*y(:,1)

    a(:,1) = a1
    a(:,2) = a2
    a(:,3) = a3
    a(:,4) = a4
    b(:,1) = b1
    b(:,2) = b2
    b(:,3) = b3
    b(:,4) = b4
    c(:,1) = c1
    c(:,2) = c2
    c(:,3) = c3
    c(:,4) = c4
    d(:,1) = d1
    d(:,2) = d2
    d(:,3) = d3
    d(:,4) = d4

    a_start(:,1) = a1
    a_start(:,2) = a1
    a_start(:,3) = a1
    a_start(:,4) = a2
    a_start(:,5) = a4
    a_start(:,6) = a3

    a_end(:,1) = a2
    a_end(:,2) = a3
    a_end(:,3) = a4
    a_end(:,4) = a3
    a_end(:,5) = a2
    a_end(:,6) = a4

    b_start(:,1) = b1
    b_start(:,2) = b1
    b_start(:,3) = b1
    b_start(:,4) = b2
    b_start(:,5) = b4
    b_start(:,6) = b3

    b_end(:,1) = b2
    b_end(:,2) = b3
    b_end(:,3) = b4
    b_end(:,4) = b3
    b_end(:,5) = b2
    b_end(:,6) = b4

    c_start(:,1) = c1
    c_start(:,2) = c1
    c_start(:,3) = c1
    c_start(:,4) = c2
    c_start(:,5) = c4
    c_start(:,6) = c3

    c_end(:,1) = c2
    c_end(:,2) = c3
    c_end(:,3) = c4
    c_end(:,4) = c3
    c_end(:,5) = c2
    c_end(:,6) = c4

    d_start(:,1) = d1
    d_start(:,2) = d1
    d_start(:,3) = d1
    d_start(:,4) = d2
    d_start(:,5) = d4
    d_start(:,6) = d3

    d_end(:,1) = d2
    d_end(:,2) = d3
    d_end(:,3) = d4
    d_end(:,4) = d3
    d_end(:,5) = d2
    d_end(:,6) = d4

    ! deallocate local arrays
    if (allocated(a1)) deallocate(a1, a2, a3, a4, b1, b2, b3, b4, &
         c1, c2, c3, c4, d1, d2, d3, d4)
    
  end subroutine calc_abcd

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine to calculate the volume of each element
  !---------------------------------------------------------------------
  subroutine calc_vol (M, x, y, z, ve)

    ! INPUT
    ! total number of elements
    integer, intent(in) :: M
    real(kind=dp), dimension (:,:), intent(in) :: x, y, z ! (1:M,1:4)

    ! OUTPUT
    real(kind=dp), allocatable, dimension(:), intent(out) :: ve

    ! LOCAL variables
    integer :: allo_stat
    !-------------------------------------------------------------------
    allocate (ve(M), stat = allo_stat)
    call allocheck &
            (log_unit, allo_stat, "calc_vol: error allocating array ve")

    ve = abs( (x(:,1)-x(:,4)) * ((y(:,2)-y(:,4))*(z(:,3)-z(:,4)) &
            - (y(:,3)-y(:,4))*(z(:,2)-z(:,4))) &
            - (x(:,2)-x(:,4)) * ((y(:,1)-y(:,4))*(z(:,3)-z(:,4)) &
            - (y(:,3)-y(:,4))*(z(:,1)-z(:,4))) &
            + (x(:,3)-x(:,4)) * ((y(:,1)-y(:,4))*(z(:,2)-z(:,4)) &
            - (y(:,2)-y(:,4))*(z(:,1)-z(:,4))) ) / 6.0_dp

  end subroutine calc_vol

  !---------------------------------------------------------------------
  !> @brief
  !> Subroutine to calculate the edge lengths
  !---------------------------------------------------------------------
  subroutine calc_edge_length (M, E, nd, ed2nd, x, y, z, el2edl, edl)

    ! INPUT
    ! total number of elements
    integer, intent(in) :: M,E
    real(kind=dp), dimension(:,:), intent(in) :: x, y, z ! (1:M,1:4)
    real(kind=dp), dimension(:,:), intent(in) :: nd 
    integer, dimension(:,:), intent(in) :: ed2nd

    ! OUTPUT
    real(kind=dp), allocatable, dimension(:,:), intent(out) :: el2edl
    real(kind=dp), allocatable, dimension(:), intent(out) :: edl

    ! LOCAL variables
    integer :: allo_stat
    integer :: i
    !-------------------------------------------------------------------
    allocate (el2edl(M,6), stat = allo_stat)
    call allocheck &
    (log_unit, allo_stat, "calc_edge_length: error allocating el2edl")

    el2edl = 0.0_dp

    do i = 1,M

       el2edl(i,1) = sqrt((x(i,1)-x(i,2))**2 &
                        + (y(i,1)-y(i,2))**2 &
                        + (z(i,1)-z(i,2))**2)
       el2edl(i,2) = sqrt((x(i,1)-x(i,3))**2 &
                        + (y(i,1)-y(i,3))**2 &
                        + (z(i,1)-z(i,3))**2)
       el2edl(i,3) = sqrt((x(i,1)-x(i,4))**2 &
                        + (y(i,1)-y(i,4))**2 &
                        + (z(i,1)-z(i,4))**2)
       el2edl(i,4) = sqrt((x(i,2)-x(i,3))**2 &
                        + (y(i,2)-y(i,3))**2 &
                        + (z(i,2)-z(i,3))**2)
       el2edl(i,5) = sqrt((x(i,4)-x(i,2))**2 &
                        + (y(i,4)-y(i,2))**2 &
                        + (z(i,4)-z(i,2))**2)
       el2edl(i,6) = sqrt((x(i,3)-x(i,4))**2 &
                        + (y(i,3)-y(i,4))**2 &
                        + (z(i,3)-z(i,4))**2)

    end do

    allocate (edl(E), stat = allo_stat)
    call allocheck &
    (log_unit, allo_stat, "calc_edge_length: error allocating edl")

    edl = 0.0_dp

    do i = 1,E

      edl(i) = sqrt((nd(ed2nd(i,1),1) - nd(ed2nd(i,2),1))**2 &
                  + (nd(ed2nd(i,1),2) - nd(ed2nd(i,2),2))**2 &
                  + (nd(ed2nd(i,1),3) - nd(ed2nd(i,2),3))**2)

    end do

  end subroutine calc_edge_length
  !---------------------------------------------------------------------

end module interp_functions
