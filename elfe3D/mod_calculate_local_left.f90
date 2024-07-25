!> @brief
!> Module of elfe3D containing subroutines to calculate local parts of
!> system matrix
!!
!> written by Paula Rulff, 28/08/2018
!!
!> Copyright (C) Paula Rulff 2020
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

module calculate_local_left

  use mod_util

  implicit none

contains

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for calculating stiffness matrix Kij
  !---------------------------------------------------------------------

  subroutine calc_stiffness_matrix (l, el2edl, Ve, ed_sign, &
       b_start, b_end, c_start, c_end, d_start, d_end, K)

    ! INPUT
    ! element index
    integer, intent(in) :: l
    real(kind=dp), dimension(:,:), intent(in) :: el2edl ! (M,6)
    real(kind=dp), dimension(:), intent(in) :: Ve ! (M)
    ! array of local edge signs
    real(kind=dp), dimension(:,:), intent(in) :: ed_sign
    real(kind=dp), dimension(:,:), intent(in) :: b_start, b_end, &
                                                 c_start, c_end, &
                                                 d_start, d_end ! (M,6)

    ! OUTPUT
    real(kind=dp), dimension(:,:), intent(inout) :: K

    ! LOCAL variables
    integer :: i,j

    !-------------------------------------------------------------------

    K = 0.0_dp

    do concurrent (j = 1:6, i = 1:6)

          K(i,j) = ((4.0_dp*el2edl(l,i)*el2edl(l,j)*ve(l))/&
                    ((6.0_dp*ve(l))**4.0_dp))*&
               (&
               (c_start(l,i)*d_end(l,i)-d_start(l,i)*c_end(l,i))*&
               (c_start(l,j)*d_end(l,j)-d_start(l,j)*c_end(l,j))+&
               
               (d_start(l,i)*b_end(l,i)-b_start(l,i)*d_end(l,i))*&
               (d_start(l,j)*b_end(l,j)-b_start(l,j)*d_end(l,j))+&
               
               (b_start(l,i)*c_end(l,i)-c_start(l,i)*b_end(l,i))*&
               (b_start(l,j)*c_end(l,j)-c_start(l,j)*b_end(l,j))&
               ) * &
               ed_sign(l,i)*ed_sign(l,j)
    end do


  end subroutine calc_stiffness_matrix

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for calculating mass matrix Mij
  !---------------------------------------------------------------------

  subroutine calc_mass_matrix (l, el2edl, Ve, ed_sign, &
       b, c, d, MM)

    ! INPUT
    ! element index
    integer, intent(in) :: l
    real(kind=dp), dimension(:,:), intent(in) :: el2edl ! (M,6)
    real(kind=dp), dimension(:), intent(in) :: ve ! (M)
    real(kind=dp), dimension(:,:), intent(in) :: b,c,d ! (M,4)
    ! array of local edge signs
    real(kind=dp), dimension(:,:), intent(in) :: ed_sign

    ! OUTPUT
    real(kind=dp), dimension(:,:), intent(inout) :: MM

    ! LOCAL variables
    integer :: allo_stat
    integer :: i,j
    real(kind=dp), allocatable, dimension(:,:) :: f

    !-------------------------------------------------------------------
    ! calculating fij
    allocate (f(4,4), stat = allo_stat)
    call allocheck(log_unit, allo_stat, "error allocating array f")

    f = 0.0_dp

    do j = 1,4
       do i = 1,4

          f(i,j) = b(l,i)*b(l,j) + c(l,i)*c(l,j) + d(l,i)*d(l,j)

       end do
    end do

    ! calculating MM_ij for one element
    MM = 0.0_dp

    MM(1,1) = (((el2edl(l,1))*(el2edl(l,1)))/ &
                (360.0_dp*Ve(l)))*(f(2,2)-f(1,2)+f(1,1))

    MM(1,2) = (((el2edl(l,1))*(el2edl(l,2)))/&
                (720.0_dp*Ve(l)))*(2*f(2,3)-f(2,1)-f(1,3)+f(1,1))

    MM(1,3) = (((el2edl(l,1))*(el2edl(l,3)))/ &
                (720.0_dp*Ve(l)))*(2*f(2,4)-f(2,1)-f(1,4)+f(1,1))

    MM(1,4) = (((el2edl(l,1))*(el2edl(l,4)))/ &
                (720.0_dp*Ve(l)))*(f(2,3)-f(2,2)-2*f(1,3)+f(1,2))

    MM(1,5) = (((el2edl(l,1))*(el2edl(l,5)))/ &
                (720.0_dp*Ve(l)))*(f(2,2)-f(2,4)-f(1,2)+2*f(1,4))

    MM(1,6) = (((el2edl(l,1))*(el2edl(l,6)))/ &
                (720.0_dp*Ve(l)))*(f(2,4)-f(2,3)-f(1,4)+f(1,3))

    MM(2,2) = (((el2edl(l,2))*(el2edl(l,2)))/ &
                (360.0_dp*Ve(l)))*(f(3,3)-f(1,3)+f(1,1))
 !PR CORRECTION f(3,1) INSTEAD OF f(1,3) - wrong in Jin 
    MM(2,3) = (((el2edl(l,2))*(el2edl(l,3)))/ &
                (720.0_dp*Ve(l)))*(2*f(3,4)-f(3,1)-f(1,4)+f(1,1))

    MM(2,4) = (((el2edl(l,2))*(el2edl(l,4)))/ &
                (720.0_dp*Ve(l)))*(f(3,3)-f(2,3)-f(1,3)+2*f(1,2))

    MM(2,5) = (((el2edl(l,2))*(el2edl(l,5)))/ &
                (720.0_dp*Ve(l)))*(f(2,3)-f(3,4)-f(1,2)+f(1,4))

    MM(2,6) = (((el2edl(l,2))*(el2edl(l,6)))/ &
                (720.0_dp*Ve(l)))*(f(1,3)-f(3,3)-2*f(1,4)+f(3,4))

    MM(3,3) = (((el2edl(l,3))*(el2edl(l,3)))/ &
                (360.0_dp*Ve(l)))*(f(4,4)-f(1,4)+f(1,1))

    MM(3,4) = (((el2edl(l,3))*(el2edl(l,4)))/ &
                (720.0_dp*Ve(l)))*(f(3,4)-f(2,4)-f(1,3)+f(1,2))

    MM(3,5) = (((el2edl(l,3))*(el2edl(l,5)))/ &
                (720.0_dp*Ve(l)))*(f(2,4)-f(4,4)-2*f(1,2)+f(1,4))

    MM(3,6) = (((el2edl(l,3))*(el2edl(l,6)))/ &
                (720.0_dp*Ve(l)))*(f(4,4)-f(3,4)-f(1,4)+2*f(1,3))

    MM(4,4) = (((el2edl(l,4))*(el2edl(l,4)))/ &
                (360.0_dp*Ve(l)))*(f(3,3)-f(2,3)+f(2,2))

    MM(4,5) = (((el2edl(l,4))*(el2edl(l,5)))/ &
                (720.0_dp*Ve(l)))*(f(2,3)-2*f(3,4)-f(2,2)+f(2,4))

    MM(4,6) = (((el2edl(l,4))*(el2edl(l,6)))/ &
                (720.0_dp*Ve(l)))*(f(3,4)-f(3,3)-2*f(2,4)+f(2,3))

    MM(5,5) = (((el2edl(l,5))*(el2edl(l,5)))/ &
                (360.0_dp*Ve(l)))*(f(2,2)-f(2,4)+f(4,4))

    MM(5,6) = (((el2edl(l,5))*(el2edl(l,6)))/ &
                (720.0_dp*Ve(l)))*(f(2,4)-2*f(2,3)-f(4,4)+f(3,4))

    MM(6,6) = (((el2edl(l,6))*(el2edl(l,6)))/ &
                (360.0_dp*Ve(l)))*(f(4,4)-f(3,4)+f(3,3))

    MM(2,1) = MM(1,2)
    MM(3,1) = MM(1,3)
    MM(3,2) = MM(2,3)
    MM(4,1) = MM(1,4)
    MM(4,2) = MM(2,4)
    MM(4,3) = MM(3,4)
    MM(5,1) = MM(1,5)
    MM(5,2) = MM(2,5)
    MM(5,3) = MM(3,5)
    MM(5,4) = MM(4,5)
    MM(6,1) = MM(1,6)
    MM(6,2) = MM(2,6)
    MM(6,3) = MM(3,6)
    MM(6,4) = MM(4,6)
    MM(6,5) = MM(5,6)

    do concurrent (i = 1:6, j = 1:6)
       MM(i,j) = MM(i,j) * ed_sign(l,i)*ed_sign(l,j)
    end do

   deallocate(f)

  end subroutine calc_mass_matrix

  !---------------------------------------------------------------------
end module calculate_local_left
