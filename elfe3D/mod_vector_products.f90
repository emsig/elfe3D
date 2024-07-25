!> @brief
!> Module of elfe3d containing functions to calculate vector products
!!
!> written by Paula Rulff, 29/01/2020
!!
!> Last change: March 2024
!!
!> Copyright (C) Paula Rulff 2020
!!
!> original functions (for real vectors) can be found at
!> https://rosettacode.org/wiki/Vector_products#Fortran
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

module vector_products

  use mod_util

  implicit none

  contains

  !---------------------------------------------------------------------
  !> @brief
  !> function for calculating cross product of two real vectors [3]
  ! A x B = (a2b3  -   a3b2,     a3b1   -   a1b3,     a1b2   -   a2b1) 
  !---------------------------------------------------------------------
    function cross_product_real(a, b)
      real(kind=dp), dimension(3) :: cross_product_real
      real(kind=dp), dimension(3), intent(in) :: a, b
   
      cross_product_real(1) = a(2)*b(3) - a(3)*b(2)
      cross_product_real(2) = a(3)*b(1) - a(1)*b(3)
      cross_product_real(3) = a(1)*b(2) - b(1)*a(2)
    end function cross_product_real


  !---------------------------------------------------------------------
  !> @brief
  !> function for calculating cross product of two complex vectors [3]
  ! A x B = (a2b3  -   a3b2,     a3b1   -   a1b3,     a1b2   -   a2b1) 
  !---------------------------------------------------------------------
    function cross_product(a, b)
      complex(kind=dp), dimension(3) :: cross_product
      complex(kind=dp), dimension(3), intent(in) :: a, b
   
      cross_product(1) = a(2)*b(3) - a(3)*b(2)
      cross_product(2) = a(3)*b(1) - a(1)*b(3)
      cross_product(3) = a(1)*b(2) - b(1)*a(2)
    end function cross_product

  end module vector_products