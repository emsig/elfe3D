!> @brief
!> Module of elfe3D containing subroutines to calculate the field
!> components
!!
!> written by Paula Rulff, 24/06/2019
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
module calculate_tf

  use mod_util
  use define_model

  implicit none

contains

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for calculating electric and magnetic fields at a
  ! certain receiver location
  !---------------------------------------------------------------------
  subroutine calculate_fields (u1, v1, w1, rec1_el, el2ed, S, &
                               a_start, a_end, b_start, b_end, &
                               c_start, c_end, d_start, d_end, &
                               el2edl, ed_sign, Ve, w, mu, &
                               E_rec1, H_rec1)
  
    ! INPUT
    ! Cartesian coordinates of current receiver
    real(kind=dp), intent(in) :: u1,v1,w1 
    integer, intent(in) :: rec1_el ! current receiver element
    integer, dimension(:,:), intent(in) :: el2ed
    complex (kind=dp), dimension(:), intent(in) :: S ! solution vector
    real(kind=dp), dimension(:,:), intent(in) :: a_start, a_end, &
                                                 b_start, b_end, &
                                                 c_start, c_end, &
                                                 d_start, d_end
    real(kind=dp), dimension(:,:), intent(in) :: el2edl
    real(kind=dp), dimension(:,:), intent(in) :: ed_sign
    real(kind=dp), dimension(:), intent(in) :: Ve
    real(kind=dp), intent(in) :: w
    real(kind=dp), dimension(:), intent(in) :: mu

    ! OUTPUT
    complex(kind=dp), dimension(3), intent(out) :: E_rec1, H_rec1
    

    ! LOCAL variables
    integer :: l
    ! edges of the element containing receiver
    integer, dimension(6) :: rec1_ed 
    real(kind=dp), dimension(3) :: grad_Lstart ! grad Lstart
    real(kind=dp), dimension(3) :: grad_Lend ! grad Lend
    ! Nedelec basis functions of edges containing receiver (vector)
    real(kind=dp), dimension(3) :: N_rec1 
    complex(kind=dp) :: factor   ! factor -1/(iwmu)
    ! curl of Nedelec basis function for one edge
    real(kind=dp), dimension(3) :: curl_N !
    !-------------------------------------------------------------------
    
     ! Find edges of receiver earth element
     rec1_ed = el2ed(rec1_el,:)

     ! Calculate electric fields:
     ! call Write_Message (log_unit, &
     !                    'Edge numbers of receiver earth element are:')
     ! do i = 1,6
     !    write(*,*) rec1_ed(i)
     ! end do

     E_rec1 = (/(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp)/)
     N_rec1 = (/(0.0_dp),(0.0_dp),(0.0_dp)/)

     do l = 1,6

        ! Calculate grad Lstart and grad Lend vectors
        grad_Lstart = (/ b_start(rec1_el,l), &
                         c_start(rec1_el,l), &
                         d_start(rec1_el,l) /)

        grad_Lend = (/ b_end(rec1_el,l), &
                       c_end(rec1_el,l), &
                       d_end(rec1_el,l) /)

        ! Calculate Nedelec basis function for edge l
        N_rec1 = (((1.0_dp/(6.0_dp*Ve(rec1_el)))**2.0_dp)* &
                  (a_start(rec1_el,l) + b_start(rec1_el,l)*u1 &
                                      + c_start(rec1_el,l)*v1 &
                                      + d_start(rec1_el,l)*w1) &
                  *grad_Lend &
             
             - &
             
             ((1.0_dp/(6.0_dp*Ve(rec1_el)))**2.0_dp)* &
              (a_end(rec1_el,l) + b_end(rec1_el,l)*u1 &
                                + c_end(rec1_el,l)*v1 &
                                + d_end(rec1_el,l)*w1)*grad_Lstart) &
             
             * el2edl(rec1_el,l) &

             * ed_sign(rec1_el,l)


        ! Calculate E_rec1 contatining Ex, Ey and Ez 
        ! as the sum over all N times S
        ! of all edges of the element contatining the receiver
        E_rec1 =  E_rec1 + cmplx(N_rec1,kind=dp)*S(rec1_ed(l))

     end do

     ! Calculate magnetic fields:
     H_rec1 = (/(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp)/)
     curl_N = (/(0.0_dp),(0.0_dp),(0.0_dp)/)
     factor = cmplx(0.0,(1/(w*mu(rec1_el))), kind=dp)

     do l = 1,6

        ! Calculate curl_N for every edge of receiver element
        curl_N = ((2.0_dp*el2edl(rec1_el,l))/ &
                  (6.0_dp*Ve(rec1_el))**2.0_dp)* &
                 (/(c_start(rec1_el,l)*d_end(rec1_el,l) &
                  - d_start(rec1_el,l)*c_end(rec1_el,l)), &
                   (d_start(rec1_el,l)*b_end(rec1_el,l) &
                  - b_start(rec1_el,l)*d_end(rec1_el,l)),&
                   (b_start(rec1_el,l)*c_end(rec1_el,l) &
                  - c_start(rec1_el,l)*b_end(rec1_el,l))/)&

                  * ed_sign(rec1_el,l)

        ! Calculate H_rec1 = -1/(iwmu)sum(Ej curl (Nj) and 
        ! add for all 6 edges
        H_rec1 =  H_rec1 + cmplx(S(rec1_ed(l))*curl_N, kind=dp)

     end do

     H_rec1 =  H_rec1*factor

   !--------------------------------------------------------------------
   end subroutine calculate_fields

  end module calculate_tf
