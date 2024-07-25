!> @brief
!> Module of elfe3d containing subroutines for error estimation
!!
!> written by Paula Rulff & Laura Maria Buntin, 25/10/2019
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

module error_estimates

  use mod_util
  use define_model
  use vector_products

  implicit none

contains

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for calculating elemental residuals res_s, res_w
  !---------------------------------------------------------------------
  subroutine calculate_elemental_residuals (M, w, el2ed, &
                                            S, Wg, Cg, Bg, &
                                            rho, Ve, el2edl, ed_sign, &
                                            a_start, a_end, &
                                            b_start, b_end, &
                                            c_start, c_end, &
                                            d_start, d_end, &
                                            el2nd, nd, res_s, res_w)
  
    ! INPUT
    integer, intent(in) :: M
    real(kind=dp), intent(in) :: w
    integer, dimension(:,:), intent(in) :: el2ed, el2nd
    complex (kind=dp), dimension(:), intent(in) :: S, Wg, Cg, Bg
    real(kind=dp), dimension(:), intent(in) :: rho, Ve
    real(kind=dp), dimension(:,:), intent(in) :: el2edl, ed_sign, nd
    real(kind=dp), dimension(:,:), intent(in) :: a_start, a_end, &
                                                 b_start, b_end, &
                                                 c_start, c_end, &
                                                 d_start, d_end

    ! OUTPUT
    real(kind=dp), allocatable, dimension(:), intent(out) :: res_s, &
                                                             res_w
    
    
    ! LOCAL variables
    integer :: mi, i, e
    integer :: allo_stat
    real(kind=dp) :: abs_res_s, abs_res_w
    complex (kind=dp) :: factor, ae, aw
    real(kind=dp), dimension(3) :: grad_Lstart
    real(kind=dp), dimension(3) :: grad_Lend
    real(kind=dp), dimension(3) :: N, rS, rW
    real(kind=dp), dimension(3) :: Gausspoint
    complex(kind=dp), dimension(3) :: Gauss_S, Gauss_Wg, &
                                      Gauss_Bg, Gauss_Cg, &
                                      Gauss_ae, Gauss_aw
    ! for Gauss order 2 tests
    integer :: GaussOrder, ngp

   !--------------------------------------------------------------------
   ! allocate space for res_s and res_w: size = number of elements
   ! will be deallocated in main program
   allocate (res_s(M),stat = allo_stat)
     call allocheck(log_unit, allo_stat, &
      "calculate_elemental_residuals: error allocating residual")
      
   allocate (res_w(M),stat = allo_stat)
     call allocheck(log_unit, allo_stat, &
      "calculate_elemental_residuals: error allocating residual")
   
   ! Initialise res_w and res_s to zero
   res_w = 0.0_dp
   res_s = 0.0_dp

   GaussOrder = 2
   

   !$OMP PARALLEL DO
   ! Loop over all elements
    do mi = 1,M

     select case(GaussOrder)
     case(1)
       !print*, 'Gauss Order 1'
       ! calculate Gauss point of current element
       Gausspoint = nd(el2nd(mi,1),:) * 0.25_dp + nd(el2nd(mi,2),:) &
                                      * 0.25_dp + nd(el2nd(mi,3),:) &
                                      * 0.25_dp + nd(el2nd(mi,4),:) &
                                      * 0.25_dp
       ! Gauss weight = Ve(mi) for 1st order

       ! factor for residual equation
       factor = cmplx(0.0_dp, (w/rho(mi)), kind=dp) * (-1.0_dp)

       ! initialise
       abs_res_s = 0.0_dp
       abs_res_w = 0.0_dp
       rS = 0.0_dp
       rW = 0.0_dp

       Gauss_S = ZEROW
       Gauss_Wg = ZEROW
       Gauss_Bg = ZEROW
       Gauss_Cg = ZEROW
       Gauss_ae = ZEROW
       Gauss_aw = ZEROW

       ! Loop over edges in element
       do i = 1,6
       
         ! element to edge:
         ! e = global edge ID of edge i of element mi
         e = el2ed(mi,i)
     
         ! S(e) = E field at current edge, Wg(e) = W at current edge
         ! Cg(e) = rhs of dual problem, Bg(e) = rhs of primal problem

         ! weighting facors to get relative-to-local-field-amplitude
         ae = 1.0_dp / (S(e) + cmplx(10E-9, kind=dp))
         aw = 1.0_dp / (Wg(e) + cmplx(10E-9, kind=dp))

         ! Calculate grad Lstart and grad Lend vectors
         grad_Lstart = (/ b_start(mi,i), c_start(mi,i), d_start(mi,i) /)

         grad_Lend = (/ b_end(mi,i), c_end(mi,i), d_end(mi,i) /)

         ! Calculate Nedelec basis function for edge e
         N = (((1/(6*Ve(mi)))**2)*(a_start(mi,i) + &
                                   b_start(mi,i)*Gausspoint(1) + &
                                   c_start(mi,i)*Gausspoint(2) + &
                                   d_start(mi,i)*Gausspoint(3)) &
                                 *grad_Lend &
            
                   - &
            
                   ((1/(6*Ve(mi)))**2)*(a_end(mi,i) + &
                                        b_end(mi,i)*Gausspoint(1) + &
                                        c_end(mi,i)*Gausspoint(2) + &
                                        d_end(mi,i)*Gausspoint(3)) &
                                      *grad_Lstart) &
            
                   * el2edl(mi,i) &

                   * ed_sign(mi,i)

         ! fields, rhs and weights (vectors) at Gauss point
         Gauss_S = Gauss_S + cmplx(N,kind=dp) * S(e)
         Gauss_Wg = Gauss_Wg + cmplx(N,kind=dp) * Wg(e)
         Gauss_Bg = Gauss_Bg + cmplx(N,kind=dp) * Bg(e)
         Gauss_Cg = Gauss_Cg + cmplx(N,kind=dp) * Cg(e)
         Gauss_ae = Gauss_ae + cmplx(N,kind=dp) * ae
         Gauss_aw = Gauss_aw + cmplx(N,kind=dp) * aw


       ! End loop over edges
       end do

       ! residuals at Gauss point
       rS =  abs(Gauss_ae*(factor * Gauss_S - Gauss_Bg))
       rW =  abs(Gauss_aw*(factor * Gauss_Wg - Gauss_Cg))


       ! absolute residuals at Gauss point
       abs_res_s = sqrt((abs(rS(1)))**2.0_dp + &
                        (abs(rS(2)))**2.0_dp + &
                        (abs(rS(3)))**2.0_dp)
       abs_res_w = sqrt((abs(rW(1)))**2.0_dp + &
                        (abs(rW(2)))**2.0_dp + &
                        (abs(rW(3)))**2.0_dp)

       ! squared residual for L2 norm times gauss weight for 
       ! the current element
       res_s(mi) = abs(abs_res_s * abs_res_s)
       res_w(mi) = abs(abs_res_w * abs_res_w)


     case(2)
      !print*, 'Gauss Order 2'
      do ngp = 1,4
       ! calculate Gauss point of current element
       select case (ngp)
       case(1)
        Gausspoint = nd(el2nd(mi,1),:) * 0.5585410197_dp + &
                     nd(el2nd(mi,2),:) * 0.1381966011_dp + &
                     nd(el2nd(mi,3),:) * 0.1381966011_dp + &
                     nd(el2nd(mi,4),:) * 0.1381966011_dp
       case(2)
        Gausspoint = nd(el2nd(mi,1),:) * 0.1381966011_dp + &
                     nd(el2nd(mi,2),:) * 0.5585410197_dp + &
                     nd(el2nd(mi,3),:) * 0.1381966011_dp + &
                     nd(el2nd(mi,4),:) * 0.1381966011_dp
       case(3)
        Gausspoint = nd(el2nd(mi,1),:) * 0.1381966011_dp + &
                     nd(el2nd(mi,2),:) * 0.1381966011_dp + &
                     nd(el2nd(mi,3),:) * 0.5585410197_dp + &
                     nd(el2nd(mi,4),:) * 0.1381966011_dp
       case(4)
        Gausspoint = nd(el2nd(mi,1),:) * 0.1381966011_dp + &
                     nd(el2nd(mi,2),:) * 0.1381966011_dp + &
                     nd(el2nd(mi,3),:) * 0.1381966011_dp + &
                     nd(el2nd(mi,4),:) * 0.5585410197_dp
       end select
       ! Gauss weight = Ve(mi)*0.25 for 2nd order
       ! factor for residual equation
       factor = cmplx(0.0_dp, (w/rho(mi)), kind=dp) * (-1.0_dp)

       ! initialise
       abs_res_s = 0.0_dp
       abs_res_w = 0.0_dp
       rS = 0.0_dp
       rW = 0.0_dp

       Gauss_S = ZEROW
       Gauss_Wg = ZEROW
       Gauss_Bg = ZEROW
       Gauss_Cg = ZEROW
       Gauss_ae = ZEROW
       Gauss_aw = ZEROW

       ! Loop over edges in element
       do i = 1,6
       
         ! element to edge:
         ! e = global edge ID of edge i of element mi
         e = el2ed(mi,i)
       
         ! S(e) = E field at current edge, Wg(e) = W at current edge
         ! Cg(e) = rhs of dual problem, Bg(e) = rhs of primal problem

         ! weighting facors to get relative-to-local-field-amplitude
         ae = 1.0_dp / (S(e) + cmplx(10E-9, kind=dp))
         aw = 1.0_dp / (Wg(e) + cmplx(10E-9, kind=dp))

         ! Calculate grad Lstart and grad Lend vectors
         grad_Lstart = (/ b_start(mi,i), c_start(mi,i), d_start(mi,i) /)

         grad_Lend = (/ b_end(mi,i), c_end(mi,i), d_end(mi,i) /)

         ! Calculate Nedelec basis function for edge e
         N = (((1/(6*Ve(mi)))**2)*(a_start(mi,i) + &
                                   b_start(mi,i)*Gausspoint(1) + &
                                   c_start(mi,i)*Gausspoint(2) + &
                                   d_start(mi,i)*Gausspoint(3)) &
                                 *grad_Lend &
            
                   - &
            
                   ((1/(6*Ve(mi)))**2)*(a_end(mi,i) + &
                                        b_end(mi,i)*Gausspoint(1) + &
                                        c_end(mi,i)*Gausspoint(2) + &
                                        d_end(mi,i)*Gausspoint(3))&
                                      *grad_Lstart) &
            
                   * el2edl(mi,i) &

                   * ed_sign(mi,i)

         ! fields, rhs and weights (vectors) at Gauss point
         Gauss_S = Gauss_S + cmplx(N,kind=dp) * S(e)
         Gauss_Wg = Gauss_Wg + cmplx(N,kind=dp) * Wg(e)
         Gauss_Bg = Gauss_Bg + cmplx(N,kind=dp) * Bg(e)
         Gauss_Cg = Gauss_Cg + cmplx(N,kind=dp) * Cg(e)
         Gauss_ae = Gauss_ae + cmplx(N,kind=dp) * ae
         Gauss_aw = Gauss_aw + cmplx(N,kind=dp) * aw


       ! End loop over edges
       end do

       ! residuals at Gauss point
       rS =  abs(Gauss_ae*(factor * Gauss_S - Gauss_Bg))
       rW =  abs(Gauss_aw*(factor * Gauss_Wg - Gauss_Cg))


       ! absolute residuals at Gauss point
       abs_res_s = sqrt((abs(rS(1)))**2.0_dp + &
                        (abs(rS(2)))**2.0_dp + &
                        (abs(rS(3)))**2.0_dp)
       abs_res_w = sqrt((abs(rW(1)))**2.0_dp + &
                        (abs(rW(2)))**2.0_dp + &
                        (abs(rW(3)))**2.0_dp)

       ! squared residual for L2 norm times gauss weight
       ! for the current element
       ! summed over all four Gauss points
       res_s(mi) = res_s(mi) + abs(abs_res_s * abs_res_s) * 0.25_dp
       res_w(mi) = res_w(mi) + abs(abs_res_w * abs_res_w) * 0.25_dp
      end do ! gauss points
     end select


   ! End loop over elements
    end do
   !$OMP END PARALLEL DO
  !---------------------------------------------------------------------
  end subroutine calculate_elemental_residuals

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for calculating face jumps
  !---------------------------------------------------------------------
  subroutine calculate_face_jumps (M, w, el2ed, el2nd, el2neigh, &
                                   el2edl, ed_sign, nd, &
                                   a_start, a_end, b_start, b_end, &
                                   c_start, c_end, d_start, d_end, &
                                   rho, mu, Ve, S, Wg, &
                                   fjJ_s, fjJ_w, fjH_s, fjH_w)
  
  
    ! INPUT
    integer, intent(in) :: M
    real(kind=dp), intent(in) :: w
    integer, dimension(:,:), intent(in) :: el2ed, el2nd, el2neigh
    real(kind=dp), dimension(:,:), intent(in) :: el2edl, ed_sign, nd
    real(kind=dp), dimension(:,:), intent(in) :: a_start, a_end, &
                                                 b_start, b_end, &
                                                 c_start, c_end, &
                                                 d_start, d_end
    real(kind=dp), dimension(:), intent(in) :: rho, mu, Ve
    complex (kind=dp), dimension(:), intent(in) :: S, Wg

    ! OUTPUT
    ! |face jump²| for each element primal and dual problem:
    real(kind=dp), allocatable, dimension(:), intent(out) :: fjJ_s, &
                                                             fjJ_w, &
                                                             fjH_s, &
                                                             fjH_w 

    ! LOCAL variables
    integer :: i, j, l
    integer :: allo_stat
    integer :: elem1, elem2, no_ele
    real(kind=dp), dimension(3) :: p1, p2, p3 ! three nodes of triangle
    real(kind=dp) :: px, py, pz ! midpoint coordinates x,y,z of triangle
    real (kind=dp), dimension(3) :: face_normal
    real (kind=dp) :: area, diameter
    integer, dimension(6) :: elem1_ed, elem2_ed 
    real(kind=dp), dimension(3) :: grad_Lstart ! grad Lstart
    real(kind=dp), dimension(3) :: grad_Lend ! grad Lend
    real(kind=dp), dimension(3) :: N, curl_N ! Nedelec basis functions 
    complex(kind=dp), dimension(3) :: S1, W1, S2, W2, &
                                      curl_S1, curl_W1, curl_S2, curl_W2
    complex(kind=dp) :: fjJ_tmp_s, fjJ_tmp_w
    complex(kind=dp), dimension(3) :: fjH_tmp_s, fjH_tmp_w
    real(kind=dp) :: abs_S1, abs_S2, abs_W1, abs_W2
    real(kind=dp) :: abs_curl_S1, abs_curl_S2, abs_curl_W1, abs_curl_W2
    complex(kind=dp) :: be, bw, ce, cw

    !-------------------------------------------------------------------
    ! allocate space for fjJ_s(M), fjJ_w(M), fjH_s(M), fjH_w(M)
    ! will be deallocated in main program
    allocate (fjJ_s(M), fjJ_w(M), fjH_s(M), fjH_w(M),stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
                    "calculate_face_jumps: error allocating face jumps")

    ! Initialise fjJ_w and fjJ_s  fjH_w and fjH_s to zero
    fjJ_s = 0.0_dp
    fjJ_w = 0.0_dp
    fjH_s = 0.0_dp
    fjH_w = 0.0_dp
  
    elem1 = 0
    elem2 = 0

    ! set value for not-existing neighbour-element at boundary 
    ! (as in .neigh file)
    no_ele = -1

    !$OMP PARALLEL DO
     ! Start loop over elements:
     do i = 1,M

      elem1 = i ! current element

      ! Start loop over 4 neighbouring elements of elem1:
      do j = 1,4

        elem2 = el2neigh(i,j) ! current neighbouring element of elem1

        ! if the neighbouring element is existing - 
        ! elem1 is no boundary element, then calculate face jump:
        if (elem2 .ne. no_ele) then

          ! Detect common 3 nodes (p1, p2, p3) of elem1 & elem2
          call find_common_nodes(elem1, elem2, el2nd, nd, p1, p2, p3)

          ! Calculate midpoint of the common face
          call triangle_midpoint(p1, p2, p3, px, py, pz)

          ! Calculate triangle-area of the common face
          call triangle_area(p1, p2, p3, area)

          ! Calculate triangle-diameter of the common face
          call triangle_diameter(p1, p2, p3, diameter) 

          ! Calculate face normal elem1 to elem2
          call calculate_face_normal(p1, p2, p3, face_normal)

          ! Calculate solution of primal (S1) and dual (W1) problem 
          ! at face midpoint (px, py, pz) for elem1
          elem1_ed = el2ed(elem1,:) ! Find edges of elem1

          ! variables for face jumps J
          S1 = (/(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp)/)
          W1 = (/(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp)/)
          N = (/(0.0_dp),(0.0_dp),(0.0_dp)/)

          ! variables for face jumps B
          curl_S1 = (/(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp)/)
          curl_W1 = (/(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp)/)
          curl_N = (/(0.0_dp),(0.0_dp),(0.0_dp)/)

          ! Loop over 6 edges of elem1
          do l = 1,6

            ! Lstart and grad Lend vectors
            grad_Lstart = (/ b_start(elem1,l), &
                             c_start(elem1,l), &
                             d_start(elem1,l) /)

            grad_Lend = (/ b_end(elem1,l), &
                           c_end(elem1,l), &
                           d_end(elem1,l) /)

            ! Calculate Nedelec basis function for edge l
            N = (((1.0_dp/(6.0_dp*Ve(elem1)))**2.0_dp)* &
                 (a_start(elem1,l) + b_start(elem1,l)*px &
                                   + c_start(elem1,l)*py &
                                   + d_start(elem1,l)*pz)&
                 *grad_Lend &
               
                      - &
               
                      ((1.0_dp/(6.0_dp*Ve(elem1)))**2.0_dp)* &
                      (a_end(elem1,l) + b_end(elem1,l)*px &
                                      + c_end(elem1,l)*py &
                                      + d_end(elem1,l)*pz) &
                      *grad_Lstart) &
               
                      * el2edl(elem1,l) &

                      * ed_sign(elem1,l)

            ! Calculate curl_N of edge current edge
            curl_N = ((2.0_dp*el2edl(elem1,l))/ &
                      (6.0_dp*Ve(elem1))**2.0_dp)* &
                     (/(c_start(elem1,l)*d_end(elem1,l) - &
                        d_start(elem1,l)*c_end(elem1,l)), &
                       (d_start(elem1,l)*b_end(elem1,l) - &
                        b_start(elem1,l)*d_end(elem1,l)), &
                       (b_start(elem1,l)*c_end(elem1,l) - &
                        c_start(elem1,l)*b_end(elem1,l))/)&
                  * ed_sign(elem1,l)


            ! for face jumps J
            S1 =  S1 + cmplx(N,kind=dp) * S(elem1_ed(l))
            W1 =  W1 + cmplx(N,kind=dp) * Wg(elem1_ed(l))
            ! for face jumps B
            curl_S1 =  curl_S1 + &
                       cmplx(0.0_dp, 1.0_dp/(w * mu(elem1)), kind=dp) &
                       * cmplx(curl_N, kind=dp) * S(elem1_ed(l))
            curl_W1 =  curl_W1 + &
                       cmplx(0.0_dp, 1.0_dp/(w * mu(elem1)), kind=dp) &
                       * cmplx(curl_N, kind=dp) * Wg(elem1_ed(l))


          end do


          ! Calculate solution of primal (S2) and dual (W2) problem
          ! at face midpoint (px, py, pz) for elem2
          elem2_ed = el2ed(elem2,:) ! Find edges of elem1

          ! variables for face jumps J
          S2 = (/(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp)/)
          W2 = (/(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp)/)
          N = (/(0.0_dp),(0.0_dp),(0.0_dp)/)

          ! variables for face jumps B
          curl_S2 = (/(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp)/)
          curl_W2 = (/(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp)/)
          curl_N = (/(0.0_dp),(0.0_dp),(0.0_dp)/)

          ! Loop over 6 edges of elem2
          do l = 1,6

            ! Calculate grad Lstart and grad Lend vectors
            grad_Lstart = (/ b_start(elem2,l), &
                             c_start(elem2,l), &
                             d_start(elem2,l) /)

            grad_Lend = (/ b_end(elem2,l), &
                           c_end(elem2,l), &
                           d_end(elem2,l) /)

            ! Calculate Nedelec basis function for edge l
            N = (((1.0_dp/(6.0_dp*Ve(elem2)))**2.0_dp)* &
                 (a_start(elem2,l) + b_start(elem2,l)*px &
                                   + c_start(elem2,l)*py &
                                   + d_start(elem2,l)*pz)&
                 *grad_Lend &
               
                      - &
               
                      ((1.0_dp/(6.0_dp*Ve(elem2)))**2.0_dp)* &
                      (a_end(elem2,l) + b_end(elem2,l)*px &
                                      + c_end(elem2,l)*py &
                                      + d_end(elem2,l)*pz) &
                      *grad_Lstart) &
               
                      * el2edl(elem2,l) &

                      * ed_sign(elem2,l)

            ! Calculate curl_N of edge current edge
            curl_N = ((2.0_dp*el2edl(elem2,l))/ &
                      (6.0_dp*Ve(elem2))**2.0_dp)* &
                     (/(c_start(elem2,l)*d_end(elem2,l) - &
                        d_start(elem2,l)*c_end(elem2,l)), &
                       (d_start(elem2,l)*b_end(elem2,l) - &
                        b_start(elem2,l)*d_end(elem2,l)), &
                       (b_start(elem2,l)*c_end(elem2,l) - &
                        c_start(elem2,l)*b_end(elem2,l))/)&
                     * ed_sign(elem2,l)


            ! for face jumps J
            S2 =  S2 + cmplx(N,kind=dp) * S(elem2_ed(l))
            W2 =  W2 + cmplx(N,kind=dp) * Wg(elem2_ed(l))
            ! for face jumps H
            curl_S2 =  curl_S2 + &
                       cmplx(0.0_dp, 1.0_dp/(w * mu(elem2)), kind=dp) &
                     * cmplx(curl_N, kind=dp) * S(elem2_ed(l))
            curl_W2 =  curl_W2 + &
                       cmplx(0.0_dp, 1.0_dp/(w * mu(elem2)), kind=dp) &
                     * cmplx(curl_N, kind=dp) * Wg(elem2_ed(l))

          end do


          ! Calculate face jumps J with dot product for current face
          !factor = cmplx(0.0_dp, w, kind=dp) 
          fjJ_tmp_s = dot_product(cmplx(face_normal,kind=dp), &
                                  cmplx(S1 * (1.0_dp/rho(elem1)) - &
                                        S2 * (1.0_dp/rho(elem2)), &
                                        kind=dp))      
          fjJ_tmp_w = dot_product(cmplx(face_normal,kind=dp), &
                                  cmplx(W1 * (1.0_dp/rho(elem1)) - &
                                        W2 * (1.0_dp/rho(elem2)), &
                                        kind=dp))


          ! Calculate face jumps H with cross product for current face
          fjH_tmp_s = cross_product(cmplx(face_normal,kind=dp), &
                                    cmplx(curl_S1 - curl_S2, kind=dp))
          fjH_tmp_w = cross_product(cmplx(face_normal,kind=dp), &
                                    cmplx(curl_W1 - curl_W2, kind=dp))


          ! absolute values of electric fields at faces
          abs_S1 = sqrt((abs(S1(1)))**2.0_dp + &
                        (abs(S1(2)))**2.0_dp + &
                        (abs(S1(3)))**2.0_dp)
          abs_S2 = sqrt((abs(S2(1)))**2.0_dp + &
                        (abs(S2(2)))**2.0_dp + &
                        (abs(S2(3)))**2.0_dp)
          abs_W1 = sqrt((abs(W1(1)))**2.0_dp + &
                        (abs(W1(2)))**2.0_dp + &
                        (abs(W1(3)))**2.0_dp)
          abs_W2 = sqrt((abs(W2(1)))**2.0_dp + &
                        (abs(W2(2)))**2.0_dp + &
                        (abs(W2(3)))**2.0_dp)

          ! absolute values of magnetic fields at faces
          abs_curl_S1 = sqrt((abs(curl_S1(1)))**2.0_dp + &
                             (abs(curl_S1(2)))**2.0_dp + &
                             (abs(curl_S1(3)))**2.0_dp)
          abs_curl_S2 = sqrt((abs(curl_S2(1)))**2.0_dp + &
                             (abs(curl_S2(2)))**2.0_dp + &
                             (abs(curl_S2(3)))**2.0_dp)
          abs_curl_W1 = sqrt((abs(curl_W1(1)))**2.0_dp + &
                             (abs(curl_W1(2)))**2.0_dp + &
                             (abs(curl_W1(3)))**2.0_dp)
          abs_curl_W2 = sqrt((abs(curl_W2(1)))**2.0_dp + &
                             (abs(curl_W2(2)))**2.0_dp + &
                             (abs(curl_W2(3)))**2.0_dp)

          ! Weighting factors for relative field amplitudes:
          be = ((abs((abs_S1 + abs_S2)/2)) + 10E-9)**2.0_dp
          bw = ((abs((abs_W1 + abs_W2)/2)) + 10E-9)**2.0_dp
          ce = ((abs((abs_curl_S1 + abs_curl_S2)/2)) + 10E-9)**2.0_dp
          cw = ((abs((abs_curl_W1 + abs_curl_W2)/2)) + 10E-9)**2.0_dp


          ! Sum the weighted norm of the face jumps up over all 4 
          ! neighbouring faces of current element i scaled with 0.5
          ! divide by abs((abs_S1**2 + abs_S2**2)/2 + 10E-7**2) to get 
          ! a relative face jump (relative to field amplitude)
          fjJ_s(i) = fjJ_s(i) + &
                     (0.5_dp * abs(fjJ_tmp_s * fjJ_tmp_s)/ abs(be)) &
                     * area
          fjJ_w(i) = fjJ_w(i) + &
                     (0.5_dp * abs(fjJ_tmp_w * fjJ_tmp_w)/ abs(bw)) &
                     * area
          fjH_s(i) = fjH_s(i) + &
                     (0.5_dp * &
                      abs((dot_product(fjH_tmp_s, fjH_tmp_s)))/abs(ce))&
                     * area
          fjH_w(i) = fjH_w(i) + &
                     (0.5_dp * &
                      abs((dot_product(fjH_tmp_w, fjH_tmp_w)))/abs(cw))&
                      * area

        ! end if of elem2 is existing
        end if

      ! end loop over 4 neighbouring elements
      end do

     ! end loop over elements
     end do
    !$OMP END PARALLEL DO
  !---------------------------------------------------------------------
  end subroutine calculate_face_jumps

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for computing error estimates for each element using 
  !> residuals and phase jumps J and H
  !---------------------------------------------------------------------

  subroutine compute_elemental_error_estimates (errorEst_method, M, &
                                                res_s, res_w, &
                                                fjJ_s, fjJ_w, &
                                                fjH_s, fjH_w, &
                                                eta_L, Sum_eta_L, &
                                                eta_s, eta_w )
  
    ! INPUT
    integer, intent(in) :: errorEst_method, M
    real(kind=dp), dimension(:), intent(in) :: res_s, res_w, &
                                               fjJ_s, fjJ_w, &
                                               fjH_s, fjH_w

    ! OUTPUT
    ! error estimator
    real(kind=dp), allocatable, dimension(:), intent(out) :: eta_L 
    ! Sum of all element error estimators
    real(kind=dp), intent(out) :: Sum_eta_L 
    !error erstimators for primal and dual problem
    real(kind=dp), allocatable, dimension(:), intent(out) :: eta_s, &
                                                             eta_w 

    ! LOCAL variables
    integer :: i, allo_stat
    !-------------------------------------------------------------------
    ! allocate space for eta_L(M), eta_s(M), eta_w(M)
    ! eta_L(M) will be deallocated in main program
    allocate (eta_L(M), eta_s(M), eta_w(M), stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
      "compute_elemental_error_estimates: error allocating arrays")

    ! Initialise res_w and res_s to zero
    eta_L = 0.0_dp
    eta_s = 0.0_dp
    eta_w = 0.0_dp

    ! test outputs
    ! print *, 'Max primal elemental residual:', maxval(res_s)
    ! print *, 'Max dual elemental residual:', maxval(res_w)
    ! print *, 'Max primal face jump J:', maxval(fjJ_s)
    ! print *, 'Max dual face jump J:', maxval(fjJ_w)
    ! print *, 'Max primal face jump H:', maxval(fjH_s)
    ! print *, 'Max dual face jump H:', maxval(fjH_w)

    select case(errorEst_method)

      case (1)
        ! Element error estimator with residuals
        ! Start loop over elements:
        do i = 1,M
          eta_s(i) = sqrt(abs(res_s(i)))
          eta_w(i) = sqrt(abs(res_w(i)))
          eta_L(i) = eta_s(i) * eta_w(i)
        end do

      case (2)
        ! Element error estimator with residuals, face jumps J
        ! Start loop over elements:
        do i = 1,M
          eta_s(i) = sqrt(abs(res_s(i)) + abs(fjJ_s(i)))
          eta_w(i) = sqrt(abs(res_w(i)) + abs(fjJ_w(i)))
          eta_L(i) = eta_s(i) * eta_w(i)
        end do

      case (3)
        ! Element error estimator with residuals, face jumps J and H
        ! Start loop over elements:
        do i = 1,M
          eta_s(i) = sqrt(abs(res_s(i)) + abs(fjJ_s(i)) + abs(fjH_s(i)))
          eta_w(i) = sqrt(abs(res_w(i)) + abs(fjJ_w(i)) + abs(fjH_w(i)))
          eta_L(i) = eta_s(i) * eta_w(i)
        end do

      case (4)
        ! Element error estimator with face jumps J
        ! Start loop over elements:
        do i = 1,M
          eta_s(i) = sqrt(abs(fjJ_s(i)))
          eta_w(i) = sqrt(abs(fjJ_w(i)))
          eta_L(i) = eta_s(i) * eta_w(i)
        end do

      case (5)
        ! Element error estimator with face jumps H
        ! Start loop over elements:
        do i = 1,M
          eta_s(i) = sqrt(abs(fjH_s(i)))
          eta_w(i) = sqrt(abs(fjH_w(i)))
          eta_L(i) = eta_s(i) * eta_w(i)
        end do

      case (6)
        ! Element error estimator with face jumps J and face jumps H
        ! Start loop over elements:
        do i = 1,M
          eta_s(i) = sqrt(abs(fjJ_s(i)) + abs(fjH_s(i)))
          eta_w(i) = sqrt(abs(fjJ_w(i)) + abs(fjH_w(i)))
          eta_L(i) = eta_s(i) * eta_w(i)
        end do

    end select


    ! Sum of element error estimates over all elements
    Sum_eta_L = 0.0_dp
    Sum_eta_L = sum(eta_L)
    !-------------------------------------------------------------------
  end subroutine compute_elemental_error_estimates


    
  !---------------------------------------------------------------------
  !> @brief
  !> TRIANGLE GEOMETRY
  !> subroutine for calculating the unit normal on the triangle with 
  !> corner nodes p1, p2 and p3
  !---------------------------------------------------------------------
  subroutine calculate_face_normal(p1, p2, p3, face_normal)
  
    ! INPUT
    real(kind=dp), dimension(3), intent(in) :: p1, p2, p3 

    ! OUTPUT
    real (kind=dp), dimension(3), intent(out) :: face_normal
    
    ! LOCAL variables
    integer :: i
    real (kind=dp) :: l
    real (kind=dp), dimension(3) :: u, v

   !--------------------------------------------------------------------
   ! Calculate vector v: from p1 to p2
   !              and u: from p1 to p3
   
    do i = 1,3
        u(i) = p2(i) - p1(i) 
        v(i) = p3(i) - p1(i) 
    end do
        
    ! Cross product of u and v is the face normal
    face_normal(1) = u(2) * v(3) - u(3) * v(2)
    face_normal(2) = u(3) * v(1) - u(1) * v(3)
    face_normal(3) = u(1) * v(2) - u(2) * v(1)
        
    ! Length l of the face normal vector
    l = sqrt(face_normal(1)**2.0_dp + &
             face_normal(2)**2.0_dp + &
             face_normal(3)**2.0_dp)
    
    ! Scale face normal by its length to get unit length
    do i = 1,3
        face_normal(i) = face_normal(i) * (1.0_dp/l)
    end do

   !--------------------------------------------------------------------
   end subroutine calculate_face_normal
   
    !-------------------------------------------------------------------
    !> @brief
    !> TRIANGLE GEOMETRY
    !> subroutine for calculating the area of a triangle: 
    !> Needs coordinates p1, p2, p3 and returns area
    !-------------------------------------------------------------------
    subroutine triangle_area(p1, p2, p3, area) 
    
    ! INPUT
     real(kind=dp), dimension(3), intent(in) :: p1, p2, p3 

    ! OUTPUT
     real (kind=dp), intent(out) :: area
    
    ! LOCAL variables
     real (kind=dp) :: a, b, c, s
     
    !-------------------------------------------------------------------   
    ! Calculate edge lengths for triangle:
    ! a: p1 -> p2, b: p2 -> p3, c: p3 -> p1
    a = sqrt( (p2(1) - p1(1))**2.0_dp + &
              (p2(2) - p1(2))**2.0_dp + &
              (p2(3) - p1(3))**2.0_dp)
    b = sqrt( (p3(1) - p2(1))**2.0_dp + &
              (p3(2) - p2(2))**2.0_dp + &
              (p3(3) - p2(3))**2.0_dp)
    c = sqrt( (p1(1) - p3(1))**2.0_dp + &
              (p1(2) - p3(2))**2.0_dp + &
              (p1(3) - p3(3))**2.0_dp)
        
    ! Semiperimeter of triangle
    s = 0.5_dp * (a + b + c)
        
    ! Area (using Heron's formula: 
    ! http://www.mathopenref.com/heronsformula.html
    area = sqrt(s * (s - a) * (s - b) * (s - c))
    
    !-------------------------------------------------------------------
    end subroutine triangle_area 
    
    !-------------------------------------------------------------------
    !> @brief
    !> TRIANGLE GEOMETRY
    !> subroutine for calculating the diameter of the triangle with
    !> corner nodes p1, p2 and p3
    !-------------------------------------------------------------------
    subroutine triangle_diameter(p1, p2, p3, diameter)
  
    ! INPUT
     real(kind=dp), dimension(3), intent(in) :: p1, p2, p3 

    ! OUTPUT
     real (kind=dp), intent(out) :: diameter
    
    ! LOCAL variables
     real (kind=dp) :: a, b, c

   !--------------------------------------------------------------------
   ! Calculate edge lengths:
   ! a: p1 -> p2, b: p2 -> p3, c: p3 -> p1
    a = sqrt( (p2(1) - p1(1))**2.0_dp + &
              (p2(2) - p1(2))**2.0_dp + &
              (p2(3) - p1(3))**2.0_dp)
    b = sqrt( (p3(1) - p2(1))**2.0_dp + &
              (p3(2) - p2(2))**2.0_dp + &
              (p3(3) - p2(3))**2.0_dp)
    c = sqrt( (p1(1) - p3(1))**2.0_dp + &
              (p1(2) - p3(2))**2.0_dp + &
              (p1(3) - p3(3))**2.0_dp)
        
    ! Diameter = length of longest edge
    ! Which is the longest edge?
    if (a > b) then
        if (a > c) then
            diameter = a
        else
            diameter = c
        endif
    else if (b > c) then
        diameter = b
    else
        diameter = c
    end if

  !---------------------------------------------------------------------
  end subroutine triangle_diameter


  !---------------------------------------------------------------------
  !> TRIANGLE GEOMETRY
  !> subroutine for calculating the midpoint of a triangle with corner 
  !> nodes p1, p2 and p3
  !---------------------------------------------------------------------
  subroutine triangle_midpoint(p1, p2, p3, px, py, pz)
  
    ! INPUT
    real(kind=dp), dimension(3), intent(in) :: p1, p2, p3 

    ! OUTPUT
    real (kind=dp), intent(out) :: px, py, pz

   !--------------------------------------------------------------------
    px = (p1(1)+p2(1)+p3(1))/3.0_dp
    py = (p1(2)+p2(2)+p3(2))/3.0_dp
    pz = (p1(3)+p2(3)+p3(3))/3.0_dp
   !--------------------------------------------------------------------
   end subroutine triangle_midpoint


  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for calculating finding common nodes p1, p2 and p3 of 
  !> 2 adjacent elements
  !---------------------------------------------------------------------
  subroutine find_common_nodes(elem1, elem2, el2nd, nd, p1, p2, p3)
  
    ! INPUT
    integer, intent(in) :: elem1, elem2
    integer, dimension(:,:), intent(in) :: el2nd
    real(kind=dp), dimension(:,:), intent(in) :: nd

    ! OUTPUT
    real(kind=dp), dimension(3), intent(out) :: p1, p2, p3 
    
    ! LOCAL variables
    integer :: i, count
    integer, dimension(3) :: common_nd

   !--------------------------------------------------------------------

     count = 1

     do i = 1,4
      if(any(el2nd(elem1,i) == el2nd(elem2,:))) then
        common_nd(count) = el2nd(elem1,i)
        count = count+1
        if(count == 5) then 
          call Write_Message &
               (log_unit, 'Error in neighbouring elements')
        end if    
      end if
     end do

     p1 = nd(common_nd(1),:)
     p2 = nd(common_nd(2),:)
     p3 = nd(common_nd(3),:)

  !---------------------------------------------------------------------
  end subroutine find_common_nodes


  !---------------------------------------------------------------------
  !> @brief
  !> TETRAHEDRON GEOMETRY
  !> subroutine for calculating the average height of a tetrahedron 
  !> using its volume
  !---------------------------------------------------------------------
  subroutine tetr_height(Ve_tet, h_tet)
  
    ! INPUT
    real(kind=dp), intent(in) :: Ve_tet 

    ! OUTPUT
    real(kind=dp), intent(out) :: h_tet

   !--------------------------------------------------------------------
    h_tet = (sqrt(6.0_dp)/3.0_dp) * &
            ((12.0_dp/sqrt(2.0_dp)) * Ve_tet)**(1.0_dp/3.0_dp)
   !--------------------------------------------------------------------
   end subroutine tetr_height


  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for calculating element influence sources at a
  !> certain receiver location
  !---------------------------------------------------------------------
  subroutine calculate_elemental_influence_source (u1, v1, w1, &
                                                   rec1_el, el2ed, &
                                                   a_start, a_end, &
                                                   b_start, b_end, &
                                                   c_start, c_end, &
                                                   d_start, d_end, &
                                                   el2edl, ed_sign, &
                                                   Ve, w, p, &
                                                   midp_source, &
                                                   rec1_ed, infl_source)
  
    ! INPUT
    real(kind=dp), intent(in) :: u1,v1,w1
    integer, intent(in) :: rec1_el
    integer, dimension(:,:), intent(in) :: el2ed
    real(kind=dp), dimension(:,:), intent(in) :: a_start, a_end, &
                                                 b_start, b_end, &
                                                 c_start, c_end, &
                                                 d_start, d_end
    real(kind=dp), dimension(:,:), intent(in) :: el2edl
    real(kind=dp), dimension(:,:), intent(in) :: ed_sign
    real(kind=dp), dimension(:), intent(in) :: Ve
    real(kind=dp), intent(in) :: w,p
    real(kind=dp), dimension(:), intent(in) :: midp_source

    ! OUTPUT
    ! edges of the element containing receiver
    integer, dimension(6), intent(out) :: rec1_ed 
    complex(kind=dp), dimension(6), intent(out) :: infl_source

    ! LOCAL variables
    integer :: i, l
    real(kind=dp), dimension(3) :: grad_Lstart ! grad Lstart
    real(kind=dp), dimension(3) :: grad_Lend ! grad Lend
    real(kind=dp), dimension(3) :: N_rec1 ! Nedelec basis functions
    real(kind=dp), dimension(3) :: II
    real(kind=dp) :: dist ! distance of current receiver to source midp.
    !-------------------------------------------------------------------
  
     ! Find edges of receiver earth element
     rec1_ed = el2ed(rec1_el,:)

     ! initialise
     infl_source = (0.0_dp,0.0_dp)
     N_rec1 = 0.0_dp
     grad_Lstart = 0.0_dp
     grad_Lend = 0.0_dp

     ! Calculate the distance of the current receiver to source midpoint
     dist = sqrt((midp_source(1)-u1)**2.0_dp + &
                 (midp_source(2)-v1)**2.0_dp + &
                 (midp_source(3)-w1)**2.0_dp)

     ! Define identity vector (weighted)
     II(:) = 1.0_dp * w * p * dist**3.0_dp

     ! Edge numbers of current receiver element
     ! call Write_Message &
     !      (log_unit, 'Edge numbers of receiver earth element are:')
     ! do i = 1,6
     !    write(*,*) rec1_ed(i)
     ! end do

     ! Start loop over 6 edges of current element
     do l = 1,6

        ! Lstart and grad Lend vectors
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
                                   + d_end(rec1_el,l)*w1) &
                 *grad_Lstart) &
             
             * el2edl(rec1_el,l) &
             * ed_sign(rec1_el,l)



        ! Calculate influence source for the current edge of
        ! the current element
        infl_source(l) =  dot_product(cmplx(N_rec1,kind=dp),&
                                      cmplx(0.0_dp,II,kind=dp))


     end do
  !---------------------------------------------------------------------
  end subroutine calculate_elemental_influence_source


  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for writing .vtk files with error estimates for
  !> viewing in paraview
  !---------------------------------------------------------------------
  subroutine write_vtk (M, refStep, StringStep, StringEnding, &
                        MeshFileName, errorEst, eta_s, eta_w, &
                        Ve, fjJ_s, fjJ_w, fjH_s, fjH_w)
  
    ! INPUT
    integer, intent(in) :: M, refStep
    real(kind=dp), dimension(:), intent(in) :: errorEst, eta_s, eta_w, &
                                               Ve, fjJ_s, fjJ_w, &
                                                   fjH_s, fjH_w


    ! LOCAL variables
    integer :: i, mi
    integer :: opening, length, io
    character(len = 255) :: vtkFile
    character(len = 50) :: StringStep, StringEnding, MeshFileName

   !--------------------------------------------------------------------
    ! initialise length to zero
    length = 0

    ! define name of vtk file
    write(StringStep , *) refStep
    StringEnding = ".vtk"
    vtkFile = trim(adjustl(MeshFileName))// &
              trim(adjustl(StringStep))//trim(adjustl(StringEnding))
    print*, vtkFile

    ! detect file length:
    ! open file
    open (999, file = trim(vtkFile), status='old', iostat = opening)
    ! was opening successful?
    if (opening /= 0) then
      call Write_Error_Message(log_unit, &
              ' file ' // trim(vtkFile) // ' could not be opened')
    else
        length = 0 !40829
        ! loop all lines
        do
            read(999, *, iostat = io)
            ! if end of file, exit loop
            if (io/=0) exit
            ! increment line counter
            length = length + 1
        end do
        close(unit=999)
    end if


    ! Write error estimates in file
    ! open file
    open (999, file = trim(vtkFile), status='old', iostat = opening)
    ! was opening successful?
    if (opening /= 0) then
       call Write_Error_Message(log_unit, &
              ' file ' // trim(vtkFile) // ' could not be opened')
    else

       do i = 1,length
           read(999,*)
       end do

       ! write header lines in the end of vtk file
       write(999,*)'SCALARS ErrorEst_data double 1'
       write(999,*)'LOOKUP_TABLE default'

       ! write errorEst array below vtk file
       do mi = 1,M
           write(999,*) errorEst(mi)
       end do

       ! write header lines in the end of vtk file
       write(999,*)' '
       write(999,*)'SCALARS eta_s_data double 1'
       write(999,*)'LOOKUP_TABLE default'

       ! write eta_s array below vtk file
       do mi = 1,M
           write(999,*) eta_s(mi)
       end do

       ! write header lines in the end of vtk file
       write(999,*)' '
       write(999,*)'SCALARS eta_w_data double 1'
       write(999,*)'LOOKUP_TABLE default'

       ! write eta_w array below vtk file
       do mi = 1,M
           write(999,*) eta_w(mi)
       end do

       ! write header lines in the end of vtk file
       write(999,*)' '
       write(999,*)'SCALARS element_volume double 1'
       write(999,*)'LOOKUP_TABLE default'

       ! write eta_w array below vtk file
       do mi = 1,M
           write(999,*) Ve(mi)
       end do

       ! write header lines in the end of vtk file
       write(999,*)' '
       write(999,*)'SCALARS element_fjJ_s double 1'
       write(999,*)'LOOKUP_TABLE default'

       ! write fjJ_s array below vtk file
       do mi = 1,M
           write(999,*) fjJ_s(mi)
       end do

       ! write header lines in the end of vtk file
       write(999,*)' '
       write(999,*)'SCALARS element_fjJ_w double 1'
       write(999,*)'LOOKUP_TABLE default'

       ! write fjJ_w array below vtk file
       do mi = 1,M
           write(999,*) fjJ_w(mi)
       end do

       ! write header lines in the end of vtk file
       write(999,*)' '
       write(999,*)'SCALARS element_fjH_s double 1'
       write(999,*)'LOOKUP_TABLE default'

       ! write fjH_s array below vtk file
       do mi = 1,M
           write(999,*) fjH_s(mi)
       end do

       ! write header lines in the end of vtk file
       write(999,*)' '
       write(999,*)'SCALARS element_fjH_w double 1'
       write(999,*)'LOOKUP_TABLE default'

       ! write fjH_w array below vtk file
       do mi = 1,M
           write(999,*) fjH_w(mi)
       end do


       close(unit=999)
    end if
  !---------------------------------------------------------------------
  end subroutine write_vtk

end module error_estimates
