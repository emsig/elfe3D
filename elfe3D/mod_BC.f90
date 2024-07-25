!> @brief
!> Module of elfe3d containing subroutines for boundary
!>  conditions-related calculations
!!
!>  written by Paula Rulff, 24/06/2019
!!
!> So far only of Dirichlet-type boundary conditions are implemented.
!> Also includes perfect electric conductor options.
!!
!>  Copyright (C) Paula Rulff 2020
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

module BC

  use mod_util
  use define_model

  implicit none

contains

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for detecting surface edges
  !---------------------------------------------------------------------
  subroutine detect_surface_edges (M, el2ed, el2neigh, &
                                       x_start, x_end, &
                                       y_start, y_end, &
                                       z_start, z_end, &
                                       s_edges, num_s_edges)

    ! INPUT
    integer, intent(in) :: M
    integer, dimension(:,:), intent(in) :: el2ed, el2neigh
    real(kind=dp), dimension(:,:) :: x_start, x_end, &
                                     y_start, y_end, &
                                     z_start, z_end

    ! OUTPUT
    integer, allocatable, dimension(:), intent(out) :: s_edges
    integer, intent(out) :: num_s_edges

    ! LOCAL variables
    ! model dimension coordinates
    real(kind=dp) :: x_min, x_max, y_min, y_max, z_min, z_max
    integer :: l,i,j, ineigh
    integer :: allo_stat, no_ele
    integer, dimension(40000) :: dummy_s_edges


    !------------------------------------------------------------------- 
    call define_model_size(x_min, x_max, y_min, y_max, z_min, z_max)

    dummy_s_edges = 0

    ! set value for not-existing neighbour-element at boundary
    ! (as in .neigh file)
    no_ele = -1

    ! Detect surface edge numbers for given mesh:

    l = 1
      do i = 1,M
        do ineigh = 1,4 ! loop over element faces
          ! if one face of the current element is a surface face
          if (el2neigh(i,ineigh) .eq. no_ele) then
              do j = 1,6    
                  if (( &

                      ! edges of constant y
                      abs(y_start(i,j) - y_min) <= eps_dp .and. &
                          x_start(i,j) .ge. x_min .and. &
                          x_start(i,j) .le. x_max .and. &
                          z_start(i,j) .ge. z_min .and. &
                          z_start(i,j) .le. z_max .and. &
                      abs(y_end(i,j) - y_min) <= eps_dp .and. &
                          x_end(i,j) .ge. x_min .and. &
                          x_end(i,j) .le. x_max .and. &
                          z_end(i,j) .ge. z_min .and. &
                          z_end(i,j) .le. z_max .or. &
                      
                      abs(y_start(i,j) - y_max) <= eps_dp .and. &
                          x_start(i,j) .ge. x_min .and. &
                          x_start(i,j) .le. x_max .and. &
                          z_start(i,j) .ge. z_min .and. &
                          z_start(i,j) .le. z_max .and. &
                      abs(y_end(i,j) - y_max) <= eps_dp .and. &
                          x_end(i,j) .ge. x_min .and. &
                          x_end(i,j) .le. x_max .and. &
                          z_end(i,j) .ge. z_min .and. &
                          z_end(i,j) .le. z_max .or. &

                      ! edges of constant x
                      abs(x_start(i,j) - x_min) <= eps_dp .and. &
                          y_start(i,j) .ge. y_min .and. &
                          y_start(i,j) .le. y_max .and. &
                          z_start(i,j) .ge. z_min .and. &
                          z_start(i,j) .le. z_max .and. &
                      abs(x_end(i,j) - x_min) <= eps_dp .and. &
                          y_end(i,j) .ge. y_min .and. &
                          y_end(i,j) .le. y_max .and. &
                          z_end(i,j) .ge. z_min .and. &
                          z_end(i,j) .le. z_max .or. &
                      
                      abs(x_start(i,j) - x_max) <= eps_dp .and. &
                          y_start(i,j) .ge. y_min .and. &
                          y_start(i,j) .le. y_max .and. &
                          z_start(i,j) .ge. z_min .and. &
                          z_start(i,j) .le. z_max .and. &
                      abs(x_end(i,j) - x_max) <= eps_dp .and. &
                          y_end(i,j) .ge. y_min .and. &
                          y_end(i,j) .le. y_max .and. &
                          z_end(i,j) .ge. z_min .and. &
                          z_end(i,j) .le. z_max .or. &
                      
                      ! edges of constant z
                      abs(z_start(i,j) - z_min) <= eps_dp .and. &
                          x_start(i,j) .ge. x_min .and. &
                          x_start(i,j) .le. x_max .and. &
                          y_start(i,j) .ge. y_min .and. &
                          y_start(i,j) .le. y_max .and. &
                      abs(z_end(i,j) - z_min) <= eps_dp .and. &
                          x_end(i,j) .ge. x_min .and. &
                          x_end(i,j) .le. x_max .and. &
                          y_end(i,j) .ge. y_min .and. &
                          y_end(i,j) .le. y_max .or. &
                      
                      abs(z_start(i,j) - z_max) <= eps_dp .and. &
                          x_start(i,j) .ge. x_min .and. &
                          x_start(i,j) .le. x_max .and. &
                          y_start(i,j) .ge. y_min .and. &
                          y_start(i,j) .le. y_max .and. &
                      abs(z_end(i,j) - z_max) <= eps_dp .and. &
                          x_end(i,j) .ge. x_min .and. &
                          x_end(i,j) .le. x_max .and. &
                          y_end(i,j) .ge. y_min .and. &
                          y_end(i,j) .le. y_max ) .and. &   
                      

                      ! write in array if it is not yey there
                      ((any(dummy_s_edges .eq. el2ed(i,j))) .eqv. &
                        .false.)) then

                          dummy_s_edges(l) = el2ed(i,j)
                          l = l+1

                  end if ! if edge is a surface edge

            end do ! edge loop
          end if ! if element has a surface face
        end do ! face loop
      end do ! element loop


    ! call Write_Message (log_unit, 'Number of surface edges:')
    ! print *, l-1

    num_s_edges = l-1

    ! Allocate array of length num_s_edges with surface edges for 
    ! storing detected surface edges
    allocate (s_edges(num_s_edges), stat = allo_stat) 
    call allocheck(log_unit, allo_stat, "error allocating array")

    s_edges = dummy_s_edges(1:num_s_edges)

  !---------------------------------------------------------------------
  end subroutine detect_surface_edges

  

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for applying boundary conditions
  !---------------------------------------------------------------------
  subroutine apply_BC(num_s_edges, s_edges, NNZ, Agcoo, Agcol, Agrow)

    ! INPUT
    integer, intent(in) :: num_s_edges
    integer, dimension(:), intent(in) :: s_edges
  
    ! IN & OUTPUT
    integer, intent(inout) :: NNZ
    complex(kind=dp), dimension(:), intent(inout) :: Agcoo
    integer, dimension(:), intent(inout) :: Agcol, Agrow
    
    ! LOCAL variables
    complex(kind=dp), allocatable, dimension(:) :: Agcoo_dummy
    integer, allocatable, dimension(:) ::  jAgcoo_dummy, iAgcoo_dummy
    ! model dimension coordinates
    integer :: l,t, num_diag_entry
    integer :: allo_stat

    !-------------------------------------------------------------------
    ! call Write_Message (log_unit, 'Implementing Dirichlet BC')
     ! l = index of boundary (surface) edges
     ! t = index of edges in Agcoo

      do l = 1,num_s_edges

        num_diag_entry = 0

        do t = 1,NNZ
           ! detect diagonal entries
           if (Agcol(t) .eq. s_edges(l) .and. &
               Agrow(t) .eq. s_edges(l)) then
              if (num_diag_entry .eq. 0) then
                ! set diagonal element to (1.0,0.0)
                Agcoo(t) = cmplx(1.0_dp, 0.0_dp, kind=dp)
                !counter = counter + 1
                num_diag_entry = 1
              else if (num_diag_entry .gt. 0) then
                ! set duplicate diagonal element to (0.0,0.0)
                Agcoo(t) = cmplx(0.0_dp, 0.0_dp, kind=dp)
              end if

           else if (Agrow(t) .eq. s_edges(l) .and. &
                    Agcol(t) .ne. s_edges(l)) then
              ! set row to zero/delete the entries
              Agcoo(t) = cmplx(0.0_dp, 0.0_dp, kind=dp)

           else if (Agrow(t) .ne. s_edges(l) .and. &
                    Agcol(t) .eq. s_edges(l)) then
              ! set column to zero/delete the entries
              Agcoo(t) = cmplx(0.0_dp, 0.0_dp, kind=dp)
           end if

        end do

     end do

     ! Delete zero entries in Agcoo

     allocate (Agcoo_dummy(NNZ), jAgcoo_dummy(NNZ), iAgcoo_dummy(NNZ), &
               stat = allo_stat)
     call allocheck(log_unit, allo_stat, "error allocating arrays")
     Agcoo_dummy = (0.0_dp,0.0_dp) 
     jAgcoo_dummy = 0 
     iAgcoo_dummy = 0 

     ! Use dummy arrays to create Agcoo, Agrow, Agcol 
     ! without zero entries:
     t = 1
     do l = 1,NNZ
      if (abs(Agcoo(l) - cmplx(0.0_dp, 0.0_dp, kind=dp)) &
         .gt. eps_dp) then

        Agcoo_dummy(t) = Agcoo(l)
        jAgcoo_dummy(t) = Agcol(l)
        iAgcoo_dummy(t) = Agrow(l)

        t = t+1

      end if

     end do

     NNZ = t-1
     Agcoo = (0.0_dp,0.0_dp)
     Agrow = 0
     Agcol = 0 

     Agcoo(1:NNZ) = Agcoo_dummy(1:NNZ)
     Agrow(1:NNZ) = iAgcoo_dummy(1:NNZ)
     Agcol(1:NNZ) = jAgcoo_dummy(1:NNZ)

     deallocate (Agcoo_dummy, jAgcoo_dummy, iAgcoo_dummy)
    !-------------------------------------------------------------------
  end subroutine apply_BC


  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for detecting perfect electric conductor (PEC) edges
  !> so far only exactly vertical boreholes, z negative downwards
  !--------------------------------------------------------------------
  subroutine detect_PEC_edges (E, nd, ed2nd, num_PEC, &
                               x_start, y_start, z_start, &
                               x_end, y_end, z_end, &
                               PEC_edges, num_PEC_edges)

    ! INPUT
    integer, intent(in) :: E, num_PEC
    real(kind=dp), dimension(:,:), intent(in) :: nd 
    integer, dimension(:,:), intent(in) :: ed2nd
    real(kind=dp),dimension(:) :: x_start, x_end, &
                                  y_start, y_end, &
                                  z_start, z_end

    ! OUTPUT
    integer, allocatable, dimension(:), intent(out) :: PEC_edges
    integer, intent(out) :: num_PEC_edges

    ! LOCAL variables
    integer :: l,i,p
    integer :: allo_stat
    integer, dimension(E) :: dummy_PEC_edges
    !-------------------------------------------------------------------
    dummy_PEC_edges = 0
    l = 0

    num_PEC_loop: do p = 1,num_PEC
      global_edge_loop: do i = 1,E  

        if(abs(nd(ed2nd(i,1),2) - y_start(p)) <= eps_dp .and. & ! y
           abs(nd(ed2nd(i,1),1) - x_start(p)) <= eps_dp .and. & ! x
           ! node 1 on PEC edge
           ! if z negative downwards
           nd(ed2nd(i,1),3) .le. z_start(p) .and. &
           nd(ed2nd(i,1),3) .ge. z_end(p) .and. &! z

           abs(nd(ed2nd(i,2),2) - y_start(p)) <= eps_dp .and. & ! y
           abs(nd(ed2nd(i,2),1) - x_start(p)) <= eps_dp .and. & ! x
           ! node 2 on PEC edge
           ! if z negative downwards
           nd(ed2nd(i,2),3) .le. z_start(p) .and. &
           nd(ed2nd(i,2),3) .ge. z_end(p)) then !z

              ! write edge number in dummy array
              l = l+1
              dummy_PEC_edges(l) = i
              
        end if
      end do global_edge_loop 
    end do num_PEC_loop


    call Write_Message (log_unit, 'Number of PEC edges:')
    print *, l
    call Write_Message (log_unit, 'PEC coordinates:')
    print*, x_start, x_end, y_start, y_end, z_start, z_end

    num_PEC_edges = l

    ! Allocate array of length num_s_edges with surface edges for 
    ! storing detected PEC edges
    allocate (PEC_edges(num_PEC_edges), stat = allo_stat) 
    call allocheck(log_unit, allo_stat, &
                  "error allocating array PEC_edges")

    PEC_edges = dummy_PEC_edges(1:num_PEC_edges)
    !print*, PEC_edges
  !---------------------------------------------------------------------
  end subroutine detect_PEC_edges


  !---------------------------------------------------------------------
  !> @brief
  !> subroutine for applying PEC boundary conditions (E = 0 on PEC edges)
  !---------------------------------------------------------------------
  subroutine apply_PEC_BC(num_PEC_edges, PEC_edges, NNZ, &
                          Agcoo, Agcol, Agrow)

    ! INPUT
    integer, intent(in) :: num_PEC_edges
    integer, dimension(:), intent(in) :: PEC_edges
  
    ! IN & OUTPUT
    integer, intent(inout) :: NNZ
    complex(kind=dp), dimension(:), intent(inout) :: Agcoo
    integer, dimension(:), intent(inout) :: Agcol, Agrow
    
    ! LOCAL variables
    complex(kind=dp), allocatable, dimension(:) :: Agcoo_dummy
    integer, allocatable, dimension(:) ::  jAgcoo_dummy, iAgcoo_dummy
    integer :: l,t, num_diag_entry 
    integer :: allo_stat

    !-------------------------------------------------------------------
    call Write_Message (log_unit, 'Implementing PECs in system matrix')
    call Write_Message (log_unit, 'Number of PEC edges:')
    print *, num_PEC_edges

     ! l = index of PEC edges
     ! t = index of edges in Agcoo

      do l = 1,num_PEC_edges

        num_diag_entry = 0

        do t = 1,NNZ
           ! detect diagonal entries
           if (Agcol(t) .eq. PEC_edges(l) .and. &
               Agrow(t) .eq. PEC_edges(l)) then
              if (num_diag_entry .eq. 0) then
                ! set diagonal element to (1.0,0.0)
                Agcoo(t) = cmplx(1.0_dp, 0.0_dp, kind=dp)
                !counter = counter + 1
                num_diag_entry = 1
              else if (num_diag_entry .gt. 0) then
                ! set duplicate diagonal element to (0.0,0.0)
                Agcoo(t) = cmplx(0.0_dp, 0.0_dp, kind=dp)
              end if

           else if (Agrow(t) .eq. PEC_edges(l) .and. &
                    Agcol(t) .ne. PEC_edges(l)) then
              ! set row to zero/delete the entries
              Agcoo(t) = cmplx(0.0_dp, 0.0_dp, kind=dp)

           else if (Agrow(t) .ne. PEC_edges(l) .and. &
                    Agcol(t) .eq. PEC_edges(l)) then
              ! set column to zero/delete the entries
              Agcoo(t) = cmplx(0.0_dp, 0.0_dp, kind=dp)
           end if

        end do

     end do



     ! Delete zero entries in Agcoo
     allocate (Agcoo_dummy(NNZ), jAgcoo_dummy(NNZ), iAgcoo_dummy(NNZ), &
               stat = allo_stat)
     call allocheck(log_unit, allo_stat, "error allocating arrays")
     Agcoo_dummy = (0.0_dp,0.0_dp) 
     jAgcoo_dummy = 0 
     iAgcoo_dummy = 0 


     ! Use dummy arrays to create Agcoo, Agrow, Agcol 
     ! without zero entries:
     t = 1
     do l = 1,NNZ
      if (abs(Agcoo(l) - cmplx(0.0_dp, 0.0_dp, kind=dp)) &
         .gt. eps_dp) then

        Agcoo_dummy(t) = Agcoo(l)
        jAgcoo_dummy(t) = Agcol(l)
        iAgcoo_dummy(t) = Agrow(l)

        t = t+1

      end if

     end do

     NNZ = t-1
     Agcoo = (0.0_dp,0.0_dp)
     Agrow = 0
     Agcol = 0 

     Agcoo(1:NNZ) = Agcoo_dummy(1:NNZ)
     Agrow(1:NNZ) = iAgcoo_dummy(1:NNZ)
     Agcol(1:NNZ) = jAgcoo_dummy(1:NNZ)

     deallocate (Agcoo_dummy, jAgcoo_dummy, iAgcoo_dummy)
    !-------------------------------------------------------------------
  end subroutine apply_PEC_BC


end module BC
