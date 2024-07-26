!> @brief
!> Module of elfe3D containing subroutineto read in/regionparameters.txt
!> and assign model parameters 
!> (resistivities and magnetic permeabilities) for each element
!!
!> written by Paula Rulff, 27/08/2018
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

module model_parameters

  use mod_util
  use mod_constant
  use define_model

  implicit none

contains
  !---------------------------------------------------------------------
  !> @brief 
  !> subroutine for assigning model parameters to elements based on 
  !> their element attributes
  !> Provide the following in input file: in/regionparameters.txt
  !> This input file must have the following 
  !> structure (example):
  !> # eleattr
  !> 3
  !> # eleattr rho mu_r epsilon_r
  !> 1 100000000.0 1.0 0.0
  !> 2 100 1.0 0.0
  !> 3 10 1.0 0.0
  !>
  !> epsilon_r is not used.
  !---------------------------------------------------------------------

  subroutine read_model_param (attr, M, rho, mu)

    ! INPUT
    integer, dimension(:), intent(in) :: attr
    ! number of elements
    integer, intent(in) :: M

    ! OUTPUT
    real(kind=dp), allocatable, dimension(:), intent(out) :: rho,mu

    ! LOCAL variables
    character(len=100) :: ModParaFileName
    integer :: num_attr
    integer, allocatable, dimension(:) :: region_attr
    real(kind=dp), allocatable, dimension(:) :: region_rho, &
                                                region_mu_r, &
                                                region_epsilon_r
    integer :: i, j, allo_stat
    !-------------------------------------------------------------------
    ! Allocate resistivity and magnetic permeability arrays
    allocate (rho(M), mu(M), stat = allo_stat)
    call allocheck(log_unit, allo_stat, &
                  "read_model_param: error allocating array rho and mu")
    ! initialise
    rho = 0.0_dp
    mu = 0.0_dp

    ! read region parameters from regionparameters.txt input file
    ModParaFileName = 'in/regionparameters.txt'
    ! open the file
    open (in_unit, file = trim(ModParaFileName), status='old', &
                   action = 'read', iostat = opening)

    ! was opening successful?
    if (opening /= 0) then
         write(*,*) 'file ' // trim(ModParaFileName) // &
                    ' could not be opened'
    else
      ! skip # line
      read (in_unit, *)
      ! read number of element attributes from the second line
      read (in_unit, *) num_attr
      ! skip # line
      read (in_unit, *)
      ! allocate region_attr, region_rho, region_mu_r, region_epsilon_r
      allocate (region_attr(num_attr), region_rho(num_attr), &
                region_mu_r(num_attr), region_epsilon_r(num_attr), &
                stat = allo_stat)  
      ! read all following lines
      do j = 1, num_attr
         read (in_unit, *) region_attr(j), &
                           region_rho(j), &
                           region_mu_r(j), &
                           region_epsilon_r(j) 
      end do

      ! print *, 'region_attribute', region_attr
      ! print *, 'region_rho', region_rho
      ! print *, 'region_mu_r', region_mu_r
      ! print *, 'region_epsilon_r', region_epsilon_r

      ! close modelpara file
      close (unit = in_unit)

      ! assign element resistivities and permeabilities
      do i = 1, M
        do j = 1, num_attr
          if (attr(i) .eq. region_attr(j)) then
            rho(i) = region_rho(j)
            mu(i) = mu_0 * region_mu_r(j)
          end if
        end do
      end do

    end if

    ! Check if arrays contain NaN elements or are zero
    do i = 1,M
      if (rho(i) .ne. rho(i) .or. mu(i) .ne. mu(i)) then
        call Write_Message (log_unit, &
                       'model parameter arrays contain NaN elements!!!')
      else if (rho(i) .eq. 0.0_dp .or. mu(i) .eq. 0.0_dp) then
        call Write_Message (log_unit, &
             'model parameter arrays contain elements equal to zero!!!')
      end if
    end do

    ! deallocate locally allocated arrays
    if(allocated(region_attr)) deallocate(region_attr)
    if(allocated(region_rho)) deallocate(region_rho)
    if(allocated(region_mu_r)) deallocate(region_mu_r)
    if(allocated(region_epsilon_r)) deallocate(region_epsilon_r)

  end subroutine read_model_param
  !---------------------------------------------------------------------
end module model_parameters
