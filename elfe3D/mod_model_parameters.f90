!> Module of eelfe3D containing subroutines to read in attr.txt file
!> and put numbers into model parameters
!!
!> written by Paula Rulff, 27/08/2018
!!
!> Copyright (C) Paula Rulff 2020
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
  ! subroutine for assigning model parameters to elements based on 
  ! their elemental attributes
  ! (schauen, wie das in EMILIA gehandhabt wird und fuer komplexere 
  ! modelle umschreiben)
  !---------------------------------------------------------------------

  subroutine read_model_param (attr, M, rho, mu)

    ! INPUT
    integer, dimension(:), intent(in) :: attr
    ! number of elements
    integer, intent(in) :: M

    ! OUTPUT
    real(kind=dp), allocatable, dimension(:), intent(out) :: rho,mu

    ! LOCAL variables
    real(kind=dp) :: res_2, res_3, res_4, res_5, res_6, res_7, res_8, res_9 ! electrical resistivity subsurface
    real(kind=dp) :: mu_r_2, mu_r_3, mu_r_4, mu_r_5, mu_r_6, mu_r_7, mu_r_8, mu_r_9 ! relative magnetic permeability
    real(kind=dp), parameter :: mu_zero = mu_0!0.000001256637061_dp
    integer :: i, allo_stat

    !-------------------------------------------------------------------

    ! Reading model parameters from mod_define_model
    ! call Write_Message (log_unit, 'Reading model parameters')
    ! uncomment, this routine not used for custEM model
    call define_model_parameters(res_2, mu_r_2, res_3, mu_r_3, res_4, mu_r_4, res_5, mu_r_5, res_6, mu_r_6, res_7, mu_r_7, res_8, mu_r_8, res_9, mu_r_9)
    
    allocate (rho(M), mu(M), stat = allo_stat)
    call allocheck(log_unit, allo_stat, "read_model_param: error allocating array rho and mu")

    rho = 0.0_dp
    mu = 0.0_dp


    ! Region attribute numbers in model parameters umwandeln:
    !call Write_Message (log_unit, 'Assigning model parameters')
!KT: besser als where statement
    do i = 1,M
       if (attr(i) .eq. 1) then
          rho(i) = 100000000.0_dp ! Ohmm (air)
          ! change back, to air, this is for custEM model:
       !    rho(i) = 8.0_dp
         mu(i) = mu_zero
       ! elseif (attr(i) .ge. 2) then
       !    rho(i) = 0.0000002_dp
       !    mu(i) = mu_zero
       ! uncomment again, this is for custEM model:
       elseif (attr(i) .eq. 2) then
          rho(i) = res_2 ! Ohmm (earth)
          mu(i) = mu_zero * mu_r_2
       elseif (attr(i) .eq. 3) then
          rho(i) = res_3 ! Ohmm (earth)
          mu(i) = mu_zero * mu_r_3
       elseif (attr(i) .eq. 4) then
          rho(i) = res_4 ! Ohmm (earth)
          mu(i) = mu_zero * mu_r_4
       elseif (attr(i) .eq. 5) then
          rho(i) = res_5 ! Ohmm (earth)
          mu(i) = mu_zero * mu_r_5
       elseif (attr(i) .eq. 6) then
          rho(i) = res_6 ! Ohmm (earth)
          mu(i) = mu_zero * mu_r_6
       elseif (attr(i) .eq. 7) then
          rho(i) = res_7 ! Ohmm (earth)
          mu(i) = mu_zero * mu_r_7
       elseif (attr(i) .eq. 8 ) then
          rho(i) = res_8 ! Ohmm (earth)
          mu(i) = mu_zero * mu_r_8
       elseif (attr(i) .eq. 9 ) then
          rho(i) = res_9 ! Ohmm (earth)
          mu(i) = mu_zero * mu_r_9
       ! elseif (attr(i) .eq. 0 ) then
       !    rho(i) = 80000.0_dp ! Ohmm (air) for Um-custEM model
       !    mu(i) = mu_zero
       end if
    end do

    ! Check if rrays contain NaN elements or are zero
    do i = 1,M
      if (rho(i) .ne. rho(i) .or. mu(i) .ne. mu(i)) then
        call Write_Message (log_unit, &
                       'model parameter arrays contain NaN elements!!!')
      else if (rho(i) .eq. 0.0_dp .or. mu(i) .eq. 0.0_dp) then
        call Write_Message (log_unit, &
             'model parameter arrays contain elements equal to zero!!!')
      end if
    end do

  end subroutine read_model_param
  !---------------------------------------------------------------------

end module model_parameters
