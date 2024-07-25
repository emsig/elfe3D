!> @brief
!> Module of elfe3d containing definitions of constants
!!
!> written by Paula Rulff and Thomas Kalscheuer, 
!> 23/07/2019, extended 04/2020
!!
!> Last change: March 2024
!!
!> Copyright (C) Paula Rulff & Thomas Kalscheuer 2024
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

module mod_constant

  implicit none

  ! loog-file unit-number
  integer, parameter, public :: log_unit = 999
  ! save screen output to log file?
  logical, parameter :: LOGFILE_SCREEN = .true.

  ! kind type parameters
  integer, parameter, public :: ssp = kind(1.0), &
       dp = selected_real_kind(2*precision(1.0_ssp))

  ! double parameters
  real(kind=dp), parameter, public :: D0 = 0.0_dp
  real(kind=dp), parameter, public :: D1 = 1.0_dp
  real(kind=dp), parameter, public :: D2 = 2.0_dp
  real(kind=dp), parameter, public :: D10 = 10.0_dp
  complex(kind=dp), parameter, public :: ZEROW = (D0, D0)

  ! mathematical and physical constants
  real(kind=dp), parameter, public :: pi = D2*asin(D1)   ! 3.14159265359
  ! natural log of 10
  real(kind=dp), parameter, public :: ln_10 = log(D10)
  real(kind=dp), parameter, public :: eps_dp = abs(epsilon(1.0_dp))
  ! resistivity of air
  real(kind=dp), parameter, public :: rho_0 = 100000000.0_dp
  ! electric permittivity of free space
  real(kind=dp), parameter, public :: epsilon_0 = 8.854E-12_dp
  ! magnetic permeability of free space
  real(kind=dp), parameter, public :: mu_0 = 4.0E-07_dp*pi

end module mod_constant
