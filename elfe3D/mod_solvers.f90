!> @brief
!> Module of elfe3D containing subroutines to call solvers 
!> currently, only MUMPS
!!
!> last change by Paula Rulff, 23/07/2024
!!
!> Copyright (C) Paula Rulff 2024
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

module solvers

  use mod_util
  use omp_lib

  implicit none

contains
!---------------------------------------------------------------------
!> @brief'
!> subroutine Pardiso_solving performs PARDISO operations including
!> matrix check and solving the system of equations
!> does nothing, solver currently not available
!---------------------------------------------------------------------
subroutine PARDISO_solving(n, a, ja, ia, b, x)

    ! INPUT
    ! length of matrix and vector (E)
    integer, intent(in) :: n
    !integer, intent(in) :: NNZ         

    ! declaration of compressed stiffness matrix Agcsr (values)
    complex(kind=dp), dimension(:), intent(in) :: a !(NNZ)

    ! declaration of compressed stiffness matrix Agcsr 
    ! (column and pointer)
    integer, dimension(:), intent(in) :: ja !(NNZ)
    integer, dimension(:), intent(in) :: ia !(n+1)

    ! declaration of source term vector Bg
    complex(kind=dp), dimension(:,:), intent(in) :: b !(n)

    ! OUTPUT
    ! Solution vector S
    complex(kind=dp), dimension(:,:), intent(inout) :: x !(n)
    ! real(kind=dp) :: mem_pardiso

    ! LOCAL variables
    integer(kind=dp), dimension(64) :: pt
    integer, dimension(64) :: iparm
    integer :: maxfct, mnum, mtype, phase, nrhs, error, msglvl, solver
    real(kind=dp), dimension(64) :: dparm 
    integer :: idum !, i 
    real(kind=dp) :: ddum
    ! number of active OpenMP thread
    integer :: omp_threads

!     !-----------------------------------------------------------------
      call Write_Error_Message(log_unit, "PARDISO is not available.")


!     ! test output
!     ! print *, 'a:'
!     ! do i = 1,9
!     !   write(*,*) a(i)
!     ! end do

!     ! print *, 'ja:'
!     ! do i = 1,9
!     !   write(*,*) ja(i)
!     ! end do

!     ! print *, 'ia:'
!     ! do i = 1,6
!     !   write(*,*) ia(i)
!     ! end do

!     omp_threads = min(omp_get_num_procs( ),omp_get_max_threads( ))
!     print  *, "OMP_NUM_THREADS = ", omp_threads

!     ! set variables
!     nrhs = 2
!     maxfct = 1
!     mnum = 1
!     mtype = 6  ! complex symmetric 6 (unsymmetric = 13)
!     solver =  0!0  ! use sparse direct method
!     !iparm(1) = 0  ! entries to default values exept iparm(3)
!     !iparm(27)= 1 ! check matrices
!     !dparm(35) = 1000
!     mem_pardiso = 0.0_dp


!     ! following lines are from developer site 
!     ! (fortran code cmplex, non-symmetric system)

!     ! Setup Pardiso control parameters und initialize the solvers     
!     ! internal adress pointers. This is only necessary for the FIRST   
!     ! call of the PARDISO solver.   

!     !  matrix
! #ifdef IFORT
!     !call pardisoinit(pt, mtype, iparm)
! #else
!     !call pardisoinit(pt, mtype, solver, iparm, dparm, error)
!     if (error .ne. 0) then
!        if (error.eq.-10 ) write(*,*) 'No license file found'
!        if (error.eq.-11 ) write(*,*) 'License is expired'
!        if (error.eq.-12 ) write(*,*) 'Wrong username or hostname'
!        stop
!     else
!        write(*,*) '[PARDISO]: License check was successful ... '
!     end if
! #endif

!     ! Numbers of Processors ( value of OMP_NUM_THREADS )
!     iparm(3) = omp_threads

!     ! Return inverse elements in full symmetric compressed CSR format 
!     ! (1-index based)  
!     ! iparm(37) = 1  

! #ifndef IFORT
!     print *, 'Check matrix'
!     ! pardiso_chk_matrix(...)
!     ! Checks the consistency of the given matrix.
!     ! Use this functionality only for debugging purposes
!     !call pardiso_chkmatrix_z  (mtype, n, a, ia, ja, error);
!     if (error .ne. 0) then
!        write(*,*) 'The following ERROR was detected: ', error
!        stop
!     end if

!     print *, 'Check vector'
!     ! pardiso_chkvec(...)
!     ! Checks the given vectors for infinite and NaN values
!     ! Input parameters (see PARDISO user manual for a description):
!     ! Use this functionality only for debugging purposes
!     !call pardiso_chkvec_z (n, nrhs, b, error);
!     if (error .ne. 0) then
!        write(*,*) 'The following ERROR was detected: ', error
!        stop
!     end if

!     ! pardiso_printstats(...) 
!     ! prints information on the matrix to STDOUT.
!     ! Use this functionality only for debugging purposes
!     !call pardiso_printstats_z (mtype, n, a, ia, ja, nrhs, b, error);
!     if (error .ne. 0) then
!        write(*,*) 'The following ERROR was detected: ', error
!        stop
!     end if
! #endif

!     ! Reordering and Symbolic Factorization, This step also allocates
!     ! all memory that is necessary for the factorization 

!     phase = 11      ! only reordering and symbolic factorization
!     msglvl = 1       ! with statistical information
! #ifdef IFORT
!     !call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
!          !idum, nrhs, iparm, msglvl, ddum, ddum, error)
! #else
!     !call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
!          !idum, nrhs, iparm, msglvl, ddum, ddum, error, dparm)
! #endif
    
!     write(*,*) 'Reordering completed ... '
!     if (error .ne. 0) then
!        write(*,*) 'The following ERROR was detected: ', error
!        stop
!     end if

!     write(*,*) 'Number of nonzeros in factors   = ',iparm(18)
!     write(*,*) 'Number of factorization MFLOPS  = ',iparm(19)

!     ! Factorization.
!     phase = 22  ! only factorization
! #ifdef IFORT
!     !call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, & 
!     !     idum, nrhs, iparm, msglvl, ddum, ddum, error) 
! #else
!     !call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, & 
!     !     idum, nrhs, iparm, msglvl, ddum, ddum, error) 
! #endif

!     write(*,*) 'Factorization completed ...  '
!     if (error .ne. 0) then
!        write(*,*) 'The following ERROR was detected: ', error
!        stop
!     end if

!     ! Back substitution and iterative refinement
!     phase = 33  ! only solve
!     iparm(8)  = 20   ! max numbers of iterative refinement steps
!     dparm(2) = 1e-12 ! relative residual runtersetzen, default ist 1e-6

! #ifdef IFORT
!     !call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, & 
!     !     idum, nrhs, iparm, msglvl, b, x, error) 
! #else
!     !call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, & 
!     !     idum, nrhs, iparm, msglvl, b, x, error, dparm) 
! #endif
!     write(*,*) 'Solve completed ... '

! !!$    write(*,*) 'The solution of the system is '
! !!$    do i = 1, n
! !!$       write(*,*) ' x(',i,') = ', x(i)
! !!$    end do
! ! memory consuption
! write(*,*) 'peak memory in KBytes symbolic factorisation within Pardiso: ', iparm(15)
! write(*,*) 'Permanent memory symbolic factorization within Pardiso: ', iparm(16)
! write(*,*) 'Memory numerical factorization and solution within Pardiso: ', iparm(17)
! ! total peak memory consumtion
! write(*,*) 'total peak memory consumption within Pardiso: ', max(iparm(15), iparm(16)+iparm(17))
! !mem_pardiso = max(iparm(15), iparm(16)+iparm(17))
! ! Termination and release of memory
! !phase = -1           ! release internal memory

! #ifdef IFORT
!     !call pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum,&
!     !     idum, nrhs, iparm, msglvl, ddum, ddum, error)
! #else
!     !call pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum,&
!     !     idum, nrhs, iparm, msglvl, ddum, ddum, error, dparm)
! #endif


!     write(*,*) 'Number of iterative refinement steps within Pardiso: ', iparm(7)



   end subroutine PARDISO_solving


  !---------------------------------------------------------------------
  !> @brief
  !> subroutine MUMPS_solving solves the system of equations
  !---------------------------------------------------------------------
  subroutine MUMPS_solving(VAL_CSR, J_CSR, I_CSR, RHSin, &
                           Solution_out, nrhs)

    implicit none

    include 'mpif.h' 
    include 'zmumps_struc.h'


    ! INPUT
    ! number of RHS
    integer, intent(in) :: nrhs  
    ! Storage for row and column indices
    integer, dimension(:), intent(in) :: I_CSR, J_CSR
    ! Storage for corresponding entry value
    complex(kind=dp), dimension(:), intent(in) :: VAL_CSR
    ! RHS
    complex(kind=dp), dimension(:), intent(in) :: RHSin

    ! OUTPUT
    ! Solution
    complex(kind=dp), dimension(:), intent(inout) :: Solution_out

    ! Local variables
    integer :: omp_threads

    ! Variables for MUMPS and MPI
    integer :: ierr
    logical :: flag

    type(zmumps_struc) mumps_par

    ! Check if mpi is initialised
    flag = .false.
    call mpi_initialized(flag, ierr)

    ! If it is not initialised, initialise
    if (flag .eqv. .false.) then
      ! Initialise mpi
      call mpi_init(ierr)
    end if

    ! ! Define communicator
    mumps_par%comm = mpi_comm_world

    ! Initialise an instance of the package for LU factorisation
    ! (sym = 0, with working host)
    mumps_par%JOB = -1
    mumps_par%SYM = 2
    mumps_par%PAR = 1

    ! Set MUMPS ordering scheme
    mumps_par%icntl(7) = 5
    mumps_par%icntl(22) = 0 ! 1 --> out of core 
    mumps_par%icntl(28) = 0 ! Parallel analysis tool
    mumps_par%icntl(29) = 0 ! 0: automatic choice, 
                            ! 1: Parmetis ordering tool
    mumps_par%icntl(11) = 2 ! Error analysis, 
                            ! set to 0 for multiple RHS, otherwise 2

    ! Set OpenMP threads
    ! Set number of threads used for solving procedures
    omp_threads = min(omp_get_num_procs(),omp_get_max_threads())
    mumps_par%icntl(16) = omp_threads
    print*,'OMP num threads', omp_threads

    call zmumps(mumps_par)

    ! Define problem on the host (processor 0)
    if ( mumps_par%MYID .eq. 0 ) then
       mumps_par%N = size(RHSin)        ! Dimension of matrix
       mumps_par%nrhs = nrhs            ! Number of right-hand-sides
       mumps_par%NNZ = size(VAL_CSR)    ! Number of non-zeros

       ! Allocate space
       allocate(mumps_par%IRN (mumps_par%NNZ) )
       allocate(mumps_par%JCN (mumps_par%NNZ) )
       allocate(mumps_par%A (mumps_par%NNZ))
       allocate(mumps_par%RHS (mumps_par%N))


       ! Set indices and values for allocated arrays
       mumps_par%IRN = I_CSR(1:mumps_par%NNZ)
       mumps_par%JCN = J_CSR(1:mumps_par%NNZ)
       mumps_par%A = VAL_CSR(1:mumps_par%NNZ)
       mumps_par%RHS = RHSin

    end if

    ! Specify element entry
    mumps_par%icntl(5) = 0  ! Centralized assembled matrix
    mumps_par%icntl(18) = 0

    mumps_par%icntl(7) = 3  ! 0 -> AMD, 2 -> AMF, 3 -> SCOTCH, 4 -> PORD
                            ! 5 -> METIS

    mumps_par%icntl(22) = 0 ! 1 --> out of core 

    mumps_par%icntl(28) = 0 ! Sequential or parallel analysis tool
                            ! 0 -> automatic choice, 1 -> sequential
                            ! 2 -> parallel

    mumps_par%icntl(29) = 0 ! 0 -> automatic choice, 1 -> PTSCOTCH,
                            ! 2 -> PARMETIS

    mumps_par%icntl(11) = 2 ! Error analysis

    ! Set OpenMP threads
    ! Set number of threads used for solving procedures
    omp_threads = min(omp_get_num_procs(),omp_get_max_threads())
    mumps_par%icntl(16) = omp_threads

    ! Call package for solution (analyse, factorise and solve)
    mumps_par%JOB = 6
    ! Allow for higher memory consumption (20% is default)
    mumps_par%icntl(14) = 20

    call zmumps(mumps_par)

    if (mumps_par%INFOG(1) .lt. 0) then
       write(6,'(A,a,I6,A,I9)') ' Error ', &
            '  mumps_par%INFOG(1)= ', mumps_par%INFOG(1), &
            '  mumps_par%INFOG(1)= ', mumps_par%INFOG(2)

       call MPI_FINALIZE(ierr)
       stop

    end if


    ! Solution has been assembled on host
    ! Write to solution vector E
    !call Allocate_cmplx1d(E,size(RHSin))
    !E = cmplx(0.0_dp,0.0_dp,kind=dp)
    if (mumps_par%MYID .eq. 0 ) then

       Solution_out = mumps_par%rhs

       deallocate(mumps_par%IRN,mumps_par%JCN,mumps_par%A,mumps_par%RHS)

    end if

    ! Destroy internal structure of MUMPS
    mumps_par%JOB = -2
    call zmumps(mumps_par)
    !call MPI_FINALIZE(ierr)

  end subroutine MUMPS_solving

  !---------------------------------------------------------------------
  !> @brief
  !> subroutine MUMPS_solving_multiple_RHS solves the system of equations
  !> for multiple RHS
  !---------------------------------------------------------------------
  subroutine MUMPS_solving_multiple_RHS(VAL_CSR, J_CSR, I_CSR, RHSin, &
                                        Solution_out, nrhs)

      implicit none

      include 'mpif.h' 
      include 'zmumps_struc.h'


      ! INPUT
      ! number of RHS
      integer, intent(in) :: nrhs  
      ! Storage for row and column indices
      integer, dimension(:), intent(in) :: I_CSR, J_CSR
      ! Storage for corresponding entry value
      complex(kind=dp), dimension(:), intent(in) :: VAL_CSR
      ! RHS
      complex(kind=dp), dimension(:,:), intent(in) :: RHSin

      ! OUTPUT
      ! Solution
      complex(kind=dp), dimension(:,:), intent(inout) :: Solution_out

      ! Local variables
      integer :: omp_threads

      ! Variables for MUMPS and MPI
      integer :: ierr, irhs
      logical :: flag
      integer :: array_start, array_end

      type(zmumps_struc) mumps_par

      ! Check if mpi is initialised
      flag = .false.
      call mpi_initialized(flag, ierr)

      ! If it is not initialised, initialise
      if (flag .eqv. .false.) then
        ! Initialise mpi
        call mpi_init(ierr)
      end if

      ! ! Define communicator
      mumps_par%comm = mpi_comm_world

      ! Initialise an instance of the package for LU factorisation
      ! (sym = 0, with working host)
      mumps_par%JOB = -1
      mumps_par%SYM = 2
      mumps_par%PAR = 1

      ! Set MUMPS ordering scheme
      mumps_par%icntl(7) = 5
      mumps_par%icntl(22) = 0 ! 1 --> out of core 
      mumps_par%icntl(28) = 0 ! Parallel analysis tool
      mumps_par%icntl(29) = 0 ! 0: automatic choice, 
                              ! 1: Parmetis ordering tool
      mumps_par%icntl(11) = 0 ! Error analysis, 
                              ! set to 0 for multiple RHS, otherwise 2

      ! Set OpenMP threads
      ! Set number of threads used for solving procedures
      omp_threads = min(omp_get_num_procs(),omp_get_max_threads())
      mumps_par%icntl(16) = omp_threads
      print*,'OMP num threads', omp_threads

      call zmumps(mumps_par)

      ! Define problem on the host (processor 0)
      if ( mumps_par%MYID .eq. 0 ) then
         mumps_par%N = size(RHSin(:,1))     ! Dimension of matrix
         mumps_par%NRHS = nrhs              ! Number of right-hand-sides
         mumps_par%NNZ = size(VAL_CSR)      ! Number of non-zeros
         mumps_par%LRHS = size(RHSin(:,1))  ! leading dimension of RHS

         ! Allocate space
         allocate(mumps_par%IRN (mumps_par%NNZ) )
         allocate(mumps_par%JCN (mumps_par%NNZ) )
         allocate(mumps_par%A (mumps_par%NNZ))
         allocate(mumps_par%RHS (mumps_par%LRHS * mumps_par%NRHS))


         ! Set indices and values for allocated arrays
         mumps_par%IRN = I_CSR(1:mumps_par%NNZ)
         mumps_par%JCN = J_CSR(1:mumps_par%NNZ)
         mumps_par%A = VAL_CSR(1:mumps_par%NNZ)

         ! Change RHS from multi-dimensional array 
         ! to one-dimensional array
         array_start = 1
         array_end = mumps_par%LRHS
         do irhs = 1, mumps_par%NRHS
             mumps_par%RHS(array_start:array_end) = RHSin(:,irhs)
             ! update array indices for next loop step
             array_start = array_start + mumps_par%LRHS
          array_end = array_end + mumps_par%LRHS
         end do

      end if

      ! Specify element entry
      mumps_par%icntl(5) = 0  ! Centralized assembled matrix
      mumps_par%icntl(18) = 0

      mumps_par%icntl(7) = 3  ! 0 -> AMD, 2 ->AMF, 3 -> COTCH, 4 ->PORD
                              ! 5 -> METIS

      mumps_par%icntl(22) = 0 ! 1 --> out of core 

      mumps_par%icntl(28) = 0 ! Sequential or parallel analysis tool
                              ! 0 -> automatic choice, 1 -> sequential
                              ! 2 -> parallel

      mumps_par%icntl(29) = 0 ! 0 -> automatic choice, 1 -> PTSCOTCH,
                              ! 2 -> PARMETIS

      mumps_par%icntl(11) = 0 ! Error analysis
      mumps_par%icntl(24) = 1 ! Null pivot row detection
      


      ! Set OpenMP threads
      ! Set number of threads used for solving procedures
      omp_threads = min(omp_get_num_procs(),omp_get_max_threads())
      mumps_par%icntl(16) = omp_threads

      ! Call package for solution (analyse, factorise and solve)
      mumps_par%JOB = 6
      ! Allow for higher memory consumption (20% is default)
      mumps_par%icntl(14) = 20

      call zmumps(mumps_par)

      if (mumps_par%INFOG(1) .lt. 0) then
         write(6,'(A,a,I6,A,I9)') ' Error ', &
              '  mumps_par%INFOG(1)= ', mumps_par%INFOG(1), &
              '  mumps_par%INFOG(1)= ', mumps_par%INFOG(2)

         call MPI_FINALIZE(ierr)
         stop

      end if


      ! Solution has been assembled on host
      if (mumps_par%MYID .eq. 0 ) then

         ! Change solution from one-dimensional array to
         ! multi-dimensional array
         array_start = 1
         array_end = mumps_par%LRHS
         do irhs = 1, mumps_par%NRHS
             Solution_out(:,irhs) = mumps_par%RHS(array_start:array_end)
             ! update array indices for next loop step
             array_start = array_start + mumps_par%LRHS
          array_end = array_end + mumps_par%LRHS
         end do

         deallocate(mumps_par%IRN,mumps_par%JCN, &
                    mumps_par%A,mumps_par%RHS)

      end if

      ! Destroy internal structure of MUMPS
      mumps_par%JOB = -2
      call zmumps(mumps_par)
      !call MPI_FINALIZE(ierr)

  end subroutine MUMPS_solving_multiple_RHS
  !---------------------------------------------------------------------

end module solvers
