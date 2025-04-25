!> @brief
!> Module of elfe3d containing subroutines mostly adapted from
!> SPARSEKIT package (Saad, 2004)
!!
!> written by Paula Rulff, 12/2018
!!
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

module sparse_matrix_operations

  use mod_util

  implicit none

contains

  !---------------------------------------------------------------------
  !> @brief
  !> getu extracts the upper triangular matrix.
  !> Changed to compgetu to work with complex matrices by 
  !> Paula Rulff, Jan 2019.
  !---------------------------------------------------------------------
  subroutine compgetu ( n, NNZ, a, ja, ia, ao, jao, iao )
    !
    !  GETU extracts the upper triangular part of a matrix.
    !
    !  Discussion:
    !
    !    The routine writes the result ao, jao, iao. 
    !
    !    The routine is in place in that ao, jao, iao can be the same 
    !    as a, ja, ia if desired.
    !
    !    The diagonal element is the last element in each row.
    !    i.e. in  a(ia(i+1)-1 )
    !    ao, jao, iao may be the same as a, ja, ia on entry -- 
    !    in which case getu will overwrite the result on a, ja, ia.
    !
    !  Modified:
    !
    !    07 January 2004
    !
    !   Author:
    !
    !     Youcef Saad
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the order of the matrix.
    ! 
    !    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), 
    !    the matrix in CSR Compressed Sparse Row format.
    !
    !    Input, real AO(*), JAO(*), IAO(N+1), the upper triangular
    !    part of the input matrix, in CSR Compressed Sparse Row format.
    !

    ! INPUT
    integer, intent(in) :: n
    integer, intent(in) :: NNZ
    complex(kind=dp), dimension(:), intent(in) :: a
    integer, dimension(:), intent(in) :: ia
    integer, dimension(:), intent(in) :: ja


    ! OUTPUT
    complex(kind=dp), dimension(:), intent(inout) :: ao
    integer, dimension(:), intent(inout) :: iao
    integer, dimension(:), intent(inout) :: jao

    ! LOCAL variables
    integer :: i
    integer :: k
    integer :: kdiag
    integer :: kfirst
    integer :: ko
    complex(kind=dp) :: t

    !-------------------------------------------------------------------
    ! check array dimensions
    if (size(ao) .ne. NNZ) &
      call Write_Error_Message(log_unit, "compgetu: ao has wrong size")
    if (size(iao) .ne. n+1) &
      call Write_Error_Message(log_unit, "compgetu: iao has wrong size")
    if (size(jao) .ne. NNZ) &
      call Write_Error_Message(log_unit, "compgetu: jao has wrong size")


    ko = 0

    do i = 1, n

       kfirst = ko + 1
       kdiag = 0

       do k = ia(i), ia(i+1)-1

          if ( i <= ja(k) ) then
             ko = ko + 1
             ao(ko) = a(k)
             jao(ko) = ja(k)
             if ( ja(k) == i ) then
                kdiag = ko
             end if
          end if

       end do
       !
       !  Exchange.
       !
       if ( kdiag /= 0 .and. kdiag /= kfirst ) then

          t = ao(kdiag)
          ao(kdiag) = ao(kfirst)
          ao(kfirst) = t

          k = jao(kdiag)
          jao(kdiag) = jao(kfirst)
          jao(kfirst) = k

       end if

       iao(i) = kfirst

    end do
    !
    !  Redefine IAO(N+1).
    !
    iao(n+1) = ko + 1

  end subroutine compgetu

  !----------------------------------------------------------------------
  !> @brief
  !> coicsr converts coordinate format to CSR format (SPARSEKIT).
  !> Changed to compcoicsr to work with complex matrices by Paula Rulff
  !> July 2019.
  !----------------------------------------------------------------------
  subroutine compcoicsr (n,nnz,job,a,ja,ia)
    
    ! INPUT
    integer, intent(in) :: n, nnz, job

    ! IN/OUTPUT
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer, dimension(:), intent(inout) :: ja, ia

    ! LOCAL variables
    integer, dimension(n+1) :: iwk
    complex(kind=dp) :: t, tnext
    logical :: values, run_init_loop
    integer :: i, inext, init, j, jnext, k, ipos

    !-------------------------------------------------------------------
    ! IN-PLACE coo-csr conversion routine.
    !-------------------------------------------------------------------
    ! this subroutine converts a matrix stored in coordinate format into 
    ! the csr format. The conversion is done in place in that the arrays 
    ! a,ja,ia of the result are overwritten onto the original arrays.
    !-------------------------------------------------------------------
    ! on entry:
    !--------- 
    ! n = integer. row dimension of A.
    ! nnz = integer. number of nonzero elements in A.
    ! job   = integer. Job indicator. when job=1, the real values in a are
    !         filled. Otherwise a is not touched and the structure of the
    !         array only (i.e. ja, ia)  is obtained.
    ! a = real array of size nnz (number of nonzero elements in A)
    !         containing the nonzero elements 
    ! ja  = integer array of length nnz containing the column positions
    !     of the corresponding elements in a.
    ! ia  = integer array of length max(nnz,n+1) containing the row 
    !     positions of the corresponding elements in a.
    ! iwk = integer work array of length n+1 
    ! on return:
    !----------
    ! a
    ! ja 
    ! ia  = contains the compressed sparse row data structure for the 
    !         resulting matrix. 
    ! Note: 
    !-------
    !         the entries of the output matrix are not sorted (the column
    !         indices in each are not in increasing order) use coocsr
    !         if you want them sorted. Note also that ia has to have at
    !         least n+1 locations 
    !-------------------------------------------------------------------
    !  Coded by Y. Saad, Sep. 26 1989                                      
    !-------------------------------------------------------------------
    !-------------------------------------------------------------------
    t = (0.0_dp, 0.0_dp)
    tnext = (0.0_dp, 0.0_dp)
    values = (job .eq. 1) 

    
    ! find pointer array for resulting matrix. 
    do i = 1,n+1
       iwk(i) = 0
    end do
    do k = 1,nnz
       i = ia(k)
       iwk(i+1) = iwk(i+1)+1
    end do
    !-------------------------------------------------------------------
    iwk(1) = 1 
    do i = 2,n
       iwk(i) = iwk(i-1) + iwk(i)
    end do
    !
    !     loop for a cycle in chasing process. 
    !
    init = 1
    k = 0
    chasing_loop: do
       if (values) t = a(init)
       i = ia(init)
       j = ja(init)
       ia(init) = -1

       !----------------------------------------------------------------
       inner_loop: do while (k .lt. nnz)
          k = k+1                 
          !     current row number is i.  determine  where to go. 
          ipos = iwk(i)
          !     save the chased element. 
          if (values) tnext = a(ipos)
          inext = ia(ipos)
          jnext = ja(ipos)
          !     then occupy its location.
          if (values) a(ipos) = t
          ja(ipos) = j
          ! update pointer information for next element to come in row i 
          iwk(i) = ipos+1
          !     determine  next element to be chased,
          if (ia(ipos) .lt. 0) then
             run_init_loop = .true.
             exit inner_loop
          else
             run_init_loop = .false.
          end if
          t = tnext
          i = inext
          j = jnext 
          ia(ipos) = -1
       end do inner_loop

       select case (run_init_loop)
       case (.false.)
          exit chasing_loop
       case (.true.)
          init_loop: do
             init = init+1
             if (init .gt. nnz) exit chasing_loop
             if (ia(init) .ge. 0) exit init_loop
          end do init_loop
       end select
       !     restart chasing --  
    end do chasing_loop

    do i = 1,n
       ia(i+1) = iwk(i)
    end do
    ia(1) = 1

  end subroutine compcoicsr
 
  

  !---------------------------------------------------------------------
  !> @brief
  !> csort routine from SPARSEKIT
  !> modified to work with complex matrices by Paula Rulff, 2020
  !---------------------------------------------------------------------
  subroutine csort (n, values, a, ja, ia) 

  ! INPUT
  integer, intent(in) :: n
  logical, intent(in) :: values

  ! IN/OUTPUT
  complex(kind=dp), dimension(:), intent(inout) :: a
  integer, dimension(:), intent(inout) :: ja, ia

  ! LOCAL variables
  integer :: row, i, k, j
  complex(kind=dp) :: rj

  !---------------------------------------------------------------------
  ! This routine sorts the elements of  a matrix (stored in Compressed
  ! Sparse Row Format in increasing order of their column indices within 
  ! each row. It uses insertion sort.
  !---------------------------------------------------------------------
  ! on entry:
  !--------- 
  ! n     = the row dimension of the matrix
  ! a     = the matrix A in compressed sparse row format.
  ! ja    = the array of column indices of the elements in array a.
  ! ia    = the array of pointers to the rows. 
  ! values= logical indicating whether or not the real values a(*) must 
  !         also be permuted. if (.not. values) then the array a is not
  !         touched by csort and can be a dummy array. 
  ! 
  ! on return:
  !----------
  ! the matrix stored in the structure a, ja, ia is permuted in such a
  ! way that the column indices are in increasing order within each row.
  ! iwork(1:nnz) contains the permutation used to rearrange the elements
  !---------------------------------------------------------------------
  ! Y. Saad - recoded Dec. 20th -  2017 - 
  !---------------------------------------------------------------------

  rj = (0.0_dp,0.0_dp)

  ! for each row

        do row=1, n

  ! scan row and do an insertion sort
   
           do k=ia(row)+1, ia(row+1)-1 
              j = ja(k)
              if (values) rj = a(k)
              i = k-1;
              do while ((i .ge. ia(row)) .and. (j < ja(i)) ) 
                 ja(i+1) = ja(i)
                 if (values) a(i+1) = a(i)
                 i = i-1
                 if (i .eq. 0) exit
              enddo
              ja(i+1) = j
              if (values) a(i+1) = rj
           enddo
        enddo
        return 

  end subroutine csort

  !---------------------------------------------------------------------
  !> @brief
  !> csumdup taken from stackoverflow example
  !> modified to work with complex matrices by Paula Rulff, 2020
  !---------------------------------------------------------------------

  subroutine csumdup (Aj, Ap, Ax)

  ! Sum together duplicate column entries in each row of CSR matrix A
  ! The column indicies within each row must be in sorted order.
  ! Explicit zeros are retained.
  ! Ap, Aj, and Ax will be modified *inplace*
  !---------------------------------------------------------------------

  ! INPUT/Output
  integer, intent(inout) :: Ap(:), Aj(:)
  complex(kind=dp), intent(inout) :: Ax(:)

  !Local variables
  integer :: nnz, r1, r2, i, j, jj
  complex(kind=dp) :: x

  !---------------------------------------------------------------------


  nnz = 1
  r2 = 1
  do i = 1, size(Ap) - 1
      r1 = r2
      r2 = Ap(i+1)
      jj = r1
      do while (jj < r2)
          j = Aj(jj)
          x = Ax(jj)
          jj = jj + 1
          do while (jj < r2)
              if (Aj(jj) == j) then
                  x = x + Ax(jj)
                  jj = jj + 1
              else
                  exit
              end if
          end do
          Aj(nnz) = j
          Ax(nnz) = x
          nnz = nnz + 1
      end do
      Ap(i+1) = nnz
  end do

  end subroutine csumdup

  !---------------------------------------------------------------------
  !> @brief
  !> csrcoo from SPARSEKIT multiplies a CSR matrix A times a vector
  !> Changed to compcsrcoo to work with complex matrices by Paula Rulff,
  !> June 2022.
  !----------------------------------------------------------------------
  subroutine compcsrcoo (nrow,job,nzmax,a,ja,ia,nnz,ao,ir,jc,ierr)
    
    ! INPUT
    integer, intent(in) :: nrow, job, nzmax
    complex(kind=dp), dimension(:), intent(inout) :: a
    integer, dimension(:), intent(inout) :: ja, ia
    ! OUTPUT
    complex(kind=dp), dimension(:), intent(out) :: ao
    integer, dimension(:), intent(out) :: ir, jc
    integer, intent(out) ::  nnz, ierr
    ! LOCAL variables
    integer :: i, k, k1, k2

    !
    ! on entry: 
    !---------
    ! nrow  = dimension of the matrix.
    ! job   = integer serving as a job indicator. 
    !         if job = 1 fill in only the array ir, ignore jc, and ao.
    !         if job = 2 fill in ir, and jc but not ao 
    !         if job = 3 fill in everything.
    !     The reason why these options are provided is that on return 
    !     ao and jc are the same as a, ja. So when job = 3, a and ja are
    !     simply copied into ao, jc.  When job=2, only jc and ir are
    !     returned. With job=1 only the array ir is returned. Moreover,
    !     the algorithm is in place:
    !     call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr) 
    !     will write the output matrix in coordinate format on a, ja,ia.
    !
    ! a,
    ! ja,
    ! ia    = matrix in compressed sparse row format.
    ! nzmax = length of space available in ao, ir, jc.
    !         the code will stop immediatly if the number of
    !         nonzero elements found in input matrix exceeds nzmax.
    ! 
    ! on return:
    !----------- 
    ! ao, ir, jc = matrix in coordinate format.
    !
    ! nnz        = number of nonzero elements in matrix.
    ! ierr       = integer error indicator.
    !         ierr .eq. 0 means normal retur
    !         ierr .eq. 1 means that the the code stopped 
    !         because there was no space in ao, ir, jc 
    !         (according to the value of  nzmax).
    ! 
    ! NOTES: 1)This routine is PARTIALLY in place: csrcoo can be called 
    !         with ao being the same array as as a, and jc the same  
    !         array as ja. but ir CANNOT be the same as ia. 
    !         2) note the order in the output arrays, 
    !-------------------------------------------------------------------
    ierr = 0
    nnz = ia(nrow+1)-1
    ! print *,'nnz in iAgcsr_o', nnz
    if (nnz .gt. nzmax) then
       ierr = 1
       return
    end if
    !-------------------------------------------------------------------
    select case (job)
    case (1)
       !
       !     copy backward to allow for in-place processing. 
       !
       do i = nrow,1,-1
          k1 = ia(i+1)-1
          k2 = ia(i)
          do k=k1,k2,-1
             ir(k) = i
          end do
       end do
    case (2)
       do k = 1,nnz
          jc(k) = ja(k)
       end do
    case (3)
       do k = 1,nnz
          ao(k) = a(k)
       end do
    end select

  end subroutine compcsrcoo


end module sparse_matrix_operations
