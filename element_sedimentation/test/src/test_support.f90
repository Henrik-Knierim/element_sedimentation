! ***********************************************************************
!
!   Copyright (C) 2025-2026  Henrik Knierim
!
!   This file is part of the element_sedimentation extension for MESA.
!   It is distributed under the MIT License; see the LICENSE file in the
!   root of the element_sedimentation repository for the full text.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
! ***********************************************************************

! Minimal self-checking test helpers for the element_sedimentation unit tests.

module test_support

   use const_def, only: dp, i8

   implicit none

   private

   public :: section, assert_true, assert_close, assert_equal_int, print_test_summary, next_random
   public :: num_pass, num_fail

   integer :: num_pass = 0   ! checks passed (all runners)
   integer :: num_fail = 0   ! checks failed (all runners)

contains

   subroutine section(title)
      character(len=*), intent(in) :: title
      write(*,'(a)') ''
      write(*,'(a,a)') ' === ', trim(title)
   end subroutine section


   subroutine assert_true(name, cond, nfail)
      character(len=*), intent(in) :: name
      logical, intent(in) :: cond
      integer, intent(inout) :: nfail
      if (cond) then
         num_pass = num_pass + 1
         write(*,'(a,a)') '    PASS  ', trim(name)
      else
         num_fail = num_fail + 1
         nfail = nfail + 1
         write(*,'(a,a)') '    FAIL  ', trim(name)
      end if
   end subroutine assert_true


   subroutine assert_close(name, a, b, rtol, atol, nfail)
      ! passes when |a - b| <= atol + rtol*max(|a|, |b|) (and neither is NaN)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: a, b, rtol, atol
      integer, intent(inout) :: nfail
      real(dp) :: diff, tol
      logical :: ok
      diff = abs(a - b)
      tol = atol + rtol*max(abs(a), abs(b))
      ok = (diff <= tol) .and. .not. (ieee_is_nan(a) .or. ieee_is_nan(b))
      if (ok) then
         num_pass = num_pass + 1
         write(*,'(a,a)') '    PASS  ', trim(name)
      else
         num_fail = num_fail + 1
         nfail = nfail + 1
         write(*,'(a,a,3(a,1pd26.16))') '    FAIL  ', trim(name), '  a = ', a, '  b = ', b, '  |a-b| = ', diff
      end if
   end subroutine assert_close


   subroutine assert_equal_int(name, a, b, nfail)
      character(len=*), intent(in) :: name
      integer, intent(in) :: a, b
      integer, intent(inout) :: nfail
      if (a == b) then
         num_pass = num_pass + 1
         write(*,'(a,a)') '    PASS  ', trim(name)
      else
         num_fail = num_fail + 1
         nfail = nfail + 1
         write(*,'(a,a,2(a,i12))') '    FAIL  ', trim(name), '  a = ', a, '  b = ', b
      end if
   end subroutine assert_equal_int


   subroutine print_test_summary(nfail)
      integer, intent(in) :: nfail
      write(*,'(a)') ''
      write(*,'(a,i8)') ' checks passed: ', num_pass
      write(*,'(a,i8)') ' checks failed: ', num_fail
      if (nfail == 0 .and. num_fail == 0) then
         write(*,'(a)') ' element_sedimentation tests: ALL PASSED'
      else
         write(*,'(a)') ' element_sedimentation tests: FAILED'
      end if
   end subroutine print_test_summary


   real(dp) function next_random(state)
      ! deterministic uniform pseudo-random number in (0,1) (Park-Miller minimal standard)
      integer(i8), intent(inout) :: state
      integer(i8), parameter :: a = 16807_i8, m = 2147483647_i8
      if (state <= 0_i8 .or. state >= m) state = 12345_i8
      state = mod(a*state, m)
      next_random = real(state, dp)/real(m, dp)
   end function next_random

end module test_support
