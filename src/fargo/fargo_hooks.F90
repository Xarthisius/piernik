! $Id$
!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.h"
#include "macros.h"
!>
!! Various subroutine that are necessary for a fast eulerian transport algorithm for differentially rotating disks (Masset 2000)
!<
module fargo_hooks
! pulled by ANY
   implicit none
   private
   public :: avg_ang_vel

   interface
      subroutine average_angular_vel(cg, avgvphi)

         use constants,        only: xdim
         use grid_cont,        only: grid_container

         implicit none

         type(grid_container), pointer, intent(in) :: cg
         real, dimension(:), allocatable, intent(inout) :: avgvphi
      end subroutine average_angular_vel
   end interface

   procedure(average_angular_vel), pointer :: avg_ang_vel => dummy_avg_ang_vel

contains

   subroutine dummy_avg_ang_vel(cg, avgvphi)

      use constants,        only: xdim
      use dataio_pub,       only: msg, die
      use grid_cont,        only: grid_container

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      real, dimension(:), allocatable, intent(inout) :: avgvphi

      write(msg, '(a,I5,a)') "[fargo_hooks:dummy_avg_ang_vel] This subroutine should return average angular velocity in vector of ", cg%n_(xdim), " cells"
      call die(msg)

      if (allocated(avgvphi)) then
         write(msg, '(a,a)') "[fargo_hooks:dummy_avg_ang_vel] If you are reading this, today is the day you should buy a lottery ticket...\n", &
            "Cause man! Chances of this are close to nil!"
         call die(msg)
      endif

   end subroutine dummy_avg_ang_vel

end module fargo_hooks
