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

module initproblem

   use constants, only: dsetnamelen, ndims, LO, HI

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   ! namelist parameters
   real                   :: sigma0      !< equilibrium density
   real, dimension(ndims) :: pulse_vel   !< uniform velocity components
   real                   :: pulse_amp   !< relative amplitude of the density wave

   namelist /PROBLEM_CONTROL/  sigma0, pulse_vel, pulse_amp

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use constants,  only: I_ONE, xdim, zdim
      use dataio_pub, only: nh      ! QA_WARN required for diff_nml
      use dataio_pub, only: warn, die
      use domain,     only: dom
      use fluidindex, only: flind
      use global,     only: smalld, smallei
      use mpisetup,   only: rbuff, ibuff, lbuff, master, slave, proc, have_mpi, LAST, piernik_MPI_Bcast
      use refinement, only: set_n_updAMR, n_updAMR
      use user_hooks, only: problem_refine_derefine

      implicit none

      ! namelist default parameter values
      sigma0        = 6e-4
      pulse_vel(:)  = 0.0
      pulse_amp     = 1e-2

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1)   = sigma0
         rbuff(2)   = pulse_amp
         rbuff(2+xdim:2+zdim) = pulse_vel(:)

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         sigma0     = rbuff(1)
         pulse_amp  = rbuff(2)
         pulse_vel  = rbuff(2+xdim:2+zdim)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------
!>
!! http://www.astro.virginia.edu/VITA/ATHENA/linear_waves.html
!! http://www.astro.virginia.edu/VITA/ATHENA/lw/wave_modes.html
!< 
   subroutine problem_initial_conditions

      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use constants,        only: xdim, ydim, zdim, GEO_XYZ, GEO_RPZ
      use dataio_pub,       only: die
      use domain,           only: dom
      use fluidindex,       only: flind
      use global,           only: smallei, t
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      integer :: i, j, k

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%b = 0.0

         do j = cg%js, cg%je
            cg%u(flind%neu%idn, cg%is:cg%ie, j, cg%ks:cg%ke) = sigma0 * (1.0 + pulse_amp * sin(cg%y(j)))
         enddo
         
         if (dom%geometry_type == GEO_RPZ) then
            cg%u(flind%neu%imx, :, :, :) = pulse_vel(xdim) * cg%u(flind%neu%idn, :, :, :)
            cg%u(flind%neu%imy, :, :, :) = pulse_vel(ydim) * cg%u(flind%neu%idn, :, :, :)
            cg%u(flind%neu%imz, :, :, :) = pulse_vel(zdim) * cg%u(flind%neu%idn, :, :, :)
         else
            call die("[initproblem:problem_initial_conditions] only cylindrical geometry is supported")
         end if

         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions

!-----------------------------------------------------------------------------

end module initproblem
