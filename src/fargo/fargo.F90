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
!! Implementation of a fast eulerian transport algorithm for differentially rotating disks (Masset 2000)
!!
!! See also:
!!   1. Masset, F. "FARGO: A fast eulerian transport algorithm for differentially rotating disks" (2000) A&A, 141:165-173, arXiv:astro-ph/9910390
!!   2. Kley, W., Bitsch, B., Klahr, H. "Planet migration in three-dimensional radiative discs" (2009) A&A, 506:971-987, arXiv:0908.1863
!!
!<
module fargo
! pulled by ANY
   implicit none
   real,    dimension(:, :),     allocatable :: vphi_mean
   real,    dimension(:, :, :),  allocatable :: vphi_cr
   integer, dimension(:, :),     allocatable :: nshift

   private
   public :: init_fargo, vphi_mean, vphi_cr, nshift, subtract_mean

contains

   subroutine init_fargo

      use constants,    only: GEO_RPZ
      use dataio_pub,   only: die
      use domain,       only: dom
      use global,       only: use_fargo
      use cg_leaves,    only: leaves

      implicit none

      if (.not. use_fargo) return

      if (dom%geometry_type /= GEO_RPZ) call die("[fargo:init_fargo] FARGO works only for cylindrical geometry")

   end subroutine init_fargo

   subroutine subtract_mean(dt)

      use constants,        only: xdim, ydim, zdim, LO, HI
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use grid_cont,        only: grid_container
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      
      implicit none

      real, intent(in) :: dt

      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      integer :: ifl, icg, i
      class(component_fluid), pointer :: pfl

      cgl => leaves%first
      cg => cgl%cg
      if (.not. allocated(vphi_mean)) allocate(vphi_mean(leaves%cnt, cg%lhn(xdim, LO):cg%lhn(xdim, HI)))
      if (.not. allocated(vphi_cr)) allocate(vphi_cr(leaves%cnt, flind%fluids, cg%lhn(xdim, LO):cg%lhn(xdim, HI)))
      if (.not. allocated(nshift)) allocate(nshift(leaves%cnt, cg%lhn(xdim, LO):cg%lhn(xdim, HI)))

      icg = 1
      do while (associated(cgl))
         cg => cgl%cg
         do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)
            do ifl = 1, flind%fluids
               pfl => flind%all_fluids(ifl)%fl
               vphi_mean(icg, i) = vphi_mean(icg, i) + sum(cg%u(pfl%imy, i, :, :) / cg%u(pfl%idn, i, :, :)) 
            enddo
         enddo
         vphi_mean(icg, :) = vphi_mean(icg, :) / (flind%fluids * cg%n_(ydim) * cg%n_(zdim))
         nshift(icg, :) = nint(vphi_mean(icg, :) * dt / (cg%x(:) * cg%dl(ydim)))
         do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)
            do ifl = 1, flind%fluids
               pfl => flind%all_fluids(ifl)%fl
               cg%u(pfl%imy, i, :, :) = cg%u(pfl%imy, i, :, :) - vphi_mean(icg, i) * cg%u(pfl%idn, i, :, :) 
               vphi_cr(icg, ifl, :) = vphi_mean(icg, :) - nshift(icg, :) * (cg%x(:) * cg%dl(ydim)) / dt
            enddo
         enddo
         cgl => cgl%nxt
         icg = icg + 1
      enddo

   end subroutine subtract_mean

end module fargo
