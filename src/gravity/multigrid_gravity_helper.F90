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

!> \brief Multigrid routines and parameters useful for various variants of the solver

module multigrid_gravity_helper
! pulled by MULTIGRID && SELF_GRAV

   implicit none

   private
   public :: approximate_solution, fft_solve_level, &
        &    nsmoob

   ! namelist parameters
   integer(kind=4)    :: nsmoob      !< maximum number of smoothing cycles on coarsest level when cannot use FFT.

contains

!>
!! \brief This routine has to find an approximate solution for given source field and implemented differential operator
!!
!! \details First, it initializes the solution. It means that the solution is:
!! * set to 0 for coarsest level,
!! * prolonged from coarser level otherwise.
!! Then a relaxation is called. The number of relaxation passes is determined according to the level.
!!
!! If for some reason it is undesirable to initialize the solution, call approximate_solution_relax directly.
!<

   subroutine approximate_solution(curl, src, soln)

      use cg_level_coarsest,  only: coarsest
      use cg_level_connected, only: cg_level_connected_T
      use constants,          only: BND_NEGREF, fft_none
      use multigridvars,      only: nsmool
      use multigrid_Laplace,  only: approximate_solution_relax

      implicit none

      type(cg_level_connected_T), pointer, intent(inout) :: curl !< pointer to a level for which we approximate the solution
      integer(kind=4),                     intent(in)    :: src  !< index of source in cg%q(:)
      integer(kind=4),                     intent(in)    :: soln !< index of solution in cg%q(:)

      integer(kind=4) :: nsmoo

      call curl%check_dirty(src, "approx_soln src-")

      if (associated(curl, coarsest%level) .and. curl%fft_type /= fft_none) then
         call fft_solve_level(curl, src, soln)
      else
         if (associated(curl, coarsest%level)) then
            !> \todo Implement alternative bottom-solvers
            nsmoo = nsmoob
            call coarsest%level%set_q_value(soln, 0.)
         else
            nsmoo = nsmool
            call curl%coarser%prolong_q_1var(soln, bnd_type = BND_NEGREF)
            !> \warning when this is incompatible with V-cycle or other scheme, use direct call to approximate_solution_relax
         endif

         call approximate_solution_relax(curl, src, soln, nsmoo)
      endif

      call curl%check_dirty(soln, "approx_soln soln+")

   end subroutine approximate_solution

!> \brief Solve given level if allowed using FFT

   subroutine fft_solve_level(curl, src, soln)

      use cg_level_connected,  only: cg_level_connected_T
      use cg_list,             only: cg_list_element
      use constants,           only: fft_none, fft_rcr, fft_dst
      use dataio_pub,          only: die
      use grid_cont,           only: grid_container
      use named_array,         only: p3

      implicit none

      type(cg_level_connected_T), pointer, intent(inout) :: curl !< the level on which we want the solution to be performed
      integer(kind=4),                     intent(in)    :: src  !< index of source in cg%q(:)
      integer(kind=4),                     intent(in)    :: soln !< index of solution in cg%q(:)

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (associated(curl%first)) then
         if (associated(curl%first%nxt)) call die("[multigrid_gravity_helper:fft_solve_level] multicg not possible")
      else
         return
      endif
      !> \todo Check if there is one and only one cg on app processes, die if not

      if (curl%fft_type == fft_none) call die("[multigrid_gravity_helper:fft_solve_level] FFT type not set")

      cgl => curl%first
      do while (associated(cgl))
         cg => cgl%cg

         p3 => cg%q(src)%span(cg%ijkse)
         cg%mg%src(:, :, :) = p3

         ! do the convolution in Fourier space; cg%mg%src(:,:,:) -> cg%mg%fft{r}(:,:,:)
         call dfftw_execute(cg%mg%planf)

         select case (curl%fft_type)
            case (fft_rcr)
               cg%mg%fft(:,:,:)  = cg%mg%fft(:,:,:)  * cg%mg%Green3D(:,:,:)
            case (fft_dst)
               cg%mg%fftr(:,:,:) = cg%mg%fftr(:,:,:) * cg%mg%Green3D(:,:,:)
            case default
               call die("[multigrid_gravity_helper:fft_solve_level] Unknown FFT type.")
         end select

         call dfftw_execute(cg%mg%plani) ! cg%mg%fft{r}(:,:,:) -> cg%mg%src(:,:,:)

         p3 => cg%q(soln)%span(cg%ijkse)
         p3 = cg%mg%src(:, :, :)

         cgl => cgl%nxt
      enddo

   end subroutine fft_solve_level

end module multigrid_gravity_helper
