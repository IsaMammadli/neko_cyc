! Copyright (c) 2008-2020, UCHICAGO ARGONNE, LLC.
!
! The UChicago Argonne, LLC as Operator of Argonne National
! Laboratory holds copyright in the Software. The copyright holder
! reserves all rights except those expressly granted to licensees,
! and U.S. Government license rights.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the disclaimer below.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the disclaimer (as noted below)
! in the documentation and/or other materials provided with the
! distribution.
!
! 3. Neither the name of ANL nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
! UCHICAGO ARGONNE, LLC, THE U.S. DEPARTMENT OF
! ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
! TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Additional BSD Notice
! ---------------------
! 1. This notice is required to be provided under our contract with
! the U.S. Department of Energy (DOE). This work was produced at
! Argonne National Laboratory under Contract
! No. DE-AC02-06CH11357 with the DOE.
!
! 2. Neither the United States Government nor UCHICAGO ARGONNE,
! LLC nor any of their employees, makes any warranty,
! express or implied, or assumes any liability or responsibility for the
! accuracy, completeness, or usefulness of any information, apparatus,
! product, or process disclosed, or represents that its use would not
! infringe privately-owned rights.
!
! 3. Also, reference herein to any specific commercial products, process,
! or services by trade name, trademark, manufacturer or otherwise does
! not necessarily constitute or imply its endorsement, recommendation,
! or favoring by the United States Government or UCHICAGO ARGONNE LLC.
! The views and opinions of authors expressed
! herein do not necessarily state or reflect those of the United States
! Government or UCHICAGO ARGONNE, LLC, and shall
! not be used for advertising or product endorsement purposes.
!
module fluid_volflow
  use operators, only : opgrad, cdtp, rotate_cyc
  use num_types, only : rp
  use mathops, only : opchsign
  use krylov, only : ksp_t, ksp_monitor_t
  use precon, only : pc_t
  use dofmap, only : dofmap_t
  use field, only : field_t
  use field_registry, only : neko_field_registry
  use coefs, only : coef_t
  use time_scheme_controller, only : time_scheme_controller_t
  use math, only : copy, glsc2, glmin, glmax, add2, glsum, add3
  use comm
  use neko_config, only : NEKO_BCKND_DEVICE
  use device_math, only : device_cfill, device_rzero, device_copy, &
       device_add2, device_add2s2, device_glsc2, device_add3, device_addcol3
  use device_mathops, only : device_opchsign
  use gather_scatter, only : gs_t, GS_OP_ADD
  use json_module, only : json_file
  use json_utils, only: json_get
  use scratch_registry, only : scratch_registry_t
  use bc_list, only : bc_list_t
  use ax_product, only : ax_t
  implicit none
  private

   

  !> Defines volume flow
  type, public :: fluid_volflow_t
     integer :: flow_dir !< these two should be moved to params
     logical :: avflow
     real(kind=rp) :: flow_rate
     real(kind=rp) :: dtlag = 0d0
     real(kind=rp) :: bdlag = 0d0 !< Really quite pointless since we do not vary the timestep
     type(field_t) :: u_vol, v_vol, w_vol, p_vol
     real(kind=rp) :: domain_length, base_flow

     !helical_mod
     type(field_t), pointer :: xax, yax, zax!, helix_data
     real(kind=rp) :: Rc, Rpipe, pitch !helix radius, pipe radius, pitch, atan2(Rc, phi).
     !> Manager for temporary fields
     type(scratch_registry_t) :: scratch
     ! Procedure pointers (with explicit interfaces)
     procedure(adjust_interface), pointer :: adjust => null()
     procedure(compute_interface), pointer :: compute => null()
     
   contains
     procedure, pass(this) :: init => fluid_vol_flow_init
     procedure, pass(this) :: free => fluid_vol_flow_free
     procedure, pass(this) :: makebf_str => makebf_str_helix
     procedure, pass(this) :: get_ubar_str => get_ubar_str_helix
     !procedure, pass(this) :: adjust => fluid_vol_flow
     !procedure, private, pass(this) :: compute => fluid_vol_flow_compute

  end type fluid_volflow_t

    abstract interface
      subroutine adjust_interface(this, u, v, w, p, u_res, v_res, w_res, p_res, &
          c_Xh, gs_Xh, ext_bdf, rho, mu, dt, &
          bclst_dp, bclst_du, bclst_dv, bclst_dw, bclst_vel_res, &
          Ax_vel, Ax_prs, ksp_prs, ksp_vel, pc_prs, pc_vel, prs_max_iter, vel_max_iter)
          import fluid_volflow_t 
          import field_t
          import coef_t
          import gs_t
          import time_scheme_controller_t
          import bc_list_t
          import ax_t
          import ksp_t 
          import pc_t
          import rp


          class(fluid_volflow_t), intent(inout) :: this
          type(field_t), intent(inout) :: u, v, w, p
          type(field_t), intent(inout) :: u_res, v_res, w_res, p_res
          type(coef_t), intent(inout) :: c_Xh
          type(gs_t), intent(inout) :: gs_Xh
          type(time_scheme_controller_t), intent(in) :: ext_bdf
          real(kind=rp), intent(in) :: rho, mu, dt
          type(bc_list_t), intent(inout) :: bclst_dp, bclst_du, bclst_dv, bclst_dw
          type(bc_list_t), intent(inout) :: bclst_vel_res
          class(ax_t), intent(in) :: Ax_vel
          class(ax_t), intent(in) :: Ax_prs
          class(ksp_t), intent(inout) :: ksp_prs, ksp_vel
          class(pc_t), intent(inout) :: pc_prs, pc_vel
          integer, intent(in) :: prs_max_iter, vel_max_iter
      end subroutine adjust_interface

      subroutine compute_interface(this, u_res, v_res, w_res, p_res, &
          ext_bdf, gs_Xh, c_Xh, rho, mu, bd, dt, &
          bclst_dp, bclst_du, bclst_dv, bclst_dw, bclst_vel_res, &
          Ax_vel, Ax_prs, ksp_prs, ksp_vel, pc_prs, pc_vel, prs_max_iter, vel_max_iter)
          import fluid_volflow_t 
          import field_t
          import coef_t
          import gs_t
          import time_scheme_controller_t
          import bc_list_t
          import ax_t
          import ksp_t 
          import pc_t
          import rp

          class(fluid_volflow_t), intent(inout) :: this
          type(field_t), intent(inout) :: u_res, v_res, w_res, p_res
          type(coef_t), intent(inout) :: c_Xh
          type(gs_t), intent(inout) :: gs_Xh
          type(time_scheme_controller_t), intent(in) :: ext_bdf
          type(bc_list_t), intent(inout) :: bclst_dp, bclst_du, bclst_dv, bclst_dw
          type(bc_list_t), intent(inout) :: bclst_vel_res
          class(ax_t), intent(in) :: Ax_vel
          class(ax_t), intent(in) :: Ax_prs
          class(ksp_t), intent(inout) :: ksp_prs, ksp_vel
          class(pc_t), intent(inout) :: pc_prs, pc_vel
          real(kind=rp), intent(in) :: bd
          real(kind=rp), intent(in) :: rho, mu, dt
          integer, intent(in) :: vel_max_iter, prs_max_iter
      end subroutine compute_interface
   end interface

contains

  subroutine fluid_vol_flow_init(this, dm_Xh, params)
    class(fluid_volflow_t), intent(inout) :: this
    type(dofmap_t), target, intent(in) :: dm_Xh
    type(json_file), intent(inout) :: params
    logical average
    integer :: direction
    real(kind=rp) :: rate

    call this%free()

    !Initialize vol_flow (if there is a forced volume flow)
    call json_get(params, 'case.fluid.flow_rate_force.direction', direction)
    call json_get(params, 'case.fluid.flow_rate_force.value', rate)
    call json_get(params, 'case.fluid.flow_rate_force.use_averaged_flow',&
                  average) ! logical -> True=use bulk velocity / False=use volumetric velocity given by ratea above
 
    this%flow_dir = direction
    this%avflow = average 
    this%flow_rate = rate

      
   if (this%flow_dir .le. 3) then
      this%adjust => fluid_vol_flow
      this%compute => fluid_vol_flow_compute
   else if (this%flow_dir  .eq. 4) then
      this%adjust => fluid_vol_flow_str
      this%compute => fluid_vol_flow_compute_str

      !helical_mod
      call neko_field_registry%add_field(dm_Xh, "xax")
      call neko_field_registry%add_field(dm_Xh, "yax")
      call neko_field_registry%add_field(dm_Xh, "zax")
      !call neko_field_registry%add_field(dm_Xh, "helix_data")
      this%xax => neko_field_registry%get_field("xax")
      this%yax => neko_field_registry%get_field("yax")
      this%zax => neko_field_registry%get_field("zax")
      !this%helix_data => neko_field_registry%get_field("helix_data")

      call json_get(params, 'case.fluid.flow_rate_force.pitch', this%pitch)
      call json_get(params, 'case.fluid.flow_rate_force.rp', this%Rpipe)
      call json_get(params, 'case.fluid.flow_rate_force.rc', this%Rc)
      !write(*, *) 'HELIX pitch, Rpipe, Rc', this%pitch, this%Rpipe, this%Rc
   end if

    
    if (this%flow_dir .ne. 0) then
       call this%u_vol%init(dm_Xh, 'u_vol')
       call this%v_vol%init(dm_Xh, 'v_vol')
       call this%w_vol%init(dm_Xh, 'w_vol')
       call this%p_vol%init(dm_Xh, 'p_vol')
    end if



    this%scratch = scratch_registry_t(dm_Xh, 3, 1) 

  end subroutine fluid_vol_flow_init

  subroutine fluid_vol_flow_free(this)
    class(fluid_volflow_t), intent(inout) :: this

    call this%u_vol%free()
    call this%v_vol%free()
    call this%w_vol%free()
    call this%p_vol%free()

    call this%scratch%free()

  end subroutine fluid_vol_flow_free
  !> Compute flow adjustment
  !! @brief Compute pressure and velocity using fractional step method.
  !! (Tombo splitting scheme).
!-----------------------------Base compute---------------------------------------------------
 
  subroutine fluid_vol_flow_compute(this, u_res, v_res, w_res, p_res, &
       ext_bdf, gs_Xh, c_Xh, rho, mu, bd, dt, &
       bclst_dp, bclst_du, bclst_dv, bclst_dw, bclst_vel_res, &
       Ax_vel, Ax_prs, ksp_prs, ksp_vel, pc_prs, pc_vel, prs_max_iter, vel_max_iter)
    class(fluid_volflow_t), intent(inout) :: this
    type(field_t), intent(inout) :: u_res, v_res, w_res, p_res
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    type(time_scheme_controller_t), intent(in) :: ext_bdf
    type(bc_list_t), intent(inout) :: bclst_dp, bclst_du, bclst_dv, bclst_dw
    type(bc_list_t), intent(inout) :: bclst_vel_res
    class(ax_t), intent(in) :: Ax_vel
    class(ax_t), intent(in) :: Ax_prs
    class(ksp_t), intent(inout) :: ksp_prs, ksp_vel
    class(pc_t), intent(inout) :: pc_prs, pc_vel
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: rho, mu, dt
    integer, intent(in) :: vel_max_iter, prs_max_iter
    integer :: n, i
    real(kind=rp) :: xlmin, xlmax
    real(kind=rp) :: ylmin, ylmax
    real(kind=rp) :: zlmin, zlmax
    type(ksp_monitor_t) :: ksp_results(4)
    type(field_t), pointer :: ta1, ta2, ta3
    integer :: temp_indices(3)


    call this%scratch%request_field(ta1, temp_indices(1))
    call this%scratch%request_field(ta2, temp_indices(2))
    call this%scratch%request_field(ta3, temp_indices(3))


    associate(msh => c_Xh%msh, p_vol => this%p_vol, &
         u_vol => this%u_vol, v_vol => this%v_vol, w_vol => this%w_vol)

      n = c_Xh%dof%size()
      xlmin = glmin(c_Xh%dof%x, n)
      xlmax = glmax(c_Xh%dof%x, n)
      ylmin = glmin(c_Xh%dof%y, n)          !  for Y!
      ylmax = glmax(c_Xh%dof%y, n)
      zlmin = glmin(c_Xh%dof%z, n)          !  for Z!
      zlmax = glmax(c_Xh%dof%z, n)
      if (this%flow_dir .eq. 1) then
         this%domain_length = xlmax - xlmin
      end if
      if (this%flow_dir .eq. 2) then
         this%domain_length = ylmax - ylmin
      end if
      if (this%flow_dir .eq. 3) then
         this%domain_length = zlmax - zlmin
      end if

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_cfill(c_Xh%h1_d, 1.0_rp/rho, n)
         call device_rzero(c_Xh%h2_d, n)
      else
         do i = 1, n
            c_Xh%h1(i,1,1,1) = 1.0_rp / rho
            c_Xh%h2(i,1,1,1) = 0.0_rp
         end do
      end if
      c_Xh%ifh2 = .false.

      !   Compute pressure

      if (this%flow_dir .eq. 1) then
         call cdtp(p_res%x, c_Xh%h1, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
      end if

      if (this%flow_dir .eq. 2) then
         call cdtp(p_res%x, c_Xh%h1, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
      end if

      if (this%flow_dir .eq. 3) then
         call cdtp(p_res%x, c_Xh%h1, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)
      end if

      call gs_Xh%op(p_res, GS_OP_ADD)
      call bclst_dp%apply_scalar(p_res%x, n)  !apply scalar boundary conditions
      call pc_prs%update() 
      ksp_results(1) = ksp_prs%solve(Ax_prs, p_vol, p_res%x, n, &
           c_Xh, bclst_dp, gs_Xh, prs_max_iter)

      !   Compute velocity

      call opgrad(u_res%x, v_res%x, w_res%x, p_vol%x, c_Xh)

      if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
           (NEKO_BCKND_OPENCL .eq. 1)) then
         call device_opchsign(u_res%x_d, v_res%x_d, w_res%x_d, msh%gdim, n)
         call device_copy(ta1%x_d, c_Xh%B_d, n)
         call device_copy(ta2%x_d, c_Xh%B_d, n)
         call device_copy(ta3%x_d, c_Xh%B_d, n)
      else
         call opchsign(u_res%x, v_res%x, w_res%x, msh%gdim, n)
         call copy(ta1%x, c_Xh%B, n)
         call copy(ta2%x, c_Xh%B, n)
         call copy(ta3%x, c_Xh%B, n)
      end if
      call bclst_vel_res%apply_vector(ta1%x, ta2%x, ta3%x, n)

      ! add forcing

      if (NEKO_BCKND_DEVICE .eq. 1) then
         if (this%flow_dir .eq. 1) then
            call device_add2(u_res%x_d, ta1%x_d, n)
         else if (this%flow_dir .eq. 2) then
            call device_add2(v_res%x_d, ta2%x_d, n)
         else if (this%flow_dir .eq. 3) then
            call device_add2(w_res%x_d, ta3%x_d, n)
         end if
      else
         if (this%flow_dir .eq. 1) then
            call add2(u_res%x, ta1%x, n)
         else if (this%flow_dir .eq. 2) then
            call add2(v_res%x, ta2%x, n)
         else if (this%flow_dir .eq. 3) then
            call add2(w_res%x, ta3%x, n)
         end if
      end if

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_cfill(c_Xh%h1_d, mu, n)
         call device_cfill(c_Xh%h2_d, rho * (bd / dt), n)
      else
         do i = 1, n
            c_Xh%h1(i,1,1,1) = mu
            c_Xh%h2(i,1,1,1) = rho * (bd / dt)
         end do
      end if
      c_Xh%ifh2 = .true.
      !cyclic_mod_ToBeDecided -- probably needed but only for stress_formulation
      !write(*, *) 'Rotate from fluid_volflow 1'
      call rotate_cyc(u_res%x, v_res%x, w_res%x, 1, c_Xh)
      call gs_Xh%op(u_res, GS_OP_ADD)
      call gs_Xh%op(v_res, GS_OP_ADD)
      call gs_Xh%op(w_res, GS_OP_ADD)
      call rotate_cyc(u_res%x, v_res%x, w_res%x, 0, c_Xh)

      call bclst_vel_res%apply_vector(u_res%x, v_res%x, w_res%x, n) 
      call pc_vel%update()

      ksp_results(2:4) = ksp_vel%solve_coupled(Ax_vel, &
           u_vol, v_vol, w_vol,  &
           u_res%x, v_res%x, w_res%x, &
           n, c_Xh, &
           bclst_du, bclst_dv, bclst_dw, &
           gs_Xh, vel_max_iter)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         if (this%flow_dir .eq. 1) then
            this%base_flow = &
                 device_glsc2(u_vol%x_d, c_Xh%B_d, n) / this%domain_length
         end if

         if (this%flow_dir .eq. 2) then
            this%base_flow = &
                 device_glsc2(v_vol%x_d, c_Xh%B_d, n) / this%domain_length
         end if

         if (this%flow_dir .eq. 3) then
            this%base_flow = &
                 device_glsc2(w_vol%x_d, c_Xh%B_d, n) / this%domain_length
         end if
      else
         if (this%flow_dir .eq. 1) then
            this%base_flow = glsc2(u_vol%x, c_Xh%B, n) / this%domain_length
         end if

         if (this%flow_dir .eq. 2) then
            this%base_flow = glsc2(v_vol%x, c_Xh%B, n) / this%domain_length
         end if

         if (this%flow_dir .eq. 3) then
            this%base_flow = glsc2(w_vol%x, c_Xh%B, n) / this%domain_length
         end if
      end if
    end associate

    call this%scratch%relinquish_field(temp_indices)
  end subroutine fluid_vol_flow_compute
!-----------------------------Base volflow---------------------------------------------------
  !> Adjust flow volume
  !! @brief  Adjust flow volume at end of time step to keep flow rate fixed by
  !! adding an appropriate multiple of the linear solution to the Stokes
  !! problem arising from a unit forcing in the X-direction.  This assumes
  !! that the flow rate in the X-direction is to be fixed (as opposed to Y-
  !! or Z-) *and* that the periodic boundary conditions in the X-direction
  !! occur at the extreme left and right ends of the mesh.
  !!
  !! pff 6/28/98
   subroutine fluid_vol_flow(this, u, v, w, p, u_res, v_res, w_res, p_res, &
       c_Xh, gs_Xh, ext_bdf, rho, mu, dt, &
       bclst_dp, bclst_du, bclst_dv, bclst_dw, bclst_vel_res, &
       Ax_vel, Ax_prs, ksp_prs, ksp_vel, pc_prs, pc_vel, prs_max_iter, vel_max_iter)

    class(fluid_volflow_t), intent(inout) :: this
    type(field_t), intent(inout) :: u, v, w, p
    type(field_t), intent(inout) :: u_res, v_res, w_res, p_res
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    type(time_scheme_controller_t), intent(in) :: ext_bdf
    real(kind=rp), intent(in) :: rho, mu, dt
    type(bc_list_t), intent(inout) :: bclst_dp, bclst_du, bclst_dv, bclst_dw
    type(bc_list_t), intent(inout) :: bclst_vel_res
    class(ax_t), intent(in) :: Ax_vel
    class(ax_t), intent(in) :: Ax_prs
    class(ksp_t), intent(inout) :: ksp_prs, ksp_vel
    class(pc_t), intent(inout) :: pc_prs, pc_vel
    integer, intent(in) :: prs_max_iter, vel_max_iter
    real(kind=rp) :: ifcomp, flow_rate, xsec
    real(kind=rp) :: current_flow, delta_flow, scale
    integer :: n, ierr, i

    associate(u_vol => this%u_vol, v_vol => this%v_vol, &
         w_vol => this%w_vol, p_vol => this%p_vol)

      n = c_Xh%dof%size() 

      ! If either dt or the backwards difference coefficient change,
      ! then recompute base flow solution corresponding to unit forcing:

      ifcomp = 0.0_rp

      if (dt .ne. this%dtlag .or. &
           ext_bdf%diffusion_coeffs(1) .ne. this%bdlag) then
         ifcomp = 1.0_rp
      end if

      this%dtlag = dt
      this%bdlag = ext_bdf%diffusion_coeffs(1)

      call MPI_Allreduce(MPI_IN_PLACE, ifcomp, 1, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

      if (ifcomp .gt. 0d0) then
         call this%compute(u_res, v_res, w_res, p_res, &
              ext_bdf, gs_Xh, c_Xh, rho, mu, ext_bdf%diffusion_coeffs(1), dt, &
              bclst_dp, bclst_du, bclst_dv, bclst_dw, bclst_vel_res, &
              Ax_vel, Ax_prs, ksp_prs, ksp_vel, pc_prs, pc_vel, prs_max_iter, &
              vel_max_iter)
      end if

      if (NEKO_BCKND_DEVICE .eq. 1) then
         if (this%flow_dir .eq. 1) then
            current_flow = &
                 device_glsc2(u%x_d, c_Xh%B_d, n) / this%domain_length  ! for X
         else if (this%flow_dir .eq. 2) then
            current_flow = &
                 device_glsc2(v%x_d, c_Xh%B_d, n) / this%domain_length  ! for Y
         else if (this%flow_dir .eq. 3) then
            current_flow = &
                 device_glsc2(w%x_d, c_Xh%B_d, n) / this%domain_length  ! for Z
         end if
      else
         if (this%flow_dir .eq. 1) then
            current_flow = glsc2(u%x, c_Xh%B, n) / this%domain_length  ! for X
         else if (this%flow_dir .eq. 2) then
            current_flow = glsc2(v%x, c_Xh%B, n) / this%domain_length  ! for Y
         else if (this%flow_dir .eq. 3) then
            current_flow = glsc2(w%x, c_Xh%B, n) / this%domain_length  ! for Z
         end if
      end if

      if (this%avflow) then
         xsec = c_Xh%volume / this%domain_length
         flow_rate = this%flow_rate*xsec
      else
         flow_rate = this%flow_rate
      end if

      delta_flow = flow_rate - current_flow
      scale = delta_flow / this%base_flow

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_add2s2(u%x_d, u_vol%x_d, scale, n)
         call device_add2s2(v%x_d, v_vol%x_d, scale, n)
         call device_add2s2(w%x_d, w_vol%x_d, scale, n)
         call device_add2s2(p%x_d, p_vol%x_d, scale, n)
      else
         do concurrent (i = 1: n)
            u%x(i,1,1,1) = u%x(i,1,1,1) + scale * u_vol%x(i,1,1,1)
            v%x(i,1,1,1) = v%x(i,1,1,1) + scale * v_vol%x(i,1,1,1)
            w%x(i,1,1,1) = w%x(i,1,1,1) + scale * w_vol%x(i,1,1,1)
            p%x(i,1,1,1) = p%x(i,1,1,1) + scale * p_vol%x(i,1,1,1)
         end do
      end if
    end associate

  end subroutine fluid_vol_flow

!-----------------------------STR compute------------------------------------------------------------
   !compute_vol_soln+plan4_vol
   subroutine fluid_vol_flow_compute_str(this, u_res, v_res, w_res, p_res, &
       ext_bdf, gs_Xh, c_Xh, rho, mu, bd, dt, &
       bclst_dp, bclst_du, bclst_dv, bclst_dw, bclst_vel_res, &
       Ax_vel, Ax_prs, ksp_prs, ksp_vel, pc_prs, pc_vel, prs_max_iter, vel_max_iter)

    class(fluid_volflow_t), intent(inout) :: this
    type(field_t), intent(inout) :: u_res, v_res, w_res, p_res
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    type(time_scheme_controller_t), intent(in) :: ext_bdf
    type(bc_list_t), intent(inout) :: bclst_dp, bclst_du, bclst_dv, bclst_dw
    type(bc_list_t), intent(inout) :: bclst_vel_res
    class(ax_t), intent(in) :: Ax_vel
    class(ax_t), intent(in) :: Ax_prs
    class(ksp_t), intent(inout) :: ksp_prs, ksp_vel
    class(pc_t), intent(inout) :: pc_prs, pc_vel
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: rho, mu, dt
    integer, intent(in) :: vel_max_iter, prs_max_iter
    integer :: n, i
    real(kind=rp) :: xlmin, xlmax
    real(kind=rp) :: ylmin, ylmax
    real(kind=rp) :: zlmin, zlmax
    type(ksp_monitor_t) :: ksp_results(4)
    type(field_t), pointer :: ta1, ta2, ta3
    integer :: temp_indices(3)
    real(kind=rp), allocatable :: vxc(:,:,:,:), vyc(:,:,:,:), vzc(:,:,:,:), p_res_y(:,:,:,:), p_res_z(:,:,:,:)
    !type(field_t), pointer :: fx, fy, fz
    real(kind=rp), allocatable, dimension(:,:,:,:) :: fx, fy, fz
    !type(field_t), pointer :: helix_data


    call this%scratch%request_field(ta1, temp_indices(1))
    call this%scratch%request_field(ta2, temp_indices(2))
    call this%scratch%request_field(ta3, temp_indices(3))

    associate(msh => c_Xh%msh, p_vol => this%p_vol, &
         u_vol => this%u_vol, v_vol => this%v_vol, w_vol => this%w_vol)

      n = c_Xh%dof%size()
   !-------compute_vol_soln????

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_cfill(c_Xh%h1_d, 1.0_rp/rho, n)
         call device_rzero(c_Xh%h2_d, n)
      else
         do i = 1, n
            c_Xh%h1(i,1,1,1) = 1.0_rp / rho
            c_Xh%h2(i,1,1,1) = 0.0_rp
         end do
      end if
      c_Xh%ifh2 = .false.
   !--------plan4_vol
      !   Compute pressure
      !TODOS: ortho exist . Check fluid_pnpn line 750
      allocate(vxc(c_Xh%Xh%lx, c_Xh%Xh%ly, c_Xh%Xh%lz, c_Xh%msh%nelv))
      allocate(vyc(c_Xh%Xh%lx, c_Xh%Xh%ly, c_Xh%Xh%lz, c_Xh%msh%nelv))
      allocate(vzc(c_Xh%Xh%lx, c_Xh%Xh%ly, c_Xh%Xh%lz, c_Xh%msh%nelv))
      allocate(fx(c_Xh%Xh%lx, c_Xh%Xh%ly, c_Xh%Xh%lz, c_Xh%msh%nelv))
      allocate(fy(c_Xh%Xh%lx, c_Xh%Xh%ly, c_Xh%Xh%lz, c_Xh%msh%nelv))
      allocate(fz(c_Xh%Xh%lx, c_Xh%Xh%ly, c_Xh%Xh%lz, c_Xh%msh%nelv))

      !call makebf_str(u_vol%x, v_vol%x, w_vol%x, c_Xh) ! Store Bf in uvol,vvol,wvol
      !for gs_Xh%op field_t is needed. So we cannot pass vxc but u_vol can be used
      !call gs_Xh%op(u_vol, GS_OP_ADD) !call opdssum   (vxc,vyc,vzc) ! Where is the opdssum
      !call gs_Xh%op(v_vol, GS_OP_ADD)
      !call gs_Xh%op(w_vol, GS_OP_ADD)
      call this%makebf_str(fx, fy, fz, c_Xh)
      !cyclic_mod_checked
      !write(*, *) 'Rotate from fluid_volflow 2'
      call rotate_cyc(fx, fy, fz, 1, c_Xh)
      call gs_Xh%op(fx, n, GS_OP_ADD)
      call gs_Xh%op(fy, n, GS_OP_ADD)
      call gs_Xh%op(fz, n, GS_OP_ADD)
      call rotate_cyc(fx, fy, fz, 0, c_Xh)

      !this%case%fluid%c_Xh%gs_h%op(this%veldiv,GS_OP_ADD) ! for opdssum
      if (NEKO_BCKND_DEVICE .eq. 1) then
         do i = 1, n
            !FIX! vxc(i,1,1,1) = vxc(i,1,1,1)*c_Xh%Binv_d(i,1,1,1)/rho ! /vtrans which is convective coefficient rho?
            !FIX!vyc(i,1,1,1) = vyc(i,1,1,1)*c_Xh%Binv_d(i,1,1,1)/rho 
            !FIX!vzc(i,1,1,1) = vzc(i,1,1,1)*c_Xh%Binv_d(i,1,1,1)/rho
         end do
      else
         !do i = 1, n
         do concurrent (i = 1: n)
            vxc(i,1,1,1) = fx(i,1,1,1)*c_Xh%Binv(i,1,1,1)/rho !after this we start using vxc
            vyc(i,1,1,1) = fy(i,1,1,1)*c_Xh%Binv(i,1,1,1)/rho
            vzc(i,1,1,1) = fz(i,1,1,1)*c_Xh%Binv(i,1,1,1)/rho
         end do
      end if

      !first argument stores the result
      allocate(p_res_y(c_Xh%Xh%lx, c_Xh%Xh%ly, c_Xh%Xh%lz, c_Xh%msh%nelv))
      allocate(p_res_z(c_Xh%Xh%lx, c_Xh%Xh%ly, c_Xh%Xh%lz, c_Xh%msh%nelv))
      call cdtp(p_res%x,   vxc, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh) !2nd argument c_Xh%h1 changed to vxc
      call cdtp(p_res_y,   vyc, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh) !2nd argument c_Xh%h1 changed to vycc
      call cdtp(p_res_z,   vzc, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh) !2nd argument c_Xh%h1 changed to vzc !skipped if argument for 2D

      if (NEKO_BCKND_DEVICE .eq. 1) then
         !FIX!call device_add3(p_res%x_d, p_res_y, p_res_z, n) ! instead of add2 in 2 steps, used add3
      else
         call add3(p_res%x, p_res_y, p_res_z, n)
      end if
      
      !Untouched as plan4_vol and plan4_vol_azm are same here
      call gs_Xh%op(p_res, GS_OP_ADD) !gather scatter operation
      call bclst_dp%apply_scalar(p_res%x, n)
      call pc_prs%update()
      ksp_results(1) = ksp_prs%solve(Ax_prs, p_vol, p_res%x, n, &
           c_Xh, bclst_dp, gs_Xh, prs_max_iter)

      !   Compute velocity
      !Untouched as plan4_vol and plan4_vol_azm are same here
      call opgrad(u_res%x, v_res%x, w_res%x, p_vol%x, c_Xh)

      if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
           (NEKO_BCKND_OPENCL .eq. 1)) then
         call device_opchsign(u_res%x_d, v_res%x_d, w_res%x_d, msh%gdim, n)
         call device_copy(ta1%x_d, c_Xh%B_d, n)       !!copy c_Xh%B_d to ta1%x_d
         call device_copy(ta2%x_d, c_Xh%B_d, n)
         call device_copy(ta3%x_d, c_Xh%B_d, n)
      else
         call opchsign(u_res%x, v_res%x, w_res%x, msh%gdim, n)
         call copy(ta1%x, c_Xh%B, n)
         call copy(ta2%x, c_Xh%B, n)
         call copy(ta3%x, c_Xh%B, n)
      end if
      call bclst_vel_res%apply_vector(ta1%x, ta2%x, ta3%x, n) !apply bc to vector field

      ! add forcing
      if (NEKO_BCKND_DEVICE .eq. 1) then
         do i=1, n 
            !FIX!call device_addcol3(u_res%x_d, ta1%x_d, vxc, n)
            !FIX!call device_addcol3(u_res%x_d, ta2%x_d, vyc, n)
            !FIX!call device_addcol3(u_res%x_d, ta3%x_d, vzc, n)
         end do
      else
         do concurrent (i = 1: n)
            u_res%x(i,1,1,1) = u_res%x(i,1,1,1)+ta1%x(i,1,1,1)*vxc(i,1,1,1)
            v_res%x(i,1,1,1) = v_res%x(i,1,1,1)+ta2%x(i,1,1,1)*vyc(i,1,1,1)
            w_res%x(i,1,1,1) = w_res%x(i,1,1,1)+ta3%x(i,1,1,1)*vzc(i,1,1,1)
         end do
      end if

      !TODOS: if(ifexplvis) may this be needed in case of cyclic
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_cfill(c_Xh%h1_d, mu, n)
         call device_cfill(c_Xh%h2_d, rho * (bd / dt), n)
      else
         do i = 1, n
            c_Xh%h1(i,1,1,1) = mu
            c_Xh%h2(i,1,1,1) = rho * (bd / dt)
         end do
      end if
      c_Xh%ifh2 = .true.
      !cyclic_mod_checked
      !write(*, *) 'Rotate from fluid_volflow 3'
      call rotate_cyc(u_res%x, v_res%x, w_res%x, 1, c_Xh)
      call gs_Xh%op(u_res, GS_OP_ADD)
      call gs_Xh%op(v_res, GS_OP_ADD)
      call gs_Xh%op(w_res, GS_OP_ADD)
      call rotate_cyc(u_res%x, v_res%x, w_res%x, 0, c_Xh)

      call bclst_vel_res%apply_vector(u_res%x, v_res%x, w_res%x, n)!apply bc to vector field
      call pc_vel%update()
      ksp_results(2:4) = ksp_vel%solve_coupled(Ax_vel, &
           u_vol, v_vol, w_vol, &
           !fx, fy, fz, &
           u_res%x, v_res%x, w_res%x, &
           n, c_Xh, &
           bclst_du, bclst_dv, bclst_dw, &
           gs_Xh, vel_max_iter)

      !compute_vol_soln 
      !TODOS: write device version
      if (NEKO_BCKND_DEVICE .eq. 1) then 
            !FIX!this%base_flow = get_ubar_str(u_vol,v_vol,w_vol, c_Xh) 
      else    
            this%base_flow = this%get_ubar_str(u_vol%x, v_vol%x, w_vol%x, c_Xh) 
            !this%base_flow = get_ubar_str(vxc,vyc,vzc, c_Xh) 
      end if

      !TODOS: if(ifexplvis) may this be needed in case of cyclic

    end associate
    if (allocated(vxc)) then
       deallocate(vxc)
    end if
    if (allocated(vyc)) then
       deallocate(vyc)
    end if
    if (allocated(vzc)) then
       deallocate(vzc)
    end if
    if (allocated(fx)) then
       deallocate(fx)
    end if
    if (allocated(fy)) then
       deallocate(fy)
    end if
    if (allocated(fz)) then
       deallocate(fz)
    end if
    if (allocated(p_res_y)) then
       deallocate(p_res_y)
    end if
    if (allocated(p_res_z)) then
       deallocate(p_res_z)
    end if
    
    call this%scratch%relinquish_field(temp_indices)
  end subroutine fluid_vol_flow_compute_str


!-----------------------------STR volflow------------------------------------------------------------

   subroutine fluid_vol_flow_str(this, u, v, w, p, u_res, v_res, w_res, p_res, &
       c_Xh, gs_Xh, ext_bdf, rho, mu, dt, &
       bclst_dp, bclst_du, bclst_dv, bclst_dw, bclst_vel_res, &
       Ax_vel, Ax_prs, ksp_prs, ksp_vel, pc_prs, pc_vel, prs_max_iter, vel_max_iter)

    class(fluid_volflow_t), intent(inout) :: this
    type(field_t), intent(inout) :: u, v, w, p
    type(field_t), intent(inout) :: u_res, v_res, w_res, p_res
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    type(time_scheme_controller_t), intent(in) :: ext_bdf
    real(kind=rp), intent(in) :: rho, mu, dt
    type(bc_list_t), intent(inout) :: bclst_dp, bclst_du, bclst_dv, bclst_dw
    type(bc_list_t), intent(inout) :: bclst_vel_res
    class(ax_t), intent(in) :: Ax_vel
    class(ax_t), intent(in) :: Ax_prs
    class(ksp_t), intent(inout) :: ksp_prs, ksp_vel
    class(pc_t), intent(inout) :: pc_prs, pc_vel
    integer, intent(in) :: prs_max_iter, vel_max_iter
    real(kind=rp) :: ifcomp, flow_rate, xsec, target_ubar, base_ubar
    real(kind=rp) :: current_flow, delta_flow, scale, current_ubar
    integer :: n, ierr, i

    associate(u_vol => this%u_vol, v_vol => this%v_vol, &
         w_vol => this%w_vol, p_vol => this%p_vol)

      n = c_Xh%dof%size() !ntot1

      ! If either dt or the backwards difference coefficient change,
      ! then recompute base flow solution corresponding to unit forcing:

      ifcomp = 0.0_rp

      if (dt .ne. this%dtlag .or. &
           ext_bdf%diffusion_coeffs(1) .ne. this%bdlag) then !add other conditions? : ifuservp and ifexplvis
         call this%compute(u_res, v_res, w_res, p_res, &
              ext_bdf, gs_Xh, c_Xh, rho, mu, ext_bdf%diffusion_coeffs(1), dt, &
              bclst_dp, bclst_du, bclst_dv, bclst_dw, bclst_vel_res, &
              Ax_vel, Ax_prs, ksp_prs, ksp_vel, pc_prs, pc_vel, prs_max_iter, &
              vel_max_iter)
      end if
      this%dtlag = dt
      this%bdlag = ext_bdf%diffusion_coeffs(1)

      !Should we keep this
      !call MPI_Allreduce(MPI_IN_PLACE, ifcomp, 1,MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

      !Inside the dt ne nlock
      !DELETED: if (ifcomp .gt. 0d0) then call this%compute(u_res, v_res, w_res, p_res, &
      !DELETED: if (NEKO_BCKND_DEVICE .eq. 1) then
      !DELETED: if (this%avflow) then
      current_ubar = this%get_ubar_str(u%x, v%x, w%x, c_Xh)
      target_ubar = this%flow_rate ! 1.0 !param(56) = 1.0
      delta_flow = target_ubar - current_ubar
      scale = delta_flow / this%base_flow

      if (NEKO_BCKND_DEVICE .eq. 1) then
         !need device implementation      
         call device_add2s2(u%x_d, u_vol%x_d, scale, n)
         call device_add2s2(v%x_d, v_vol%x_d, scale, n)
         call device_add2s2(w%x_d, w_vol%x_d, scale, n)
         call device_add2s2(p%x_d, p_vol%x_d, scale, n)
      else
         do concurrent (i = 1: n)
            u%x(i,1,1,1) = u%x(i,1,1,1) + scale * u_vol%x(i,1,1,1)
            v%x(i,1,1,1) = v%x(i,1,1,1) + scale * v_vol%x(i,1,1,1)
            w%x(i,1,1,1) = w%x(i,1,1,1) + scale * w_vol%x(i,1,1,1)
            p%x(i,1,1,1) = p%x(i,1,1,1) + scale * p_vol%x(i,1,1,1)
         end do
      end if
    end associate

   end subroutine fluid_vol_flow_str
  

!-----------------------------STR helpers---------------------------------------------------
   function get_ubar_str_helix(this, u, v, w, c_Xh) result(ubar)
      implicit none
      class(fluid_volflow_t), intent(inout) :: this
      real(kind=rp), intent(inout) :: u(:,:,:,:),v(:,:,:,:),w(:,:,:,:)
      type(coef_t), intent(inout) :: c_Xh
      real(kind=rp) :: Rc, pitch, delta, phi, pitch_s, pi
      real(kind=rp) :: num, den
      real(kind=rp) :: x, y, z, y0, alpha, angle_t, r, us, usr, ubar, ri
      integer :: n, i, ierr

      pi = 4*atan(1.0)

      n = c_Xh%dof%size()
      num=0.
      den=0.
      Rc = this%Rc
      pitch = this%pitch
      pitch_s = pitch/(2*pi)
      phi   = atan2(pitch_s, Rc)
      delta = Rc/(Rc**2+pitch_s**2)

      do concurrent (i = 1: n)
         x = this%xax%x(i,1,1,1) ! Corresponding toroid when pitch_s = 0 for the helix
         y = this%yax%x(i,1,1,1) 
         z = this%zax%x(i,1,1,1)
         y0 = sqrt(x**2+y**2) - Rc !y of straight pipe before deformation (when x is axial direction). Minus Rc shifts Coordinate system to center of pipe 
         alpha = atan2(y0, z)
         r = sqrt(y0**2+z**2) !r of helical coordinate system (s, r, theta)
         angle_t = atan2(x,y) !phi measured clockwise from 12 noon
         ri = 1./abs(1+delta*r*sin(alpha)) 
         us = cos(phi)*(u(i,1,1,1)*cos(angle_t)-v(i,1,1,1)*sin(angle_t))+w(i,1,1,1)*sin(phi)!dot product (u,v).(cos(phi), -sin(phi))

         usr = us*ri
         num = num + usr*c_Xh%B(i,1,1,1)
         den = den +  ri*c_Xh%B(i,1,1,1)
   
      end do

      call MPI_Allreduce(MPI_IN_PLACE, num, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, den, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

      ubar = num/den ! "1/r"-weighted volumetric average of streamwise velocity

   end function get_ubar_str_helix


   subroutine makebf_str_helix(this, bax,bay,baz, c_Xh)
      implicit none
      !c_Xh is of type coef_t and c_Xh%Xh is space_t
      !xm1,ym1,zm1 can be accessed from c_Xh%dof%x(i,1,1,1) - which is of type dofmap_t
      class(fluid_volflow_t), intent(inout) :: this
      real(kind=rp), intent(inout) :: bax(:,:,:,:), bay(:,:,:,:), baz(:,:,:,:)
      type(coef_t), intent(inout) :: c_Xh
      real(kind=rp) :: Rc, pitch, delta, phi, pitch_s
      integer :: n, i
      real(kind=rp) :: pi, x, y, z, y0, alpha, angle_t, r, dpds 
   
      pi = 4*atan(1.0)
      !required scalars
      Rc    = this%Rc
      pitch = this%pitch

      pitch_s = pitch/(2*pi)
      phi = atan2(pitch_s, Rc)
      delta = Rc/(Rc**2+pitch_s**2)
      !tau = pitch_s/(Rc**2+pitch_s**2) !not needed for calculations
      n = c_Xh%dof%size()
      !if (n .ne. size(bax)) then
      !   write(*,*) "MAKEBF_STR ERROR: Mismatch between sizes: Expected", n, " but got ", size(bax)
      !   stop
      !endif 

      do concurrent (i = 1: n)
         x = this%xax%x(i,1,1,1) ! Corresponding toroid when pitch_s = 0 for the helix
         y = this%yax%x(i,1,1,1) 
         z = this%zax%x(i,1,1,1)
         y0 = sqrt(x**2+y**2) - Rc !gives us y0 - y of straight pipe before deformation (when x is axial direction).  
         alpha = atan2(y0, z)
         r = sqrt(y0**2+z**2) !r of helical coordinate system (s, r, theta)

         dpds = 1./abs(1+delta*r*sin(alpha))   
         angle_t = atan2(x,y) !phi measured clockwise from 12 noon

         bax(i,1,1,1) = dpds*cos(phi)*cos(angle_t)*c_Xh%B(i,1,1,1)
         bay(i,1,1,1) =-dpds*cos(phi)*sin(angle_t)*c_Xh%B(i,1,1,1)
         baz(i,1,1,1) = dpds*sin(phi)*c_Xh%B(i,1,1,1)
      end do
   end subroutine makebf_str_helix

!---------------------toroid helpers---------------------------------------------------------
function get_ubar_azm(this, u,v,w,c_Xh) result(ubar)
      implicit none
      class(fluid_volflow_t), intent(inout) :: this
      real(kind=rp), intent(inout) :: u(:,:,:,:),v(:,:,:,:),w(:,:,:,:)
      type(coef_t), intent(inout) :: c_Xh
      real(kind=rp) :: num, den
      real(kind=rp) :: x, y, z, rr, ri, phi, us, usr, ubar
      integer :: n, i, ierr

      n = c_Xh%dof%size()
      num=0.
      den=0.

      do concurrent (i = 1: n)
         x = c_Xh%dof%x(i,1,1,1)
         y = c_Xh%dof%y(i,1,1,1)
         z = c_Xh%dof%z(i,1,1,1)

         rr = x*x+y*y
         ri = 1./sqrt(rr) 

         phi = atan2(x,y)
         us = u(i,1,1,1)*cos(phi)-v(i,1,1,1)*sin(phi) !dot product (u,v).(cos(phi), -sin(phi))
         usr = us*ri
         num = num + usr*c_Xh%B(i,1,1,1)
         den = den +  ri*c_Xh%B(i,1,1,1)
   
      end do
      call MPI_Allreduce(MPI_IN_PLACE, num, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, den, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      !num = glsum(num,1)
      !den = glsum(den,1)
      ubar = num/den ! "1/r"-weighted volumetric average of streamwise velocity

end function get_ubar_azm

   subroutine makebf_azm(this, bax,bay,baz,c_Xh)
      implicit none
      class(fluid_volflow_t), intent(inout) :: this
      !c_Xh is of type coef_t and c_Xh%Xh is space_t
      !xm1,ym1,zm1 can be accessed from c_Xh%dof%x(i,1,1,1) - which of type dofmap_t
      real(kind=rp), intent(inout) :: bax(:,:,:,:), bay(:,:,:,:), baz(:,:,:,:)
      type(coef_t), intent(inout) :: c_Xh
      integer :: n, i
      real(kind=rp) :: r0, x, y, z, rr, r, dpds, phi
      
      n = c_Xh%dof%size()
      !if (n .ne. size(bax)) then
      !   write(*,*) "MAKEBF_STR ERROR: Mismatch between sizes: Expected", n, " but got ", size(bax)
      !   stop
      !endif 

      r0 = glmax(c_Xh%dof%y, n)-1.0 !Major radius (minor radius is 1.0)
      do concurrent (i = 1: n)
         x = c_Xh%dof%x(i,1,1,1)
         y = c_Xh%dof%y(i,1,1,1)
         z = c_Xh%dof%z(i,1,1,1)
         rr = x*x+y*y
         r = sqrt(rr) 

         dpds = r0/r      !F/r forcing
         phi = atan2(x,y) !phi measured clockwise from 12 noon

         bax(i,1,1,1) = dpds*cos(phi)*c_Xh%B(i,1,1,1)
         bay(i,1,1,1) =-dpds*sin(phi)*c_Xh%B(i,1,1,1)
         baz(i,1,1,1) = 0.
      end do
   end subroutine makebf_azm
end module fluid_volflow


