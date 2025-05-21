! Copyright (c) 2021-2024, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Operators CPU backend
module opr_cpu
  use num_types, only : rp, dp, xp
  use space, only : space_t
  use coefs, only : coef_t
  use math, only : sub3, copy, rzero, PI
  use field, only : field_t
  use gather_scatter, only : GS_OP_ADD
  use interpolation, only : interpolator_t
  use mathops, only : opcolv
  use utils, only : linear_index, index_is_on_facet, neko_error !cyclic_mod_other
  implicit none
  private

  public :: opr_cpu_dudxyz, opr_cpu_opgrad, opr_cpu_cdtp, &
       opr_cpu_conv1, opr_cpu_curl, opr_cpu_cfl, opr_cpu_lambda2, &
       opr_cpu_convect_scalar, opr_cpu_set_convect_rst, &
       opr_cpu_rotate_cyc_r4, opr_cpu_rotate_cyc_vector 
       !opr_cpu_rotate_cyc_loc, set_periodic_idx

  interface
     module subroutine opr_cpu_dudxyz(du, u, dr, ds, dt, coef)
       type(coef_t), intent(in), target :: coef
       real(kind=rp), intent(inout), &
            dimension(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv) :: du
       real(kind=rp), intent(in), &
            dimension(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv) :: &
            u, dr, ds, dt
     end subroutine opr_cpu_dudxyz

     module subroutine opr_cpu_opgrad(ux, uy, uz, u, coef, e_start, e_end)
       type(coef_t), intent(in) :: coef
       integer, intent(in) :: e_start, e_end
       real(kind=rp), intent(inout) :: ux(coef%Xh%lxyz, e_end - e_start + 1)
       real(kind=rp), intent(inout) :: uy(coef%Xh%lxyz, e_end - e_start + 1)
       real(kind=rp), intent(inout) :: uz(coef%Xh%lxyz, e_end - e_start + 1)
       real(kind=rp), intent(in) :: u(coef%Xh%lxyz, e_end - e_start + 1)
     end subroutine opr_cpu_opgrad

     module subroutine opr_cpu_cdtp(dtx, x, dr, ds, dt, coef, e_start, e_end)
       type(coef_t), intent(in) :: coef
       integer, intent(in) :: e_start, e_end
       real(kind=rp), intent(inout) :: dtx(coef%Xh%lxyz, e_end - e_start + 1)
       real(kind=rp), intent(inout) :: x(coef%Xh%lxyz, e_end - e_start + 1)
       real(kind=rp), intent(in) :: dr(coef%Xh%lxyz, e_end - e_start + 1)
       real(kind=rp), intent(in) :: ds(coef%Xh%lxyz, e_end - e_start + 1)
       real(kind=rp), intent(in) :: dt(coef%Xh%lxyz, e_end - e_start + 1)
     end subroutine opr_cpu_cdtp

     module subroutine opr_cpu_conv1(du, u, vx, vy, vz, Xh, &
          coef, e_start, e_end)
       type(space_t), intent(in) :: Xh
       type(coef_t), intent(in) :: coef
       integer, intent(in) :: e_start, e_end
       real(kind=rp), intent(inout) :: du(Xh%lxyz, e_end - e_start + 1)
       real(kind=rp), intent(inout) :: &
            u(Xh%lx, Xh%ly, Xh%lz, e_end - e_start + 1)
       real(kind=rp), intent(inout) :: &
            vx(Xh%lx, Xh%ly, Xh%lz, e_end - e_start + 1)
       real(kind=rp), intent(inout) :: &
            vy(Xh%lx, Xh%ly, Xh%lz, e_end - e_start + 1)
       real(kind=rp), intent(inout) :: &
            vz(Xh%lx, Xh%ly, Xh%lz, e_end - e_start + 1)
     end subroutine opr_cpu_conv1

     module subroutine opr_cpu_convect_scalar(du, u, c, Xh_GLL, Xh_GL, &
                                              coef_GLL, coef_GL, GLL_to_GL)
        type(space_t), intent(in) :: Xh_GL
        type(space_t), intent(in) :: Xh_GLL
        type(coef_t), intent(in) :: coef_GLL
        type(coef_t), intent(in) :: coef_GL
        type(interpolator_t), intent(inout) :: GLL_to_GL
        real(kind=rp), intent(inout) :: &
                   du(Xh_GLL%lx, Xh_GLL%ly, Xh_GLL%lz, coef_GL%msh%nelv)
        real(kind=rp), intent(inout) :: &
                   u(Xh_GL%lx, Xh_GL%lx, Xh_GL%lx, coef_GL%msh%nelv)
        real(kind=rp), intent(inout) :: c(Xh_GL%lxyz, coef_GL%msh%nelv, 3)

      end subroutine opr_cpu_convect_scalar

      module subroutine opr_cpu_set_convect_rst(cr, cs, ct, cx, cy, cz, &
                                                Xh, coef)
         type(space_t), intent(inout) :: Xh
         type(coef_t), intent(inout) :: coef
         real(kind=rp), dimension(Xh%lxyz, coef%msh%nelv), &
                        intent(inout) :: cr, cs, ct
         real(kind=rp), dimension(Xh%lxyz, coef%msh%nelv), &
                        intent(in) :: cx, cy, cz
       end subroutine opr_cpu_set_convect_rst
  end interface

contains

  subroutine opr_cpu_curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)
    type(field_t), intent(inout) :: w1
    type(field_t), intent(inout) :: w2
    type(field_t), intent(inout) :: w3
    type(field_t), intent(in) :: u1
    type(field_t), intent(in) :: u2
    type(field_t), intent(in) :: u3
    type(field_t), intent(inout) :: work1
    type(field_t), intent(inout) :: work2
    type(coef_t), intent(in) :: c_Xh
    integer :: gdim, n

    n = w1%dof%size()
    gdim = c_Xh%msh%gdim

    !     this%work1=dw/dy ; this%work2=dv/dz
    call opr_cpu_dudxyz(work1%x, u3%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    if (gdim .eq. 3) then
       call opr_cpu_dudxyz(work2%x, u2%x, c_Xh%drdz, c_Xh%dsdz, &
                           c_Xh%dtdz, c_Xh)
       call sub3(w1%x, work1%x, work2%x, n)
    else
       call copy(w1%x, work1%x, n)
    end if
    !     this%work1=du/dz ; this%work2=dw/dx
    if (gdim .eq. 3) then
       call opr_cpu_dudxyz(work1%x, u1%x, c_Xh%drdz, c_Xh%dsdz, &
                           c_Xh%dtdz, c_Xh)
       call opr_cpu_dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, &
                           c_Xh%dtdx, c_Xh)
       call sub3(w2%x, work1%x, work2%x, n)
    else
       call rzero(work1%x, n)
       call opr_cpu_dudxyz(work2%x, u3%x, c_Xh%drdx, c_Xh%dsdx, &
                           c_Xh%dtdx, c_Xh)
       call sub3(w2%x, work1%x, work2%x, n)
    end if
    !     this%work1=dv/dx ; this%work2=du/dy
    call opr_cpu_dudxyz(work1%x, u2%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call opr_cpu_dudxyz(work2%x, u1%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call sub3(w3%x, work1%x, work2%x, n)
    !!    BC dependent, Needs to change if cyclic

    call opcolv(w1%x, w2%x, w3%x, c_Xh%B, gdim, n)
    !cyclic_mod_checked - this is opdssum from plan4.f/op_curl
    !write(*, *) 'Rotate from opr_cpu_curl'
    call opr_cpu_rotate_cyc_r4(w1%x, w2%x, w3%x, 1, c_Xh)
    call c_Xh%gs_h%op(w1, GS_OP_ADD)
    call c_Xh%gs_h%op(w2, GS_OP_ADD)
    call c_Xh%gs_h%op(w3, GS_OP_ADD)
    call opr_cpu_rotate_cyc_r4(w1%x, w2%x, w3%x, 0, c_Xh)
    call opcolv(w1%x, w2%x, w3%x, c_Xh%Binv, gdim, n)

  end subroutine opr_cpu_curl

  function opr_cpu_cfl(dt, u, v, w, Xh, coef, nelv, gdim) result(cfl)
    type(space_t) :: Xh
    type(coef_t) :: coef
    integer :: nelv, gdim
    real(kind=rp) :: dt
    real(kind=rp), dimension(Xh%lx, Xh%ly, Xh%lz, nelv) :: u, v, w
    real(kind=rp) :: cflr, cfls, cflt, cflm
    real(kind=rp) :: ur, us, ut
    real(kind=rp) :: cfl
    integer :: i, j, k, e
    cfl = 0d0
    if (gdim .eq. 3) then
       do e = 1, nelv
          do k = 1, Xh%lz
             do j = 1, Xh%ly
                do i = 1, Xh%lx
                   ur = ( u(i,j,k,e)*coef%drdx(i,j,k,e) &
                        + v(i,j,k,e)*coef%drdy(i,j,k,e) &
                        + w(i,j,k,e)*coef%drdz(i,j,k,e) ) &
                        * coef%jacinv(i,j,k,e)
                   us = ( u(i,j,k,e)*coef%dsdx(i,j,k,e) &
                        + v(i,j,k,e)*coef%dsdy(i,j,k,e) &
                        + w(i,j,k,e)*coef%dsdz(i,j,k,e) ) &
                        * coef%jacinv(i,j,k,e)
                   ut = ( u(i,j,k,e)*coef%dtdx(i,j,k,e) &
                        + v(i,j,k,e)*coef%dtdy(i,j,k,e) &
                        + w(i,j,k,e)*coef%dtdz(i,j,k,e) ) &
                        * coef%jacinv(i,j,k,e)

                   cflr = abs(dt*ur*Xh%dr_inv(i))
                   cfls = abs(dt*us*Xh%ds_inv(j))
                   cflt = abs(dt*ut*Xh%dt_inv(k))

                   cflm = cflr + cfls + cflt
                   cfl = max(cfl, cflm)
                end do
             end do
          end do
       end do
    else
       do e = 1, nelv
          do j = 1, Xh%ly
             do i = 1, Xh%lx
                ur = ( u(i,j,1,e)*coef%drdx(i,j,1,e) &
                     + v(i,j,1,e)*coef%drdy(i,j,1,e) ) * coef%jacinv(i,j,1,e)
                us = ( u(i,j,1,e)*coef%dsdx(i,j,1,e) &
                     + v(i,j,1,e)*coef%dsdy(i,j,1,e) ) * coef%jacinv(i,j,1,e)

                cflr = abs(dt*ur*Xh%dr_inv(i))
                cfls = abs(dt*us*Xh%ds_inv(j))

                cflm = cflr + cfls
                cfl = max(cfl, cflm)

             end do
          end do
       end do
    end if
  end function opr_cpu_cfl

  subroutine opr_cpu_lambda2(lambda2, u, v, w, coef)
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: lambda2
    type(field_t), intent(in) :: u, v, w
    real(kind=rp) :: grad(coef%Xh%lxyz,3,3)
    integer :: e, i
    real(kind=xp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=xp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=xp) :: a11, a22, a33, a12, a13, a23
    real(kind=xp) :: msk1, msk2, msk3

    do e = 1, coef%msh%nelv
       call opr_cpu_opgrad(grad(1,1,1), grad(1,1,2), grad(1,1,3), &
            u%x(1,1,1,e), coef,e,e)
       call opr_cpu_opgrad(grad(1,2,1), grad(1,2,2), grad(1,2,3), &
            v%x(1,1,1,e), coef,e,e)
       call opr_cpu_opgrad(grad(1,3,1), grad(1,3,2), grad(1,3,3), &
            w%x(1,1,1,e), coef,e,e)

       do i = 1, coef%Xh%lxyz
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)


          s12 = 0.5_xp*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5_xp*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5_xp*(grad(i,2,3) + grad(i,3,2))

          o12 = 0.5_xp*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5_xp*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5_xp*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12 + s12 * s22 + s13 * s23 - o13 * o23
          a13 = s11 * s13 + s12 * s23 + s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
               - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0_xp * a12 * a13 * a23 - a11 * a23*a23 &
               - a22 * a13*a13 - a33 * a12*a12 + a11 * a22 * a33)


          q = (3.0_xp * C - B*B) / 9.0_xp
          r = (9.0_xp * C * B - 27.0_xp * D - 2.0_xp * B*B*B) / 54.0_xp
          theta = acos( r / sqrt(-q*q*q) )

          eigen(1) = 2.0_xp * sqrt(-q) * cos(theta / 3.0_xp) - B / 3.0_xp
          eigen(2) = 2.0_xp * sqrt(-q) * cos((theta + 2.0_xp * pi) / 3.0_xp) - B / 3.0_xp
          eigen(3) = 2.0_xp * sqrt(-q) * cos((theta + 4.0_xp * pi) / 3.0_xp) - B / 3.0_xp
          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
               .and. eigen(1) .le. eigen(3) .or. eigen(3) &
               .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
               .and. eigen(2) .le. eigen(3) .or. eigen(3) &
               .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
               .and. eigen(3) .le. eigen(2) .or. eigen(2) &
               .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)
          lambda2%x(i,1,1,e) = l2/(real(coef%B(i,1,1,e)**2,xp))
       end do
    end do

  end subroutine opr_cpu_lambda2


  subroutine opr_cpu_rotate_cyc_vector(rx, ry, rz, idir, coef)
    use iso_c_binding
    real(kind=rp), dimension(:), intent(inout) :: rx, ry, rz
    integer, intent(in) :: idir
    type(coef_t), intent(in) :: coef
    integer :: i, ntot
    real(kind=rp) :: un(3), x, y, z, theta, dotprod, cost, sint, rnor, rtn1,tol
    !integer :: nel, lx, ly, lz, pf, pe, j, k, lidx, nidx
    !nel =  coef%msh%nelv
    !lx = coef%Xh%lx
    !ly = coef%Xh%ly
    !lz = coef%Xh%lz
    tol = 1.0e-10_rp
    ntot = coef%dof%size()

    do i = 1, ntot
          theta = coef%cyc_angle(i, 1, 1, 1)
          if(abs(theta) .gt. tol) then
               x = coef%dof%x(i, 1, 1, 1)
               y = coef%dof%y(i, 1, 1, 1)
               z = coef%dof%z(i, 1, 1, 1)
               cost = cos(theta)
               sint = sin(theta)
               dotprod = cost*y-sint*x

               !do forward (idir=1) or backward rotation (idir=0)
               if (idir.eq.1) then
                    rnor = rx(i)*cost+ry(i)*sint
                    rtn1 =-rx(i)*sint+ry(i)*cost
               else if(idir.eq.0) then
                    rnor = rx(i)*cost-ry(i)*sint
                    rtn1 = rx(i)*sint+ry(i)*cost
               else
                    write(*, *) "idir must be 0 or 1 in rotate_cyc"
                    !stop !must be removed if do concurrent is used.
               end if

               !invert direction if vector comes out of plane (i.e. dotprod<0)
               if (dotprod.ge.0) then
                    rx(i) = rnor
                    ry(i) = rtn1
               else
                    rx(i) =-rnor
                    ry(i) =-rtn1
               end if
          end if
    end do

  end subroutine opr_cpu_rotate_cyc_vector


  subroutine opr_cpu_rotate_cyc_r4(rx, ry, rz, idir, coef)
    use iso_c_binding
    real(kind=rp), dimension(:,:,:,:), intent(inout) :: rx, ry, rz
    integer, intent(in) ::  idir
    type(coef_t), intent(in) :: coef
    integer :: i, ntot
    real(kind=rp) :: un(3), x, y, z, theta, dotprod, cost, sint, rnor, rtn1, tol
    !integer :: nel, lx, ly, lz, j, k, lidx, nidx, pf, pe, nel
    !nel =  coef%msh%nelv
    !lx = coef%Xh%lx
    !ly = coef%Xh%ly
    !lz = coef%Xh%lz
    tol = 1.0e-15_rp
    ntot = coef%dof%size()

    do i = 1, ntot
          theta = coef%cyc_angle(i, 1, 1, 1)
          if(abs(theta) .gt. tol) then
               x = coef%dof%x(i, 1, 1, 1)
               y = coef%dof%y(i, 1, 1, 1)
               z = coef%dof%z(i, 1, 1, 1)
               cost = cos(theta)
               sint = sin(theta)
               dotprod = cost*y-sint*x
               !can avoid this if block using idir=1 and -1 mutliplied by sint term
               if (idir.eq.1) then
                    rnor = rx(i, 1, 1, 1)*cost+ry(i, 1, 1, 1)*sint
                    rtn1 =-rx(i, 1, 1, 1)*sint+ry(i, 1, 1, 1)*cost
               else if(idir.eq.0) then
                    rnor = rx(i, 1, 1, 1)*cost-ry(i, 1, 1, 1)*sint
                    rtn1 = rx(i, 1, 1, 1)*sint+ry(i, 1, 1, 1)*cost
               else
                    write(*, *) "idir must be 0 or 1 in rotate_cyc"
                    !stop !must be removed if do concurrent is used.
               end if

               !invert direction if vector comes out of plane (i.e. dotprod<0)
               if (dotprod.ge.0) then
                    rx(i, 1, 1, 1) = rnor
                    ry(i, 1, 1, 1) = rtn1
               else
                    rx(i, 1, 1, 1) =-rnor
                    ry(i, 1, 1, 1) =-rtn1
               end if
          end if 
     end do

  end subroutine opr_cpu_rotate_cyc_r4

!   subroutine opr_cpu_rotate_cyc_loc(rx, ry, rz, nsize, idir, c_Xh)
!     use iso_c_binding
!     real(kind=rp), intent(inout) :: rx(nsize), ry(nsize), rz(nsize)
!     integer, intent(in) :: nsize, idir
!     type(coef_t), intent(in) :: c_Xh
!     integer :: np, n, lx, ly, lz, pf, pe, i, j, k, lidx, nidx
!     real(kind=rp) :: un(3), x, y, z, dotprod, length, cost, sint, rnor, rtn1

!     np =  c_Xh%msh%periodic%size
!     lx = c_Xh%Xh%lx
!     ly = c_Xh%Xh%ly
!     lz = c_Xh%Xh%lz

!     do n = 1, np
!         do k = 1, lz
!         do j = 1, ly
!         do i = 1, lx
!           !pf = c_Xh%msh%periodic%p_facet_el(n)%x(1) !doesnt work parallel
!           !pe = c_Xh%msh%periodic%p_facet_el(n)%x(2) !doesnt work parallel
!           pf = c_Xh%msh%periodic%facet_el(n)%x(1)
!           pe = c_Xh%msh%periodic%facet_el(n)%x(2)
!           if (index_is_on_facet(i, j, k, lx, ly, lz, pf)) then
!                lidx = linear_index(i, j, k, pe, lx, ly, lz)
!                select case (pf)
!                case(1,2)
!                     nidx = j+(k-1)*ly+(pf-1)*ly*lz+(pe-1)*ly*lz*6
!                case(3,4)
!                     nidx = i+(k-1)*lx+(pf-1)*lx*lz+(pe-1)*lx*lz*6
!                case(5,6)
!                     nidx = i+(j-1)*lx+(pf-1)*lx*ly+(pe-1)*lx*ly*6
!                end select

!                un(1) = c_Xh%nx(nidx,1,1,1)
!                un(2) = c_Xh%ny(nidx,1,1,1)
!                x = c_Xh%dof%x(lidx, 1, 1, 1)
!                y = c_Xh%dof%y(lidx, 1, 1, 1)
!                z = c_Xh%dof%z(lidx, 1, 1, 1)
!                dotprod = un(1)*y-un(2)*x

!                length = sqrt(un(1)*un(1)+un(2)*un(2))
!                cost = un(1)/length
!                sint = un(2)/length

!                !do forward (idir=1) or backward rotation (idir=0)
!                if (idir.eq.1) then
!                     rnor = rx(lidx)*cost+ry(lidx)*sint
!                     rtn1 =-rx(lidx)*sint+ry(lidx)*cost
!                else if(idir.eq.0) then
!                     rnor = rx(lidx)*cost-ry(lidx)*sint
!                     rtn1 = rx(lidx)*sint+ry(lidx)*cost
!                else
!                     write(*, *) "idir must be 0 or 1 in rotate_cyc"
!                     !stop !must be removed if do concurrent is used.
!                end if

!                !invert direction if vector comes out of plane (i.e. dotprod<0)
!                if (dotprod.ge.0) then
!                     rx(lidx) = rnor
!                     ry(lidx) = rtn1
!                else
!                     rx(lidx) =-rnor
!                     ry(lidx) =-rtn1
!                end if
!           end if

!         end do
!         end do
!         end do
!     end do

!   end subroutine opr_cpu_rotate_cyc_loc




!   subroutine opr_cpu_rotate_cyc_vector(rx, ry, rz, idir, c_Xh)
!     real(kind=rp), dimension(:), intent(inout) :: rx, ry, rz
!     integer, intent(in) :: idir
!     type(coef_t), intent(in) :: c_Xh
!     real(kind=rp), allocatable :: periodic_idx(:), periodic_normal_idx(:)
!     integer :: np, n, lx, ly, lz, pf, pe, i, j, k, lidx, nidx
!     real(kind=rp) :: un(3), x, y, z, dotprod, length, cost, sint, rnor, rtn1

!     np =  c_Xh%msh%periodic%size
!     lx = c_Xh%Xh%lx
!     ly = c_Xh%Xh%ly
!     lz = c_Xh%Xh%lz
!     call set_periodic_idx(c_Xh, periodic_idx, periodic_normal_idx)

!     do n = 1, np*lx*lx
!           nidx = periodic_normal_idx(n)
!           un(1) = c_Xh%nx(nidx,1,1,1)
!           un(2) = c_Xh%ny(nidx,1,1,1)

!           lidx = periodic_idx(n)
!           x = c_Xh%dof%x(lidx, 1, 1, 1)
!           y = c_Xh%dof%y(lidx, 1, 1, 1)
!           z = c_Xh%dof%z(lidx, 1, 1, 1)

!           dotprod = un(1)*y-un(2)*x
!           length = sqrt(un(1)*un(1)+un(2)*un(2))
!           cost = un(1)/length
!           sint = un(2)/length

!           !do forward (idir=1) or backward rotation (idir=0)
!           if (idir.eq.1) then
!                rnor = rx(lidx)*cost+ry(lidx)*sint
!                rtn1 =-rx(lidx)*sint+ry(lidx)*cost
!           else if(idir.eq.0) then
!                rnor = rx(lidx)*cost-ry(lidx)*sint
!                rtn1 = rx(lidx)*sint+ry(lidx)*cost
!           else
!                write(*, *) "idir must be 0 or 1 in rotate_cyc"
!                !stop !must be removed if do concurrent is used.
!           end if

!           !invert direction if vector comes out of plane (i.e. dotprod<0)
!           if (dotprod.ge.0) then
!                rx(lidx) = rnor
!                ry(lidx) = rtn1
!           else
!                rx(lidx) =-rnor
!                ry(lidx) =-rtn1
!           end if
!     end do

!     end subroutine opr_cpu_rotate_cyc_vector

!     subroutine opr_cpu_rotate_cyc_r4(rx, ry, rz, idir, c_Xh)
!     real(kind=rp), dimension(:,:,:,:), intent(inout) :: rx, ry, rz
!     integer, intent(in) :: idir
!     type(coef_t), intent(in) :: c_Xh
!     real(kind=rp), allocatable :: periodic_idx(:), periodic_normal_idx(:)
!     integer :: np, n, lx, ly, lz, pf, pe, i, j, k, pi, pj, pk, lidx, nidx
!     real(kind=rp) :: un(3), x, y, z, dotprod, length, cost, sint, rnor, rtn1
!     real(kind=rp) :: periodic_idx_

!     np =  c_Xh%msh%periodic%size
!     lx = c_Xh%Xh%lx
!     ly = c_Xh%Xh%ly
!     lz = c_Xh%Xh%lz
!     call set_periodic_idx(c_Xh, periodic_idx, periodic_normal_idx)

!     do n = 1, np*lx*lx
!           nidx = periodic_normal_idx(n)
!           un(1) = c_Xh%nx(nidx,1,1,1)
!           un(2) = c_Xh%ny(nidx,1,1,1)

!           lidx = periodic_idx(n)
!           x = c_Xh%dof%x(lidx, 1, 1, 1)
!           y = c_Xh%dof%y(lidx, 1, 1, 1)
!           z = c_Xh%dof%z(lidx, 1, 1, 1)

!           dotprod = un(1)*y-un(2)*x
!           length = sqrt(un(1)*un(1)+un(2)*un(2))
!           cost = un(1)/length
!           sint = un(2)/length
!           !do forward (idir=1) or backward rotation (idir=0)
          ! if (idir.eq.1) then
          !      rnor = rx(lidx, 1, 1, 1)*cost+ry(lidx, 1, 1, 1)*sint
          !      rtn1 =-rx(lidx, 1, 1, 1)*sint+ry(lidx, 1, 1, 1)*cost
          ! else if(idir.eq.0) then
          !      rnor = rx(lidx, 1, 1, 1)*cost-ry(lidx, 1, 1, 1)*sint
          !      rtn1 = rx(lidx, 1, 1, 1)*sint+ry(lidx, 1, 1, 1)*cost
          ! else
          !      write(*, *) "idir must be 0 or 1 in rotate_cyc"
          !      !stop !must be removed if do concurrent is used.
          ! end if

          ! !invert direction if vector comes out of plane (i.e. dotprod<0)
          ! if (dotprod.ge.0) then
          !      rx(lidx, 1, 1, 1) = rnor
          !      ry(lidx, 1, 1, 1) = rtn1
          ! else
          !      rx(lidx, 1, 1, 1) =-rnor
          !      ry(lidx, 1, 1, 1) =-rtn1
          ! end if
!     end do

!    end subroutine opr_cpu_rotate_cyc_r4

!   !cyclic_gpu_mod
!   subroutine set_periodic_idx(coef, periodic_idx, periodic_normal_idx)
!    type(coef_t), intent(in) :: coef
!    real(kind=rp), allocatable, intent(out) :: periodic_idx(:), periodic_normal_idx(:)
!    integer :: np, lx,ly, lz, n, i, j, k, m, pf, pe, nidx
!    np =  coef%msh%periodic%size
!    lx = coef%Xh%lx
!    ly = coef%Xh%ly
!    lz = coef%Xh%lz

!    allocate(periodic_idx(np*lx*lx ))
!    allocate(periodic_normal_idx(np*lx*lx ))
!    m = 0

!    do n = 1, np
!       pf = coef%msh%periodic%facet_el(n)%x(1)
!       pe = coef%msh%periodic%facet_el(n)%x(2)
!       do k = 1, lz
!       do j = 1, ly
!       do i = 1, lx   
!          if (index_is_on_facet(i, j, k, lx, ly, lz, pf)) then
!             select case (pf)
!             case(1,2)
!                   nidx = j+(k-1)*ly+(pf-1)*ly*lz+(pe-1)*ly*lz*6
!             case(3,4)
!                   nidx = i+(k-1)*lx+(pf-1)*lx*lz+(pe-1)*lx*lz*6
!             case(5,6)
!                   nidx = i+(j-1)*lx+(pf-1)*lx*ly+(pe-1)*lx*ly*6
!             end select
!             m = m+1
!             periodic_idx(m) = linear_index(i, j, k, pe, lx, ly, lz)
!             periodic_normal_idx(m) = nidx
!          end if
!       end do
!       end do 
!       end do
!    end do

!    if ( m .ne. np*lx*lx ) then
!       call neko_error('Periodic index lengths are not consistent.')
!       !this assumes lx=ly=lz. if the assumption fails, dry run can be done 
!       !to obtain m first, allocate arrays and then come back to actual iterations.
!    end if

!   end subroutine set_periodic_idx

end module opr_cpu
