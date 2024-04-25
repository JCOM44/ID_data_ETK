! evolve_utils.F90: Misc stuff that is needed for evolution
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine JordanFBSSN_remove_trA( CCTK_ARGUMENTS )

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL                hh(3,3), aa(3,3), hu(3,3), tra, dethh
  CCTK_INT                 i, j, k, m, n


  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE(i, j, k, m, n,  &
  !$OMP hh, aa, hu, tra, dethh)
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           hh(1,1) = hxx(i,j,k)
           hh(1,2) = hxy(i,j,k)
           hh(1,3) = hxz(i,j,k)
           hh(2,2) = hyy(i,j,k)
           hh(2,3) = hyz(i,j,k)
           hh(3,3) = hzz(i,j,k)
           hh(2,1) = hh(1,2)
           hh(3,1) = hh(1,3)
           hh(3,2) = hh(2,3)

           aa(1,1) = axx(i,j,k)
           aa(1,2) = axy(i,j,k)
           aa(1,3) = axz(i,j,k)
           aa(2,2) = ayy(i,j,k)
           aa(2,3) = ayz(i,j,k)
           aa(3,3) = azz(i,j,k)
           aa(2,1) = aa(1,2)
           aa(3,1) = aa(1,3)
           aa(3,2) = aa(2,3)

           dethh =    hh(1,1) * hh(2,2) * hh(3,3)                            &
                + 2 * hh(1,2) * hh(1,3) * hh(2,3)                            &
                -     hh(1,1) * hh(2,3) ** 2                                 &
                -     hh(2,2) * hh(1,3) ** 2                                 &
                -     hh(3,3) * hh(1,2) ** 2
           hu(1,1) = (hh(2,2) * hh(3,3) - hh(2,3) ** 2     ) / dethh
           hu(2,2) = (hh(1,1) * hh(3,3) - hh(1,3) ** 2     ) / dethh
           hu(3,3) = (hh(1,1) * hh(2,2) - hh(1,2) ** 2     ) / dethh
           hu(1,2) = (hh(1,3) * hh(2,3) - hh(1,2) * hh(3,3)) / dethh
           hu(1,3) = (hh(1,2) * hh(2,3) - hh(1,3) * hh(2,2)) / dethh
           hu(2,3) = (hh(1,3) * hh(1,2) - hh(2,3) * hh(1,1)) / dethh
           hu(2,1) = hu(1,2)
           hu(3,1) = hu(1,3)
           hu(3,2) = hu(2,3)

           tra = 0
           do m = 1, 3
              do n = 1, 3
                 tra = tra + hu(m,n) * aa(m,n)
              end do
           end do

           aa = aa - hh * tra / 3

           axx(i,j,k) = aa(1,1)
           axy(i,j,k) = aa(1,2)
           axz(i,j,k) = aa(1,3)
           ayy(i,j,k) = aa(2,2)
           ayz(i,j,k) = aa(2,3)
           azz(i,j,k) = aa(3,3)

        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine JordanFBSSN_remove_trA
!
!===========================================================================
!
subroutine JordanFBSSN_reset_detmetric( CCTK_ARGUMENTS )

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL                hh(3,3), dethh
  CCTK_INT                 i, j, k

  !$OMP PARALLEL DO COLLAPSE(3)         &
  !$OMP PRIVATE(i, j, k, hh, dethh )
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           hh(1,1) = hxx(i,j,k)
           hh(1,2) = hxy(i,j,k)
           hh(1,3) = hxz(i,j,k)
           hh(2,2) = hyy(i,j,k)
           hh(2,3) = hyz(i,j,k)
           hh(3,3) = hzz(i,j,k)
           hh(2,1) = hh(1,2)
           hh(3,1) = hh(1,3)
           hh(3,2) = hh(2,3)

           dethh =       hh(1,1) * hh(2,2) * hh(3,3)                            &
                   + 2 * hh(1,2) * hh(1,3) * hh(2,3)                            &
                   -     hh(1,1) * hh(2,3) ** 2                                 &
                   -     hh(2,2) * hh(1,3) ** 2                                 &
                   -     hh(3,3) * hh(1,2) ** 2

           hh = hh / dethh**(1.0d0/3.0d0)

           hxx(i,j,k) = hh(1,1)
           hxy(i,j,k) = hh(1,2)
           hxz(i,j,k) = hh(1,3)
           hyy(i,j,k) = hh(2,2)
           hyz(i,j,k) = hh(2,3)
           hzz(i,j,k) = hh(3,3)

        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine JordanFBSSN_reset_detmetric
!
!===========================================================================
!
subroutine JordanFBSSN_impose_conf_fac_floor( CCTK_ARGUMENTS )

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT                 i, j, k

  !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i, j, k)
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           if( conf_fac(i,j,k) < conf_fac_floor ) then
               conf_fac(i,j,k) = conf_fac_floor
           end if

        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine JordanFBSSN_impose_conf_fac_floor
!
!==========================================================================
!
subroutine JordanFBSSN_apply_jacobian(dvar, jac)
  implicit none

  CCTK_REAL, intent(inout) :: dvar(3)
  CCTK_REAL, intent(in)    :: jac(3,3)
  CCTK_REAL                :: xdvar(3)
  CCTK_INT                 :: a, b

  xdvar = 0
  do a = 1, 3
     do b = 1, 3
        xdvar(a) = xdvar(a) + dvar(b) * jac(b,a)
     end do
  end do

  dvar = xdvar

end subroutine JordanFBSSN_apply_jacobian
!
!==========================================================================
!
subroutine JordanFBSSN_apply_jacobian2(dvar, ddvar, jac, hes)
  implicit none

  CCTK_REAL, intent(inout) :: ddvar(3,3), dvar(3)
  CCTK_REAL, intent(in)    :: jac(3,3), hes(3,3,3)
  CCTK_REAL                :: xddvar(3,3), xdvar(3)
  CCTK_INT                 :: a, b, c, d

  xdvar = 0
  do a = 1, 3
     do b = 1, 3
        xdvar(a) = xdvar(a) + dvar(b) * jac(b,a)
     end do
  end do

  xddvar = 0
  do a = 1, 3
     do b = 1, 3
        do c = 1, 3
           xddvar(a,b) = xddvar(a,b) + dvar(c) * hes(c,a,b)
           do d = 1, 3
              xddvar(a,b) = xddvar(a,b) + ddvar(c,d) * jac(c,a) * jac(d,b)
           end do
        end do
     end do
  end do

  dvar  = xdvar
  ddvar = xddvar

end subroutine JordanFBSSN_apply_jacobian2
