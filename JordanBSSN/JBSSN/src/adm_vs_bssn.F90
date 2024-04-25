! adm_vs_bssn.F90 : Functions for converting between ADM and BSSN variables
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine JordanFBSSN_adm2bssn( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr

  CCTK_REAL gg(3,3), kk(3,3), alph, beta(3)
  CCTK_REAL gu(3,3), detg
  CCTK_REAL ww, hh(3,3), hu(3,3), trk, aa(3,3), gammat(3)
  CCTK_REAL d1_hh11(3), d1_hh12(3), d1_hh13(3), d1_hh22(3), d1_hh23(3), d1_hh33(3)
  CCTK_REAL d1_hh(3,3,3)
  CCTK_REAL cf1(3,3,3), cf2(3,3,3)
  CCTK_REAL dx(3), deltax12, deltay12, deltaz12
  CCTK_REAL, parameter :: zero = 0
  CCTK_INT  a, b, c, m, n
  CCTK_INT  i, j, k

  ! Jacobian
  CCTK_REAL jac(3,3)

  integer                  istat
  logical                  use_jacobian
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ11, lJ12, lJ13
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ21, lJ22, lJ23
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ31, lJ32, lJ33

  CCTK_POINTER             lJ11_ptr, lJ12_ptr, lJ13_ptr
  CCTK_POINTER             lJ21_ptr, lJ22_ptr, lJ23_ptr
  CCTK_POINTER             lJ31_ptr, lJ32_ptr, lJ33_ptr

  pointer (lJ11_ptr, lJ11), (lJ12_ptr, lJ12), (lJ13_ptr, lJ13)
  pointer (lJ21_ptr, lJ21), (lJ22_ptr, lJ22), (lJ23_ptr, lJ23)
  pointer (lJ31_ptr, lJ31), (lJ32_ptr, lJ32), (lJ33_ptr, lJ33)

  ! TODO: can this be active but with a cartesian mapping choice?
  call CCTK_IsFunctionAliased(istat, "MultiPatch_GetDomainSpecification")
  if (istat == 0) then
     use_jacobian = .false.
  else
     use_jacobian = .true.
  end if

  if (use_jacobian) then
     call CCTK_VarDataPtr(lJ11_ptr, cctkGH, 0, "Coordinates::J11")
     call CCTK_VarDataPtr(lJ12_ptr, cctkGH, 0, "Coordinates::J12")
     call CCTK_VarDataPtr(lJ13_ptr, cctkGH, 0, "Coordinates::J13")
     call CCTK_VarDataPtr(lJ21_ptr, cctkGH, 0, "Coordinates::J21")
     call CCTK_VarDataPtr(lJ22_ptr, cctkGH, 0, "Coordinates::J22")
     call CCTK_VarDataPtr(lJ23_ptr, cctkGH, 0, "Coordinates::J23")
     call CCTK_VarDataPtr(lJ31_ptr, cctkGH, 0, "Coordinates::J31")
     call CCTK_VarDataPtr(lJ32_ptr, cctkGH, 0, "Coordinates::J32")
     call CCTK_VarDataPtr(lJ33_ptr, cctkGH, 0, "Coordinates::J33")
  end if


  dx(:) = CCTK_DELTA_SPACE(:)
  deltax12 = 12 * dx(1)
  deltay12 = 12 * dx(2)
  deltaz12 = 12 * dx(3)

  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE( i, j, k, &
  !$OMP gg, kk, detg, gu, &
  !$OMP alph, beta, &
  !$OMP trk, ww, hh, aa, &
  !$OMP a, b)
  do k = 1, cctk_lsh(3)
  do j = 1, cctk_lsh(2)
  do i = 1, cctk_lsh(1)

    !------------ Get local vars ---------------
    gg(1,1) = gxx(i,j,k)
    gg(1,2) = gxy(i,j,k)
    gg(1,3) = gxz(i,j,k)
    gg(2,2) = gyy(i,j,k)
    gg(2,3) = gyz(i,j,k)
    gg(3,3) = gzz(i,j,k)
    gg(2,1) = gg(1,2)
    gg(3,1) = gg(1,3)
    gg(3,2) = gg(2,3)

    kk(1,1) = kxx(i,j,k)
    kk(1,2) = kxy(i,j,k)
    kk(1,3) = kxz(i,j,k)
    kk(2,2) = kyy(i,j,k)
    kk(2,3) = kyz(i,j,k)
    kk(3,3) = kzz(i,j,k)
    kk(2,1) = kk(1,2)
    kk(3,1) = kk(1,3)
    kk(3,2) = kk(2,3)

    alph    = alp(i,j,k)

    beta(1) = betax(i,j,k)
    beta(2) = betay(i,j,k)
    beta(3) = betaz(i,j,k)
    

    
    !-------------------------------------------


    !------------- Calculate detg --------------
    detg    =       gg(1,1) * gg(2,2) * gg(3,3)                            &
              + 2 * gg(1,2) * gg(1,3) * gg(2,3)                            &
              -     gg(1,1) * gg(2,3) ** 2                                 &
              -     gg(2,2) * gg(1,3) ** 2                                 &
              -     gg(3,3) * gg(1,2) ** 2
    !-------------------------------------------


    !------------- Invert metric ---------------
    gu(1,1) = (gg(2,2) * gg(3,3) - gg(2,3) ** 2     ) / detg
    gu(2,2) = (gg(1,1) * gg(3,3) - gg(1,3) ** 2     ) / detg
    gu(3,3) = (gg(1,1) * gg(2,2) - gg(1,2) ** 2     ) / detg
    gu(1,2) = (gg(1,3) * gg(2,3) - gg(1,2) * gg(3,3)) / detg
    gu(1,3) = (gg(1,2) * gg(2,3) - gg(1,3) * gg(2,2)) / detg
    gu(2,3) = (gg(1,3) * gg(1,2) - gg(2,3) * gg(1,1)) / detg
    gu(2,1) = gu(1,2)
    gu(3,1) = gu(1,3)
    gu(3,2) = gu(2,3)
    !-------------------------------------------


    !------------ Convert to BSSN --------------
    trk = 0
    do a = 1, 3
      do b = 1, 3
        trk = trk + gu(a,b) * kk(a,b)
      end do
    end do

    ww = detg**(-1.0d0/6.0d0)
    hh = ww*ww * gg
    aa = ww*ww * (kk - trk / 3 * gg)
    !-------------------------------------------


    !------------ Write to grid functions ------
    conf_fac(i,j,k)     = ww
    if( impose_conf_fac_floor_at_initial /= 0 ) then
      if( conf_fac(i,j,k) < conf_fac_floor ) conf_fac(i,j,k) = conf_fac_floor
    end if

    if( precollapsed_lapse /= 0 ) then
      alp(i,j,k)   = ww
    end if

    if( rescale_shift_initial /= 0 ) then
      betax(i,j,k)    = ww * betax(i,j,k)
      betay(i,j,k)    = ww * betay(i,j,k)
      betaz(i,j,k)    = ww * betaz(i,j,k)
    end if


    hxx(i,j,k)     = hh(1,1)
    hxy(i,j,k)     = hh(1,2)
    hxz(i,j,k)     = hh(1,3)
    hyy(i,j,k)     = hh(2,2)
    hyz(i,j,k)     = hh(2,3)
    hzz(i,j,k)     = hh(3,3)

    tracek(i,j,k)  = trk

    axx(i,j,k)     = aa(1,1)
    axy(i,j,k)     = aa(1,2)
    axz(i,j,k)     = aa(1,3)
    ayy(i,j,k)     = aa(2,2)
    ayz(i,j,k)     = aa(2,3)
    azz(i,j,k)     = aa(3,3)

    ! Just to make sure gammat has no uninitialized values anywhere
    gammatx(i,j,k)     = 0
    gammaty(i,j,k)     = 0
    gammatz(i,j,k)     = 0
    !-------------------------------------------

  end do
  end do
  end do


  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE( i, j, k, &
  !$OMP hh, jac, hu, &
  !$OMP d1_hh11, d1_hh12, d1_hh13, d1_hh22, d1_hh23, d1_hh33, d1_hh, &
  !$OMP cf1, cf2, gammat, &
  !$OMP a, b, c, m, n)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

    !------------ Get local variables ----------
    hh(1,1) = hxx(i,j,k)
    hh(1,2) = hxy(i,j,k)
    hh(1,3) = hxz(i,j,k)
    hh(2,2) = hyy(i,j,k)
    hh(2,3) = hyz(i,j,k)
    hh(3,3) = hzz(i,j,k)
    hh(2,1) = hh(1,2)
    hh(3,1) = hh(1,3)
    hh(3,2) = hh(2,3)
    !-------------------------------------------

    if (use_jacobian) then
       jac(1,1) = lJ11(i,j,k)
       jac(1,2) = lJ12(i,j,k)
       jac(1,3) = lJ13(i,j,k)
       jac(2,1) = lJ21(i,j,k)
       jac(2,2) = lJ22(i,j,k)
       jac(2,3) = lJ23(i,j,k)
       jac(3,1) = lJ31(i,j,k)
       jac(3,2) = lJ32(i,j,k)
       jac(3,3) = lJ33(i,j,k)
    else
       jac      = 0.0d0
       jac(1,1) = 1.0d0
       jac(2,2) = 1.0d0
       jac(3,3) = 1.0d0
    end if

    !------------- Invert metric ---------------
    ! NOTE: deth = 1 by construction
    hu(1,1) = (hh(2,2) * hh(3,3) - hh(2,3) ** 2     )
    hu(2,2) = (hh(1,1) * hh(3,3) - hh(1,3) ** 2     )
    hu(3,3) = (hh(1,1) * hh(2,2) - hh(1,2) ** 2     )
    hu(1,2) = (hh(1,3) * hh(2,3) - hh(1,2) * hh(3,3))
    hu(1,3) = (hh(1,2) * hh(2,3) - hh(1,3) * hh(2,2))
    hu(2,3) = (hh(1,3) * hh(1,2) - hh(2,3) * hh(1,1))
    hu(2,1) = hu(1,2)
    hu(3,1) = hu(1,3)
    hu(3,2) = hu(2,3)
    !-------------------------------------------


    !----------- Centered derivatives ----------
    d1_hh11(1)  = (   -hxx(i+2,j,k) + 8*hxx(i+1,j,k)                      &
                    - 8*hxx(i-1,j,k) +   hxx(i-2,j,k) ) / deltax12
    d1_hh12(1)  = (   -hxy(i+2,j,k) + 8*hxy(i+1,j,k)                      &
                    - 8*hxy(i-1,j,k) +   hxy(i-2,j,k) ) / deltax12
    d1_hh13(1)  = (   -hxz(i+2,j,k) + 8*hxz(i+1,j,k)                      &
                    - 8*hxz(i-1,j,k) +   hxz(i-2,j,k) ) / deltax12
    d1_hh22(1)  = (   -hyy(i+2,j,k) + 8*hyy(i+1,j,k)                      &
                    - 8*hyy(i-1,j,k) +   hyy(i-2,j,k) ) / deltax12
    d1_hh23(1)  = (  -hyz(i+2,j,k) + 8*hyz(i+1,j,k)                      &
                    - 8*hyz(i-1,j,k) +   hyz(i-2,j,k) ) / deltax12
    d1_hh33(1)  = (   -hzz(i+2,j,k) + 8*hzz(i+1,j,k)                      &
                    - 8*hzz(i-1,j,k) +   hzz(i-2,j,k) ) / deltax12

    d1_hh11(2) = (   -hxx(i,j+2,k) + 8*hxx(i,j+1,k)                      &
                    - 8*hxx(i,j-1,k) +   hxx(i,j-2,k) ) / deltay12
    d1_hh12(2) = (   -hxy(i,j+2,k) + 8*hxy(i,j+1,k)                      &
                    - 8*hxy(i,j-1,k) +   hxy(i,j-2,k) ) / deltay12
    d1_hh13(2) = (   -hxz(i,j+2,k) + 8*hxz(i,j+1,k)                      &
                    - 8*hxz(i,j-1,k) +   hxz(i,j-2,k) ) / deltay12
    d1_hh22(2) = (   -hyy(i,j+2,k) + 8*hyy(i,j+1,k)                      &
                    - 8*hyy(i,j-1,k) +   hyy(i,j-2,k) ) / deltay12
    d1_hh23(2) = (   -hyz(i,j+2,k) + 8*hyz(i,j+1,k)                      &
                    - 8*hyz(i,j-1,k) +   hyz(i,j-2,k) ) / deltay12
    d1_hh33(2) = (   -hzz(i,j+2,k) + 8*hzz(i,j+1,k)                      &
                    - 8*hzz(i,j-1,k) +   hzz(i,j-2,k) ) / deltay12

    d1_hh11(3) = (   -hxx(i,j,k+2) + 8*hxx(i,j,k+1)                      &
                    - 8*hxx(i,j,k-1) +   hxx(i,j,k-2) ) / deltaz12
    d1_hh12(3) = (   -hxy(i,j,k+2) + 8*hxy(i,j,k+1)                      &
                    - 8*hxy(i,j,k-1) +   hxy(i,j,k-2) ) / deltaz12
    d1_hh13(3) = (   -hxz(i,j,k+2) + 8*hxz(i,j,k+1)                      &
                    - 8*hxz(i,j,k-1) +   hxz(i,j,k-2) ) / deltaz12
    d1_hh22(3) = (   -hyy(i,j,k+2) + 8*hyy(i,j,k+1)                      &
                    - 8*hyy(i,j,k-1) +   hyy(i,j,k-2) ) / deltaz12
    d1_hh23(3) = (   -hyz(i,j,k+2) + 8*hyz(i,j,k+1)                      &
                    - 8*hyz(i,j,k-1) +   hyz(i,j,k-2) ) / deltaz12
    d1_hh33(3) = (   -hzz(i,j,k+2) + 8*hzz(i,j,k+1)                      &
                    - 8*hzz(i,j,k-1) +   hzz(i,j,k-2) ) / deltaz12

   if (use_jacobian) then
      call JordanFBSSN_apply_jacobian(d1_hh11, jac)
      call JordanFBSSN_apply_jacobian(d1_hh12, jac)
      call JordanFBSSN_apply_jacobian(d1_hh13, jac)
      call JordanFBSSN_apply_jacobian(d1_hh22, jac)
      call JordanFBSSN_apply_jacobian(d1_hh23, jac)
      call JordanFBSSN_apply_jacobian(d1_hh33, jac)
   end if


    d1_hh(1,1,:) = d1_hh11(:)
    d1_hh(1,2,:) = d1_hh12(:)
    d1_hh(1,3,:) = d1_hh13(:)
    d1_hh(2,2,:) = d1_hh22(:)
    d1_hh(2,3,:) = d1_hh23(:)
    d1_hh(3,3,:) = d1_hh33(:)
    d1_hh(2,1,:) = d1_hh(1,2,:)
    d1_hh(3,1,:) = d1_hh(1,3,:)
    d1_hh(3,2,:) = d1_hh(2,3,:)


    !-------------------------------------------


    !------------ Connections ------------------
    cf1 = 0
    do a = 1, 3
      do b = 1, 3
        do c = 1, 3
          cf1(a,b,c) = 0.5d0 * (d1_hh(a,b,c) + d1_hh(a,c,b) - d1_hh(b,c,a))
        end do
      end do
    end do

    cf2 = 0
    do a = 1, 3
      do b = 1, 3
        do c = 1, 3
          do m = 1, 3
            cf2(a,b,c) = cf2(a,b,c) + hu(a,m) * cf1(m,b,c)
          end do
        end do
      end do
    end do
    !-------------------------------------------


    !------------ Gammat ------------------------
    gammat = 0
    do a = 1, 3
      do m = 1, 3
        do n = 1, 3
          gammat(a) = gammat(a) + hu(m,n) * cf2(a,m,n)
        end do
      end do
    end do
    !-------------------------------------------


    !------------ Write to gf ------------------
    gammatx(i,j,k) = gammat(1)
    gammaty(i,j,k) = gammat(2)
    gammatz(i,j,k) = gammat(3)
    !-------------------------------------------

  end do
  end do
  end do

  ! now extrapolate Gamma on the outer boundary
  ierr = ExtrapolateGammas(cctkGH, gammatx)
  ierr = ExtrapolateGammas(cctkGH, gammaty)
  ierr = ExtrapolateGammas(cctkGH, gammatz)

end subroutine JordanFBSSN_adm2bssn
!
!===========================================================================
!
subroutine JordanFBSSN_bssn2adm( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  gxx = hxx / (conf_fac * conf_fac)
  gxy = hxy / (conf_fac * conf_fac)
  gxz = hxz / (conf_fac * conf_fac)
  gyy = hyy / (conf_fac * conf_fac)
  gyz = hyz / (conf_fac * conf_fac)
  gzz = hzz / (conf_fac * conf_fac)

  kxx = axx / (conf_fac * conf_fac) + tracek * gxx / 3
  kxy = axy / (conf_fac * conf_fac) + tracek * gxy / 3
  kxz = axz / (conf_fac * conf_fac) + tracek * gxz / 3
  kyy = ayy / (conf_fac * conf_fac) + tracek * gyy / 3
  kyz = ayz / (conf_fac * conf_fac) + tracek * gyz / 3
  kzz = azz / (conf_fac * conf_fac) + tracek * gzz / 3



  if (CCTK_EQUALS(lapse_evolution_method, "JordanFBSSNMoL")) then
     dtalp   = rhs_alp
  end if

  if (CCTK_EQUALS(shift_evolution_method, "JordanFBSSNMoL")) then
     dtbetax = rhs_betax
     dtbetay = rhs_betay
     dtbetaz = rhs_betaz
  end if

end subroutine JordanFBSSN_bssn2adm
