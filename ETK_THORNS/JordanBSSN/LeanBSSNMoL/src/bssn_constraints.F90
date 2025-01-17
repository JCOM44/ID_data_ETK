! bssn_constraints.F90
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine LeanBSSN_bssn_constraints( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ! Fundamental variables
  CCTK_REAL                alph, beta(3)
  CCTK_REAL                ww, hh(3,3), hu(3,3), trk, aa(3,3), gammat(3),   &
                           dethh, Tab(4,4)

  ! First derivatives
  CCTK_REAL                d1_hh11(3), d1_hh12(3), d1_hh13(3), d1_hh22(3), d1_hh23(3), d1_hh33(3)
  CCTK_REAL                d1_aa11(3), d1_aa12(3), d1_aa13(3), d1_aa22(3), d1_aa23(3), d1_aa33(3)
  CCTK_REAL                d1_gammat1(3), d1_gammat2(3), d1_gammat3(3)
  CCTK_REAL                d1_ww(3), d1_hh(3,3,3), d1_trk(3), d1_aa(3,3,3),&
                           d1_gammat(3,3)

  ! Second derivatives
  CCTK_REAL                d2_hh11(3,3), d2_hh12(3,3), d2_hh13(3,3), d2_hh22(3,3), d2_hh23(3,3), d2_hh33(3,3)
  CCTK_REAL                d2_ww(3,3), d2_hh(3,3,3,3)


  ! Covariant derivatives
  CCTK_REAL                cd2_ww(3,3), cd1_aa(3,3,3)

  ! Constraints
  CCTK_REAL                ham, mom(3)

  ! Ricci tensor
  CCTK_REAL                cf1(3,3,3), cf2(3,3,3), c_ri(3,3)
  CCTK_REAL                c_ri_ww(3,3), c_ri_hh(3,3), ri_1(3,3), ri_2(3,3),&
                           ri_3(3,3), sq_aa, a2(3,3), trr
  CCTK_REAL                tr_cd2_ww, tr_dww_dww

  ! Matter variables
  CCTK_REAL                srcE, srcS, srcjdi(3)

  ! Misc variables
  CCTK_REAL                dx12, dy12, dz12, dxsq12, dysq12, dzsq12,         &
                           dxdy144, dxdz144, dydz144
  CCTK_REAL                odx60, ody60, odz60, odxsq180, odysq180, odzsq180,&
                           odxdy3600, odxdz3600, odydz3600
  CCTK_REAL, parameter ::  one  = 1
  CCTK_REAL, parameter ::  pi   = acos(-one)
  CCTK_REAL, parameter ::  pi8  = 8*pi
  CCTK_REAL, parameter ::  pi16 = 16*pi
  CCTK_INT                 i, j, k
  CCTK_INT                 a, b, c, l, m, n, p, q

  ! Jacobian
  CCTK_REAL                jac(3,3), hes(3,3,3)

  integer                  istat
  logical                  use_jacobian
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ11, lJ12, lJ13
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ21, lJ22, lJ23
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ31, lJ32, lJ33
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: ldJ111, ldJ112, ldJ113, ldJ122, ldJ123, ldJ133
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: ldJ211, ldJ212, ldJ213, ldJ222, ldJ223, ldJ233
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: ldJ311, ldJ312, ldJ313, ldJ322, ldJ323, ldJ333

  CCTK_POINTER             lJ11_ptr, lJ12_ptr, lJ13_ptr
  CCTK_POINTER             lJ21_ptr, lJ22_ptr, lJ23_ptr
  CCTK_POINTER             lJ31_ptr, lJ32_ptr, lJ33_ptr
  CCTK_POINTER             ldJ111_ptr, ldJ112_ptr, ldJ113_ptr, ldJ122_ptr, ldJ123_ptr, ldJ133_ptr
  CCTK_POINTER             ldJ211_ptr, ldJ212_ptr, ldJ213_ptr, ldJ222_ptr, ldJ223_ptr, ldJ233_ptr
  CCTK_POINTER             ldJ311_ptr, ldJ312_ptr, ldJ313_ptr, ldJ322_ptr, ldJ323_ptr, ldJ333_ptr

  pointer (lJ11_ptr, lJ11), (lJ12_ptr, lJ12), (lJ13_ptr, lJ13)
  pointer (lJ21_ptr, lJ21), (lJ22_ptr, lJ22), (lJ23_ptr, lJ23)
  pointer (lJ31_ptr, lJ31), (lJ32_ptr, lJ32), (lJ33_ptr, lJ33)

  pointer (ldJ111_ptr, ldJ111), (ldJ112_ptr, ldJ112), (ldJ113_ptr, ldJ113), (ldJ122_ptr, ldJ122), (ldJ123_ptr, ldJ123), (ldJ133_ptr, ldJ133)
  pointer (ldJ211_ptr, ldJ211), (ldJ212_ptr, ldJ212), (ldJ213_ptr, ldJ213), (ldJ222_ptr, ldJ222), (ldJ223_ptr, ldJ223), (ldJ233_ptr, ldJ233)
  pointer (ldJ311_ptr, ldJ311), (ldJ312_ptr, ldJ312), (ldJ313_ptr, ldJ313), (ldJ322_ptr, ldJ322), (ldJ323_ptr, ldJ323), (ldJ333_ptr, ldJ333)

   if (calculate_constraints_every .le. 0) then
      return
   end if

   if (MOD(cctk_iteration, calculate_constraints_every) .ne. 0 ) then
      return
   endif

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

     call CCTK_VarDataPtr(ldJ111_ptr, cctkGH, 0, "Coordinates::dJ111")
     call CCTK_VarDataPtr(ldJ112_ptr, cctkGH, 0, "Coordinates::dJ112")
     call CCTK_VarDataPtr(ldJ113_ptr, cctkGH, 0, "Coordinates::dJ113")
     call CCTK_VarDataPtr(ldJ122_ptr, cctkGH, 0, "Coordinates::dJ122")
     call CCTK_VarDataPtr(ldJ123_ptr, cctkGH, 0, "Coordinates::dJ123")
     call CCTK_VarDataPtr(ldJ133_ptr, cctkGH, 0, "Coordinates::dJ133")

     call CCTK_VarDataPtr(ldJ211_ptr, cctkGH, 0, "Coordinates::dJ211")
     call CCTK_VarDataPtr(ldJ212_ptr, cctkGH, 0, "Coordinates::dJ212")
     call CCTK_VarDataPtr(ldJ213_ptr, cctkGH, 0, "Coordinates::dJ213")
     call CCTK_VarDataPtr(ldJ222_ptr, cctkGH, 0, "Coordinates::dJ222")
     call CCTK_VarDataPtr(ldJ223_ptr, cctkGH, 0, "Coordinates::dJ223")
     call CCTK_VarDataPtr(ldJ233_ptr, cctkGH, 0, "Coordinates::dJ233")

     call CCTK_VarDataPtr(ldJ311_ptr, cctkGH, 0, "Coordinates::dJ311")
     call CCTK_VarDataPtr(ldJ312_ptr, cctkGH, 0, "Coordinates::dJ312")
     call CCTK_VarDataPtr(ldJ313_ptr, cctkGH, 0, "Coordinates::dJ313")
     call CCTK_VarDataPtr(ldJ322_ptr, cctkGH, 0, "Coordinates::dJ322")
     call CCTK_VarDataPtr(ldJ323_ptr, cctkGH, 0, "Coordinates::dJ323")
     call CCTK_VarDataPtr(ldJ333_ptr, cctkGH, 0, "Coordinates::dJ333")
  end if

  dx12 = 12*CCTK_DELTA_SPACE(1)
  dy12 = 12*CCTK_DELTA_SPACE(2)
  dz12 = 12*CCTK_DELTA_SPACE(3)

  dxsq12 = 12*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(1)
  dysq12 = 12*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(2)
  dzsq12 = 12*CCTK_DELTA_SPACE(3)*CCTK_DELTA_SPACE(3)

  dxdy144 = 144*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2)
  dxdz144 = 144*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(3)
  dydz144 = 144*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(3)


  odx60 = 1 / (60 * CCTK_DELTA_SPACE(1))
  ody60 = 1 / (60 * CCTK_DELTA_SPACE(2))
  odz60 = 1 / (60 * CCTK_DELTA_SPACE(3))

  odxsq180 = 1 / (180*CCTK_DELTA_SPACE(1)**2)
  odysq180 = 1 / (180*CCTK_DELTA_SPACE(2)**2)
  odzsq180 = 1 / (180*CCTK_DELTA_SPACE(3)**2)

  odxdy3600 = 1 / (3600*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2))
  odxdz3600 = 1 / (3600*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(3))
  odydz3600 = 1 / (3600*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(3))

  ! make sure there are no uninitialised values anywhere
  hc    = 0
  mcx   = 0
  mcy   = 0
  mcz   = 0

  !$OMP PARALLEL DO COLLAPSE(3)                                 &
  !$OMP PRIVATE( alph, beta,                                    &
  !$OMP ww, hh, hu, trk, aa, gammat, dethh, Tab,                &
  !$OMP d1_hh11, d1_hh12, d1_hh13, d1_hh22, d1_hh23, d1_hh33,   &
  !$OMP d1_aa11, d1_aa12, d1_aa13, d1_aa22, d1_aa23, d1_aa33,   &
  !$OMP d1_gammat1, d1_gammat2, d1_gammat3,                     &
  !$OMP d1_ww, d1_hh, d1_trk, d1_aa, d1_gammat,                 &
  !$OMP d2_ww, d2_hh, cd2_ww, cd1_aa, ham, mom,                 &
  !$OMP cf1, cf2, c_ri, c_ri_ww, c_ri_hh, ri_1,                 &
  !$OMP ri_2, ri_3, sq_aa, a2, trr,                             &
  !$OMP tr_cd2_ww, tr_dww_dww, srcE, srcS, srcjdi,              &
  !$OMP i, j, k, a, b, c, l, m, n, p, q, jac, hes)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

    !------------ Get local variables ----------
    ww        = conf_fac(i,j,k)

    hh(1,1)   = hxx(i,j,k)
    hh(1,2)   = hxy(i,j,k)
    hh(1,3)   = hxz(i,j,k)
    hh(2,2)   = hyy(i,j,k)
    hh(2,3)   = hyz(i,j,k)
    hh(3,3)   = hzz(i,j,k)
    hh(2,1)   = hh(1,2)
    hh(3,1)   = hh(1,3)
    hh(3,2)   = hh(2,3)

    trk       = tracek(i,j,k)

    aa(1,1)   = axx(i,j,k)
    aa(1,2)   = axy(i,j,k)
    aa(1,3)   = axz(i,j,k)
    aa(2,2)   = ayy(i,j,k)
    aa(2,3)   = ayz(i,j,k)
    aa(3,3)   = azz(i,j,k)
    aa(2,1)   = aa(1,2)
    aa(3,1)   = aa(1,3)
    aa(3,2)   = aa(2,3)

    gammat(1) = gammatx(i,j,k)
    gammat(2) = gammaty(i,j,k)
    gammat(3) = gammatz(i,j,k)

    alph      = alp(i,j,k)

    beta(1)   = betax(i,j,k)
    beta(2)   = betay(i,j,k)
    beta(3)   = betaz(i,j,k)

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

       hes(1,1,1) = ldJ111(i,j,k)
       hes(1,1,2) = ldJ112(i,j,k)
       hes(1,1,3) = ldJ113(i,j,k)
       hes(1,2,1) = ldJ112(i,j,k)
       hes(1,2,2) = ldJ122(i,j,k)
       hes(1,2,3) = ldJ123(i,j,k)
       hes(1,3,1) = ldJ113(i,j,k)
       hes(1,3,2) = ldJ123(i,j,k)
       hes(1,3,3) = ldJ133(i,j,k)

       hes(2,1,1) = ldJ211(i,j,k)
       hes(2,1,2) = ldJ212(i,j,k)
       hes(2,1,3) = ldJ213(i,j,k)
       hes(2,2,1) = ldJ212(i,j,k)
       hes(2,2,2) = ldJ222(i,j,k)
       hes(2,2,3) = ldJ223(i,j,k)
       hes(2,3,1) = ldJ213(i,j,k)
       hes(2,3,2) = ldJ223(i,j,k)
       hes(2,3,3) = ldJ233(i,j,k)

       hes(3,1,1) = ldJ311(i,j,k)
       hes(3,1,2) = ldJ312(i,j,k)
       hes(3,1,3) = ldJ313(i,j,k)
       hes(3,2,1) = ldJ312(i,j,k)
       hes(3,2,2) = ldJ322(i,j,k)
       hes(3,2,3) = ldJ323(i,j,k)
       hes(3,3,1) = ldJ313(i,j,k)
       hes(3,3,2) = ldJ323(i,j,k)
       hes(3,3,3) = ldJ333(i,j,k)
    else
       jac      = 0.0
       jac(1,1) = 1.0
       jac(2,2) = 1.0
       jac(3,3) = 1.0
       hes      = 0.0
    end if

    !------------ Invert metric ----------------
    ! NOTE: deth = 1 by construction, but that is not satisfied numerically
    dethh =       hh(1,1) * hh(2,2) * hh(3,3)                              &
            + 2 * hh(1,2) * hh(1,3) * hh(2,3)                              &
            -     hh(1,1) * hh(2,3) ** 2                                   &
            -     hh(2,2) * hh(1,3) ** 2                                   &
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
    !-------------------------------------------


    if (derivs_order == 6) then

      !------------ Centered 1st derivatives -----
      ! d1_ww(3)
      d1_ww(1) = (  conf_fac(i+3,j,k) - 9*conf_fac(i+2,j,k) + 45*conf_fac(i+1,j,k)          &
                  - conf_fac(i-3,j,k) + 9*conf_fac(i-2,j,k) - 45*conf_fac(i-1,j,k) ) * odx60

      d1_ww(2) = (  conf_fac(i,j+3,k) - 9*conf_fac(i,j+2,k) + 45*conf_fac(i,j+1,k)          &
                  - conf_fac(i,j-3,k) + 9*conf_fac(i,j-2,k) - 45*conf_fac(i,j-1,k) ) * ody60

      d1_ww(3) = (  conf_fac(i,j,k+3) - 9*conf_fac(i,j,k+2) + 45*conf_fac(i,j,k+1)          &
                  - conf_fac(i,j,k-3) + 9*conf_fac(i,j,k-2) - 45*conf_fac(i,j,k-1) ) * odz60

      ! d1_hh(3,3,3)
      d1_hh11(1) = (  hxx(i+3,j,k) - 9*hxx(i+2,j,k) + 45*hxx(i+1,j,k)      &
                    - hxx(i-3,j,k) + 9*hxx(i-2,j,k) - 45*hxx(i-1,j,k) ) * odx60
      d1_hh12(1) = (  hxy(i+3,j,k) - 9*hxy(i+2,j,k) + 45*hxy(i+1,j,k)      &
                    - hxy(i-3,j,k) + 9*hxy(i-2,j,k) - 45*hxy(i-1,j,k) ) * odx60
      d1_hh13(1) = (  hxz(i+3,j,k) - 9*hxz(i+2,j,k) + 45*hxz(i+1,j,k)      &
                    - hxz(i-3,j,k) + 9*hxz(i-2,j,k) - 45*hxz(i-1,j,k) ) * odx60
      d1_hh22(1) = (  hyy(i+3,j,k) - 9*hyy(i+2,j,k) + 45*hyy(i+1,j,k)      &
                    - hyy(i-3,j,k) + 9*hyy(i-2,j,k) - 45*hyy(i-1,j,k) ) * odx60
      d1_hh23(1) = (  hyz(i+3,j,k) - 9*hyz(i+2,j,k) + 45*hyz(i+1,j,k)      &
                    - hyz(i-3,j,k) + 9*hyz(i-2,j,k) - 45*hyz(i-1,j,k) ) * odx60
      d1_hh33(1) = (  hzz(i+3,j,k) - 9*hzz(i+2,j,k) + 45*hzz(i+1,j,k)      &
                    - hzz(i-3,j,k) + 9*hzz(i-2,j,k) - 45*hzz(i-1,j,k) ) * odx60

      d1_hh11(2) = (  hxx(i,j+3,k) - 9*hxx(i,j+2,k) + 45*hxx(i,j+1,k)      &
                    - hxx(i,j-3,k) + 9*hxx(i,j-2,k) - 45*hxx(i,j-1,k) ) * ody60
      d1_hh12(2) = (  hxy(i,j+3,k) - 9*hxy(i,j+2,k) + 45*hxy(i,j+1,k)      &
                    - hxy(i,j-3,k) + 9*hxy(i,j-2,k) - 45*hxy(i,j-1,k) ) * ody60
      d1_hh13(2) = (  hxz(i,j+3,k) - 9*hxz(i,j+2,k) + 45*hxz(i,j+1,k)      &
                    - hxz(i,j-3,k) + 9*hxz(i,j-2,k) - 45*hxz(i,j-1,k) ) * ody60
      d1_hh22(2) = (  hyy(i,j+3,k) - 9*hyy(i,j+2,k) + 45*hyy(i,j+1,k)      &
                    - hyy(i,j-3,k) + 9*hyy(i,j-2,k) - 45*hyy(i,j-1,k) ) * ody60
      d1_hh23(2) = (  hyz(i,j+3,k) - 9*hyz(i,j+2,k) + 45*hyz(i,j+1,k)      &
                    - hyz(i,j-3,k) + 9*hyz(i,j-2,k) - 45*hyz(i,j-1,k) ) * ody60
      d1_hh33(2) = (  hzz(i,j+3,k) - 9*hzz(i,j+2,k) + 45*hzz(i,j+1,k)      &
                    - hzz(i,j-3,k) + 9*hzz(i,j-2,k) - 45*hzz(i,j-1,k) ) * ody60

      d1_hh11(3) = (  hxx(i,j,k+3) - 9*hxx(i,j,k+2) + 45*hxx(i,j,k+1)      &
                    - hxx(i,j,k-3) + 9*hxx(i,j,k-2) - 45*hxx(i,j,k-1) ) * odz60
      d1_hh12(3) = (  hxy(i,j,k+3) - 9*hxy(i,j,k+2) + 45*hxy(i,j,k+1)      &
                    - hxy(i,j,k-3) + 9*hxy(i,j,k-2) - 45*hxy(i,j,k-1) ) * odz60
      d1_hh13(3) = (  hxz(i,j,k+3) - 9*hxz(i,j,k+2) + 45*hxz(i,j,k+1)      &
                    - hxz(i,j,k-3) + 9*hxz(i,j,k-2) - 45*hxz(i,j,k-1) ) * odz60
      d1_hh22(3) = (  hyy(i,j,k+3) - 9*hyy(i,j,k+2) + 45*hyy(i,j,k+1)      &
                    - hyy(i,j,k-3) + 9*hyy(i,j,k-2) - 45*hyy(i,j,k-1) ) * odz60
      d1_hh23(3) = (  hyz(i,j,k+3) - 9*hyz(i,j,k+2) + 45*hyz(i,j,k+1)      &
                    - hyz(i,j,k-3) + 9*hyz(i,j,k-2) - 45*hyz(i,j,k-1) ) * odz60
      d1_hh33(3) = (  hzz(i,j,k+3) - 9*hzz(i,j,k+2) + 45*hzz(i,j,k+1)      &
                    - hzz(i,j,k-3) + 9*hzz(i,j,k-2) - 45*hzz(i,j,k-1) ) * odz60

      ! d1_trk(3)
      d1_trk(1) = (  tracek(i+3,j,k) - 9*tracek(i+2,j,k) + 45*tracek(i+1,j,k)          &
                   - tracek(i-3,j,k) + 9*tracek(i-2,j,k) - 45*tracek(i-1,j,k) ) * odx60

      d1_trk(2) = (  tracek(i,j+3,k) - 9*tracek(i,j+2,k) + 45*tracek(i,j+1,k)          &
                   - tracek(i,j-3,k) + 9*tracek(i,j-2,k) - 45*tracek(i,j-1,k) ) * ody60

      d1_trk(3) = (  tracek(i,j,k+3) - 9*tracek(i,j,k+2) + 45*tracek(i,j,k+1)          &
                   - tracek(i,j,k-3) + 9*tracek(i,j,k-2) - 45*tracek(i,j,k-1) ) * odz60


      ! d1_aa(3,3,3)
      d1_aa11(1) = (  axx(i+3,j,k) - 9*axx(i+2,j,k) + 45*axx(i+1,j,k)      &
                    - axx(i-3,j,k) + 9*axx(i-2,j,k) - 45*axx(i-1,j,k) ) * odx60
      d1_aa12(1) = (  axy(i+3,j,k) - 9*axy(i+2,j,k) + 45*axy(i+1,j,k)      &
                    - axy(i-3,j,k) + 9*axy(i-2,j,k) - 45*axy(i-1,j,k) ) * odx60
      d1_aa13(1) = (  axz(i+3,j,k) - 9*axz(i+2,j,k) + 45*axz(i+1,j,k)      &
                    - axz(i-3,j,k) + 9*axz(i-2,j,k) - 45*axz(i-1,j,k) ) * odx60
      d1_aa22(1) = (  ayy(i+3,j,k) - 9*ayy(i+2,j,k) + 45*ayy(i+1,j,k)      &
                    - ayy(i-3,j,k) + 9*ayy(i-2,j,k) - 45*ayy(i-1,j,k) ) * odx60
      d1_aa23(1) = (  ayz(i+3,j,k) - 9*ayz(i+2,j,k) + 45*ayz(i+1,j,k)      &
                    - ayz(i-3,j,k) + 9*ayz(i-2,j,k) - 45*ayz(i-1,j,k) ) * odx60
      d1_aa33(1) = (  azz(i+3,j,k) - 9*azz(i+2,j,k) + 45*azz(i+1,j,k)      &
                    - azz(i-3,j,k) + 9*azz(i-2,j,k) - 45*azz(i-1,j,k) ) * odx60

      d1_aa11(2) = (  axx(i,j+3,k) - 9*axx(i,j+2,k) + 45*axx(i,j+1,k)      &
                    - axx(i,j-3,k) + 9*axx(i,j-2,k) - 45*axx(i,j-1,k) ) * ody60
      d1_aa12(2) = (  axy(i,j+3,k) - 9*axy(i,j+2,k) + 45*axy(i,j+1,k)      &
                    - axy(i,j-3,k) + 9*axy(i,j-2,k) - 45*axy(i,j-1,k) ) * ody60
      d1_aa13(2) = (  axz(i,j+3,k) - 9*axz(i,j+2,k) + 45*axz(i,j+1,k)      &
                    - axz(i,j-3,k) + 9*axz(i,j-2,k) - 45*axz(i,j-1,k) ) * ody60
      d1_aa22(2) = (  ayy(i,j+3,k) - 9*ayy(i,j+2,k) + 45*ayy(i,j+1,k)      &
                    - ayy(i,j-3,k) + 9*ayy(i,j-2,k) - 45*ayy(i,j-1,k) ) * ody60
      d1_aa23(2) = (  ayz(i,j+3,k) - 9*ayz(i,j+2,k) + 45*ayz(i,j+1,k)      &
                    - ayz(i,j-3,k) + 9*ayz(i,j-2,k) - 45*ayz(i,j-1,k) ) * ody60
      d1_aa33(2) = (  azz(i,j+3,k) - 9*azz(i,j+2,k) + 45*azz(i,j+1,k)      &
                    - azz(i,j-3,k) + 9*azz(i,j-2,k) - 45*azz(i,j-1,k) ) * ody60

      d1_aa11(3) = (  axx(i,j,k+3) - 9*axx(i,j,k+2) + 45*axx(i,j,k+1)      &
                    - axx(i,j,k-3) + 9*axx(i,j,k-2) - 45*axx(i,j,k-1) ) * odz60
      d1_aa12(3) = (  axy(i,j,k+3) - 9*axy(i,j,k+2) + 45*axy(i,j,k+1)      &
                    - axy(i,j,k-3) + 9*axy(i,j,k-2) - 45*axy(i,j,k-1) ) * odz60
      d1_aa13(3) = (  axz(i,j,k+3) - 9*axz(i,j,k+2) + 45*axz(i,j,k+1)      &
                    - axz(i,j,k-3) + 9*axz(i,j,k-2) - 45*axz(i,j,k-1) ) * odz60
      d1_aa22(3) = (  ayy(i,j,k+3) - 9*ayy(i,j,k+2) + 45*ayy(i,j,k+1)      &
                    - ayy(i,j,k-3) + 9*ayy(i,j,k-2) - 45*ayy(i,j,k-1) ) * odz60
      d1_aa23(3) = (  ayz(i,j,k+3) - 9*ayz(i,j,k+2) + 45*ayz(i,j,k+1)      &
                    - ayz(i,j,k-3) + 9*ayz(i,j,k-2) - 45*ayz(i,j,k-1) ) * odz60
      d1_aa33(3) = (  azz(i,j,k+3) - 9*azz(i,j,k+2) + 45*azz(i,j,k+1)      &
                    - azz(i,j,k-3) + 9*azz(i,j,k-2) - 45*azz(i,j,k-1) ) * odz60

      ! d1_gammat(3,3)
      d1_gammat1(1) = (  gammatx(i+3,j,k) - 9*gammatx(i+2,j,k) + 45*gammatx(i+1,j,k) &
                       - gammatx(i-3,j,k) + 9*gammatx(i-2,j,k) - 45*gammatx(i-1,j,k) ) * odx60
      d1_gammat2(1) = (  gammaty(i+3,j,k) - 9*gammaty(i+2,j,k) + 45*gammaty(i+1,j,k) &
                       - gammaty(i-3,j,k) + 9*gammaty(i-2,j,k) - 45*gammaty(i-1,j,k) ) * odx60
      d1_gammat3(1) = (  gammatz(i+3,j,k) - 9*gammatz(i+2,j,k) + 45*gammatz(i+1,j,k) &
                       - gammatz(i-3,j,k) + 9*gammatz(i-2,j,k) - 45*gammatz(i-1,j,k) ) * odx60

      d1_gammat1(2) = (  gammatx(i,j+3,k) - 9*gammatx(i,j+2,k) + 45*gammatx(i,j+1,k) &
                       - gammatx(i,j-3,k) + 9*gammatx(i,j-2,k) - 45*gammatx(i,j-1,k) ) * ody60
      d1_gammat2(2) = (  gammaty(i,j+3,k) - 9*gammaty(i,j+2,k) + 45*gammaty(i,j+1,k) &
                       - gammaty(i,j-3,k) + 9*gammaty(i,j-2,k) - 45*gammaty(i,j-1,k) ) * ody60
      d1_gammat3(2) = (  gammatz(i,j+3,k) - 9*gammatz(i,j+2,k) + 45*gammatz(i,j+1,k) &
                       - gammatz(i,j-3,k) + 9*gammatz(i,j-2,k) - 45*gammatz(i,j-1,k) ) * ody60

      d1_gammat1(3) = (  gammatx(i,j,k+3) - 9*gammatx(i,j,k+2) + 45*gammatx(i,j,k+1) &
                       - gammatx(i,j,k-3) + 9*gammatx(i,j,k-2) - 45*gammatx(i,j,k-1) ) * odz60
      d1_gammat2(3) = (  gammaty(i,j,k+3) - 9*gammaty(i,j,k+2) + 45*gammaty(i,j,k+1) &
                       - gammaty(i,j,k-3) + 9*gammaty(i,j,k-2) - 45*gammaty(i,j,k-1) ) * odz60
      d1_gammat3(3) = (  gammatz(i,j,k+3) - 9*gammatz(i,j,k+2) + 45*gammatz(i,j,k+1) &
                       - gammatz(i,j,k-3) + 9*gammatz(i,j,k-2) - 45*gammatz(i,j,k-1) ) * odz60

      !--------------------------------------------------


      !------------- Centered 2nd derivatives -----------

      ! d2_ww(3,3)
      d2_ww(1,1) = (  2*conf_fac(i+3,j,k) - 27*conf_fac(i+2,j,k) + 270*conf_fac(i+1,j,k) - 490*conf_fac(i,j,k)&
                    + 2*conf_fac(i-3,j,k) - 27*conf_fac(i-2,j,k) + 270*conf_fac(i-1,j,k) ) * odxsq180

      d2_ww(2,2) = (  2*conf_fac(i,j+3,k) - 27*conf_fac(i,j+2,k) + 270*conf_fac(i,j+1,k) - 490*conf_fac(i,j,k)&
                    + 2*conf_fac(i,j-3,k) - 27*conf_fac(i,j-2,k) + 270*conf_fac(i,j-1,k) ) * odysq180

      d2_ww(3,3) = (  2*conf_fac(i,j,k+3) - 27*conf_fac(i,j,k+2) + 270*conf_fac(i,j,k+1) - 490*conf_fac(i,j,k)&
                    + 2*conf_fac(i,j,k-3) - 27*conf_fac(i,j,k-2) + 270*conf_fac(i,j,k-1) ) * odzsq180

      d2_ww(1,2) = (    -conf_fac(i-3,j+3,k) +   9*conf_fac(i-2,j+3,k) -   45*conf_fac(i-1,j+3,k) +   45*conf_fac(i+1,j+3,k) -   9*conf_fac(i+2,j+3,k) +    conf_fac(i+3,j+3,k) &
                    +  9*conf_fac(i-3,j+2,k) -  81*conf_fac(i-2,j+2,k) +  405*conf_fac(i-1,j+2,k) -  405*conf_fac(i+1,j+2,k) +  81*conf_fac(i+2,j+2,k) -  9*conf_fac(i+3,j+2,k) &
                    - 45*conf_fac(i-3,j+1,k) + 405*conf_fac(i-2,j+1,k) - 2025*conf_fac(i-1,j+1,k) + 2025*conf_fac(i+1,j+1,k) - 405*conf_fac(i+2,j+1,k) + 45*conf_fac(i+3,j+1,k) &
                    + 45*conf_fac(i-3,j-1,k) - 405*conf_fac(i-2,j-1,k) + 2025*conf_fac(i-1,j-1,k) - 2025*conf_fac(i+1,j-1,k) + 405*conf_fac(i+2,j-1,k) - 45*conf_fac(i+3,j-1,k) &
                    -  9*conf_fac(i-3,j-2,k) +  81*conf_fac(i-2,j-2,k) -  405*conf_fac(i-1,j-2,k) +  405*conf_fac(i+1,j-2,k) -  81*conf_fac(i+2,j-2,k) +  9*conf_fac(i+3,j-2,k) &
                    +    conf_fac(i-3,j-3,k) -   9*conf_fac(i-2,j-3,k) +   45*conf_fac(i-1,j-3,k) -   45*conf_fac(i+1,j-3,k) +   9*conf_fac(i+2,j-3,k) -    conf_fac(i+3,j-3,k) ) * odxdy3600

      d2_ww(1,3) = (    -conf_fac(i-3,j,k+3) +   9*conf_fac(i-2,j,k+3) -   45*conf_fac(i-1,j,k+3) +   45*conf_fac(i+1,j,k+3) -   9*conf_fac(i+2,j,k+3) +    conf_fac(i+3,j,k+3) &
                    +  9*conf_fac(i-3,j,k+2) -  81*conf_fac(i-2,j,k+2) +  405*conf_fac(i-1,j,k+2) -  405*conf_fac(i+1,j,k+2) +  81*conf_fac(i+2,j,k+2) -  9*conf_fac(i+3,j,k+2) &
                    - 45*conf_fac(i-3,j,k+1) + 405*conf_fac(i-2,j,k+1) - 2025*conf_fac(i-1,j,k+1) + 2025*conf_fac(i+1,j,k+1) - 405*conf_fac(i+2,j,k+1) + 45*conf_fac(i+3,j,k+1) &
                    + 45*conf_fac(i-3,j,k-1) - 405*conf_fac(i-2,j,k-1) + 2025*conf_fac(i-1,j,k-1) - 2025*conf_fac(i+1,j,k-1) + 405*conf_fac(i+2,j,k-1) - 45*conf_fac(i+3,j,k-1) &
                    -  9*conf_fac(i-3,j,k-2) +  81*conf_fac(i-2,j,k-2) -  405*conf_fac(i-1,j,k-2) +  405*conf_fac(i+1,j,k-2) -  81*conf_fac(i+2,j,k-2) +  9*conf_fac(i+3,j,k-2) &
                    +    conf_fac(i-3,j,k-3) -   9*conf_fac(i-2,j,k-3) +   45*conf_fac(i-1,j,k-3) -   45*conf_fac(i+1,j,k-3) +   9*conf_fac(i+2,j,k-3) -    conf_fac(i+3,j,k-3) ) * odxdz3600

      d2_ww(2,3) = (    -conf_fac(i,j-3,k+3) +   9*conf_fac(i,j-2,k+3) -   45*conf_fac(i,j-1,k+3) +   45*conf_fac(i,j+1,k+3) -   9*conf_fac(i,j+2,k+3) +    conf_fac(i,j+3,k+3) &
                    +  9*conf_fac(i,j-3,k+2) -  81*conf_fac(i,j-2,k+2) +  405*conf_fac(i,j-1,k+2) -  405*conf_fac(i,j+1,k+2) +  81*conf_fac(i,j+2,k+2) -  9*conf_fac(i,j+3,k+2) &
                    - 45*conf_fac(i,j-3,k+1) + 405*conf_fac(i,j-2,k+1) - 2025*conf_fac(i,j-1,k+1) + 2025*conf_fac(i,j+1,k+1) - 405*conf_fac(i,j+2,k+1) + 45*conf_fac(i,j+3,k+1) &
                    + 45*conf_fac(i,j-3,k-1) - 405*conf_fac(i,j-2,k-1) + 2025*conf_fac(i,j-1,k-1) - 2025*conf_fac(i,j+1,k-1) + 405*conf_fac(i,j+2,k-1) - 45*conf_fac(i,j+3,k-1) &
                    -  9*conf_fac(i,j-3,k-2) +  81*conf_fac(i,j-2,k-2) -  405*conf_fac(i,j-1,k-2) +  405*conf_fac(i,j+1,k-2) -  81*conf_fac(i,j+2,k-2) +  9*conf_fac(i,j+3,k-2) &
                    +    conf_fac(i,j-3,k-3) -   9*conf_fac(i,j-2,k-3) +   45*conf_fac(i,j-1,k-3) -   45*conf_fac(i,j+1,k-3) +   9*conf_fac(i,j+2,k-3) -    conf_fac(i,j+3,k-3) ) * odydz3600

      d2_ww(2,1) = d2_ww(1,2)
      d2_ww(3,1) = d2_ww(1,3)
      d2_ww(3,2) = d2_ww(2,3)


      ! d2_hh(3,3,3,3)
      d2_hh11(1,1) = (  2*hxx(i+3,j,k) - 27*hxx(i+2,j,k) + 270*hxx(i+1,j,k) - 490*hxx(i,j,k)&
                      + 2*hxx(i-3,j,k) - 27*hxx(i-2,j,k) + 270*hxx(i-1,j,k) ) * odxsq180
      d2_hh12(1,1) = (  2*hxy(i+3,j,k) - 27*hxy(i+2,j,k) + 270*hxy(i+1,j,k) - 490*hxy(i,j,k)&
                      + 2*hxy(i-3,j,k) - 27*hxy(i-2,j,k) + 270*hxy(i-1,j,k) ) * odxsq180
      d2_hh13(1,1) = (  2*hxz(i+3,j,k) - 27*hxz(i+2,j,k) + 270*hxz(i+1,j,k) - 490*hxz(i,j,k)&
                      + 2*hxz(i-3,j,k) - 27*hxz(i-2,j,k) + 270*hxz(i-1,j,k) ) * odxsq180
      d2_hh22(1,1) = (  2*hyy(i+3,j,k) - 27*hyy(i+2,j,k) + 270*hyy(i+1,j,k) - 490*hyy(i,j,k)&
                      + 2*hyy(i-3,j,k) - 27*hyy(i-2,j,k) + 270*hyy(i-1,j,k) ) * odxsq180
      d2_hh23(1,1) = (  2*hyz(i+3,j,k) - 27*hyz(i+2,j,k) + 270*hyz(i+1,j,k) - 490*hyz(i,j,k)&
                      + 2*hyz(i-3,j,k) - 27*hyz(i-2,j,k) + 270*hyz(i-1,j,k) ) * odxsq180
      d2_hh33(1,1) = (  2*hzz(i+3,j,k) - 27*hzz(i+2,j,k) + 270*hzz(i+1,j,k) - 490*hzz(i,j,k)&
                      + 2*hzz(i-3,j,k) - 27*hzz(i-2,j,k) + 270*hzz(i-1,j,k) ) * odxsq180

      d2_hh11(2,2) = (  2*hxx(i,j+3,k) - 27*hxx(i,j+2,k) + 270*hxx(i,j+1,k) - 490*hxx(i,j,k)&
                      + 2*hxx(i,j-3,k) - 27*hxx(i,j-2,k) + 270*hxx(i,j-1,k) ) * odysq180
      d2_hh12(2,2) = (  2*hxy(i,j+3,k) - 27*hxy(i,j+2,k) + 270*hxy(i,j+1,k) - 490*hxy(i,j,k)&
                      + 2*hxy(i,j-3,k) - 27*hxy(i,j-2,k) + 270*hxy(i,j-1,k) ) * odysq180
      d2_hh13(2,2) = (  2*hxz(i,j+3,k) - 27*hxz(i,j+2,k) + 270*hxz(i,j+1,k) - 490*hxz(i,j,k)&
                      + 2*hxz(i,j-3,k) - 27*hxz(i,j-2,k) + 270*hxz(i,j-1,k) ) * odysq180
      d2_hh22(2,2) = (  2*hyy(i,j+3,k) - 27*hyy(i,j+2,k) + 270*hyy(i,j+1,k) - 490*hyy(i,j,k)&
                      + 2*hyy(i,j-3,k) - 27*hyy(i,j-2,k) + 270*hyy(i,j-1,k) ) * odysq180
      d2_hh23(2,2) = (  2*hyz(i,j+3,k) - 27*hyz(i,j+2,k) + 270*hyz(i,j+1,k) - 490*hyz(i,j,k)&
                      + 2*hyz(i,j-3,k) - 27*hyz(i,j-2,k) + 270*hyz(i,j-1,k) ) * odysq180
      d2_hh33(2,2) = (  2*hzz(i,j+3,k) - 27*hzz(i,j+2,k) + 270*hzz(i,j+1,k) - 490*hzz(i,j,k)&
                      + 2*hzz(i,j-3,k) - 27*hzz(i,j-2,k) + 270*hzz(i,j-1,k) ) * odysq180

      d2_hh11(3,3) = (  2*hxx(i,j,k+3) - 27*hxx(i,j,k+2) + 270*hxx(i,j,k+1) - 490*hxx(i,j,k)&
                      + 2*hxx(i,j,k-3) - 27*hxx(i,j,k-2) + 270*hxx(i,j,k-1) ) * odzsq180
      d2_hh12(3,3) = (  2*hxy(i,j,k+3) - 27*hxy(i,j,k+2) + 270*hxy(i,j,k+1) - 490*hxy(i,j,k)&
                      + 2*hxy(i,j,k-3) - 27*hxy(i,j,k-2) + 270*hxy(i,j,k-1) ) * odzsq180
      d2_hh13(3,3) = (  2*hxz(i,j,k+3) - 27*hxz(i,j,k+2) + 270*hxz(i,j,k+1) - 490*hxz(i,j,k)&
                      + 2*hxz(i,j,k-3) - 27*hxz(i,j,k-2) + 270*hxz(i,j,k-1) ) * odzsq180
      d2_hh22(3,3) = (  2*hyy(i,j,k+3) - 27*hyy(i,j,k+2) + 270*hyy(i,j,k+1) - 490*hyy(i,j,k)&
                      + 2*hyy(i,j,k-3) - 27*hyy(i,j,k-2) + 270*hyy(i,j,k-1) ) * odzsq180
      d2_hh23(3,3) = (  2*hyz(i,j,k+3) - 27*hyz(i,j,k+2) + 270*hyz(i,j,k+1) - 490*hyz(i,j,k)&
                      + 2*hyz(i,j,k-3) - 27*hyz(i,j,k-2) + 270*hyz(i,j,k-1) ) * odzsq180
      d2_hh33(3,3) = (  2*hzz(i,j,k+3) - 27*hzz(i,j,k+2) + 270*hzz(i,j,k+1) - 490*hzz(i,j,k)&
                      + 2*hzz(i,j,k-3) - 27*hzz(i,j,k-2) + 270*hzz(i,j,k-1) ) * odzsq180

      d2_hh11(1,2) = (    -hxx(i-3,j+3,k) +   9*hxx(i-2,j+3,k) -   45*hxx(i-1,j+3,k) +   45*hxx(i+1,j+3,k) -   9*hxx(i+2,j+3,k) +    hxx(i+3,j+3,k) &
                      +  9*hxx(i-3,j+2,k) -  81*hxx(i-2,j+2,k) +  405*hxx(i-1,j+2,k) -  405*hxx(i+1,j+2,k) +  81*hxx(i+2,j+2,k) -  9*hxx(i+3,j+2,k) &
                      - 45*hxx(i-3,j+1,k) + 405*hxx(i-2,j+1,k) - 2025*hxx(i-1,j+1,k) + 2025*hxx(i+1,j+1,k) - 405*hxx(i+2,j+1,k) + 45*hxx(i+3,j+1,k) &
                      + 45*hxx(i-3,j-1,k) - 405*hxx(i-2,j-1,k) + 2025*hxx(i-1,j-1,k) - 2025*hxx(i+1,j-1,k) + 405*hxx(i+2,j-1,k) - 45*hxx(i+3,j-1,k) &
                      -  9*hxx(i-3,j-2,k) +  81*hxx(i-2,j-2,k) -  405*hxx(i-1,j-2,k) +  405*hxx(i+1,j-2,k) -  81*hxx(i+2,j-2,k) +  9*hxx(i+3,j-2,k) &
                      +    hxx(i-3,j-3,k) -   9*hxx(i-2,j-3,k) +   45*hxx(i-1,j-3,k) -   45*hxx(i+1,j-3,k) +   9*hxx(i+2,j-3,k) -    hxx(i+3,j-3,k) ) * odxdy3600
      d2_hh12(1,2) = (    -hxy(i-3,j+3,k) +   9*hxy(i-2,j+3,k) -   45*hxy(i-1,j+3,k) +   45*hxy(i+1,j+3,k) -   9*hxy(i+2,j+3,k) +    hxy(i+3,j+3,k) &
                      +  9*hxy(i-3,j+2,k) -  81*hxy(i-2,j+2,k) +  405*hxy(i-1,j+2,k) -  405*hxy(i+1,j+2,k) +  81*hxy(i+2,j+2,k) -  9*hxy(i+3,j+2,k) &
                      - 45*hxy(i-3,j+1,k) + 405*hxy(i-2,j+1,k) - 2025*hxy(i-1,j+1,k) + 2025*hxy(i+1,j+1,k) - 405*hxy(i+2,j+1,k) + 45*hxy(i+3,j+1,k) &
                      + 45*hxy(i-3,j-1,k) - 405*hxy(i-2,j-1,k) + 2025*hxy(i-1,j-1,k) - 2025*hxy(i+1,j-1,k) + 405*hxy(i+2,j-1,k) - 45*hxy(i+3,j-1,k) &
                      -  9*hxy(i-3,j-2,k) +  81*hxy(i-2,j-2,k) -  405*hxy(i-1,j-2,k) +  405*hxy(i+1,j-2,k) -  81*hxy(i+2,j-2,k) +  9*hxy(i+3,j-2,k) &
                      +    hxy(i-3,j-3,k) -   9*hxy(i-2,j-3,k) +   45*hxy(i-1,j-3,k) -   45*hxy(i+1,j-3,k) +   9*hxy(i+2,j-3,k) -    hxy(i+3,j-3,k) ) * odxdy3600
      d2_hh13(1,2) = (    -hxz(i-3,j+3,k) +   9*hxz(i-2,j+3,k) -   45*hxz(i-1,j+3,k) +   45*hxz(i+1,j+3,k) -   9*hxz(i+2,j+3,k) +    hxz(i+3,j+3,k) &
                      +  9*hxz(i-3,j+2,k) -  81*hxz(i-2,j+2,k) +  405*hxz(i-1,j+2,k) -  405*hxz(i+1,j+2,k) +  81*hxz(i+2,j+2,k) -  9*hxz(i+3,j+2,k) &
                      - 45*hxz(i-3,j+1,k) + 405*hxz(i-2,j+1,k) - 2025*hxz(i-1,j+1,k) + 2025*hxz(i+1,j+1,k) - 405*hxz(i+2,j+1,k) + 45*hxz(i+3,j+1,k) &
                      + 45*hxz(i-3,j-1,k) - 405*hxz(i-2,j-1,k) + 2025*hxz(i-1,j-1,k) - 2025*hxz(i+1,j-1,k) + 405*hxz(i+2,j-1,k) - 45*hxz(i+3,j-1,k) &
                      -  9*hxz(i-3,j-2,k) +  81*hxz(i-2,j-2,k) -  405*hxz(i-1,j-2,k) +  405*hxz(i+1,j-2,k) -  81*hxz(i+2,j-2,k) +  9*hxz(i+3,j-2,k) &
                      +    hxz(i-3,j-3,k) -   9*hxz(i-2,j-3,k) +   45*hxz(i-1,j-3,k) -   45*hxz(i+1,j-3,k) +   9*hxz(i+2,j-3,k) -    hxz(i+3,j-3,k) ) * odxdy3600
      d2_hh22(1,2) = (    -hyy(i-3,j+3,k) +   9*hyy(i-2,j+3,k) -   45*hyy(i-1,j+3,k) +   45*hyy(i+1,j+3,k) -   9*hyy(i+2,j+3,k) +    hyy(i+3,j+3,k) &
                      +  9*hyy(i-3,j+2,k) -  81*hyy(i-2,j+2,k) +  405*hyy(i-1,j+2,k) -  405*hyy(i+1,j+2,k) +  81*hyy(i+2,j+2,k) -  9*hyy(i+3,j+2,k) &
                      - 45*hyy(i-3,j+1,k) + 405*hyy(i-2,j+1,k) - 2025*hyy(i-1,j+1,k) + 2025*hyy(i+1,j+1,k) - 405*hyy(i+2,j+1,k) + 45*hyy(i+3,j+1,k) &
                      + 45*hyy(i-3,j-1,k) - 405*hyy(i-2,j-1,k) + 2025*hyy(i-1,j-1,k) - 2025*hyy(i+1,j-1,k) + 405*hyy(i+2,j-1,k) - 45*hyy(i+3,j-1,k) &
                      -  9*hyy(i-3,j-2,k) +  81*hyy(i-2,j-2,k) -  405*hyy(i-1,j-2,k) +  405*hyy(i+1,j-2,k) -  81*hyy(i+2,j-2,k) +  9*hyy(i+3,j-2,k) &
                      +    hyy(i-3,j-3,k) -   9*hyy(i-2,j-3,k) +   45*hyy(i-1,j-3,k) -   45*hyy(i+1,j-3,k) +   9*hyy(i+2,j-3,k) -    hyy(i+3,j-3,k) ) * odxdy3600
      d2_hh23(1,2) = (    -hyz(i-3,j+3,k) +   9*hyz(i-2,j+3,k) -   45*hyz(i-1,j+3,k) +   45*hyz(i+1,j+3,k) -   9*hyz(i+2,j+3,k) +    hyz(i+3,j+3,k) &
                      +  9*hyz(i-3,j+2,k) -  81*hyz(i-2,j+2,k) +  405*hyz(i-1,j+2,k) -  405*hyz(i+1,j+2,k) +  81*hyz(i+2,j+2,k) -  9*hyz(i+3,j+2,k) &
                      - 45*hyz(i-3,j+1,k) + 405*hyz(i-2,j+1,k) - 2025*hyz(i-1,j+1,k) + 2025*hyz(i+1,j+1,k) - 405*hyz(i+2,j+1,k) + 45*hyz(i+3,j+1,k) &
                      + 45*hyz(i-3,j-1,k) - 405*hyz(i-2,j-1,k) + 2025*hyz(i-1,j-1,k) - 2025*hyz(i+1,j-1,k) + 405*hyz(i+2,j-1,k) - 45*hyz(i+3,j-1,k) &
                      -  9*hyz(i-3,j-2,k) +  81*hyz(i-2,j-2,k) -  405*hyz(i-1,j-2,k) +  405*hyz(i+1,j-2,k) -  81*hyz(i+2,j-2,k) +  9*hyz(i+3,j-2,k) &
                      +    hyz(i-3,j-3,k) -   9*hyz(i-2,j-3,k) +   45*hyz(i-1,j-3,k) -   45*hyz(i+1,j-3,k) +   9*hyz(i+2,j-3,k) -    hyz(i+3,j-3,k) ) * odxdy3600
      d2_hh33(1,2) = (    -hzz(i-3,j+3,k) +   9*hzz(i-2,j+3,k) -   45*hzz(i-1,j+3,k) +   45*hzz(i+1,j+3,k) -   9*hzz(i+2,j+3,k) +    hzz(i+3,j+3,k) &
                      +  9*hzz(i-3,j+2,k) -  81*hzz(i-2,j+2,k) +  405*hzz(i-1,j+2,k) -  405*hzz(i+1,j+2,k) +  81*hzz(i+2,j+2,k) -  9*hzz(i+3,j+2,k) &
                      - 45*hzz(i-3,j+1,k) + 405*hzz(i-2,j+1,k) - 2025*hzz(i-1,j+1,k) + 2025*hzz(i+1,j+1,k) - 405*hzz(i+2,j+1,k) + 45*hzz(i+3,j+1,k) &
                      + 45*hzz(i-3,j-1,k) - 405*hzz(i-2,j-1,k) + 2025*hzz(i-1,j-1,k) - 2025*hzz(i+1,j-1,k) + 405*hzz(i+2,j-1,k) - 45*hzz(i+3,j-1,k) &
                      -  9*hzz(i-3,j-2,k) +  81*hzz(i-2,j-2,k) -  405*hzz(i-1,j-2,k) +  405*hzz(i+1,j-2,k) -  81*hzz(i+2,j-2,k) +  9*hzz(i+3,j-2,k) &
                      +    hzz(i-3,j-3,k) -   9*hzz(i-2,j-3,k) +   45*hzz(i-1,j-3,k) -   45*hzz(i+1,j-3,k) +   9*hzz(i+2,j-3,k) -    hzz(i+3,j-3,k) ) * odxdy3600

      d2_hh11(1,3) = (    -hxx(i-3,j,k+3) +   9*hxx(i-2,j,k+3) -   45*hxx(i-1,j,k+3) +   45*hxx(i+1,j,k+3) -   9*hxx(i+2,j,k+3) +    hxx(i+3,j,k+3) &
                      +  9*hxx(i-3,j,k+2) -  81*hxx(i-2,j,k+2) +  405*hxx(i-1,j,k+2) -  405*hxx(i+1,j,k+2) +  81*hxx(i+2,j,k+2) -  9*hxx(i+3,j,k+2) &
                      - 45*hxx(i-3,j,k+1) + 405*hxx(i-2,j,k+1) - 2025*hxx(i-1,j,k+1) + 2025*hxx(i+1,j,k+1) - 405*hxx(i+2,j,k+1) + 45*hxx(i+3,j,k+1) &
                      + 45*hxx(i-3,j,k-1) - 405*hxx(i-2,j,k-1) + 2025*hxx(i-1,j,k-1) - 2025*hxx(i+1,j,k-1) + 405*hxx(i+2,j,k-1) - 45*hxx(i+3,j,k-1) &
                      -  9*hxx(i-3,j,k-2) +  81*hxx(i-2,j,k-2) -  405*hxx(i-1,j,k-2) +  405*hxx(i+1,j,k-2) -  81*hxx(i+2,j,k-2) +  9*hxx(i+3,j,k-2) &
                      +    hxx(i-3,j,k-3) -   9*hxx(i-2,j,k-3) +   45*hxx(i-1,j,k-3) -   45*hxx(i+1,j,k-3) +   9*hxx(i+2,j,k-3) -    hxx(i+3,j,k-3) ) * odxdz3600
      d2_hh12(1,3) = (    -hxy(i-3,j,k+3) +   9*hxy(i-2,j,k+3) -   45*hxy(i-1,j,k+3) +   45*hxy(i+1,j,k+3) -   9*hxy(i+2,j,k+3) +    hxy(i+3,j,k+3) &
                      +  9*hxy(i-3,j,k+2) -  81*hxy(i-2,j,k+2) +  405*hxy(i-1,j,k+2) -  405*hxy(i+1,j,k+2) +  81*hxy(i+2,j,k+2) -  9*hxy(i+3,j,k+2) &
                        - 45*hxy(i-3,j,k+1) + 405*hxy(i-2,j,k+1) - 2025*hxy(i-1,j,k+1) + 2025*hxy(i+1,j,k+1) - 405*hxy(i+2,j,k+1) + 45*hxy(i+3,j,k+1) &
                      + 45*hxy(i-3,j,k-1) - 405*hxy(i-2,j,k-1) + 2025*hxy(i-1,j,k-1) - 2025*hxy(i+1,j,k-1) + 405*hxy(i+2,j,k-1) - 45*hxy(i+3,j,k-1) &
                      -  9*hxy(i-3,j,k-2) +  81*hxy(i-2,j,k-2) -  405*hxy(i-1,j,k-2) +  405*hxy(i+1,j,k-2) -  81*hxy(i+2,j,k-2) +  9*hxy(i+3,j,k-2) &
                      +    hxy(i-3,j,k-3) -   9*hxy(i-2,j,k-3) +   45*hxy(i-1,j,k-3) -   45*hxy(i+1,j,k-3) +   9*hxy(i+2,j,k-3) -    hxy(i+3,j,k-3) ) * odxdz3600
      d2_hh13(1,3) = (    -hxz(i-3,j,k+3) +   9*hxz(i-2,j,k+3) -   45*hxz(i-1,j,k+3) +   45*hxz(i+1,j,k+3) -   9*hxz(i+2,j,k+3) +    hxz(i+3,j,k+3) &
                      +  9*hxz(i-3,j,k+2) -  81*hxz(i-2,j,k+2) +  405*hxz(i-1,j,k+2) -  405*hxz(i+1,j,k+2) +  81*hxz(i+2,j,k+2) -  9*hxz(i+3,j,k+2) &
                      - 45*hxz(i-3,j,k+1) + 405*hxz(i-2,j,k+1) - 2025*hxz(i-1,j,k+1) + 2025*hxz(i+1,j,k+1) - 405*hxz(i+2,j,k+1) + 45*hxz(i+3,j,k+1) &
                      + 45*hxz(i-3,j,k-1) - 405*hxz(i-2,j,k-1) + 2025*hxz(i-1,j,k-1) - 2025*hxz(i+1,j,k-1) + 405*hxz(i+2,j,k-1) - 45*hxz(i+3,j,k-1) &
                      -  9*hxz(i-3,j,k-2) +  81*hxz(i-2,j,k-2) -  405*hxz(i-1,j,k-2) +  405*hxz(i+1,j,k-2) -  81*hxz(i+2,j,k-2) +  9*hxz(i+3,j,k-2) &
                      +    hxz(i-3,j,k-3) -   9*hxz(i-2,j,k-3) +   45*hxz(i-1,j,k-3) -   45*hxz(i+1,j,k-3) +   9*hxz(i+2,j,k-3) -    hxz(i+3,j,k-3) ) * odxdz3600
      d2_hh22(1,3) = (    -hyy(i-3,j,k+3) +   9*hyy(i-2,j,k+3) -   45*hyy(i-1,j,k+3) +   45*hyy(i+1,j,k+3) -   9*hyy(i+2,j,k+3) +    hyy(i+3,j,k+3) &
                      +  9*hyy(i-3,j,k+2) -  81*hyy(i-2,j,k+2) +  405*hyy(i-1,j,k+2) -  405*hyy(i+1,j,k+2) +  81*hyy(i+2,j,k+2) -  9*hyy(i+3,j,k+2) &
                      - 45*hyy(i-3,j,k+1) + 405*hyy(i-2,j,k+1) - 2025*hyy(i-1,j,k+1) + 2025*hyy(i+1,j,k+1) - 405*hyy(i+2,j,k+1) + 45*hyy(i+3,j,k+1) &
                      + 45*hyy(i-3,j,k-1) - 405*hyy(i-2,j,k-1) + 2025*hyy(i-1,j,k-1) - 2025*hyy(i+1,j,k-1) + 405*hyy(i+2,j,k-1) - 45*hyy(i+3,j,k-1) &
                      -  9*hyy(i-3,j,k-2) +  81*hyy(i-2,j,k-2) -  405*hyy(i-1,j,k-2) +  405*hyy(i+1,j,k-2) -  81*hyy(i+2,j,k-2) +  9*hyy(i+3,j,k-2) &
                      +    hyy(i-3,j,k-3) -   9*hyy(i-2,j,k-3) +   45*hyy(i-1,j,k-3) -   45*hyy(i+1,j,k-3) +   9*hyy(i+2,j,k-3) -    hyy(i+3,j,k-3) ) * odxdz3600
      d2_hh23(1,3) = (    -hyz(i-3,j,k+3) +   9*hyz(i-2,j,k+3) -   45*hyz(i-1,j,k+3) +   45*hyz(i+1,j,k+3) -   9*hyz(i+2,j,k+3) +    hyz(i+3,j,k+3) &
                      +  9*hyz(i-3,j,k+2) -  81*hyz(i-2,j,k+2) +  405*hyz(i-1,j,k+2) -  405*hyz(i+1,j,k+2) +  81*hyz(i+2,j,k+2) -  9*hyz(i+3,j,k+2) &
                      - 45*hyz(i-3,j,k+1) + 405*hyz(i-2,j,k+1) - 2025*hyz(i-1,j,k+1) + 2025*hyz(i+1,j,k+1) - 405*hyz(i+2,j,k+1) + 45*hyz(i+3,j,k+1) &
                      + 45*hyz(i-3,j,k-1) - 405*hyz(i-2,j,k-1) + 2025*hyz(i-1,j,k-1) - 2025*hyz(i+1,j,k-1) + 405*hyz(i+2,j,k-1) - 45*hyz(i+3,j,k-1) &
                      -  9*hyz(i-3,j,k-2) +  81*hyz(i-2,j,k-2) -  405*hyz(i-1,j,k-2) +  405*hyz(i+1,j,k-2) -  81*hyz(i+2,j,k-2) +  9*hyz(i+3,j,k-2) &
                      +    hyz(i-3,j,k-3) -   9*hyz(i-2,j,k-3) +   45*hyz(i-1,j,k-3) -   45*hyz(i+1,j,k-3) +   9*hyz(i+2,j,k-3) -    hyz(i+3,j,k-3) ) * odxdz3600
      d2_hh33(1,3) = (    -hzz(i-3,j,k+3) +   9*hzz(i-2,j,k+3) -   45*hzz(i-1,j,k+3) +   45*hzz(i+1,j,k+3) -   9*hzz(i+2,j,k+3) +    hzz(i+3,j,k+3) &
                      +  9*hzz(i-3,j,k+2) -  81*hzz(i-2,j,k+2) +  405*hzz(i-1,j,k+2) -  405*hzz(i+1,j,k+2) +  81*hzz(i+2,j,k+2) -  9*hzz(i+3,j,k+2) &
                      - 45*hzz(i-3,j,k+1) + 405*hzz(i-2,j,k+1) - 2025*hzz(i-1,j,k+1) + 2025*hzz(i+1,j,k+1) - 405*hzz(i+2,j,k+1) + 45*hzz(i+3,j,k+1) &
                      + 45*hzz(i-3,j,k-1) - 405*hzz(i-2,j,k-1) + 2025*hzz(i-1,j,k-1) - 2025*hzz(i+1,j,k-1) + 405*hzz(i+2,j,k-1) - 45*hzz(i+3,j,k-1) &
                      -  9*hzz(i-3,j,k-2) +  81*hzz(i-2,j,k-2) -  405*hzz(i-1,j,k-2) +  405*hzz(i+1,j,k-2) -  81*hzz(i+2,j,k-2) +  9*hzz(i+3,j,k-2) &
                      +    hzz(i-3,j,k-3) -   9*hzz(i-2,j,k-3) +   45*hzz(i-1,j,k-3) -   45*hzz(i+1,j,k-3) +   9*hzz(i+2,j,k-3) -    hzz(i+3,j,k-3) ) * odxdz3600

      d2_hh11(2,3) = (    -hxx(i,j-3,k+3) +   9*hxx(i,j-2,k+3) -   45*hxx(i,j-1,k+3) +   45*hxx(i,j+1,k+3) -   9*hxx(i,j+2,k+3) +    hxx(i,j+3,k+3) &
                      +  9*hxx(i,j-3,k+2) -  81*hxx(i,j-2,k+2) +  405*hxx(i,j-1,k+2) -  405*hxx(i,j+1,k+2) +  81*hxx(i,j+2,k+2) -  9*hxx(i,j+3,k+2) &
                      - 45*hxx(i,j-3,k+1) + 405*hxx(i,j-2,k+1) - 2025*hxx(i,j-1,k+1) + 2025*hxx(i,j+1,k+1) - 405*hxx(i,j+2,k+1) + 45*hxx(i,j+3,k+1) &
                      + 45*hxx(i,j-3,k-1) - 405*hxx(i,j-2,k-1) + 2025*hxx(i,j-1,k-1) - 2025*hxx(i,j+1,k-1) + 405*hxx(i,j+2,k-1) - 45*hxx(i,j+3,k-1) &
                      -  9*hxx(i,j-3,k-2) +  81*hxx(i,j-2,k-2) -  405*hxx(i,j-1,k-2) +  405*hxx(i,j+1,k-2) -  81*hxx(i,j+2,k-2) +  9*hxx(i,j+3,k-2) &
                      +    hxx(i,j-3,k-3) -   9*hxx(i,j-2,k-3) +   45*hxx(i,j-1,k-3) -   45*hxx(i,j+1,k-3) +   9*hxx(i,j+2,k-3) -    hxx(i,j+3,k-3) ) * odydz3600
      d2_hh12(2,3) = (    -hxy(i,j-3,k+3) +   9*hxy(i,j-2,k+3) -   45*hxy(i,j-1,k+3) +   45*hxy(i,j+1,k+3) -   9*hxy(i,j+2,k+3) +    hxy(i,j+3,k+3) &
                      +  9*hxy(i,j-3,k+2) -  81*hxy(i,j-2,k+2) +  405*hxy(i,j-1,k+2) -  405*hxy(i,j+1,k+2) +  81*hxy(i,j+2,k+2) -  9*hxy(i,j+3,k+2) &
                      - 45*hxy(i,j-3,k+1) + 405*hxy(i,j-2,k+1) - 2025*hxy(i,j-1,k+1) + 2025*hxy(i,j+1,k+1) - 405*hxy(i,j+2,k+1) + 45*hxy(i,j+3,k+1) &
                      + 45*hxy(i,j-3,k-1) - 405*hxy(i,j-2,k-1) + 2025*hxy(i,j-1,k-1) - 2025*hxy(i,j+1,k-1) + 405*hxy(i,j+2,k-1) - 45*hxy(i,j+3,k-1) &
                      -  9*hxy(i,j-3,k-2) +  81*hxy(i,j-2,k-2) -  405*hxy(i,j-1,k-2) +  405*hxy(i,j+1,k-2) -  81*hxy(i,j+2,k-2) +  9*hxy(i,j+3,k-2) &
                      +    hxy(i,j-3,k-3) -   9*hxy(i,j-2,k-3) +   45*hxy(i,j-1,k-3) -   45*hxy(i,j+1,k-3) +   9*hxy(i,j+2,k-3) -    hxy(i,j+3,k-3) ) * odydz3600
      d2_hh13(2,3) = (    -hxz(i,j-3,k+3) +   9*hxz(i,j-2,k+3) -   45*hxz(i,j-1,k+3) +   45*hxz(i,j+1,k+3) -   9*hxz(i,j+2,k+3) +    hxz(i,j+3,k+3) &
                      +  9*hxz(i,j-3,k+2) -  81*hxz(i,j-2,k+2) +  405*hxz(i,j-1,k+2) -  405*hxz(i,j+1,k+2) +  81*hxz(i,j+2,k+2) -  9*hxz(i,j+3,k+2) &
                      - 45*hxz(i,j-3,k+1) + 405*hxz(i,j-2,k+1) - 2025*hxz(i,j-1,k+1) + 2025*hxz(i,j+1,k+1) - 405*hxz(i,j+2,k+1) + 45*hxz(i,j+3,k+1) &
                      + 45*hxz(i,j-3,k-1) - 405*hxz(i,j-2,k-1) + 2025*hxz(i,j-1,k-1) - 2025*hxz(i,j+1,k-1) + 405*hxz(i,j+2,k-1) - 45*hxz(i,j+3,k-1) &
                      -  9*hxz(i,j-3,k-2) +  81*hxz(i,j-2,k-2) -  405*hxz(i,j-1,k-2) +  405*hxz(i,j+1,k-2) -  81*hxz(i,j+2,k-2) +  9*hxz(i,j+3,k-2) &
                      +    hxz(i,j-3,k-3) -   9*hxz(i,j-2,k-3) +   45*hxz(i,j-1,k-3) -   45*hxz(i,j+1,k-3) +   9*hxz(i,j+2,k-3) -    hxz(i,j+3,k-3) ) * odydz3600
      d2_hh22(2,3) = (    -hyy(i,j-3,k+3) +   9*hyy(i,j-2,k+3) -   45*hyy(i,j-1,k+3) +   45*hyy(i,j+1,k+3) -   9*hyy(i,j+2,k+3) +    hyy(i,j+3,k+3) &
                      +  9*hyy(i,j-3,k+2) -  81*hyy(i,j-2,k+2) +  405*hyy(i,j-1,k+2) -  405*hyy(i,j+1,k+2) +  81*hyy(i,j+2,k+2) -  9*hyy(i,j+3,k+2) &
                      - 45*hyy(i,j-3,k+1) + 405*hyy(i,j-2,k+1) - 2025*hyy(i,j-1,k+1) + 2025*hyy(i,j+1,k+1) - 405*hyy(i,j+2,k+1) + 45*hyy(i,j+3,k+1) &
                      + 45*hyy(i,j-3,k-1) - 405*hyy(i,j-2,k-1) + 2025*hyy(i,j-1,k-1) - 2025*hyy(i,j+1,k-1) + 405*hyy(i,j+2,k-1) - 45*hyy(i,j+3,k-1) &
                      -  9*hyy(i,j-3,k-2) +  81*hyy(i,j-2,k-2) -  405*hyy(i,j-1,k-2) +  405*hyy(i,j+1,k-2) -  81*hyy(i,j+2,k-2) +  9*hyy(i,j+3,k-2) &
                      +    hyy(i,j-3,k-3) -   9*hyy(i,j-2,k-3) +   45*hyy(i,j-1,k-3) -   45*hyy(i,j+1,k-3) +   9*hyy(i,j+2,k-3) -    hyy(i,j+3,k-3) ) * odydz3600
      d2_hh23(2,3) = (    -hyz(i,j-3,k+3) +   9*hyz(i,j-2,k+3) -   45*hyz(i,j-1,k+3) +   45*hyz(i,j+1,k+3) -   9*hyz(i,j+2,k+3) +    hyz(i,j+3,k+3) &
                      +  9*hyz(i,j-3,k+2) -  81*hyz(i,j-2,k+2) +  405*hyz(i,j-1,k+2) -  405*hyz(i,j+1,k+2) +  81*hyz(i,j+2,k+2) -  9*hyz(i,j+3,k+2) &
                      - 45*hyz(i,j-3,k+1) + 405*hyz(i,j-2,k+1) - 2025*hyz(i,j-1,k+1) + 2025*hyz(i,j+1,k+1) - 405*hyz(i,j+2,k+1) + 45*hyz(i,j+3,k+1) &
                      + 45*hyz(i,j-3,k-1) - 405*hyz(i,j-2,k-1) + 2025*hyz(i,j-1,k-1) - 2025*hyz(i,j+1,k-1) + 405*hyz(i,j+2,k-1) - 45*hyz(i,j+3,k-1) &
                      -  9*hyz(i,j-3,k-2) +  81*hyz(i,j-2,k-2) -  405*hyz(i,j-1,k-2) +  405*hyz(i,j+1,k-2) -  81*hyz(i,j+2,k-2) +  9*hyz(i,j+3,k-2) &
                      +    hyz(i,j-3,k-3) -   9*hyz(i,j-2,k-3) +   45*hyz(i,j-1,k-3) -   45*hyz(i,j+1,k-3) +   9*hyz(i,j+2,k-3) -    hyz(i,j+3,k-3) ) * odydz3600
      d2_hh33(2,3) = (    -hzz(i,j-3,k+3) +   9*hzz(i,j-2,k+3) -   45*hzz(i,j-1,k+3) +   45*hzz(i,j+1,k+3) -   9*hzz(i,j+2,k+3) +    hzz(i,j+3,k+3) &
                        +  9*hzz(i,j-3,k+2) -  81*hzz(i,j-2,k+2) +  405*hzz(i,j-1,k+2) -  405*hzz(i,j+1,k+2) +  81*hzz(i,j+2,k+2) -  9*hzz(i,j+3,k+2) &
                        - 45*hzz(i,j-3,k+1) + 405*hzz(i,j-2,k+1) - 2025*hzz(i,j-1,k+1) + 2025*hzz(i,j+1,k+1) - 405*hzz(i,j+2,k+1) + 45*hzz(i,j+3,k+1) &
                        + 45*hzz(i,j-3,k-1) - 405*hzz(i,j-2,k-1) + 2025*hzz(i,j-1,k-1) - 2025*hzz(i,j+1,k-1) + 405*hzz(i,j+2,k-1) - 45*hzz(i,j+3,k-1) &
                        -  9*hzz(i,j-3,k-2) +  81*hzz(i,j-2,k-2) -  405*hzz(i,j-1,k-2) +  405*hzz(i,j+1,k-2) -  81*hzz(i,j+2,k-2) +  9*hzz(i,j+3,k-2) &
                        +    hzz(i,j-3,k-3) -   9*hzz(i,j-2,k-3) +   45*hzz(i,j-1,k-3) -   45*hzz(i,j+1,k-3) +   9*hzz(i,j+2,k-3) -    hzz(i,j+3,k-3) ) * odydz3600

      d2_hh11(2,1) = d2_hh11(1,2)
      d2_hh12(2,1) = d2_hh12(1,2)
      d2_hh13(2,1) = d2_hh13(1,2)
      d2_hh22(2,1) = d2_hh22(1,2)
      d2_hh23(2,1) = d2_hh23(1,2)
      d2_hh33(2,1) = d2_hh33(1,2)

      d2_hh11(3,1) = d2_hh11(1,3)
      d2_hh12(3,1) = d2_hh12(1,3)
      d2_hh13(3,1) = d2_hh13(1,3)
      d2_hh22(3,1) = d2_hh22(1,3)
      d2_hh23(3,1) = d2_hh23(1,3)
      d2_hh33(3,1) = d2_hh33(1,3)

      d2_hh11(3,2) = d2_hh11(2,3)
      d2_hh12(3,2) = d2_hh12(2,3)
      d2_hh13(3,2) = d2_hh13(2,3)
      d2_hh22(3,2) = d2_hh22(2,3)
      d2_hh23(3,2) = d2_hh23(2,3)
      d2_hh33(3,2) = d2_hh33(2,3)

    else ! if derivs_order == 4

      !------------ Centered 1st derivatives -----
      ! d1_ww(3)
      d1_ww(1) = (   -conf_fac(i+2,j,k) + 8*conf_fac(i+1,j,k)                          &
                  - 8*conf_fac(i-1,j,k) +   conf_fac(i-2,j,k) ) / dx12

      d1_ww(2) = (   -conf_fac(i,j+2,k) + 8*conf_fac(i,j+1,k)                          &
                  - 8*conf_fac(i,j-1,k) +   conf_fac(i,j-2,k) ) / dy12

      d1_ww(3) = (   -conf_fac(i,j,k+2) + 8*conf_fac(i,j,k+1)                          &
                  - 8*conf_fac(i,j,k-1) +   conf_fac(i,j,k-2) ) / dz12

      ! d1_hh(3,3,3)
      d1_hh11(1) = (   -hxx(i+2,j,k) + 8*hxx(i+1,j,k)                      &
                    - 8*hxx(i-1,j,k) +   hxx(i-2,j,k) ) / dx12
      d1_hh12(1) = (   -hxy(i+2,j,k) + 8*hxy(i+1,j,k)                      &
                    - 8*hxy(i-1,j,k) +   hxy(i-2,j,k) ) / dx12
      d1_hh13(1) = (   -hxz(i+2,j,k) + 8*hxz(i+1,j,k)                      &
                    - 8*hxz(i-1,j,k) +   hxz(i-2,j,k) ) / dx12
      d1_hh22(1) = (   -hyy(i+2,j,k) + 8*hyy(i+1,j,k)                      &
                    - 8*hyy(i-1,j,k) +   hyy(i-2,j,k) ) / dx12
      d1_hh23(1) = (   -hyz(i+2,j,k) + 8*hyz(i+1,j,k)                      &
                    - 8*hyz(i-1,j,k) +   hyz(i-2,j,k) ) / dx12
      d1_hh33(1) = (   -hzz(i+2,j,k) + 8*hzz(i+1,j,k)                      &
                    - 8*hzz(i-1,j,k) +   hzz(i-2,j,k) ) / dx12

      d1_hh11(2) = (   -hxx(i,j+2,k) + 8*hxx(i,j+1,k)                      &
                    - 8*hxx(i,j-1,k) +   hxx(i,j-2,k) ) / dy12
      d1_hh12(2) = (   -hxy(i,j+2,k) + 8*hxy(i,j+1,k)                      &
                    - 8*hxy(i,j-1,k) +   hxy(i,j-2,k) ) / dy12
      d1_hh13(2) = (   -hxz(i,j+2,k) + 8*hxz(i,j+1,k)                      &
                    - 8*hxz(i,j-1,k) +   hxz(i,j-2,k) ) / dy12
      d1_hh22(2) = (   -hyy(i,j+2,k) + 8*hyy(i,j+1,k)                      &
                    - 8*hyy(i,j-1,k) +   hyy(i,j-2,k) ) / dy12
      d1_hh23(2) = (   -hyz(i,j+2,k) + 8*hyz(i,j+1,k)                      &
                    - 8*hyz(i,j-1,k) +   hyz(i,j-2,k) ) / dy12
      d1_hh33(2) = (   -hzz(i,j+2,k) + 8*hzz(i,j+1,k)                      &
                    - 8*hzz(i,j-1,k) +   hzz(i,j-2,k) ) / dy12

      d1_hh11(3) = (   -hxx(i,j,k+2) + 8*hxx(i,j,k+1)                      &
                    - 8*hxx(i,j,k-1) +   hxx(i,j,k-2) ) / dz12
      d1_hh12(3) = (   -hxy(i,j,k+2) + 8*hxy(i,j,k+1)                      &
                    - 8*hxy(i,j,k-1) +   hxy(i,j,k-2) ) / dz12
      d1_hh13(3) = (   -hxz(i,j,k+2) + 8*hxz(i,j,k+1)                      &
                    - 8*hxz(i,j,k-1) +   hxz(i,j,k-2) ) / dz12
      d1_hh22(3) = (   -hyy(i,j,k+2) + 8*hyy(i,j,k+1)                      &
                    - 8*hyy(i,j,k-1) +   hyy(i,j,k-2) ) / dz12
      d1_hh23(3) = (   -hyz(i,j,k+2) + 8*hyz(i,j,k+1)                      &
                    - 8*hyz(i,j,k-1) +   hyz(i,j,k-2) ) / dz12
      d1_hh33(3) = (   -hzz(i,j,k+2) + 8*hzz(i,j,k+1)                      &
                    - 8*hzz(i,j,k-1) +   hzz(i,j,k-2) ) / dz12

      ! d1_trk(3)
      d1_trk(1) = (   -tracek(i+2,j,k) + 8*tracek(i+1,j,k)                   &
                   - 8*tracek(i-1,j,k) +   tracek(i-2,j,k) ) / dx12

      d1_trk(2) = (   -tracek(i,j+2,k) + 8*tracek(i,j+1,k)                   &
                   - 8*tracek(i,j-1,k) +   tracek(i,j-2,k) ) / dy12

      d1_trk(3) = (   -tracek(i,j,k+2) + 8*tracek(i,j,k+1)                   &
                   - 8*tracek(i,j,k-1) +   tracek(i,j,k-2) ) / dz12

      ! d1_aa(3,3,3)
      d1_aa11(1) = (   -axx(i+2,j,k) + 8*axx(i+1,j,k)                      &
                    - 8*axx(i-1,j,k) +   axx(i-2,j,k) ) / dx12
      d1_aa12(1) = (   -axy(i+2,j,k) + 8*axy(i+1,j,k)                      &
                    - 8*axy(i-1,j,k) +   axy(i-2,j,k) ) / dx12
      d1_aa13(1) = (   -axz(i+2,j,k) + 8*axz(i+1,j,k)                      &
                    - 8*axz(i-1,j,k) +   axz(i-2,j,k) ) / dx12
      d1_aa22(1) = (   -ayy(i+2,j,k) + 8*ayy(i+1,j,k)                      &
                    - 8*ayy(i-1,j,k) +   ayy(i-2,j,k) ) / dx12
      d1_aa23(1) = (   -ayz(i+2,j,k) + 8*ayz(i+1,j,k)                      &
                    - 8*ayz(i-1,j,k) +   ayz(i-2,j,k) ) / dx12
      d1_aa33(1) = (   -azz(i+2,j,k) + 8*azz(i+1,j,k)                      &
                    - 8*azz(i-1,j,k) +   azz(i-2,j,k) ) / dx12

      d1_aa11(2) = (   -axx(i,j+2,k) + 8*axx(i,j+1,k)                      &
                    - 8*axx(i,j-1,k) +   axx(i,j-2,k) ) / dy12
      d1_aa12(2) = (   -axy(i,j+2,k) + 8*axy(i,j+1,k)                      &
                    - 8*axy(i,j-1,k) +   axy(i,j-2,k) ) / dy12
      d1_aa13(2) = (   -axz(i,j+2,k) + 8*axz(i,j+1,k)                      &
                    - 8*axz(i,j-1,k) +   axz(i,j-2,k) ) / dy12
      d1_aa22(2) = (   -ayy(i,j+2,k) + 8*ayy(i,j+1,k)                      &
                    - 8*ayy(i,j-1,k) +   ayy(i,j-2,k) ) / dy12
      d1_aa23(2) = (   -ayz(i,j+2,k) + 8*ayz(i,j+1,k)                      &
                    - 8*ayz(i,j-1,k) +   ayz(i,j-2,k) ) / dy12
      d1_aa33(2) = (   -azz(i,j+2,k) + 8*azz(i,j+1,k)                      &
                    - 8*azz(i,j-1,k) +   azz(i,j-2,k) ) / dy12

      d1_aa11(3) = (   -axx(i,j,k+2) + 8*axx(i,j,k+1)                      &
                    - 8*axx(i,j,k-1) +   axx(i,j,k-2) ) / dz12
      d1_aa12(3) = (   -axy(i,j,k+2) + 8*axy(i,j,k+1)                      &
                    - 8*axy(i,j,k-1) +   axy(i,j,k-2) ) / dz12
      d1_aa13(3) = (   -axz(i,j,k+2) + 8*axz(i,j,k+1)                      &
                    - 8*axz(i,j,k-1) +   axz(i,j,k-2) ) / dz12
      d1_aa22(3) = (   -ayy(i,j,k+2) + 8*ayy(i,j,k+1)                      &
                    - 8*ayy(i,j,k-1) +   ayy(i,j,k-2) ) / dz12
      d1_aa23(3) = (   -ayz(i,j,k+2) + 8*ayz(i,j,k+1)                      &
                    - 8*ayz(i,j,k-1) +   ayz(i,j,k-2) ) / dz12
      d1_aa33(3) = (   -azz(i,j,k+2) + 8*azz(i,j,k+1)                      &
                    - 8*azz(i,j,k-1) +   azz(i,j,k-2) ) / dz12

      ! d1_gammat(3,3)
      d1_gammat1(1) = (   -gammatx(i+2,j,k) + 8*gammatx(i+1,j,k)               &
                       - 8*gammatx(i-1,j,k) +   gammatx(i-2,j,k) ) / dx12
      d1_gammat2(1) = (   -gammaty(i+2,j,k) + 8*gammaty(i+1,j,k)               &
                       - 8*gammaty(i-1,j,k) +   gammaty(i-2,j,k) ) / dx12
      d1_gammat3(1) = (   -gammatz(i+2,j,k) + 8*gammatz(i+1,j,k)               &
                       - 8*gammatz(i-1,j,k) +   gammatz(i-2,j,k) ) / dx12

      d1_gammat1(2) = (   -gammatx(i,j+2,k) + 8*gammatx(i,j+1,k)               &
                       - 8*gammatx(i,j-1,k) +   gammatx(i,j-2,k) ) / dy12
      d1_gammat2(2) = (   -gammaty(i,j+2,k) + 8*gammaty(i,j+1,k)               &
                       - 8*gammaty(i,j-1,k) +   gammaty(i,j-2,k) ) / dy12
      d1_gammat3(2) = (   -gammatz(i,j+2,k) + 8*gammatz(i,j+1,k)               &
                       - 8*gammatz(i,j-1,k) +   gammatz(i,j-2,k) ) / dy12

      d1_gammat1(3) = (   -gammatx(i,j,k+2) + 8*gammatx(i,j,k+1)               &
                       - 8*gammatx(i,j,k-1) +   gammatx(i,j,k-2) ) / dz12
      d1_gammat2(3) = (   -gammaty(i,j,k+2) + 8*gammaty(i,j,k+1)               &
                       - 8*gammaty(i,j,k-1) +   gammaty(i,j,k-2) ) / dz12
      d1_gammat3(3) = (   -gammatz(i,j,k+2) + 8*gammatz(i,j,k+1)               &
                       - 8*gammatz(i,j,k-1) +   gammatz(i,j,k-2) ) / dz12

      !--------------------------------------------------


      !------------- Centered 2nd derivatives -----------

      ! d2_ww(3,3)
      d2_ww(1,1) = (   -conf_fac(i+2,j,k) + 16*conf_fac(i+1,j,k) - 30*conf_fac(i,j,k)       &
                   + 16*conf_fac(i-1,j,k) -    conf_fac(i-2,j,k) ) / dxsq12

      d2_ww(2,2) = (   -conf_fac(i,j+2,k) + 16*conf_fac(i,j+1,k) - 30*conf_fac(i,j,k)       &
                   + 16*conf_fac(i,j-1,k) -    conf_fac(i,j-2,k) ) / dysq12

      d2_ww(3,3) = (   -conf_fac(i,j,k+2) + 16*conf_fac(i,j,k+1) - 30*conf_fac(i,j,k)       &
                   + 16*conf_fac(i,j,k-1) -    conf_fac(i,j,k-2) ) / dzsq12

      d2_ww(1,2) = (   -conf_fac(i-2,j+2,k) +  8*conf_fac(i-1,j+2,k) -  8*conf_fac(i+1,j+2,k) +   conf_fac(i+2,j+2,k)   &
                    + 8*conf_fac(i-2,j+1,k) - 64*conf_fac(i-1,j+1,k) + 64*conf_fac(i+1,j+1,k) - 8*conf_fac(i+2,j+1,k)   &
                    - 8*conf_fac(i-2,j-1,k) + 64*conf_fac(i-1,j-1,k) - 64*conf_fac(i+1,j-1,k) + 8*conf_fac(i+2,j-1,k)   &
                    +   conf_fac(i-2,j-2,k) -  8*conf_fac(i-1,j-2,k) +  8*conf_fac(i+1,j-2,k) -   conf_fac(i+2,j-2,k) ) / dxdy144

      d2_ww(1,3) = (   -conf_fac(i-2,j,k+2) +  8*conf_fac(i-1,j,k+2) -  8*conf_fac(i+1,j,k+2) +   conf_fac(i+2,j,k+2)   &
                    + 8*conf_fac(i-2,j,k+1) - 64*conf_fac(i-1,j,k+1) + 64*conf_fac(i+1,j,k+1) - 8*conf_fac(i+2,j,k+1)   &
                    - 8*conf_fac(i-2,j,k-1) + 64*conf_fac(i-1,j,k-1) - 64*conf_fac(i+1,j,k-1) + 8*conf_fac(i+2,j,k-1)   &
                    +   conf_fac(i-2,j,k-2) -  8*conf_fac(i-1,j,k-2) +  8*conf_fac(i+1,j,k-2) -   conf_fac(i+2,j,k-2) ) / dxdz144

      d2_ww(2,3) = (   -conf_fac(i,j-2,k+2) +  8*conf_fac(i,j-1,k+2) -  8*conf_fac(i,j+1,k+2) +   conf_fac(i,j+2,k+2)   &
                    + 8*conf_fac(i,j-2,k+1) - 64*conf_fac(i,j-1,k+1) + 64*conf_fac(i,j+1,k+1) - 8*conf_fac(i,j+2,k+1)   &
                    - 8*conf_fac(i,j-2,k-1) + 64*conf_fac(i,j-1,k-1) - 64*conf_fac(i,j+1,k-1) + 8*conf_fac(i,j+2,k-1)   &
                    +   conf_fac(i,j-2,k-2) -  8*conf_fac(i,j-1,k-2) +  8*conf_fac(i,j+1,k-2) -   conf_fac(i,j+2,k-2) ) / dydz144

      d2_ww(2,1) = d2_ww(1,2)
      d2_ww(3,1) = d2_ww(1,3)
      d2_ww(3,2) = d2_ww(2,3)

      ! d2_hh(3,3,3,3)
      d2_hh11(1,1) = (   -hxx(i+2,j,k) + 16*hxx(i+1,j,k) - 30*hxx(i,j,k)   &
                     + 16*hxx(i-1,j,k) -    hxx(i-2,j,k) ) / dxsq12
      d2_hh12(1,1) = (   -hxy(i+2,j,k) + 16*hxy(i+1,j,k) - 30*hxy(i,j,k)   &
                     + 16*hxy(i-1,j,k) -    hxy(i-2,j,k) ) / dxsq12
      d2_hh13(1,1) = (   -hxz(i+2,j,k) + 16*hxz(i+1,j,k) - 30*hxz(i,j,k)   &
                     + 16*hxz(i-1,j,k) -    hxz(i-2,j,k) ) / dxsq12
      d2_hh22(1,1) = (   -hyy(i+2,j,k) + 16*hyy(i+1,j,k) - 30*hyy(i,j,k)   &
                     + 16*hyy(i-1,j,k) -    hyy(i-2,j,k) ) / dxsq12
      d2_hh23(1,1) = (   -hyz(i+2,j,k) + 16*hyz(i+1,j,k) - 30*hyz(i,j,k)   &
                     + 16*hyz(i-1,j,k) -    hyz(i-2,j,k) ) / dxsq12
      d2_hh33(1,1) = (   -hzz(i+2,j,k) + 16*hzz(i+1,j,k) - 30*hzz(i,j,k)   &
                     + 16*hzz(i-1,j,k) -    hzz(i-2,j,k) ) / dxsq12

      d2_hh11(2,2) = (   -hxx(i,j+2,k) + 16*hxx(i,j+1,k) - 30*hxx(i,j,k)   &
                     + 16*hxx(i,j-1,k) -    hxx(i,j-2,k) ) / dysq12
      d2_hh12(2,2) = (   -hxy(i,j+2,k) + 16*hxy(i,j+1,k) - 30*hxy(i,j,k)   &
                     + 16*hxy(i,j-1,k) -    hxy(i,j-2,k) ) / dysq12
      d2_hh13(2,2) = (   -hxz(i,j+2,k) + 16*hxz(i,j+1,k) - 30*hxz(i,j,k)   &
                     + 16*hxz(i,j-1,k) -    hxz(i,j-2,k) ) / dysq12
      d2_hh22(2,2) = (   -hyy(i,j+2,k) + 16*hyy(i,j+1,k) - 30*hyy(i,j,k)   &
                     + 16*hyy(i,j-1,k) -    hyy(i,j-2,k) ) / dysq12
      d2_hh23(2,2) = (   -hyz(i,j+2,k) + 16*hyz(i,j+1,k) - 30*hyz(i,j,k)   &
                     + 16*hyz(i,j-1,k) -    hyz(i,j-2,k) ) / dysq12
      d2_hh33(2,2) = (   -hzz(i,j+2,k) + 16*hzz(i,j+1,k) - 30*hzz(i,j,k)   &
                     + 16*hzz(i,j-1,k) -    hzz(i,j-2,k) ) / dysq12

      d2_hh11(3,3) = (   -hxx(i,j,k+2) + 16*hxx(i,j,k+1) - 30*hxx(i,j,k)   &
                     + 16*hxx(i,j,k-1) -    hxx(i,j,k-2) ) / dzsq12
      d2_hh12(3,3) = (   -hxy(i,j,k+2) + 16*hxy(i,j,k+1) - 30*hxy(i,j,k)   &
                     + 16*hxy(i,j,k-1) -    hxy(i,j,k-2) ) / dzsq12
      d2_hh13(3,3) = (   -hxz(i,j,k+2) + 16*hxz(i,j,k+1) - 30*hxz(i,j,k)   &
                     + 16*hxz(i,j,k-1) -    hxz(i,j,k-2) ) / dzsq12
      d2_hh22(3,3) = (   -hyy(i,j,k+2) + 16*hyy(i,j,k+1) - 30*hyy(i,j,k)   &
                     + 16*hyy(i,j,k-1) -    hyy(i,j,k-2) ) / dzsq12
      d2_hh23(3,3) = (   -hyz(i,j,k+2) + 16*hyz(i,j,k+1) - 30*hyz(i,j,k)   &
                     + 16*hyz(i,j,k-1) -    hyz(i,j,k-2) ) / dzsq12
      d2_hh33(3,3) = (   -hzz(i,j,k+2) + 16*hzz(i,j,k+1) - 30*hzz(i,j,k)   &
                     + 16*hzz(i,j,k-1) -    hzz(i,j,k-2) ) / dzsq12

      d2_hh11(1,2) = (   -hxx(i-2,j+2,k) +  8*hxx(i-1,j+2,k) -  8*hxx(i+1,j+2,k) +   hxx(i+2,j+2,k)   &
                      + 8*hxx(i-2,j+1,k) - 64*hxx(i-1,j+1,k) + 64*hxx(i+1,j+1,k) - 8*hxx(i+2,j+1,k)   &
                      - 8*hxx(i-2,j-1,k) + 64*hxx(i-1,j-1,k) - 64*hxx(i+1,j-1,k) + 8*hxx(i+2,j-1,k)   &
                      +   hxx(i-2,j-2,k) -  8*hxx(i-1,j-2,k) +  8*hxx(i+1,j-2,k) -   hxx(i+2,j-2,k) ) / dxdy144
      d2_hh12(1,2) = (   -hxy(i-2,j+2,k) +  8*hxy(i-1,j+2,k) -  8*hxy(i+1,j+2,k) +   hxy(i+2,j+2,k)   &
                      + 8*hxy(i-2,j+1,k) - 64*hxy(i-1,j+1,k) + 64*hxy(i+1,j+1,k) - 8*hxy(i+2,j+1,k)   &
                      - 8*hxy(i-2,j-1,k) + 64*hxy(i-1,j-1,k) - 64*hxy(i+1,j-1,k) + 8*hxy(i+2,j-1,k)   &
                      +   hxy(i-2,j-2,k) -  8*hxy(i-1,j-2,k) +  8*hxy(i+1,j-2,k) -   hxy(i+2,j-2,k) ) / dxdy144
      d2_hh13(1,2) = (   -hxz(i-2,j+2,k) +  8*hxz(i-1,j+2,k) -  8*hxz(i+1,j+2,k) +   hxz(i+2,j+2,k)   &
                      + 8*hxz(i-2,j+1,k) - 64*hxz(i-1,j+1,k) + 64*hxz(i+1,j+1,k) - 8*hxz(i+2,j+1,k)   &
                      - 8*hxz(i-2,j-1,k) + 64*hxz(i-1,j-1,k) - 64*hxz(i+1,j-1,k) + 8*hxz(i+2,j-1,k)   &
                      +   hxz(i-2,j-2,k) -  8*hxz(i-1,j-2,k) +  8*hxz(i+1,j-2,k) -   hxz(i+2,j-2,k) ) / dxdy144
      d2_hh22(1,2) = (   -hyy(i-2,j+2,k) +  8*hyy(i-1,j+2,k) -  8*hyy(i+1,j+2,k) +   hyy(i+2,j+2,k)   &
                      + 8*hyy(i-2,j+1,k) - 64*hyy(i-1,j+1,k) + 64*hyy(i+1,j+1,k) - 8*hyy(i+2,j+1,k)   &
                      - 8*hyy(i-2,j-1,k) + 64*hyy(i-1,j-1,k) - 64*hyy(i+1,j-1,k) + 8*hyy(i+2,j-1,k)   &
                      +   hyy(i-2,j-2,k) -  8*hyy(i-1,j-2,k) +  8*hyy(i+1,j-2,k) -   hyy(i+2,j-2,k) ) / dxdy144
      d2_hh23(1,2) = (   -hyz(i-2,j+2,k) +  8*hyz(i-1,j+2,k) -  8*hyz(i+1,j+2,k) +   hyz(i+2,j+2,k)   &
                      + 8*hyz(i-2,j+1,k) - 64*hyz(i-1,j+1,k) + 64*hyz(i+1,j+1,k) - 8*hyz(i+2,j+1,k)   &
                      - 8*hyz(i-2,j-1,k) + 64*hyz(i-1,j-1,k) - 64*hyz(i+1,j-1,k) + 8*hyz(i+2,j-1,k)   &
                      +   hyz(i-2,j-2,k) -  8*hyz(i-1,j-2,k) +  8*hyz(i+1,j-2,k) -   hyz(i+2,j-2,k) ) / dxdy144
      d2_hh33(1,2) = (   -hzz(i-2,j+2,k) +  8*hzz(i-1,j+2,k) -  8*hzz(i+1,j+2,k) +   hzz(i+2,j+2,k)   &
                      + 8*hzz(i-2,j+1,k) - 64*hzz(i-1,j+1,k) + 64*hzz(i+1,j+1,k) - 8*hzz(i+2,j+1,k)   &
                      - 8*hzz(i-2,j-1,k) + 64*hzz(i-1,j-1,k) - 64*hzz(i+1,j-1,k) + 8*hzz(i+2,j-1,k)   &
                      +   hzz(i-2,j-2,k) -  8*hzz(i-1,j-2,k) +  8*hzz(i+1,j-2,k) -   hzz(i+2,j-2,k) ) / dxdy144

      d2_hh11(1,3) = (   -hxx(i-2,j,k+2) +  8*hxx(i-1,j,k+2) -  8*hxx(i+1,j,k+2) +   hxx(i+2,j,k+2)   &
                      + 8*hxx(i-2,j,k+1) - 64*hxx(i-1,j,k+1) + 64*hxx(i+1,j,k+1) - 8*hxx(i+2,j,k+1)   &
                      - 8*hxx(i-2,j,k-1) + 64*hxx(i-1,j,k-1) - 64*hxx(i+1,j,k-1) + 8*hxx(i+2,j,k-1)   &
                      +   hxx(i-2,j,k-2) -  8*hxx(i-1,j,k-2) +  8*hxx(i+1,j,k-2) -   hxx(i+2,j,k-2) ) / dxdz144
      d2_hh12(1,3) = (   -hxy(i-2,j,k+2) +  8*hxy(i-1,j,k+2) -  8*hxy(i+1,j,k+2) +   hxy(i+2,j,k+2)   &
                      + 8*hxy(i-2,j,k+1) - 64*hxy(i-1,j,k+1) + 64*hxy(i+1,j,k+1) - 8*hxy(i+2,j,k+1)   &
                      - 8*hxy(i-2,j,k-1) + 64*hxy(i-1,j,k-1) - 64*hxy(i+1,j,k-1) + 8*hxy(i+2,j,k-1)   &
                      +   hxy(i-2,j,k-2) -  8*hxy(i-1,j,k-2) +  8*hxy(i+1,j,k-2) -   hxy(i+2,j,k-2) ) / dxdz144
      d2_hh13(1,3) = (   -hxz(i-2,j,k+2) +  8*hxz(i-1,j,k+2) -  8*hxz(i+1,j,k+2) +   hxz(i+2,j,k+2)   &
                      + 8*hxz(i-2,j,k+1) - 64*hxz(i-1,j,k+1) + 64*hxz(i+1,j,k+1) - 8*hxz(i+2,j,k+1)   &
                      - 8*hxz(i-2,j,k-1) + 64*hxz(i-1,j,k-1) - 64*hxz(i+1,j,k-1) + 8*hxz(i+2,j,k-1)   &
                      +   hxz(i-2,j,k-2) -  8*hxz(i-1,j,k-2) +  8*hxz(i+1,j,k-2) -   hxz(i+2,j,k-2) ) / dxdz144
      d2_hh22(1,3) = (   -hyy(i-2,j,k+2) +  8*hyy(i-1,j,k+2) -  8*hyy(i+1,j,k+2) +   hyy(i+2,j,k+2)   &
                      + 8*hyy(i-2,j,k+1) - 64*hyy(i-1,j,k+1) + 64*hyy(i+1,j,k+1) - 8*hyy(i+2,j,k+1)   &
                      - 8*hyy(i-2,j,k-1) + 64*hyy(i-1,j,k-1) - 64*hyy(i+1,j,k-1) + 8*hyy(i+2,j,k-1)   &
                      +   hyy(i-2,j,k-2) -  8*hyy(i-1,j,k-2) +  8*hyy(i+1,j,k-2) -   hyy(i+2,j,k-2) ) / dxdz144
      d2_hh23(1,3) = (   -hyz(i-2,j,k+2) +  8*hyz(i-1,j,k+2) -  8*hyz(i+1,j,k+2) +   hyz(i+2,j,k+2)   &
                      + 8*hyz(i-2,j,k+1) - 64*hyz(i-1,j,k+1) + 64*hyz(i+1,j,k+1) - 8*hyz(i+2,j,k+1)   &
                      - 8*hyz(i-2,j,k-1) + 64*hyz(i-1,j,k-1) - 64*hyz(i+1,j,k-1) + 8*hyz(i+2,j,k-1)   &
                      +   hyz(i-2,j,k-2) -  8*hyz(i-1,j,k-2) +  8*hyz(i+1,j,k-2) -   hyz(i+2,j,k-2) ) / dxdz144
      d2_hh33(1,3) = (   -hzz(i-2,j,k+2) +  8*hzz(i-1,j,k+2) -  8*hzz(i+1,j,k+2) +   hzz(i+2,j,k+2)   &
                      + 8*hzz(i-2,j,k+1) - 64*hzz(i-1,j,k+1) + 64*hzz(i+1,j,k+1) - 8*hzz(i+2,j,k+1)   &
                      - 8*hzz(i-2,j,k-1) + 64*hzz(i-1,j,k-1) - 64*hzz(i+1,j,k-1) + 8*hzz(i+2,j,k-1)   &
                      +   hzz(i-2,j,k-2) -  8*hzz(i-1,j,k-2) +  8*hzz(i+1,j,k-2) -   hzz(i+2,j,k-2) ) / dxdz144

      d2_hh11(2,3) = (   -hxx(i,j-2,k+2) +  8*hxx(i,j-1,k+2) -  8*hxx(i,j+1,k+2) +   hxx(i,j+2,k+2)   &
                      + 8*hxx(i,j-2,k+1) - 64*hxx(i,j-1,k+1) + 64*hxx(i,j+1,k+1) - 8*hxx(i,j+2,k+1)   &
                      - 8*hxx(i,j-2,k-1) + 64*hxx(i,j-1,k-1) - 64*hxx(i,j+1,k-1) + 8*hxx(i,j+2,k-1)   &
                      +   hxx(i,j-2,k-2) -  8*hxx(i,j-1,k-2) +  8*hxx(i,j+1,k-2) -   hxx(i,j+2,k-2) ) / dydz144
      d2_hh12(2,3) = (   -hxy(i,j-2,k+2) +  8*hxy(i,j-1,k+2) -  8*hxy(i,j+1,k+2) +   hxy(i,j+2,k+2)   &
                      + 8*hxy(i,j-2,k+1) - 64*hxy(i,j-1,k+1) + 64*hxy(i,j+1,k+1) - 8*hxy(i,j+2,k+1)   &
                      - 8*hxy(i,j-2,k-1) + 64*hxy(i,j-1,k-1) - 64*hxy(i,j+1,k-1) + 8*hxy(i,j+2,k-1)   &
                      +   hxy(i,j-2,k-2) -  8*hxy(i,j-1,k-2) +  8*hxy(i,j+1,k-2) -   hxy(i,j+2,k-2) ) / dydz144
      d2_hh13(2,3) = (   -hxz(i,j-2,k+2) +  8*hxz(i,j-1,k+2) -  8*hxz(i,j+1,k+2) +   hxz(i,j+2,k+2)   &
                      + 8*hxz(i,j-2,k+1) - 64*hxz(i,j-1,k+1) + 64*hxz(i,j+1,k+1) - 8*hxz(i,j+2,k+1)   &
                      - 8*hxz(i,j-2,k-1) + 64*hxz(i,j-1,k-1) - 64*hxz(i,j+1,k-1) + 8*hxz(i,j+2,k-1)   &
                      +   hxz(i,j-2,k-2) -  8*hxz(i,j-1,k-2) +  8*hxz(i,j+1,k-2) -   hxz(i,j+2,k-2) ) / dydz144
      d2_hh22(2,3) = (   -hyy(i,j-2,k+2) +  8*hyy(i,j-1,k+2) -  8*hyy(i,j+1,k+2) +   hyy(i,j+2,k+2)   &
                      + 8*hyy(i,j-2,k+1) - 64*hyy(i,j-1,k+1) + 64*hyy(i,j+1,k+1) - 8*hyy(i,j+2,k+1)   &
                      - 8*hyy(i,j-2,k-1) + 64*hyy(i,j-1,k-1) - 64*hyy(i,j+1,k-1) + 8*hyy(i,j+2,k-1)   &
                      +   hyy(i,j-2,k-2) -  8*hyy(i,j-1,k-2) +  8*hyy(i,j+1,k-2) -   hyy(i,j+2,k-2) ) / dydz144
      d2_hh23(2,3) = (   -hyz(i,j-2,k+2) +  8*hyz(i,j-1,k+2) -  8*hyz(i,j+1,k+2) +   hyz(i,j+2,k+2)   &
                      + 8*hyz(i,j-2,k+1) - 64*hyz(i,j-1,k+1) + 64*hyz(i,j+1,k+1) - 8*hyz(i,j+2,k+1)   &
                      - 8*hyz(i,j-2,k-1) + 64*hyz(i,j-1,k-1) - 64*hyz(i,j+1,k-1) + 8*hyz(i,j+2,k-1)   &
                      +   hyz(i,j-2,k-2) -  8*hyz(i,j-1,k-2) +  8*hyz(i,j+1,k-2) -   hyz(i,j+2,k-2) ) / dydz144
      d2_hh33(2,3) = (   -hzz(i,j-2,k+2) +  8*hzz(i,j-1,k+2) -  8*hzz(i,j+1,k+2) +   hzz(i,j+2,k+2)   &
                      + 8*hzz(i,j-2,k+1) - 64*hzz(i,j-1,k+1) + 64*hzz(i,j+1,k+1) - 8*hzz(i,j+2,k+1)   &
                      - 8*hzz(i,j-2,k-1) + 64*hzz(i,j-1,k-1) - 64*hzz(i,j+1,k-1) + 8*hzz(i,j+2,k-1)   &
                      +   hzz(i,j-2,k-2) -  8*hzz(i,j-1,k-2) +  8*hzz(i,j+1,k-2) -   hzz(i,j+2,k-2) ) / dydz144

      d2_hh11(2,1) = d2_hh11(1,2)
      d2_hh12(2,1) = d2_hh12(1,2)
      d2_hh13(2,1) = d2_hh13(1,2)
      d2_hh22(2,1) = d2_hh22(1,2)
      d2_hh23(2,1) = d2_hh23(1,2)
      d2_hh33(2,1) = d2_hh33(1,2)

      d2_hh11(3,1) = d2_hh11(1,3)
      d2_hh12(3,1) = d2_hh12(1,3)
      d2_hh13(3,1) = d2_hh13(1,3)
      d2_hh22(3,1) = d2_hh22(1,3)
      d2_hh23(3,1) = d2_hh23(1,3)
      d2_hh33(3,1) = d2_hh33(1,3)

      d2_hh11(3,2) = d2_hh11(2,3)
      d2_hh12(3,2) = d2_hh12(2,3)
      d2_hh13(3,2) = d2_hh13(2,3)
      d2_hh22(3,2) = d2_hh22(2,3)
      d2_hh23(3,2) = d2_hh23(2,3)
      d2_hh33(3,2) = d2_hh33(2,3)

    end if
    !--------------------------------------------------

    if (use_jacobian) then
      call LeanBSSN_apply_jacobian(d1_trk, jac)
      call LeanBSSN_apply_jacobian(d1_gammat1, jac)
      call LeanBSSN_apply_jacobian(d1_gammat2, jac)
      call LeanBSSN_apply_jacobian(d1_gammat3, jac)
      call LeanBSSN_apply_jacobian(d1_aa11, jac)
      call LeanBSSN_apply_jacobian(d1_aa12, jac)
      call LeanBSSN_apply_jacobian(d1_aa13, jac)
      call LeanBSSN_apply_jacobian(d1_aa22, jac)
      call LeanBSSN_apply_jacobian(d1_aa23, jac)
      call LeanBSSN_apply_jacobian(d1_aa33, jac)

      call LeanBSSN_apply_jacobian2(d1_ww, d2_ww, jac, hes)
      call LeanBSSN_apply_jacobian2(d1_hh11, d2_hh11, jac, hes)
      call LeanBSSN_apply_jacobian2(d1_hh12, d2_hh12, jac, hes)
      call LeanBSSN_apply_jacobian2(d1_hh13, d2_hh13, jac, hes)
      call LeanBSSN_apply_jacobian2(d1_hh22, d2_hh22, jac, hes)
      call LeanBSSN_apply_jacobian2(d1_hh23, d2_hh23, jac, hes)
      call LeanBSSN_apply_jacobian2(d1_hh33, d2_hh33, jac, hes)
    end if

    d1_gammat(1,:) = d1_gammat1(:)
    d1_gammat(2,:) = d1_gammat2(:)
    d1_gammat(3,:) = d1_gammat3(:)


    d1_hh(1,1,:) = d1_hh11(:)
    d1_hh(1,2,:) = d1_hh12(:)
    d1_hh(1,3,:) = d1_hh13(:)
    d1_hh(2,2,:) = d1_hh22(:)
    d1_hh(2,3,:) = d1_hh23(:)
    d1_hh(3,3,:) = d1_hh33(:)
    d1_hh(2,1,:) = d1_hh(1,2,:)
    d1_hh(3,1,:) = d1_hh(1,3,:)
    d1_hh(3,2,:) = d1_hh(2,3,:)

    d1_aa(1,1,:) = d1_aa11(:)
    d1_aa(1,2,:) = d1_aa12(:)
    d1_aa(1,3,:) = d1_aa13(:)
    d1_aa(2,2,:) = d1_aa22(:)
    d1_aa(2,3,:) = d1_aa23(:)
    d1_aa(3,3,:) = d1_aa33(:)
    d1_aa(2,1,:) = d1_aa(1,2,:)
    d1_aa(3,1,:) = d1_aa(1,3,:)
    d1_aa(3,2,:) = d1_aa(2,3,:)

    d2_hh(1,1,:,:) = d2_hh11(:,:)
    d2_hh(1,2,:,:) = d2_hh12(:,:)
    d2_hh(1,3,:,:) = d2_hh13(:,:)
    d2_hh(2,2,:,:) = d2_hh22(:,:)
    d2_hh(2,3,:,:) = d2_hh23(:,:)
    d2_hh(3,3,:,:) = d2_hh33(:,:)
    d2_hh(2,1,:,:) = d2_hh(1,2,:,:)
    d2_hh(3,1,:,:) = d2_hh(1,3,:,:)
    d2_hh(3,2,:,:) = d2_hh(2,3,:,:)
    !------------ Christoffel symbols ----------
    cf1 = 0
    do a = 1, 3
      do b = 1, 3
        do c = b, 3
          cf1(a,b,c) = 0.5d0 * (d1_hh(a,b,c) + d1_hh(a,c,b) - d1_hh(b,c,a))
        end do
      end do
    end do
    cf1(:,2,1) = cf1(:,1,2)
    cf1(:,3,1) = cf1(:,1,3)
    cf1(:,3,2) = cf1(:,2,3)

    cf2 = 0
    do a = 1, 3
      do b = 1, 3
        do c = b, 3
          do m = 1, 3
            cf2(a,b,c) = cf2(a,b,c) + hu(a,m) * cf1(m,b,c)
          end do
        end do
      end do
    end do
    cf2(:,2,1) = cf2(:,1,2)
    cf2(:,3,1) = cf2(:,1,3)
    cf2(:,3,2) = cf2(:,2,3)
    !-------------------------------------------


    !------------ Covariant derivatives --------
    cd2_ww = d2_ww
    cd1_aa   = d1_aa
    do a = 1, 3
      do b = a, 3
        do l = 1, 3
          cd2_ww(a,b) = cd2_ww(a,b) - cf2(l,a,b) * d1_ww(l)
        do m = 1, 3
            cd1_aa(a,b,l) = cd1_aa(a,b,l) - cf2(m,a,l) * aa(b,m) - cf2(m,b,l) * aa(a,m)
        end do
          end do
        end do
      end do
    cd2_ww(2,1)   = cd2_ww(1,2)
    cd2_ww(3,1)   = cd2_ww(1,3)
    cd2_ww(3,2)   = cd2_ww(2,3)
    cd1_aa(2,1,:) = cd1_aa(1,2,:)
    cd1_aa(3,1,:) = cd1_aa(1,3,:)
    cd1_aa(3,2,:) = cd1_aa(2,3,:)
    !-------------------------------------------


    !------------ Ricci Tensor -----------------
    ! Note: we are implementing W^2 R_{ij} here
    ri_1 = 0
    ri_2 = 0
    ri_3 = 0
    c_ri_ww = 0
    c_ri_hh = 0

    tr_cd2_ww = 0
    tr_dww_dww = 0
    do l = 1, 3
      do m = 1, 3
        tr_cd2_ww  = tr_cd2_ww  + hu(l,m) * cd2_ww(l,m)
        tr_dww_dww = tr_dww_dww + hu(l,m) * d1_ww(l) * d1_ww(m)
      end do
    end do
    ! Note: we are implementing W^2 R_{ij} here
    c_ri_ww = ww * ( cd2_ww + hh * tr_cd2_ww ) - 2 * hh * tr_dww_dww


    do a = 1, 3
      do b = a, 3
        do l = 1, 3
          ri_1(a,b) = ri_1(a,b) + hh(l,a) * d1_gammat(l,b) / 2              &
                                + hh(l,b) * d1_gammat(l,a) / 2              &
                                + gammat(l) * (cf1(a,b,l) + cf1(b,a,l)) / 2
          do m = 1, 3
            ri_2(a,b) = ri_2(a,b) - hu(l,m) * d2_hh(a,b,l,m) / 2
            do n = 1, 3
              ri_3(a,b) = ri_3(a,b) + hu(l,m)                              &
                                    * ( cf2(n,l,a) * cf1(b,n,m)            &
                                    + cf2(n,l,b) * (cf1(a,n,m) + cf1(n,m,a)) )
            end do
          end do
        end do
      end do
    end do
    ! Note: we are implementing W^2 R_{ij} here
    c_ri_hh = ww*ww * (ri_1 + ri_2 + ri_3)
    c_ri    = c_ri_ww + c_ri_hh

    c_ri(2,1) = c_ri(1,2)
    c_ri(3,1) = c_ri(1,3)
    c_ri(3,2) = c_ri(2,3)

    trr = 0
    do m = 1, 3
      do n = 1, 3
        trr = trr + hu(m,n) * c_ri(m,n)
      end do
    end do
    !-------------------------------------------


    !------------ Source terms -----------------
    sq_aa = 0
    a2    = 0
    do m = 1, 3
      do n = 1, 3
        do p = 1, 3
          do q = 1, 3
            a2(m,n) = a2(m,n) + hu(p,q) * aa(m,p) * aa(n,q)
          end do
        end do
        sq_aa = sq_aa + hu(m,n) * a2(m,n)
      end do
    end do
    !-------------------------------------------


    !------------ Constraints ------------------
    ham = trr + 2 * trk**2 / 3 - sq_aa

    mom = -2 * d1_trk / 3
    do a = 1, 3
      do l = 1, 3
        do m = 1, 3
          mom(a) = mom(a) + hu(l,m) * ( cd1_aa(a,l,m) - 3 * aa(a,l) * d1_ww(m) / ww )
        end do
      end do
    end do
    !-------------------------------------------


    !------------ Matter terms -----------------
    ! n_mu = (-alph, 0, 0, 0)
    ! n^mu = (1, -betax, -betay, -betaz)/alph
    !
    ! E   = n^mu n^nu T_{mu nu}
    !     = (T_{00} - 2 beta^i T_{i0} + beta^i beta^j T_{ij})/alph^2
    !
    ! j_a = -h_a^mu n^nu T_{mu nu}
    !     = -(T_{a 0} - beta^j T_{a j})/alph
    ! stress-energy tensor variables
    Tab = 0
    if (stress_energy_state /= 0) then
       Tab(4,4) = eTtt(i,j,k)
       Tab(4,1) = eTtx(i,j,k)
       Tab(4,2) = eTty(i,j,k)
       Tab(4,3) = eTtz(i,j,k)
       Tab(1,1) = eTxx(i,j,k)
       Tab(1,2) = eTxy(i,j,k)
       Tab(1,3) = eTxz(i,j,k)
       Tab(2,2) = eTyy(i,j,k)
       Tab(2,3) = eTyz(i,j,k)
       Tab(3,3) = eTzz(i,j,k)
       Tab(1,4) = Tab(4,1)
       Tab(2,4) = Tab(4,2)
       Tab(3,4) = Tab(4,3)
       Tab(2,1) = Tab(1,2)
       Tab(3,1) = Tab(1,3)
       Tab(3,2) = Tab(2,3)

       srcE = Tab(4,4)
       do m = 1, 3
          srcE = srcE - 2 * beta(m) * Tab(m,4)
          do n = 1, 3
             srcE = srcE + beta(m) * beta(n) * Tab(m,n)
          end do
       end do
       srcE = srcE / (alph * alph)

       srcjdi = 0
       do a = 1, 3
          do m = 1, 3
             srcjdi(a) = srcjdi(a) + beta(m) * Tab(a,m)
          end do
       end do
       srcjdi = (srcjdi - Tab(1:3,4)) / alph

       !------------ Correct source terms ---------
       ham = ham - pi16 * srcE
       mom = mom - pi8  * srcjdi

    end if

    !------------ Write to grid functions ------
    hc(i,j,k)  = ham

    mcx(i,j,k) = mom(1)
    mcy(i,j,k) = mom(2)
    mcz(i,j,k) = mom(3)
    !-------------------------------------------

  end do
  end do
  end do
  !$OMP END PARALLEL DO

end subroutine LeanBSSN_bssn_constraints
