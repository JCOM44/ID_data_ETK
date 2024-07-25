! calc_bssn_rhs.F90 : Calculate the right hand side of the BSSN equations
!
! 1. Set local variables
! 2. Calculate derivatives
! 3. Calculate covariant derivatives
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


subroutine LeanBSSN_calc_bssn_rhs( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ! Fundamental variables
  CCTK_REAL                alph, beta(3), beta_l(3)
  CCTK_REAL                ww, hh(3,3), hu(3,3), trk, aa(3,3), gammat(3),  &
                           au(3,3), dethh, Tab(4,4)

  ! Scalar field
  CCTK_REAL                lphi, lKphi, Bphi, BKphi, lphi2
  CCTK_REAL                d1_lphi(3), d1_lKphi(3)
  CCTK_REAL                omega_Bphi, domega_Bphi, F_phi
  CCTK_REAL                d2_lphi(3,3), cd2_lphi(3,3), cW_dphi(3,3) 
  CCTK_REAL                trk_phi,trk_omega, aa_omega(3,3), aa_phi(3,3), gammat_omega(3), gammat_phi(3)

  ! First derivatives
  CCTK_REAL                d1_beta1(3), d1_beta2(3), d1_beta3(3)
  CCTK_REAL                d1_hh11(3), d1_hh12(3), d1_hh13(3), d1_hh22(3), d1_hh23(3), d1_hh33(3)
  CCTK_REAL                d1_aa11(3), d1_aa12(3), d1_aa13(3), d1_aa22(3), d1_aa23(3), d1_aa33(3)
  CCTK_REAL                d1_gammat1(3), d1_gammat2(3), d1_gammat3(3)
  CCTK_REAL                d1_alph(3), d1_beta(3,3)
  CCTK_REAL                d1_ww(3), d1_hh(3,3,3), d1_trk(3), d1_aa(3,3,3),&
                           d1_gammat(3,3)

  ! Second derivatives
  CCTK_REAL                d2_beta1(3,3), d2_beta2(3,3), d2_beta3(3,3)
  CCTK_REAL                d2_hh11(3,3), d2_hh12(3,3), d2_hh13(3,3), d2_hh22(3,3), d2_hh23(3,3), d2_hh33(3,3)
  CCTK_REAL                d2_alph(3,3), d2_beta(3,3,3)
  CCTK_REAL                d2_ww(3,3), d2_hh(3,3,3,3)

  ! Advection derivatives
  CCTK_REAL                ad1_alph, ad1_beta(3)
  CCTK_REAL                ad1_ww, ad1_hh(3,3), ad1_trk, ad1_aa(3,3),      &
                           ad1_gammat(3)
  CCTK_REAL                d1_f(3)   ! Place holder for the advection derivs
  CCTK_REAL                ad1_lphi, ad1_lKphi

  ! Covaraint derivatives
  CCTK_REAL                cd2_ww(3,3), cd2_alph(3,3)

  ! Auxiliary variables
  CCTK_REAL                cf1(3,3,3), cf2(3,3,3), c_ri(3,3), c_ll(3,3)
  CCTK_REAL                c_ri_ww(3,3), c_ri_hh(3,3), ri_1(3,3), ri_2(3,3), &
                           ri_3(3,3), tr_ll, sq_aa, a2(3,3), trr,            &
                           tf_c_ll(3,3), tf_c_ri(3,3), gamcon(3)
  CCTK_REAL                tr_cd2_ww, tr_dww_dww, aux
  CCTK_REAL                tr_cd2_phi, tr_dww_dphi, tr_dalp_dphi, tr_dphi_dphi, tr_cd2_phi_new
  CCTK_REAL                divbeta, aadbeta(3,3), hhdbeta(3,3)

  ! Matter variables
  CCTK_REAL                srcE, srcjdi(3), srcji(3), srcSij(3,3),    &
                           srcSijTF(3,3), srcS_ww2, src_trT

  ! Right hand sides
  CCTK_REAL                rhs_ww, rhs_hh(3,3), rhs_trk, rhs_aa(3,3),      &
                           rhs_gammat(3), rhs_beta(3), rhs_lphi, rhs_lKphi

  ! Misc variables
  CCTK_REAL                dx12, dy12, dz12, dxsq12, dysq12, dzsq12,         &
                           dxdy144, dxdz144, dydz144
  CCTK_REAL                odx60, ody60, odz60, odxsq180, odysq180, odzsq180,&
                           odxdy3600, odxdz3600, odydz3600
  CCTK_INT                 i, j, k
  CCTK_INT                 di, dj, dk
  CCTK_REAL, parameter ::  one  = 1
  CCTK_REAL, parameter ::  zero = 0
  CCTK_REAL, parameter ::  pi   = acos(-one)
  CCTK_REAL, parameter ::  pi4  = 4*pi
  CCTK_REAL, parameter ::  pi8  = 8*pi
  CCTK_REAL, parameter ::  pi16 = 16*pi
  CCTK_INT                 a, b, c, l, m, n, p, q
  CCTK_REAL                myeta, r2, f_shift

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

  logical                   evolve_alp
  logical                   evolve_beta
  logical                   evolve_scalar
  logical                   evolve_Jordan
  logical                   JordanFrame

  evolve_alp    = CCTK_EQUALS(lapse_evolution_method, "LeanBSSNMoL")
  evolve_beta   = CCTK_EQUALS(shift_evolution_method, "LeanBSSNMoL")
  evolve_scalar = CCTK_EQUALS(scalar_evolution_method, "LeanBSSNMoL")
  evolve_Jordan = CCTK_EQUALS(theory,"DEFold") .OR. CCTK_EQUALS(theory,"BDold")
  JordanFrame   = CCTK_EQUALS(theory,"DEF") .OR. CCTK_EQUALS(theory,"BD") .OR. evolve_Jordan 

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

  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE( i, j, k, di, dj, dk, &
  !$OMP ww, hh, trk, aa, gammat, alph, beta, Tab, dethh, hu, beta_l, &
  !$OMP lphi, lKphi, Bphi, BKphi, omega_Bphi, domega_Bphi, F_phi, lphi2,&
  !$OMP trk_omega, trk_phi, aa_omega, aa_phi, gammat_omega, gammat_phi,&
  !$OMP d1_lphi, d1_lKphi, d2_lphi, cd2_lphi, cW_dphi, tr_cd2_phi_new, &
  !$OMP d1_beta1, d1_beta2, d1_beta3, &
  !$OMP d1_hh11, d1_hh12, d1_hh13, d1_hh22, d1_hh23, d1_hh33, &
  !$OMP d1_aa11, d1_aa12, d1_aa13, d1_aa22, d1_aa23, d1_aa33, &
  !$OMP d1_gammat1, d1_gammat2, d1_gammat3, &
  !$OMP d1_ww, d1_hh, d1_trk, d1_aa, d1_gammat, d1_alph, d1_beta, &
  !$OMP d2_beta1, d2_beta2, d2_beta3, &
  !$OMP d2_hh11, d2_hh12, d2_hh13, d2_hh22, d2_hh23, d2_hh33, &
  !$OMP d2_ww, d2_hh, d2_alph, d2_beta, &
  !$OMP d1_f, cf1, cf2, &
  !$OMP ad1_ww, ad1_hh, ad1_trk, ad1_aa, ad1_gammat, ad1_alph, ad1_beta, ad1_lphi, ad1_lKphi,&
  !$OMP cd2_ww, cd2_alph, ri_1, ri_2, ri_3, c_ri, c_ri_ww, c_ri_hh, &
  !$OMP tr_cd2_ww, tr_dww_dww, tr_dalp_dphi, tr_cd2_phi, tr_dww_dphi, tr_dphi_dphi, &
  !$OMP c_ll, aux, f_shift, divbeta, gamcon, &
  !$OMP rhs_ww, rhs_hh, rhs_trk, rhs_aa, rhs_gammat, rhs_beta, rhs_lphi, rhs_lKphi, &
  !$OMP hhdbeta, aadbeta, tr_ll, sq_aa, a2, trr, tf_c_ll, tf_c_ri, au,  &
  !$OMP myeta, r2, &
  !$OMP srcE, srcjdi, srcji, srcSij, srcS_ww2, srcSijTF, src_trT, &
  !$OMP a, b, c, l, m, n, p, q, jac, hes)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

    !------------ Get local variables ----------
    ww        = conf_fac(i,j,k)
   
    if (evolve_scalar) then
    lphi      = phi1(i,j,k)
    lphi2     = lphi*lphi
    lKphi     = Kphi1(i,j,k)
    end if 

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
    !------------------Beta in local coordinates------------------------
    beta_l(1) = jac(1,1)*beta(1) + jac(1,2)*beta(2) + jac(1,3)*beta(3)

    beta_l(2) = jac(2,1)*beta(1) + jac(2,2)*beta(2) + jac(2,3)*beta(3)

    beta_l(3) = jac(3,1)*beta(1) + jac(3,2)*beta(2) + jac(3,3)*beta(3)

    if (derivs_order == 6) then

      !------------ Centered 1st derivatives ----- 6th order
      ! d1_ww(3)
      d1_ww(1) = (  conf_fac(i+3,j,k) - 9*conf_fac(i+2,j,k) + 45*conf_fac(i+1,j,k)          &
                  - conf_fac(i-3,j,k) + 9*conf_fac(i-2,j,k) - 45*conf_fac(i-1,j,k) ) * odx60



      d1_ww(2) = (  conf_fac(i,j+3,k) - 9*conf_fac(i,j+2,k) + 45*conf_fac(i,j+1,k)          &
                  - conf_fac(i,j-3,k) + 9*conf_fac(i,j-2,k) - 45*conf_fac(i,j-1,k) ) * ody60

      d1_ww(3) = (  conf_fac(i,j,k+3) - 9*conf_fac(i,j,k+2) + 45*conf_fac(i,j,k+1)          &
                  - conf_fac(i,j,k-3) + 9*conf_fac(i,j,k-2) - 45*conf_fac(i,j,k-1) ) * odz60

      if (evolve_scalar) then
      !d1_phi(3)
      d1_lphi(1) = ( phi1(i+3,j,k) - 9*phi1(i+2,j,k) + 45*phi1(i+1,j,k)  &
                    -phi1(i-3,j,k) + 9*phi1(i-2,j,k) - 45*phi1(i-1,j,k)  ) * odx60

      d1_lphi(2) = ( phi1(i,j+3,k) - 9*phi1(i,j+2,k) + 45*phi1(i,j+1,k)  &
                    -phi1(i,j-3,k) + 9*phi1(i,j-2,k) - 45*phi1(i,j+1,k)  ) * odx60
 
      d1_lphi(3) = ( phi1(i,j,k+3) - 9*phi1(i,j,k+2) + 45*phi1(i,j,k+1)  &
                    -phi1(i,j,k-3) + 9*phi1(i,j,k-2) - 45*phi1(i,j,k-1)  ) * odx60

      !d1_Kphi(3)
      d1_lKphi(1) = ( Kphi1(i+3,j,k) - 9*Kphi1(i+2,j,k) + 45*Kphi1(i+1,j,k)  &
                     -Kphi1(i-3,j,k) + 9*Kphi1(i-2,j,k) - 45*Kphi1(i-1,j,k)  ) * odx60

      d1_lKphi(2) = ( Kphi1(i,j+3,k) - 9*Kphi1(i,j+2,k) + 45*Kphi1(i,j+1,k)  &
                     -Kphi1(i,j-3,k) + 9*Kphi1(i,j-2,k) - 45*Kphi1(i,j+1,k)  ) * odx60
 
      d1_lKphi(3) = ( Kphi1(i,j,k+3) - 9*Kphi1(i,j,k+2) + 45*Kphi1(i,j,k+1)  &
                     -Kphi1(i,j,k-3) + 9*Kphi1(i,j,k-2) - 45*Kphi1(i,j,k-1)  ) * odx60
      end if

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


      ! d1_alph(3)
      d1_alph(1) = (  alp(i+3,j,k) - 9*alp(i+2,j,k) + 45*alp(i+1,j,k) &
                    - alp(i-3,j,k) + 9*alp(i-2,j,k) - 45*alp(i-1,j,k) ) * odx60

      d1_alph(2) = (  alp(i,j+3,k) - 9*alp(i,j+2,k) + 45*alp(i,j+1,k) &
                    - alp(i,j-3,k) + 9*alp(i,j-2,k) - 45*alp(i,j-1,k) ) * ody60

      d1_alph(3) = (  alp(i,j,k+3) - 9*alp(i,j,k+2) + 45*alp(i,j,k+1) &
                    - alp(i,j,k-3) + 9*alp(i,j,k-2) - 45*alp(i,j,k-1) ) * odz60


      ! d1_beta(3,3)
      d1_beta1(1) = (  betax(i+3,j,k) - 9*betax(i+2,j,k) + 45*betax(i+1,j,k) &
                     - betax(i-3,j,k) + 9*betax(i-2,j,k) - 45*betax(i-1,j,k) ) * odx60
      d1_beta2(1) = (  betay(i+3,j,k) - 9*betay(i+2,j,k) + 45*betay(i+1,j,k) &
                     - betay(i-3,j,k) + 9*betay(i-2,j,k) - 45*betay(i-1,j,k) ) * odx60
      d1_beta3(1) = (  betaz(i+3,j,k) - 9*betaz(i+2,j,k) + 45*betaz(i+1,j,k) &
                     - betaz(i-3,j,k) + 9*betaz(i-2,j,k) - 45*betaz(i-1,j,k) ) * odx60

      d1_beta1(2) = (  betax(i,j+3,k) - 9*betax(i,j+2,k) + 45*betax(i,j+1,k) &
                     - betax(i,j-3,k) + 9*betax(i,j-2,k) - 45*betax(i,j-1,k) ) * ody60
      d1_beta2(2) = (  betay(i,j+3,k) - 9*betay(i,j+2,k) + 45*betay(i,j+1,k) &
                     - betay(i,j-3,k) + 9*betay(i,j-2,k) - 45*betay(i,j-1,k) ) * ody60
      d1_beta3(2) = (  betaz(i,j+3,k) - 9*betaz(i,j+2,k) + 45*betaz(i,j+1,k) &
                     - betaz(i,j-3,k) + 9*betaz(i,j-2,k) - 45*betaz(i,j-1,k) ) * ody60

      d1_beta1(3) = (  betax(i,j,k+3) - 9*betax(i,j,k+2) + 45*betax(i,j,k+1) &
                     - betax(i,j,k-3) + 9*betax(i,j,k-2) - 45*betax(i,j,k-1) ) * odz60
      d1_beta2(3) = (  betay(i,j,k+3) - 9*betay(i,j,k+2) + 45*betay(i,j,k+1) &
                     - betay(i,j,k-3) + 9*betay(i,j,k-2) - 45*betay(i,j,k-1) ) * odz60
      d1_beta3(3) = (  betaz(i,j,k+3) - 9*betaz(i,j,k+2) + 45*betaz(i,j,k+1) &
                     - betaz(i,j,k-3) + 9*betaz(i,j,k-2) - 45*betaz(i,j,k-1) ) * odz60


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

      if (evolve_scalar) then
      !d2_lphi(3,3) 
      d2_lphi(1,1) = ( 2*phi1(i+3,j,k) - 27*phi1(i+2,j,k) + 270*phi1(i+1,j,k) - 490*phi1(i,j,k) &
                      +2*phi1(i-3,j,k) - 27*phi1(i-2,j,k) + 270*phi1(i-1,j,k) ) / odxsq180
      
      d2_lphi(2,2) = ( 2*phi1(i,j+3,k) - 27*phi1(i,j+2,k) + 270*phi1(i,j+1,k) - 490*phi1(i,j,k) &
                      +2*phi1(i,j-3,k) - 27*phi1(i,j-2,k) + 270*phi1(i,j-1,k) ) / odysq180
      
      d2_lphi(3,3) = ( 2*phi1(i,j,k+3) - 27*phi1(i,j,k+2) + 270*phi1(i,j,k+1) - 490*phi1(i,j,k) &
                      +2*phi1(i,j,k-3) - 27*phi1(i,j,k-2) + 270*phi1(i,j,k-1) ) / odzsq180

      d2_lphi(1,2) = (  -  phi1(i-3,j+3,k) +   9*phi1(i-2,j+3,k) -   45*phi1(i-1,j+3,k) +   45*phi1(i+1,j+3,k) -   9*phi1(i+2,j+3,k) +    phi1(i+3,j+3,k) &
                      +  9*phi1(i-3,j+2,k) -  81*phi1(i-2,j+2,k) +  405*phi1(i-1,j+2,k) -  405*phi1(i+1,j+2,k) +  81*phi1(i+2,j+2,k) -  9*phi1(i+3,j+2,k) &
                      - 45*phi1(i-3,j+1,k) + 405*phi1(i-2,j+1,k) - 2025*phi1(i-1,j+1,k) + 2025*phi1(i+1,j+1,k) - 405*phi1(i+2,j+1,k) + 45*phi1(i+3,j+1,k) &
                      + 45*phi1(i-3,j-1,k) - 405*phi1(i-2,j-1,k) + 2025*phi1(i-1,j-1,k) - 2025*phi1(i+1,j-1,k) + 405*phi1(i+2,j-1,k) - 45*phi1(i+3,j-1,k) &
                      -  9*phi1(i-3,j-2,k) +  81*phi1(i-2,j-2,k) -  405*phi1(i-1,j-2,k) +  405*phi1(i+1,j-2,k) -  81*phi1(i+2,j-2,k) +  9*phi1(i+3,j-2,k) &
                      +    phi1(i-3,j-3,k) -   9*phi1(i-2,j-3,k) +   45*phi1(i-1,j-3,k) -   45*phi1(i+1,j-3,k) +   9*phi1(i+2,j-3,k) -    phi1(i+3,j-3,k) ) * odxdy3600

      d2_lphi(1,3) = (  -  phi1(i-3,j,k+3) +   9*phi1(i-2,j,k+3) -   45*phi1(i-1,j,k+3) +   45*phi1(i+1,j,k+3) -   9*phi1(i+2,j,k+3) +    phi1(i+3,j,k+3) &
                      +  9*phi1(i-3,j,k+2) -  81*phi1(i-2,j,k+2) +  405*phi1(i-1,j,k+2) -  405*phi1(i+1,j,k+2) +  81*phi1(i+2,j,k+2) -  9*phi1(i+3,j,k+2) &
                      - 45*phi1(i-3,j,k+1) + 405*phi1(i-2,j,k+1) - 2025*phi1(i-1,j,k+1) + 2025*phi1(i+1,j,k+1) - 405*phi1(i+2,j,k+1) + 45*phi1(i+3,j,k+1) &
                      + 45*phi1(i-3,j,k-1) - 405*phi1(i-2,j,k-1) + 2025*phi1(i-1,j,k-1) - 2025*phi1(i+1,j,k-1) + 405*phi1(i+2,j,k-1) - 45*phi1(i+3,j,k-1) &
                      -  9*phi1(i-3,j,k-2) +  81*phi1(i-2,j,k-2) -  405*phi1(i-1,j,k-2) +  405*phi1(i+1,j,k-2) -  81*phi1(i+2,j,k-2) +  9*phi1(i+3,j,k-2) &
                      +    phi1(i-3,j,k-3) -   9*phi1(i-2,j,k-3) +   45*phi1(i-1,j,k-3) -   45*phi1(i+1,j,k-3) +   9*phi1(i+2,j,k-3) -    phi1(i+3,j,k-3) ) * odxdz3600

      d2_lphi(2,3) = (  -  phi1(i,j-3,k+3) +   9*phi1(i,j-2,k+3) -   45*phi1(i,j-1,k+3) +   45*phi1(i,j+1,k+3) -   9*phi1(i,j+2,k+3) +    phi1(i,j+3,k+3) &
                      +  9*phi1(i,j-3,k+2) -  81*phi1(i,j-2,k+2) +  405*phi1(i,j-1,k+2) -  405*phi1(i,j+1,k+2) +  81*phi1(i,j+2,k+2) -  9*phi1(i,j+3,k+2) &
                      - 45*phi1(i,j-3,k+1) + 405*phi1(i,j-2,k+1) - 2025*phi1(i,j-1,k+1) + 2025*phi1(i,j+1,k+1) - 405*phi1(i,j+2,k+1) + 45*phi1(i,j+3,k+1) &
                      + 45*phi1(i,j-3,k-1) - 405*phi1(i,j-2,k-1) + 2025*phi1(i,j-1,k-1) - 2025*phi1(i,j+1,k-1) + 405*phi1(i,j+2,k-1) - 45*phi1(i,j+3,k-1) &
                      -  9*phi1(i,j-3,k-2) +  81*phi1(i,j-2,k-2) -  405*phi1(i,j-1,k-2) +  405*phi1(i,j+1,k-2) -  81*phi1(i,j+2,k-2) +  9*phi1(i,j+3,k-2) &
                      +    phi1(i,j-3,k-3) -   9*phi1(i,j-2,k-3) +   45*phi1(i,j-1,k-3) -   45*phi1(i,j+1,k-3) +   9*phi1(i,j+2,k-3) -    phi1(i,j+3,k-3) ) * odydz3600

      d2_lphi(2,1) = d2_lphi(1,2)
      d2_lphi(3,1) = d2_lphi(1,3)
      d2_lphi(3,2) = d2_lphi(2,3)

      end if 

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

      ! d2_alph(3,3)
      d2_alph(1,1) = (  2*alp(i+3,j,k) - 27*alp(i+2,j,k) + 270*alp(i+1,j,k) - 490*alp(i,j,k)&
                      + 2*alp(i-3,j,k) - 27*alp(i-2,j,k) + 270*alp(i-1,j,k) ) * odxsq180

      d2_alph(2,2) = (  2*alp(i,j+3,k) - 27*alp(i,j+2,k) + 270*alp(i,j+1,k) - 490*alp(i,j,k)&
                      + 2*alp(i,j-3,k) - 27*alp(i,j-2,k) + 270*alp(i,j-1,k) ) * odysq180

      d2_alph(3,3) = (  2*alp(i,j,k+3) - 27*alp(i,j,k+2) + 270*alp(i,j,k+1) - 490*alp(i,j,k)&
                      + 2*alp(i,j,k-3) - 27*alp(i,j,k-2) + 270*alp(i,j,k-1) ) * odzsq180

      d2_alph(1,2) = (    -alp(i-3,j+3,k) +   9*alp(i-2,j+3,k) -   45*alp(i-1,j+3,k) +   45*alp(i+1,j+3,k) -   9*alp(i+2,j+3,k) +    alp(i+3,j+3,k) &
                      +  9*alp(i-3,j+2,k) -  81*alp(i-2,j+2,k) +  405*alp(i-1,j+2,k) -  405*alp(i+1,j+2,k) +  81*alp(i+2,j+2,k) -  9*alp(i+3,j+2,k) &
                      - 45*alp(i-3,j+1,k) + 405*alp(i-2,j+1,k) - 2025*alp(i-1,j+1,k) + 2025*alp(i+1,j+1,k) - 405*alp(i+2,j+1,k) + 45*alp(i+3,j+1,k) &
                      + 45*alp(i-3,j-1,k) - 405*alp(i-2,j-1,k) + 2025*alp(i-1,j-1,k) - 2025*alp(i+1,j-1,k) + 405*alp(i+2,j-1,k) - 45*alp(i+3,j-1,k) &
                      -  9*alp(i-3,j-2,k) +  81*alp(i-2,j-2,k) -  405*alp(i-1,j-2,k) +  405*alp(i+1,j-2,k) -  81*alp(i+2,j-2,k) +  9*alp(i+3,j-2,k) &
                      +    alp(i-3,j-3,k) -   9*alp(i-2,j-3,k) +   45*alp(i-1,j-3,k) -   45*alp(i+1,j-3,k) +   9*alp(i+2,j-3,k) -    alp(i+3,j-3,k) ) * odxdy3600

      d2_alph(1,3) = (    -alp(i-3,j,k+3) +   9*alp(i-2,j,k+3) -   45*alp(i-1,j,k+3) +   45*alp(i+1,j,k+3) -   9*alp(i+2,j,k+3) +    alp(i+3,j,k+3) &
                      +  9*alp(i-3,j,k+2) -  81*alp(i-2,j,k+2) +  405*alp(i-1,j,k+2) -  405*alp(i+1,j,k+2) +  81*alp(i+2,j,k+2) -  9*alp(i+3,j,k+2) &
                      - 45*alp(i-3,j,k+1) + 405*alp(i-2,j,k+1) - 2025*alp(i-1,j,k+1) + 2025*alp(i+1,j,k+1) - 405*alp(i+2,j,k+1) + 45*alp(i+3,j,k+1) &
                      + 45*alp(i-3,j,k-1) - 405*alp(i-2,j,k-1) + 2025*alp(i-1,j,k-1) - 2025*alp(i+1,j,k-1) + 405*alp(i+2,j,k-1) - 45*alp(i+3,j,k-1) &
                      -  9*alp(i-3,j,k-2) +  81*alp(i-2,j,k-2) -  405*alp(i-1,j,k-2) +  405*alp(i+1,j,k-2) -  81*alp(i+2,j,k-2) +  9*alp(i+3,j,k-2) &
                      +    alp(i-3,j,k-3) -   9*alp(i-2,j,k-3) +   45*alp(i-1,j,k-3) -   45*alp(i+1,j,k-3) +   9*alp(i+2,j,k-3) -    alp(i+3,j,k-3) ) * odxdz3600

      d2_alph(2,3) = (    -alp(i,j-3,k+3) +   9*alp(i,j-2,k+3) -   45*alp(i,j-1,k+3) +   45*alp(i,j+1,k+3) -   9*alp(i,j+2,k+3) +    alp(i,j+3,k+3) &
                      +  9*alp(i,j-3,k+2) -  81*alp(i,j-2,k+2) +  405*alp(i,j-1,k+2) -  405*alp(i,j+1,k+2) +  81*alp(i,j+2,k+2) -  9*alp(i,j+3,k+2) &
                      - 45*alp(i,j-3,k+1) + 405*alp(i,j-2,k+1) - 2025*alp(i,j-1,k+1) + 2025*alp(i,j+1,k+1) - 405*alp(i,j+2,k+1) + 45*alp(i,j+3,k+1) &
                      + 45*alp(i,j-3,k-1) - 405*alp(i,j-2,k-1) + 2025*alp(i,j-1,k-1) - 2025*alp(i,j+1,k-1) + 405*alp(i,j+2,k-1) - 45*alp(i,j+3,k-1) &
                      -  9*alp(i,j-3,k-2) +  81*alp(i,j-2,k-2) -  405*alp(i,j-1,k-2) +  405*alp(i,j+1,k-2) -  81*alp(i,j+2,k-2) +  9*alp(i,j+3,k-2) &
                      +    alp(i,j-3,k-3) -   9*alp(i,j-2,k-3) +   45*alp(i,j-1,k-3) -   45*alp(i,j+1,k-3) +   9*alp(i,j+2,k-3) -    alp(i,j+3,k-3) ) * odydz3600

      d2_alph(2,1) = d2_alph(1,2)
      d2_alph(3,1) = d2_alph(1,3)
      d2_alph(3,2) = d2_alph(2,3)


      ! d2_beta(3,3,3)
      d2_beta1(1,1) = (  2*betax(i+3,j,k) - 27*betax(i+2,j,k) + 270*betax(i+1,j,k) - 490*betax(i,j,k)&
                       + 2*betax(i-3,j,k) - 27*betax(i-2,j,k) + 270*betax(i-1,j,k) ) * odxsq180
      d2_beta2(1,1) = (  2*betay(i+3,j,k) - 27*betay(i+2,j,k) + 270*betay(i+1,j,k) - 490*betay(i,j,k)&
                       + 2*betay(i-3,j,k) - 27*betay(i-2,j,k) + 270*betay(i-1,j,k) ) * odxsq180
      d2_beta3(1,1) = (  2*betaz(i+3,j,k) - 27*betaz(i+2,j,k) + 270*betaz(i+1,j,k) - 490*betaz(i,j,k)&
                       + 2*betaz(i-3,j,k) - 27*betaz(i-2,j,k) + 270*betaz(i-1,j,k) ) * odxsq180

      d2_beta1(2,2) = (  2*betax(i,j+3,k) - 27*betax(i,j+2,k) + 270*betax(i,j+1,k) - 490*betax(i,j,k)&
                       + 2*betax(i,j-3,k) - 27*betax(i,j-2,k) + 270*betax(i,j-1,k) ) * odysq180
      d2_beta2(2,2) = (  2*betay(i,j+3,k) - 27*betay(i,j+2,k) + 270*betay(i,j+1,k) - 490*betay(i,j,k)&
                       + 2*betay(i,j-3,k) - 27*betay(i,j-2,k) + 270*betay(i,j-1,k) ) * odysq180
      d2_beta3(2,2) = (  2*betaz(i,j+3,k) - 27*betaz(i,j+2,k) + 270*betaz(i,j+1,k) - 490*betaz(i,j,k)&
                       + 2*betaz(i,j-3,k) - 27*betaz(i,j-2,k) + 270*betaz(i,j-1,k) ) * odysq180

      d2_beta1(3,3) = (  2*betax(i,j,k+3) - 27*betax(i,j,k+2) + 270*betax(i,j,k+1) - 490*betax(i,j,k)&
                       + 2*betax(i,j,k-3) - 27*betax(i,j,k-2) + 270*betax(i,j,k-1) ) * odzsq180
      d2_beta2(3,3) = (  2*betay(i,j,k+3) - 27*betay(i,j,k+2) + 270*betay(i,j,k+1) - 490*betay(i,j,k)&
                       + 2*betay(i,j,k-3) - 27*betay(i,j,k-2) + 270*betay(i,j,k-1) ) * odzsq180
      d2_beta3(3,3) = (  2*betaz(i,j,k+3) - 27*betaz(i,j,k+2) + 270*betaz(i,j,k+1) - 490*betaz(i,j,k)&
                       + 2*betaz(i,j,k-3) - 27*betaz(i,j,k-2) + 270*betaz(i,j,k-1) ) * odzsq180

      d2_beta1(1,2) = (    -betax(i-3,j+3,k) +   9*betax(i-2,j+3,k) -   45*betax(i-1,j+3,k) +   45*betax(i+1,j+3,k) -   9*betax(i+2,j+3,k) +    betax(i+3,j+3,k) &
                       +  9*betax(i-3,j+2,k) -  81*betax(i-2,j+2,k) +  405*betax(i-1,j+2,k) -  405*betax(i+1,j+2,k) +  81*betax(i+2,j+2,k) -  9*betax(i+3,j+2,k) &
                       - 45*betax(i-3,j+1,k) + 405*betax(i-2,j+1,k) - 2025*betax(i-1,j+1,k) + 2025*betax(i+1,j+1,k) - 405*betax(i+2,j+1,k) + 45*betax(i+3,j+1,k) &
                       + 45*betax(i-3,j-1,k) - 405*betax(i-2,j-1,k) + 2025*betax(i-1,j-1,k) - 2025*betax(i+1,j-1,k) + 405*betax(i+2,j-1,k) - 45*betax(i+3,j-1,k) &
                       -  9*betax(i-3,j-2,k) +  81*betax(i-2,j-2,k) -  405*betax(i-1,j-2,k) +  405*betax(i+1,j-2,k) -  81*betax(i+2,j-2,k) +  9*betax(i+3,j-2,k) &
                       +    betax(i-3,j-3,k) -   9*betax(i-2,j-3,k) +   45*betax(i-1,j-3,k) -   45*betax(i+1,j-3,k) +   9*betax(i+2,j-3,k) -    betax(i+3,j-3,k) ) * odxdy3600
      d2_beta2(1,2) = (    -betay(i-3,j+3,k) +   9*betay(i-2,j+3,k) -   45*betay(i-1,j+3,k) +   45*betay(i+1,j+3,k) -   9*betay(i+2,j+3,k) +    betay(i+3,j+3,k) &
                       +  9*betay(i-3,j+2,k) -  81*betay(i-2,j+2,k) +  405*betay(i-1,j+2,k) -  405*betay(i+1,j+2,k) +  81*betay(i+2,j+2,k) -  9*betay(i+3,j+2,k) &
                       - 45*betay(i-3,j+1,k) + 405*betay(i-2,j+1,k) - 2025*betay(i-1,j+1,k) + 2025*betay(i+1,j+1,k) - 405*betay(i+2,j+1,k) + 45*betay(i+3,j+1,k) &
                       + 45*betay(i-3,j-1,k) - 405*betay(i-2,j-1,k) + 2025*betay(i-1,j-1,k) - 2025*betay(i+1,j-1,k) + 405*betay(i+2,j-1,k) - 45*betay(i+3,j-1,k) &
                       -  9*betay(i-3,j-2,k) +  81*betay(i-2,j-2,k) -  405*betay(i-1,j-2,k) +  405*betay(i+1,j-2,k) -  81*betay(i+2,j-2,k) +  9*betay(i+3,j-2,k) &
                         +    betay(i-3,j-3,k) -   9*betay(i-2,j-3,k) +   45*betay(i-1,j-3,k) -   45*betay(i+1,j-3,k) +   9*betay(i+2,j-3,k) -    betay(i+3,j-3,k) ) * odxdy3600
      d2_beta3(1,2) = (    -betaz(i-3,j+3,k) +   9*betaz(i-2,j+3,k) -   45*betaz(i-1,j+3,k) +   45*betaz(i+1,j+3,k) -   9*betaz(i+2,j+3,k) +    betaz(i+3,j+3,k) &
                       +  9*betaz(i-3,j+2,k) -  81*betaz(i-2,j+2,k) +  405*betaz(i-1,j+2,k) -  405*betaz(i+1,j+2,k) +  81*betaz(i+2,j+2,k) -  9*betaz(i+3,j+2,k) &
                       - 45*betaz(i-3,j+1,k) + 405*betaz(i-2,j+1,k) - 2025*betaz(i-1,j+1,k) + 2025*betaz(i+1,j+1,k) - 405*betaz(i+2,j+1,k) + 45*betaz(i+3,j+1,k) &
                       + 45*betaz(i-3,j-1,k) - 405*betaz(i-2,j-1,k) + 2025*betaz(i-1,j-1,k) - 2025*betaz(i+1,j-1,k) + 405*betaz(i+2,j-1,k) - 45*betaz(i+3,j-1,k) &
                       -  9*betaz(i-3,j-2,k) +  81*betaz(i-2,j-2,k) -  405*betaz(i-1,j-2,k) +  405*betaz(i+1,j-2,k) -  81*betaz(i+2,j-2,k) +  9*betaz(i+3,j-2,k) &
                       +    betaz(i-3,j-3,k) -   9*betaz(i-2,j-3,k) +   45*betaz(i-1,j-3,k) -   45*betaz(i+1,j-3,k) +   9*betaz(i+2,j-3,k) -    betaz(i+3,j-3,k) ) * odxdy3600


      d2_beta1(1,3) = (    -betax(i-3,j,k+3) +   9*betax(i-2,j,k+3) -   45*betax(i-1,j,k+3) +   45*betax(i+1,j,k+3) -   9*betax(i+2,j,k+3) +    betax(i+3,j,k+3) &
                       +  9*betax(i-3,j,k+2) -  81*betax(i-2,j,k+2) +  405*betax(i-1,j,k+2) -  405*betax(i+1,j,k+2) +  81*betax(i+2,j,k+2) -  9*betax(i+3,j,k+2) &
                       - 45*betax(i-3,j,k+1) + 405*betax(i-2,j,k+1) - 2025*betax(i-1,j,k+1) + 2025*betax(i+1,j,k+1) - 405*betax(i+2,j,k+1) + 45*betax(i+3,j,k+1) &
                       + 45*betax(i-3,j,k-1) - 405*betax(i-2,j,k-1) + 2025*betax(i-1,j,k-1) - 2025*betax(i+1,j,k-1) + 405*betax(i+2,j,k-1) - 45*betax(i+3,j,k-1) &
                       -  9*betax(i-3,j,k-2) +  81*betax(i-2,j,k-2) -  405*betax(i-1,j,k-2) +  405*betax(i+1,j,k-2) -  81*betax(i+2,j,k-2) +  9*betax(i+3,j,k-2) &
                       +    betax(i-3,j,k-3) -   9*betax(i-2,j,k-3) +   45*betax(i-1,j,k-3) -   45*betax(i+1,j,k-3) +   9*betax(i+2,j,k-3) -    betax(i+3,j,k-3) ) * odxdz3600
      d2_beta2(1,3) = (    -betay(i-3,j,k+3) +   9*betay(i-2,j,k+3) -   45*betay(i-1,j,k+3) +   45*betay(i+1,j,k+3) -   9*betay(i+2,j,k+3) +    betay(i+3,j,k+3) &
                       +  9*betay(i-3,j,k+2) -  81*betay(i-2,j,k+2) +  405*betay(i-1,j,k+2) -  405*betay(i+1,j,k+2) +  81*betay(i+2,j,k+2) -  9*betay(i+3,j,k+2) &
                       - 45*betay(i-3,j,k+1) + 405*betay(i-2,j,k+1) - 2025*betay(i-1,j,k+1) + 2025*betay(i+1,j,k+1) - 405*betay(i+2,j,k+1) + 45*betay(i+3,j,k+1) &
                       + 45*betay(i-3,j,k-1) - 405*betay(i-2,j,k-1) + 2025*betay(i-1,j,k-1) - 2025*betay(i+1,j,k-1) + 405*betay(i+2,j,k-1) - 45*betay(i+3,j,k-1) &
                       -  9*betay(i-3,j,k-2) +  81*betay(i-2,j,k-2) -  405*betay(i-1,j,k-2) +  405*betay(i+1,j,k-2) -  81*betay(i+2,j,k-2) +  9*betay(i+3,j,k-2) &
                       +    betay(i-3,j,k-3) -   9*betay(i-2,j,k-3) +   45*betay(i-1,j,k-3) -   45*betay(i+1,j,k-3) +   9*betay(i+2,j,k-3) -    betay(i+3,j,k-3) ) * odxdz3600
      d2_beta3(1,3) = (    -betaz(i-3,j,k+3) +   9*betaz(i-2,j,k+3) -   45*betaz(i-1,j,k+3) +   45*betaz(i+1,j,k+3) -   9*betaz(i+2,j,k+3) +    betaz(i+3,j,k+3) &
                       +  9*betaz(i-3,j,k+2) -  81*betaz(i-2,j,k+2) +  405*betaz(i-1,j,k+2) -  405*betaz(i+1,j,k+2) +  81*betaz(i+2,j,k+2) -  9*betaz(i+3,j,k+2) &
                       - 45*betaz(i-3,j,k+1) + 405*betaz(i-2,j,k+1) - 2025*betaz(i-1,j,k+1) + 2025*betaz(i+1,j,k+1) - 405*betaz(i+2,j,k+1) + 45*betaz(i+3,j,k+1) &
                       + 45*betaz(i-3,j,k-1) - 405*betaz(i-2,j,k-1) + 2025*betaz(i-1,j,k-1) - 2025*betaz(i+1,j,k-1) + 405*betaz(i+2,j,k-1) - 45*betaz(i+3,j,k-1) &
                       -  9*betaz(i-3,j,k-2) +  81*betaz(i-2,j,k-2) -  405*betaz(i-1,j,k-2) +  405*betaz(i+1,j,k-2) -  81*betaz(i+2,j,k-2) +  9*betaz(i+3,j,k-2) &
                       +    betaz(i-3,j,k-3) -   9*betaz(i-2,j,k-3) +   45*betaz(i-1,j,k-3) -   45*betaz(i+1,j,k-3) +   9*betaz(i+2,j,k-3) -    betaz(i+3,j,k-3) ) * odxdz3600

      d2_beta1(2,3) = (    -betax(i,j-3,k+3) +   9*betax(i,j-2,k+3) -   45*betax(i,j-1,k+3) +   45*betax(i,j+1,k+3) -   9*betax(i,j+2,k+3) +    betax(i,j+3,k+3) &
                       +  9*betax(i,j-3,k+2) -  81*betax(i,j-2,k+2) +  405*betax(i,j-1,k+2) -  405*betax(i,j+1,k+2) +  81*betax(i,j+2,k+2) -  9*betax(i,j+3,k+2) &
                       - 45*betax(i,j-3,k+1) + 405*betax(i,j-2,k+1) - 2025*betax(i,j-1,k+1) + 2025*betax(i,j+1,k+1) - 405*betax(i,j+2,k+1) + 45*betax(i,j+3,k+1) &
                       + 45*betax(i,j-3,k-1) - 405*betax(i,j-2,k-1) + 2025*betax(i,j-1,k-1) - 2025*betax(i,j+1,k-1) + 405*betax(i,j+2,k-1) - 45*betax(i,j+3,k-1) &
                       -  9*betax(i,j-3,k-2) +  81*betax(i,j-2,k-2) -  405*betax(i,j-1,k-2) +  405*betax(i,j+1,k-2) -  81*betax(i,j+2,k-2) +  9*betax(i,j+3,k-2) &
                       +    betax(i,j-3,k-3) -   9*betax(i,j-2,k-3) +   45*betax(i,j-1,k-3) -   45*betax(i,j+1,k-3) +   9*betax(i,j+2,k-3) -    betax(i,j+3,k-3) ) * odydz3600
      d2_beta2(2,3) = (    -betay(i,j-3,k+3) +   9*betay(i,j-2,k+3) -   45*betay(i,j-1,k+3) +   45*betay(i,j+1,k+3) -   9*betay(i,j+2,k+3) +    betay(i,j+3,k+3) &
                       +  9*betay(i,j-3,k+2) -  81*betay(i,j-2,k+2) +  405*betay(i,j-1,k+2) -  405*betay(i,j+1,k+2) +  81*betay(i,j+2,k+2) -  9*betay(i,j+3,k+2) &
                       - 45*betay(i,j-3,k+1) + 405*betay(i,j-2,k+1) - 2025*betay(i,j-1,k+1) + 2025*betay(i,j+1,k+1) - 405*betay(i,j+2,k+1) + 45*betay(i,j+3,k+1) &
                       + 45*betay(i,j-3,k-1) - 405*betay(i,j-2,k-1) + 2025*betay(i,j-1,k-1) - 2025*betay(i,j+1,k-1) + 405*betay(i,j+2,k-1) - 45*betay(i,j+3,k-1) &
                       -  9*betay(i,j-3,k-2) +  81*betay(i,j-2,k-2) -  405*betay(i,j-1,k-2) +  405*betay(i,j+1,k-2) -  81*betay(i,j+2,k-2) +  9*betay(i,j+3,k-2) &
                       +    betay(i,j-3,k-3) -   9*betay(i,j-2,k-3) +   45*betay(i,j-1,k-3) -   45*betay(i,j+1,k-3) +   9*betay(i,j+2,k-3) -    betay(i,j+3,k-3) ) * odydz3600
      d2_beta3(2,3) = (    -betaz(i,j-3,k+3) +   9*betaz(i,j-2,k+3) -   45*betaz(i,j-1,k+3) +   45*betaz(i,j+1,k+3) -   9*betaz(i,j+2,k+3) +    betaz(i,j+3,k+3) &
                       +  9*betaz(i,j-3,k+2) -  81*betaz(i,j-2,k+2) +  405*betaz(i,j-1,k+2) -  405*betaz(i,j+1,k+2) +  81*betaz(i,j+2,k+2) -  9*betaz(i,j+3,k+2) &
                       - 45*betaz(i,j-3,k+1) + 405*betaz(i,j-2,k+1) - 2025*betaz(i,j-1,k+1) + 2025*betaz(i,j+1,k+1) - 405*betaz(i,j+2,k+1) + 45*betaz(i,j+3,k+1) &
                       + 45*betaz(i,j-3,k-1) - 405*betaz(i,j-2,k-1) + 2025*betaz(i,j-1,k-1) - 2025*betaz(i,j+1,k-1) + 405*betaz(i,j+2,k-1) - 45*betaz(i,j+3,k-1) &
                       -  9*betaz(i,j-3,k-2) +  81*betaz(i,j-2,k-2) -  405*betaz(i,j-1,k-2) +  405*betaz(i,j+1,k-2) -  81*betaz(i,j+2,k-2) +  9*betaz(i,j+3,k-2) &
                       +    betaz(i,j-3,k-3) -   9*betaz(i,j-2,k-3) +   45*betaz(i,j-1,k-3) -   45*betaz(i,j+1,k-3) +   9*betaz(i,j+2,k-3) -    betaz(i,j+3,k-3) ) * odydz3600

      d2_beta1(2,1) = d2_beta1(1,2)
      d2_beta2(2,1) = d2_beta2(1,2)
      d2_beta3(2,1) = d2_beta3(1,2)
      d2_beta1(3,1) = d2_beta1(1,3)
      d2_beta2(3,1) = d2_beta2(1,3)
      d2_beta3(3,1) = d2_beta3(1,3)
      d2_beta1(3,2) = d2_beta1(2,3)
      d2_beta2(3,2) = d2_beta2(2,3)
      d2_beta3(3,2) = d2_beta3(2,3)


      !------------ Advection derivatives --------
      if( use_advection_stencils /= 0 ) then

        di = int( sign( one, beta_l(1) ) )
        dj = int( sign( one, beta_l(2) ) )
        dk = int( sign( one, beta_l(3) ) )

        ! ad1_ww
        d1_f(1) = di * (   2*conf_fac(i-2*di,j,k) - 24*conf_fac(i-di,j,k) - 35*conf_fac(i,j,k) + 80*conf_fac(i+di,j,k) &
                        - 30*conf_fac(i+2*di,j,k) +  8*conf_fac(i+3*di,j,k) -  conf_fac(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*conf_fac(i,j-2*dj,k) - 24*conf_fac(i,j-dj,k) - 35*conf_fac(i,j,k) + 80*conf_fac(i,j+dj,k) &
                        - 30*conf_fac(i,j+2*dj,k) +  8*conf_fac(i,j+3*dj,k) -  conf_fac(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*conf_fac(i,j,k-2*dk) - 24*conf_fac(i,j,k-dk) - 35*conf_fac(i,j,k) + 80*conf_fac(i,j,k+dk) &
                        - 30*conf_fac(i,j,k+2*dk) +  8*conf_fac(i,j,k+3*dk) -  conf_fac(i,j,k+4*dk) ) * odz60
        ad1_ww  = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        if (evolve_scalar) then
        ! ad1_phi 
        d1_f(1) = di * ( 2*phi1(i-2*di,j,k) - 24*phi1(i-di,j,k) - 35*phi1(i,j,k) + 80*phi1(i+di,j,k) &
                        - 30*phi1(i+2*di,j,k) + 8*phi1(i+3*di,j,k) - phi1(i+4*di,j,k) ) * odx60

        d1_f(2) = dj * ( 2*phi1(i,j-2*dj,k) - 24*phi1(i,j-dj,k) - 35*phi1(i,j,k) + 80*phi1(i,j+dj,k) &
                        - 30*phi1(i,j+2*dj,k) + 8*phi1(i,j+3*dj,k) - phi1(i,j+4*dj,k) ) * ody60

        d1_f(3) = dk * ( 2*phi1(i,j,k-2*dk) - 24*phi1(i,j,k-dk) - 35*phi1(i,j,k) + 80*phi1(i,j,k+dk) &
                        - 30*phi1(i,j,k+2*dk) + 8*phi1(i,j,k+3*dk) - phi1(i,j,k+4*dk) ) * odz60

        ad1_lphi = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        ! ad1_Kphi 
        d1_f(1) = di * ( 2*Kphi1(i-2*di,j,k) - 24*Kphi1(i-di,j,k) - 35*Kphi1(i,j,k) + 80*Kphi1(i+di,j,k) &
                        - 30*Kphi1(i+2*di,j,k) + 8*Kphi1(i+3*di,j,k) - Kphi1(i+4*di,j,k) ) * odx60

        d1_f(2) = dj * ( 2*Kphi1(i,j-2*dj,k) - 24*Kphi1(i,j-dj,k) - 35*Kphi1(i,j,k) + 80*Kphi1(i,j+dj,k) &
                        - 30*Kphi1(i,j+2*dj,k) + 8*Kphi1(i,j+3*dj,k) - Kphi1(i,j+4*dj,k) ) * ody60

        d1_f(3) = dk * ( 2*Kphi1(i,j,k-2*dk) - 24*Kphi1(i,j,k-dk) - 35*Kphi1(i,j,k) + 80*Kphi1(i,j,k+dk) &
                        - 30*Kphi1(i,j,k+2*dk) + 8*Kphi1(i,j,k+3*dk) - Kphi1(i,j,k+4*dk) ) * odz60

        ad1_lKphi = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)
        end if

        ! ad1_hh(3,3)
        d1_f(1) = di * (   2*hxx(i-2*di,j,k) - 24*hxx(i-di,j,k) - 35*hxx(i,j,k) + 80*hxx(i+di,j,k) &
                        - 30*hxx(i+2*di,j,k) +  8*hxx(i+3*di,j,k) -  hxx(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*hxx(i,j-2*dj,k) - 24*hxx(i,j-dj,k) - 35*hxx(i,j,k) + 80*hxx(i,j+dj,k) &
                        - 30*hxx(i,j+2*dj,k) +  8*hxx(i,j+3*dj,k) -  hxx(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*hxx(i,j,k-2*dk) - 24*hxx(i,j,k-dk) - 35*hxx(i,j,k) + 80*hxx(i,j,k+dk) &
                        - 30*hxx(i,j,k+2*dk) +  8*hxx(i,j,k+3*dk) -  hxx(i,j,k+4*dk) ) * odz60
        ad1_hh(1,1) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * (   2*hxy(i-2*di,j,k) - 24*hxy(i-di,j,k) - 35*hxy(i,j,k) + 80*hxy(i+di,j,k) &
                        - 30*hxy(i+2*di,j,k) +  8*hxy(i+3*di,j,k) -  hxy(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*hxy(i,j-2*dj,k) - 24*hxy(i,j-dj,k) - 35*hxy(i,j,k) + 80*hxy(i,j+dj,k) &
                        - 30*hxy(i,j+2*dj,k) +  8*hxy(i,j+3*dj,k) -  hxy(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*hxy(i,j,k-2*dk) - 24*hxy(i,j,k-dk) - 35*hxy(i,j,k) + 80*hxy(i,j,k+dk) &
                        - 30*hxy(i,j,k+2*dk) +  8*hxy(i,j,k+3*dk) -  hxy(i,j,k+4*dk) ) * odz60
        ad1_hh(1,2) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * (   2*hxz(i-2*di,j,k) - 24*hxz(i-di,j,k) - 35*hxz(i,j,k) + 80*hxz(i+di,j,k) &
                        - 30*hxz(i+2*di,j,k) +  8*hxz(i+3*di,j,k) -  hxz(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*hxz(i,j-2*dj,k) - 24*hxz(i,j-dj,k) - 35*hxz(i,j,k) + 80*hxz(i,j+dj,k) &
                        - 30*hxz(i,j+2*dj,k) +  8*hxz(i,j+3*dj,k) -  hxz(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*hxz(i,j,k-2*dk) - 24*hxz(i,j,k-dk) - 35*hxz(i,j,k) + 80*hxz(i,j,k+dk) &
                        - 30*hxz(i,j,k+2*dk) +  8*hxz(i,j,k+3*dk) -  hxz(i,j,k+4*dk) ) * odz60
        ad1_hh(1,3) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * (   2*hyy(i-2*di,j,k) - 24*hyy(i-di,j,k) - 35*hyy(i,j,k) + 80*hyy(i+di,j,k) &
                        - 30*hyy(i+2*di,j,k) +  8*hyy(i+3*di,j,k) -  hyy(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*hyy(i,j-2*dj,k) - 24*hyy(i,j-dj,k) - 35*hyy(i,j,k) + 80*hyy(i,j+dj,k) &
                        - 30*hyy(i,j+2*dj,k) +  8*hyy(i,j+3*dj,k) -  hyy(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*hyy(i,j,k-2*dk) - 24*hyy(i,j,k-dk) - 35*hyy(i,j,k) + 80*hyy(i,j,k+dk) &
                        - 30*hyy(i,j,k+2*dk) +  8*hyy(i,j,k+3*dk) -  hyy(i,j,k+4*dk) ) * odz60
        ad1_hh(2,2) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * (   2*hyz(i-2*di,j,k) - 24*hyz(i-di,j,k) - 35*hyz(i,j,k) + 80*hyz(i+di,j,k) &
                        - 30*hyz(i+2*di,j,k) +  8*hyz(i+3*di,j,k) -  hyz(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*hyz(i,j-2*dj,k) - 24*hyz(i,j-dj,k) - 35*hyz(i,j,k) + 80*hyz(i,j+dj,k) &
                        - 30*hyz(i,j+2*dj,k) +  8*hyz(i,j+3*dj,k) -  hyz(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*hyz(i,j,k-2*dk) - 24*hyz(i,j,k-dk) - 35*hyz(i,j,k) + 80*hyz(i,j,k+dk) &
                        - 30*hyz(i,j,k+2*dk) +  8*hyz(i,j,k+3*dk) -  hyz(i,j,k+4*dk) ) * odz60
        ad1_hh(2,3) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * (   2*hzz(i-2*di,j,k) - 24*hzz(i-di,j,k) - 35*hzz(i,j,k) + 80*hzz(i+di,j,k) &
                        - 30*hzz(i+2*di,j,k) +  8*hzz(i+3*di,j,k) -  hzz(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*hzz(i,j-2*dj,k) - 24*hzz(i,j-dj,k) - 35*hzz(i,j,k) + 80*hzz(i,j+dj,k) &
                        - 30*hzz(i,j+2*dj,k) +  8*hzz(i,j+3*dj,k) -  hzz(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*hzz(i,j,k-2*dk) - 24*hzz(i,j,k-dk) - 35*hzz(i,j,k) + 80*hzz(i,j,k+dk) &
                        - 30*hzz(i,j,k+2*dk) +  8*hzz(i,j,k+3*dk) -  hzz(i,j,k+4*dk) ) * odz60
        ad1_hh(3,3) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        ad1_hh(2,1) = ad1_hh(1,2)
        ad1_hh(3,1) = ad1_hh(1,3)
        ad1_hh(3,2) = ad1_hh(2,3)

        ! ad1_trk
        d1_f(1) = di * (   2*tracek(i-2*di,j,k) - 24*tracek(i-di,j,k) - 35*tracek(i,j,k) + 80*tracek(i+di,j,k) &
                        - 30*tracek(i+2*di,j,k) +  8*tracek(i+3*di,j,k) -  tracek(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*tracek(i,j-2*dj,k) - 24*tracek(i,j-dj,k) - 35*tracek(i,j,k) + 80*tracek(i,j+dj,k) &
                        - 30*tracek(i,j+2*dj,k) +  8*tracek(i,j+3*dj,k) -  tracek(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*tracek(i,j,k-2*dk) - 24*tracek(i,j,k-dk) - 35*tracek(i,j,k) + 80*tracek(i,j,k+dk) &
                        - 30*tracek(i,j,k+2*dk) +  8*tracek(i,j,k+3*dk) -  tracek(i,j,k+4*dk) ) * odz60
        ad1_trk = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        ! ad1_aa(3,3)
        d1_f(1) = di * (   2*axx(i-2*di,j,k) - 24*axx(i-di,j,k) - 35*axx(i,j,k) + 80*axx(i+di,j,k) &
                        - 30*axx(i+2*di,j,k) +  8*axx(i+3*di,j,k) -  axx(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*axx(i,j-2*dj,k) - 24*axx(i,j-dj,k) - 35*axx(i,j,k) + 80*axx(i,j+dj,k) &
                        - 30*axx(i,j+2*dj,k) +  8*axx(i,j+3*dj,k) -  axx(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*axx(i,j,k-2*dk) - 24*axx(i,j,k-dk) - 35*axx(i,j,k) + 80*axx(i,j,k+dk) &
                        - 30*axx(i,j,k+2*dk) +  8*axx(i,j,k+3*dk) -  axx(i,j,k+4*dk) ) * odz60
        ad1_aa(1,1) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * (   2*axy(i-2*di,j,k) - 24*axy(i-di,j,k) - 35*axy(i,j,k) + 80*axy(i+di,j,k) &
                        - 30*axy(i+2*di,j,k) +  8*axy(i+3*di,j,k) -  axy(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*axy(i,j-2*dj,k) - 24*axy(i,j-dj,k) - 35*axy(i,j,k) + 80*axy(i,j+dj,k) &
                        - 30*axy(i,j+2*dj,k) +  8*axy(i,j+3*dj,k) -  axy(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*axy(i,j,k-2*dk) - 24*axy(i,j,k-dk) - 35*axy(i,j,k) + 80*axy(i,j,k+dk) &
                        - 30*axy(i,j,k+2*dk) +  8*axy(i,j,k+3*dk) -  axy(i,j,k+4*dk) ) * odz60
        ad1_aa(1,2) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * (   2*axz(i-2*di,j,k) - 24*axz(i-di,j,k) - 35*axz(i,j,k) + 80*axz(i+di,j,k) &
                        - 30*axz(i+2*di,j,k) +  8*axz(i+3*di,j,k) -  axz(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*axz(i,j-2*dj,k) - 24*axz(i,j-dj,k) - 35*axz(i,j,k) + 80*axz(i,j+dj,k) &
                        - 30*axz(i,j+2*dj,k) +  8*axz(i,j+3*dj,k) -  axz(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*axz(i,j,k-2*dk) - 24*axz(i,j,k-dk) - 35*axz(i,j,k) + 80*axz(i,j,k+dk) &
                        - 30*axz(i,j,k+2*dk) +  8*axz(i,j,k+3*dk) -  axz(i,j,k+4*dk) ) * odz60
        ad1_aa(1,3) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * (   2*ayy(i-2*di,j,k) - 24*ayy(i-di,j,k) - 35*ayy(i,j,k) + 80*ayy(i+di,j,k) &
                        - 30*ayy(i+2*di,j,k) +  8*ayy(i+3*di,j,k) -  ayy(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*ayy(i,j-2*dj,k) - 24*ayy(i,j-dj,k) - 35*ayy(i,j,k) + 80*ayy(i,j+dj,k) &
                        - 30*ayy(i,j+2*dj,k) +  8*ayy(i,j+3*dj,k) -  ayy(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*ayy(i,j,k-2*dk) - 24*ayy(i,j,k-dk) - 35*ayy(i,j,k) + 80*ayy(i,j,k+dk) &
                        - 30*ayy(i,j,k+2*dk) +  8*ayy(i,j,k+3*dk) -  ayy(i,j,k+4*dk) ) * odz60
        ad1_aa(2,2) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * (   2*ayz(i-2*di,j,k) - 24*ayz(i-di,j,k) - 35*ayz(i,j,k) + 80*ayz(i+di,j,k) &
                        - 30*ayz(i+2*di,j,k) +  8*ayz(i+3*di,j,k) -  ayz(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*ayz(i,j-2*dj,k) - 24*ayz(i,j-dj,k) - 35*ayz(i,j,k) + 80*ayz(i,j+dj,k) &
                        - 30*ayz(i,j+2*dj,k) +  8*ayz(i,j+3*dj,k) -  ayz(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*ayz(i,j,k-2*dk) - 24*ayz(i,j,k-dk) - 35*ayz(i,j,k) + 80*ayz(i,j,k+dk) &
                        - 30*ayz(i,j,k+2*dk) +  8*ayz(i,j,k+3*dk) -  ayz(i,j,k+4*dk) ) * odz60
        ad1_aa(2,3) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * (   2*azz(i-2*di,j,k) - 24*azz(i-di,j,k) - 35*azz(i,j,k) + 80*azz(i+di,j,k) &
                        - 30*azz(i+2*di,j,k) +  8*azz(i+3*di,j,k) -  azz(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*azz(i,j-2*dj,k) - 24*azz(i,j-dj,k) - 35*azz(i,j,k) + 80*azz(i,j+dj,k) &
                        - 30*azz(i,j+2*dj,k) +  8*azz(i,j+3*dj,k) -  azz(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*azz(i,j,k-2*dk) - 24*azz(i,j,k-dk) - 35*azz(i,j,k) + 80*azz(i,j,k+dk) &
                        - 30*azz(i,j,k+2*dk) +  8*azz(i,j,k+3*dk) -  azz(i,j,k+4*dk) ) * odz60
        ad1_aa(3,3) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        ad1_aa(2,1) = ad1_aa(1,2)
        ad1_aa(3,1) = ad1_aa(1,3)
        ad1_aa(3,2) = ad1_aa(2,3)

        ! ad1_gammat(3)
        d1_f(1) = di * (   2*gammatx(i-2*di,j,k) - 24*gammatx(i-di,j,k) - 35*gammatx(i,j,k) + 80*gammatx(i+di,j,k) &
                        - 30*gammatx(i+2*di,j,k) + 8*gammatx(i+3*di,j,k) - gammatx(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*gammatx(i,j-2*dj,k) - 24*gammatx(i,j-dj,k) - 35*gammatx(i,j,k) + 80*gammatx(i,j+dj,k) &
                        - 30*gammatx(i,j+2*dj,k) + 8*gammatx(i,j+3*dj,k) - gammatx(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*gammatx(i,j,k-2*dk) - 24*gammatx(i,j,k-dk) - 35*gammatx(i,j,k) + 80*gammatx(i,j,k+dk) &
                        - 30*gammatx(i,j,k+2*dk) + 8*gammatx(i,j,k+3*dk) - gammatx(i,j,k+4*dk) ) * odz60
        ad1_gammat(1) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * (   2*gammaty(i-2*di,j,k) - 24*gammaty(i-di,j,k) - 35*gammaty(i,j,k) + 80*gammaty(i+di,j,k) &
                        - 30*gammaty(i+2*di,j,k) +  8*gammaty(i+3*di,j,k) -  gammaty(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*gammaty(i,j-2*dj,k) - 24*gammaty(i,j-dj,k) - 35*gammaty(i,j,k) + 80*gammaty(i,j+dj,k) &
                        - 30*gammaty(i,j+2*dj,k) +  8*gammaty(i,j+3*dj,k) -  gammaty(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*gammaty(i,j,k-2*dk) - 24*gammaty(i,j,k-dk) - 35*gammaty(i,j,k) + 80*gammaty(i,j,k+dk) &
                        - 30*gammaty(i,j,k+2*dk) +  8*gammaty(i,j,k+3*dk) -  gammaty(i,j,k+4*dk) ) * odz60
        ad1_gammat(2) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * (   2*gammatz(i-2*di,j,k) - 24*gammatz(i-di,j,k) - 35*gammatz(i,j,k) + 80*gammatz(i+di,j,k) &
                        - 30*gammatz(i+2*di,j,k) +  8*gammatz(i+3*di,j,k) -  gammatz(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * (   2*gammatz(i,j-2*dj,k) - 24*gammatz(i,j-dj,k) - 35*gammatz(i,j,k) + 80*gammatz(i,j+dj,k) &
                        - 30*gammatz(i,j+2*dj,k) +  8*gammatz(i,j+3*dj,k) -  gammatz(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * (   2*gammatz(i,j,k-2*dk) - 24*gammatz(i,j,k-dk) - 35*gammatz(i,j,k) + 80*gammatz(i,j,k+dk) &
                        - 30*gammatz(i,j,k+2*dk) +  8*gammatz(i,j,k+3*dk) -  gammatz(i,j,k+4*dk) ) * odz60
        ad1_gammat(3) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        ! ad1_alph
        if (evolve_alp) then
           d1_f(1) = di * (   2*alp(i-2*di,j,k) - 24*alp(i-di,j,k) - 35*alp(i,j,k) + 80*alp(i+di,j,k) &
                - 30*alp(i+2*di,j,k) +  8*alp(i+3*di,j,k) -  alp(i+4*di,j,k) ) * odx60
           d1_f(2) = dj * (   2*alp(i,j-2*dj,k) - 24*alp(i,j-dj,k) - 35*alp(i,j,k) + 80*alp(i,j+dj,k) &
                - 30*alp(i,j+2*dj,k) +  8*alp(i,j+3*dj,k) -  alp(i,j+4*dj,k) ) * ody60
           d1_f(3) = dk * (   2*alp(i,j,k-2*dk) - 24*alp(i,j,k-dk) - 35*alp(i,j,k) + 80*alp(i,j,k+dk) &
                - 30*alp(i,j,k+2*dk) +  8*alp(i,j,k+3*dk) -  alp(i,j,k+4*dk) ) * odz60
           ad1_alph = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)
        end if

        ! ad1_beta(3)
        if (evolve_beta) then
           d1_f(1) = di * (   2*betax(i-2*di,j,k) - 24*betax(i-di,j,k) - 35*betax(i,j,k) + 80*betax(i+di,j,k) &
                - 30*betax(i+2*di,j,k) +  8*betax(i+3*di,j,k) -  betax(i+4*di,j,k) ) * odx60
           d1_f(2) = dj * (   2*betax(i,j-2*dj,k) - 24*betax(i,j-dj,k) - 35*betax(i,j,k) + 80*betax(i,j+dj,k) &
                - 30*betax(i,j+2*dj,k) +  8*betax(i,j+3*dj,k) -  betax(i,j+4*dj,k) ) * ody60
           d1_f(3) = dk * (   2*betax(i,j,k-2*dk) - 24*betax(i,j,k-dk) - 35*betax(i,j,k) + 80*betax(i,j,k+dk) &
                - 30*betax(i,j,k+2*dk) +  8*betax(i,j,k+3*dk) -  betax(i,j,k+4*dk) ) * odz60
           ad1_beta(1) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

           d1_f(1) = di * (   2*betay(i-2*di,j,k) - 24*betay(i-di,j,k) - 35*betay(i,j,k) + 80*betay(i+di,j,k) &
                - 30*betay(i+2*di,j,k) +  8*betay(i+3*di,j,k) -  betay(i+4*di,j,k) ) * odx60
           d1_f(2) = dj * (   2*betay(i,j-2*dj,k) - 24*betay(i,j-dj,k) - 35*betay(i,j,k) + 80*betay(i,j+dj,k) &
                - 30*betay(i,j+2*dj,k) +  8*betay(i,j+3*dj,k) -  betay(i,j+4*dj,k) ) * ody60
           d1_f(3) = dk * (   2*betay(i,j,k-2*dk) - 24*betay(i,j,k-dk) - 35*betay(i,j,k) + 80*betay(i,j,k+dk) &
                - 30*betay(i,j,k+2*dk) +  8*betay(i,j,k+3*dk) -  betay(i,j,k+4*dk) ) * odz60
           ad1_beta(2) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

           d1_f(1) = di * (   2*betaz(i-2*di,j,k) - 24*betaz(i-di,j,k) - 35*betaz(i,j,k) + 80*betaz(i+di,j,k) &
                - 30*betaz(i+2*di,j,k) +  8*betaz(i+3*di,j,k) -  betaz(i+4*di,j,k) ) * odx60
           d1_f(2) = dj * (   2*betaz(i,j-2*dj,k) - 24*betaz(i,j-dj,k) - 35*betaz(i,j,k) + 80*betaz(i,j+dj,k) &
                - 30*betaz(i,j+2*dj,k) +  8*betaz(i,j+3*dj,k) -  betaz(i,j+4*dj,k) ) * ody60
           d1_f(3) = dk * (   2*betaz(i,j,k-2*dk) - 24*betaz(i,j,k-dk) - 35*betaz(i,j,k) + 80*betaz(i,j,k+dk) &
                - 30*betaz(i,j,k+2*dk) +  8*betaz(i,j,k+3*dk) -  betaz(i,j,k+4*dk) ) * odz60
           ad1_beta(3) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)
        end if
      else
        ! We can use the already calculated expressions for the
        ! first derivatives.

        ! ad1_ww
        ad1_ww  = beta_l(1)*d1_ww(1) + beta_l(2)*d1_ww(2) + beta_l(3)*d1_ww(3)
        
        if (evolve_scalar) then
        ! ad1_phi
        ad1_lphi = beta_l(1)*d1_lphi(1) + beta_l(2)*d1_lphi(2) + beta_l(3)*d1_lphi(3)
       
        ! ad1_Kphi
        ad1_lKphi = beta_l(1)*d1_lKphi(1) + beta_l(2)*d1_lKphi(2) + beta_l(3)*d1_lKphi(3)
        end if

        ! ad1_hh(3,3)
        ad1_hh(1,1) = beta_l(1)*d1_hh11(1) + beta_l(2)*d1_hh11(2) + beta_l(3)*d1_hh11(3)
        ad1_hh(1,2) = beta_l(1)*d1_hh12(1) + beta_l(2)*d1_hh12(2) + beta_l(3)*d1_hh12(3)
        ad1_hh(1,3) = beta_l(1)*d1_hh13(1) + beta_l(2)*d1_hh13(2) + beta_l(3)*d1_hh13(3)
        ad1_hh(2,2) = beta_l(1)*d1_hh22(1) + beta_l(2)*d1_hh22(2) + beta_l(3)*d1_hh22(3)
        ad1_hh(2,3) = beta_l(1)*d1_hh23(1) + beta_l(2)*d1_hh23(2) + beta_l(3)*d1_hh23(3)
        ad1_hh(3,3) = beta_l(1)*d1_hh33(1) + beta_l(2)*d1_hh33(2) + beta_l(3)*d1_hh33(3)

        ad1_hh(2,1) = ad1_hh(1,2)
        ad1_hh(3,1) = ad1_hh(1,3)
        ad1_hh(3,2) = ad1_hh(2,3)

        ! ad1_trk
        ad1_trk = beta_l(1)*d1_trk(1) + beta_l(2)*d1_trk(2) + beta_l(3)*d1_trk(3)

        ! ad1_aa(3,3)
        ad1_aa(1,1) = beta_l(1)*d1_aa11(1) + beta_l(2)*d1_aa11(2) + beta_l(3)*d1_aa11(3)
        ad1_aa(1,2) = beta_l(1)*d1_aa12(1) + beta_l(2)*d1_aa12(2) + beta_l(3)*d1_aa12(3)
        ad1_aa(1,3) = beta_l(1)*d1_aa13(1) + beta_l(2)*d1_aa13(2) + beta_l(3)*d1_aa13(3)
        ad1_aa(2,2) = beta_l(1)*d1_aa22(1) + beta_l(2)*d1_aa22(2) + beta_l(3)*d1_aa22(3)
        ad1_aa(2,3) = beta_l(1)*d1_aa23(1) + beta_l(2)*d1_aa23(2) + beta_l(3)*d1_aa23(3)
        ad1_aa(3,3) = beta_l(1)*d1_aa33(1) + beta_l(2)*d1_aa33(2) + beta_l(3)*d1_aa33(3)

        ad1_aa(2,1) = ad1_aa(1,2)
        ad1_aa(3,1) = ad1_aa(1,3)
        ad1_aa(3,2) = ad1_aa(2,3)

        ! ad1_gammat(3)
        ad1_gammat(1) = beta_l(1)*d1_gammat1(1) + beta_l(2)*d1_gammat1(2) + beta_l(3)*d1_gammat1(3)
        ad1_gammat(2) = beta_l(1)*d1_gammat2(1) + beta_l(2)*d1_gammat2(2) + beta_l(3)*d1_gammat2(3)
        ad1_gammat(3) = beta_l(1)*d1_gammat3(1) + beta_l(2)*d1_gammat3(2) + beta_l(3)*d1_gammat3(3)

        ! ad1_alph
        if (evolve_alp) then
           ad1_alph = beta_l(1)*d1_alph(1) + beta_l(2)*d1_alph(2) + beta_l(3)*d1_alph(3)
        end if

        ! ad1_beta(3)
        if (evolve_beta) then
           ad1_beta(1) = beta_l(1)*d1_beta1(1) + beta_l(2)*d1_beta1(2) + beta_l(3)*d1_beta1(3)
           ad1_beta(2) = beta_l(1)*d1_beta2(1) + beta_l(2)*d1_beta2(2) + beta_l(3)*d1_beta2(3)
           ad1_beta(3) = beta_l(1)*d1_beta3(1) + beta_l(2)*d1_beta3(2) + beta_l(3)*d1_beta3(3)
        end if

        !if( abs(y(i,j,k)) < 1.0d-05 .and. abs(z(i,j,k)) < 1.0d-05 ) then
        !  write(*,*) 'i, j, k, x = ', i, j, k, x(i,j,k)
        !  write(*,*) 'ad1_ww   = ', ad1_ww
        !  write(*,*) 'ad1_hh   = ', ad1_hh
        !  write(*,*) 'ad1_trk  = ', ad1_trk
        !  write(*,*) 'ad1_aa   = ', ad1_aa
        !  write(*,*) 'ad1_alph = ', ad1_alph
        !  write(*,*) 'ad1_gam  = ', ad1_gammat
        !  call flush(6)
        !end if
      end if ! use_advection_stencils /= 0

    else ! if derivs_order == 4

      !------------ Centered 1st derivatives -----
      ! d1_ww(3)
      d1_ww(1) = (   -conf_fac(i+2,j,k) + 8*conf_fac(i+1,j,k)                          &
                  - 8*conf_fac(i-1,j,k) +   conf_fac(i-2,j,k) ) / dx12

      d1_ww(2) = (   -conf_fac(i,j+2,k) + 8*conf_fac(i,j+1,k)                          &
                  - 8*conf_fac(i,j-1,k) +   conf_fac(i,j-2,k) ) / dy12

      d1_ww(3) = (   -conf_fac(i,j,k+2) + 8*conf_fac(i,j,k+1)                          &
                  - 8*conf_fac(i,j,k-1) +   conf_fac(i,j,k-2) ) / dz12
      
      if (evolve_scalar) then 
      ! d1_phi1(3)
      d1_lphi(1) = (   -phi1(i+2,j,k) + 8*phi1(i+1,j,k)          &
                     -8*phi1(i-1,j,k) +   phi1(i-2,j,k) ) / dx12 
      
      d1_lphi(2) = (  -phi1(i,j+2,k)  + 8*phi1(i,j+1,k)          &
                     -8*phi1(i,j-1,k) +   phi1(i,j-2,k) ) / dy12

      d1_lphi(3) = ( -phi1(i,j,k+2)  + 8*phi1(i,j,k+1)           &
                    -8*phi1(i,j,k-1) +   phi1(i,j,k-2) )  / dz12  

      ! d1_Kphi1(3)
      d1_lKphi(1) = (  -Kphi1(i+2,j,k) + 8*Kphi1(i+1,j,k)          &
                     -8*Kphi1(i-1,j,k) +   Kphi1(i-2,j,k) ) / dx12 
      
      d1_lKphi(2) = (  -Kphi1(i,j+2,k) + 8*Kphi1(i,j+1,k)          &
                     -8*Kphi1(i,j-1,k) +   Kphi1(i,j-2,k) ) / dy12

      d1_lKphi(3) = ( -Kphi1(i,j,k+2) + 8*Kphi1(i,j,k+1)           &
                    -8*Kphi1(i,j,k-1) +   Kphi1(i,j,k-2) )  / dz12  
      end if 

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

      ! d1_alph(3)
      d1_alph(1) = (   -alp(i+2,j,k) + 8*alp(i+1,j,k)                        &
                    - 8*alp(i-1,j,k) +   alp(i-2,j,k) ) / dx12

      d1_alph(2) = (   -alp(i,j+2,k) + 8*alp(i,j+1,k)                        &
                    - 8*alp(i,j-1,k) +   alp(i,j-2,k) ) / dy12

      d1_alph(3) = (   -alp(i,j,k+2) + 8*alp(i,j,k+1)                        &
                    - 8*alp(i,j,k-1) +   alp(i,j,k-2) ) / dz12

      ! d1_beta(3,3)
      d1_beta1(1)  = (   -betax(i+2,j,k) + 8*betax(i+1,j,k)                 &
                      - 8*betax(i-1,j,k) +   betax(i-2,j,k) ) / dx12
      d1_beta2(1)  = (   -betay(i+2,j,k) + 8*betay(i+1,j,k)                 &
                      - 8*betay(i-1,j,k) +   betay(i-2,j,k) ) / dx12
      d1_beta3(1)  = (   -betaz(i+2,j,k) + 8*betaz(i+1,j,k)                 &
                      - 8*betaz(i-1,j,k) +   betaz(i-2,j,k) ) / dx12

      d1_beta1(2)  = (   -betax(i,j+2,k) + 8*betax(i,j+1,k)                 &
                      - 8*betax(i,j-1,k) +   betax(i,j-2,k) ) / dy12
      d1_beta2(2)  = (   -betay(i,j+2,k) + 8*betay(i,j+1,k)                 &
                      - 8*betay(i,j-1,k) +   betay(i,j-2,k) ) / dy12
      d1_beta3(2)  = (   -betaz(i,j+2,k) + 8*betaz(i,j+1,k)                 &
                      - 8*betaz(i,j-1,k) +   betaz(i,j-2,k) ) / dy12

      d1_beta1(3)  = (   -betax(i,j,k+2) + 8*betax(i,j,k+1)                 &
                      - 8*betax(i,j,k-1) +   betax(i,j,k-2) ) / dz12
      d1_beta2(3)  = (   -betay(i,j,k+2) + 8*betay(i,j,k+1)                 &
                      - 8*betay(i,j,k-1) +   betay(i,j,k-2) ) / dz12
      d1_beta3(3)  = (   -betaz(i,j,k+2) + 8*betaz(i,j,k+1)                 &
                      - 8*betaz(i,j,k-1) +   betaz(i,j,k-2) ) / dz12
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

      if (evolve_scalar) then
      !d2_phi(3,3)


      d2_lphi(1,1) = (   -phi1(i+2,j,k) + 16*phi1(i+1,j,k) - 30*phi1(i,j,k)     &
                     + 16*phi1(i-1,j,k) -    phi1(i-2,j,k) ) / dxsq12

      d2_lphi(2,2) = (   -phi1(i,j+2,k) + 16*phi1(i,j+1,k) - 30*phi1(i,j,k)     &
                     + 16*phi1(i,j-1,k) -    phi1(i,j-2,k) ) / dysq12

      d2_lphi(3,3) = (   -phi1(i,j,k+2) + 16*phi1(i,j,k+1) - 30*phi1(i,j,k)     &
                     + 16*phi1(i,j,k-1) -    phi1(i,j,k-2) ) / dzsq12

      d2_lphi(1,2) = (   -phi1(i-2,j+2,k) +  8*phi1(i-1,j+2,k) -  8*phi1(i+1,j+2,k) +   phi1(i+2,j+2,k)   &
                      + 8*phi1(i-2,j+1,k) - 64*phi1(i-1,j+1,k) + 64*phi1(i+1,j+1,k) - 8*phi1(i+2,j+1,k)   &
                      - 8*phi1(i-2,j-1,k) + 64*phi1(i-1,j-1,k) - 64*phi1(i+1,j-1,k) + 8*phi1(i+2,j-1,k)   &
                      +   phi1(i-2,j-2,k) -  8*phi1(i-1,j-2,k) +  8*phi1(i+1,j-2,k) -   phi1(i+2,j-2,k) ) / dxdy144

      d2_lphi(1,3) = (   -phi1(i-2,j,k+2) +  8*phi1(i-1,j,k+2) -  8*phi1(i+1,j,k+2) +   phi1(i+2,j,k+2)   &
                      + 8*phi1(i-2,j,k+1) - 64*phi1(i-1,j,k+1) + 64*phi1(i+1,j,k+1) - 8*phi1(i+2,j,k+1)   &
                      - 8*phi1(i-2,j,k-1) + 64*phi1(i-1,j,k-1) - 64*phi1(i+1,j,k-1) + 8*phi1(i+2,j,k-1)   &
                      +   phi1(i-2,j,k-2) -  8*phi1(i-1,j,k-2) +  8*phi1(i+1,j,k-2) -   phi1(i+2,j,k-2) ) / dxdz144

      d2_lphi(2,3) = (   -phi1(i,j-2,k+2) +  8*phi1(i,j-1,k+2) -  8*phi1(i,j+1,k+2) +   phi1(i,j+2,k+2)   &
                      + 8*phi1(i,j-2,k+1) - 64*phi1(i,j-1,k+1) + 64*phi1(i,j+1,k+1) - 8*phi1(i,j+2,k+1)   &
                      - 8*phi1(i,j-2,k-1) + 64*phi1(i,j-1,k-1) - 64*phi1(i,j+1,k-1) + 8*phi1(i,j+2,k-1)   &
                      +   phi1(i,j-2,k-2) -  8*phi1(i,j-1,k-2) +  8*phi1(i,j+1,k-2) -   phi1(i,j+2,k-2) ) / dydz144

      d2_lphi(2,1) = d2_lphi(1,2)
      d2_lphi(3,1) = d2_lphi(1,3)
      d2_lphi(3,2) = d2_lphi(2,3)
!!!
      end if

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

      ! d2_alph(3,3)
      d2_alph(1,1) = (   -alp(i+2,j,k) + 16*alp(i+1,j,k) - 30*alp(i,j,k)     &
                     + 16*alp(i-1,j,k) -    alp(i-2,j,k) ) / dxsq12

      d2_alph(2,2) = (   -alp(i,j+2,k) + 16*alp(i,j+1,k) - 30*alp(i,j,k)     &
                     + 16*alp(i,j-1,k) -    alp(i,j-2,k) ) / dysq12

      d2_alph(3,3) = (   -alp(i,j,k+2) + 16*alp(i,j,k+1) - 30*alp(i,j,k)     &
                     + 16*alp(i,j,k-1) -    alp(i,j,k-2) ) / dzsq12

      d2_alph(1,2) = (   -alp(i-2,j+2,k) +  8*alp(i-1,j+2,k) -  8*alp(i+1,j+2,k) +   alp(i+2,j+2,k)   &
                      + 8*alp(i-2,j+1,k) - 64*alp(i-1,j+1,k) + 64*alp(i+1,j+1,k) - 8*alp(i+2,j+1,k)   &
                      - 8*alp(i-2,j-1,k) + 64*alp(i-1,j-1,k) - 64*alp(i+1,j-1,k) + 8*alp(i+2,j-1,k)   &
                      +   alp(i-2,j-2,k) -  8*alp(i-1,j-2,k) +  8*alp(i+1,j-2,k) -   alp(i+2,j-2,k) ) / dxdy144

      d2_alph(1,3) = (   -alp(i-2,j,k+2) +  8*alp(i-1,j,k+2) -  8*alp(i+1,j,k+2) +   alp(i+2,j,k+2)   &
                      + 8*alp(i-2,j,k+1) - 64*alp(i-1,j,k+1) + 64*alp(i+1,j,k+1) - 8*alp(i+2,j,k+1)   &
                      - 8*alp(i-2,j,k-1) + 64*alp(i-1,j,k-1) - 64*alp(i+1,j,k-1) + 8*alp(i+2,j,k-1)   &
                      +   alp(i-2,j,k-2) -  8*alp(i-1,j,k-2) +  8*alp(i+1,j,k-2) -   alp(i+2,j,k-2) ) / dxdz144

      d2_alph(2,3) = (   -alp(i,j-2,k+2) +  8*alp(i,j-1,k+2) -  8*alp(i,j+1,k+2) +   alp(i,j+2,k+2)   &
                      + 8*alp(i,j-2,k+1) - 64*alp(i,j-1,k+1) + 64*alp(i,j+1,k+1) - 8*alp(i,j+2,k+1)   &
                      - 8*alp(i,j-2,k-1) + 64*alp(i,j-1,k-1) - 64*alp(i,j+1,k-1) + 8*alp(i,j+2,k-1)   &
                      +   alp(i,j-2,k-2) -  8*alp(i,j-1,k-2) +  8*alp(i,j+1,k-2) -   alp(i,j+2,k-2) ) / dydz144

      d2_alph(2,1) = d2_alph(1,2)
      d2_alph(3,1) = d2_alph(1,3)
      d2_alph(3,2) = d2_alph(2,3)

      ! d2_beta(3,3,3)
      d2_beta1(1,1) = (  -betax(i+2,j,k) + 16*betax(i+1,j,k) - 30*betax(i,j,k)&
                     + 16*betax(i-1,j,k) -    betax(i-2,j,k) ) / dxsq12
      d2_beta2(1,1) = (  -betay(i+2,j,k) + 16*betay(i+1,j,k) - 30*betay(i,j,k)&
                     + 16*betay(i-1,j,k) -    betay(i-2,j,k) ) / dxsq12
      d2_beta3(1,1) = (  -betaz(i+2,j,k) + 16*betaz(i+1,j,k) - 30*betaz(i,j,k)&
                     + 16*betaz(i-1,j,k) -    betaz(i-2,j,k) ) / dxsq12

      d2_beta1(2,2) = (  -betax(i,j+2,k) + 16*betax(i,j+1,k) - 30*betax(i,j,k)&
                     + 16*betax(i,j-1,k) -    betax(i,j-2,k) ) / dysq12
      d2_beta2(2,2) = (  -betay(i,j+2,k) + 16*betay(i,j+1,k) - 30*betay(i,j,k)&
                     + 16*betay(i,j-1,k) -    betay(i,j-2,k) ) / dysq12
      d2_beta3(2,2) = (  -betaz(i,j+2,k) + 16*betaz(i,j+1,k) - 30*betaz(i,j,k)&
                     + 16*betaz(i,j-1,k) -    betaz(i,j-2,k) ) / dysq12

      d2_beta1(3,3) = (  -betax(i,j,k+2) + 16*betax(i,j,k+1) - 30*betax(i,j,k)&
                     + 16*betax(i,j,k-1) -    betax(i,j,k-2) ) / dzsq12
      d2_beta2(3,3) = (  -betay(i,j,k+2) + 16*betay(i,j,k+1) - 30*betay(i,j,k)&
                     + 16*betay(i,j,k-1) -    betay(i,j,k-2) ) / dzsq12
      d2_beta3(3,3) = (  -betaz(i,j,k+2) + 16*betaz(i,j,k+1) - 30*betaz(i,j,k)&
                     + 16*betaz(i,j,k-1) -    betaz(i,j,k-2) ) / dzsq12

      d2_beta1(1,2) = (  -betax(i-2,j+2,k) +  8*betax(i-1,j+2,k) -  8*betax(i+1,j+2,k) +   betax(i+2,j+2,k)   &
                      + 8*betax(i-2,j+1,k) - 64*betax(i-1,j+1,k) + 64*betax(i+1,j+1,k) - 8*betax(i+2,j+1,k)   &
                      - 8*betax(i-2,j-1,k) + 64*betax(i-1,j-1,k) - 64*betax(i+1,j-1,k) + 8*betax(i+2,j-1,k)   &
                      +   betax(i-2,j-2,k) -  8*betax(i-1,j-2,k) +  8*betax(i+1,j-2,k) -   betax(i+2,j-2,k) ) / dxdy144
      d2_beta2(1,2) = (  -betay(i-2,j+2,k) +  8*betay(i-1,j+2,k) -  8*betay(i+1,j+2,k) +   betay(i+2,j+2,k)   &
                      + 8*betay(i-2,j+1,k) - 64*betay(i-1,j+1,k) + 64*betay(i+1,j+1,k) - 8*betay(i+2,j+1,k)   &
                      - 8*betay(i-2,j-1,k) + 64*betay(i-1,j-1,k) - 64*betay(i+1,j-1,k) + 8*betay(i+2,j-1,k)   &
                      +   betay(i-2,j-2,k) -  8*betay(i-1,j-2,k) +  8*betay(i+1,j-2,k) -   betay(i+2,j-2,k) ) / dxdy144
      d2_beta3(1,2) = (  -betaz(i-2,j+2,k) +  8*betaz(i-1,j+2,k) -  8*betaz(i+1,j+2,k) +   betaz(i+2,j+2,k)   &
                      + 8*betaz(i-2,j+1,k) - 64*betaz(i-1,j+1,k) + 64*betaz(i+1,j+1,k) - 8*betaz(i+2,j+1,k)   &
                      - 8*betaz(i-2,j-1,k) + 64*betaz(i-1,j-1,k) - 64*betaz(i+1,j-1,k) + 8*betaz(i+2,j-1,k)   &
                      +   betaz(i-2,j-2,k) -  8*betaz(i-1,j-2,k) +  8*betaz(i+1,j-2,k) -   betaz(i+2,j-2,k) ) / dxdy144

      d2_beta1(1,3) = (  -betax(i-2,j,k+2) +  8*betax(i-1,j,k+2) -  8*betax(i+1,j,k+2) +   betax(i+2,j,k+2)   &
                      + 8*betax(i-2,j,k+1) - 64*betax(i-1,j,k+1) + 64*betax(i+1,j,k+1) - 8*betax(i+2,j,k+1)   &
                      - 8*betax(i-2,j,k-1) + 64*betax(i-1,j,k-1) - 64*betax(i+1,j,k-1) + 8*betax(i+2,j,k-1)   &
                      +   betax(i-2,j,k-2) -  8*betax(i-1,j,k-2) +  8*betax(i+1,j,k-2) -   betax(i+2,j,k-2) ) / dxdz144
      d2_beta2(1,3) = (  -betay(i-2,j,k+2) +  8*betay(i-1,j,k+2) -  8*betay(i+1,j,k+2) +   betay(i+2,j,k+2)   &
                      + 8*betay(i-2,j,k+1) - 64*betay(i-1,j,k+1) + 64*betay(i+1,j,k+1) - 8*betay(i+2,j,k+1)   &
                      - 8*betay(i-2,j,k-1) + 64*betay(i-1,j,k-1) - 64*betay(i+1,j,k-1) + 8*betay(i+2,j,k-1)   &
                      +   betay(i-2,j,k-2) -  8*betay(i-1,j,k-2) +  8*betay(i+1,j,k-2) -   betay(i+2,j,k-2) ) / dxdz144
      d2_beta3(1,3) = (  -betaz(i-2,j,k+2) +  8*betaz(i-1,j,k+2) -  8*betaz(i+1,j,k+2) +   betaz(i+2,j,k+2)   &
                      + 8*betaz(i-2,j,k+1) - 64*betaz(i-1,j,k+1) + 64*betaz(i+1,j,k+1) - 8*betaz(i+2,j,k+1)   &
                      - 8*betaz(i-2,j,k-1) + 64*betaz(i-1,j,k-1) - 64*betaz(i+1,j,k-1) + 8*betaz(i+2,j,k-1)   &
                      +   betaz(i-2,j,k-2) -  8*betaz(i-1,j,k-2) +  8*betaz(i+1,j,k-2) -   betaz(i+2,j,k-2) ) / dxdz144

      d2_beta1(2,3) = (  -betax(i,j-2,k+2) +  8*betax(i,j-1,k+2) -  8*betax(i,j+1,k+2) +   betax(i,j+2,k+2)   &
                      + 8*betax(i,j-2,k+1) - 64*betax(i,j-1,k+1) + 64*betax(i,j+1,k+1) - 8*betax(i,j+2,k+1)   &
                      - 8*betax(i,j-2,k-1) + 64*betax(i,j-1,k-1) - 64*betax(i,j+1,k-1) + 8*betax(i,j+2,k-1)   &
                      +   betax(i,j-2,k-2) -  8*betax(i,j-1,k-2) +  8*betax(i,j+1,k-2) -   betax(i,j+2,k-2) ) / dydz144
      d2_beta2(2,3) = (  -betay(i,j-2,k+2) +  8*betay(i,j-1,k+2) -  8*betay(i,j+1,k+2) +   betay(i,j+2,k+2)   &
                      + 8*betay(i,j-2,k+1) - 64*betay(i,j-1,k+1) + 64*betay(i,j+1,k+1) - 8*betay(i,j+2,k+1)   &
                      - 8*betay(i,j-2,k-1) + 64*betay(i,j-1,k-1) - 64*betay(i,j+1,k-1) + 8*betay(i,j+2,k-1)   &
                      +   betay(i,j-2,k-2) -  8*betay(i,j-1,k-2) +  8*betay(i,j+1,k-2) -   betay(i,j+2,k-2) ) / dydz144
      d2_beta3(2,3) = (  -betaz(i,j-2,k+2) +  8*betaz(i,j-1,k+2) -  8*betaz(i,j+1,k+2) +   betaz(i,j+2,k+2)   &
                      + 8*betaz(i,j-2,k+1) - 64*betaz(i,j-1,k+1) + 64*betaz(i,j+1,k+1) - 8*betaz(i,j+2,k+1)   &
                      - 8*betaz(i,j-2,k-1) + 64*betaz(i,j-1,k-1) - 64*betaz(i,j+1,k-1) + 8*betaz(i,j+2,k-1)   &
                      +   betaz(i,j-2,k-2) -  8*betaz(i,j-1,k-2) +  8*betaz(i,j+1,k-2) -   betaz(i,j+2,k-2) ) / dydz144

      d2_beta1(2,1) = d2_beta1(1,2)
      d2_beta2(2,1) = d2_beta2(1,2)
      d2_beta3(2,1) = d2_beta3(1,2)
      d2_beta1(3,1) = d2_beta1(1,3)
      d2_beta2(3,1) = d2_beta2(1,3)
      d2_beta3(3,1) = d2_beta3(1,3)
      d2_beta1(3,2) = d2_beta1(2,3)
      d2_beta2(3,2) = d2_beta2(2,3)
      d2_beta3(3,2) = d2_beta3(2,3)

      !--------------------------------------------------


      !------------ Advection derivatives --------
      if( use_advection_stencils /= 0 ) then

        di = int( sign( one, beta_l(1) ) )
        dj = int( sign( one, beta_l(2) ) )
        dk = int( sign( one, beta_l(3) ) )

        ! ad1_ph
        d1_f(1) = di * ( -3*conf_fac(i-di,j,k) - 10*conf_fac(i,j,k) + 18*conf_fac(i+di,j,k)   &
                        - 6*conf_fac(i+2*di,j,k) + conf_fac(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*conf_fac(i,j-dj,k) - 10*conf_fac(i,j,k) + 18*conf_fac(i,j+dj,k)   &
                        - 6*conf_fac(i,j+2*dj,k) + conf_fac(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*conf_fac(i,j,k-dk) - 10*conf_fac(i,j,k) + 18*conf_fac(i,j,k+dk)   &
                        - 6*conf_fac(i,j,k+2*dk) + conf_fac(i,j,k+3*dk)) / dz12
        ad1_ww  = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        if (evolve_scalar) then
        ! ad1_phi
      d1_f(1) = di * ( -3*phi1(i-di,j,k) - 10*phi1(i,j,k) + 18*phi1(i+di,j,k)   &
                      - 6*phi1(i+2*di,j,k)  + phi1(i+3*di,j,k)) / dx12
      d1_f(2) = dj * ( -3*phi1(i,j-dj,k) - 10*phi1(i,j,k) + 18*phi1(i,j+dj,k)   &
                      - 6*phi1(i,j+2*dj,k)  + phi1(i,j+3*dj,k)) / dy12
      d1_f(3) = dk * ( -3*phi1(i,j,k-dk) - 10*phi1(i,j,k) + 18*phi1(i,j,k+dk)   &
                      - 6*phi1(i,j,k+2*dk)  + phi1(i,j,k+3*dk)) / dz12
      ad1_lphi = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)




        ! ad1_Kphi
        d1_f(1) = di * ( -3*Kphi1(i-di,j,k) - 10*Kphi1(i,j,k) + 18*Kphi1(i+di,j,k) &
                        - 6*Kphi1(i+2*di,j,k) + Kphi1(i+3*di,j,k) ) / dx12
        d1_f(2) = dj * ( -3*Kphi1(i,j-dj,k) - 10*Kphi1(i,j,k) + 18*Kphi1(i,j+dj,k) &
                        - 6*Kphi1(i,j+2*dj,k) + Kphi1(i,j+3*dj,k) ) / dy12
        d1_f(3) = dk * ( -3*Kphi1(i,j,k-dk) - 10*Kphi1(i,j,k) + 18*Kphi1(i,j,k+dk) &
                        - 6*Kphi1(i,j,k+2*dk) + Kphi1(i,j,k+3*dk) ) / dz12
        ad1_lKphi = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)
        end if

        ! ad1_hh(3,3)
        d1_f(1) = di * ( -3*hxx(i-di,j,k) - 10*hxx(i,j,k) + 18*hxx(i+di,j,k)   &
                        - 6*hxx(i+2*di,j,k) + hxx(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*hxx(i,j-dj,k) - 10*hxx(i,j,k) + 18*hxx(i,j+dj,k)   &
                        - 6*hxx(i,j+2*dj,k) + hxx(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*hxx(i,j,k-dk) - 10*hxx(i,j,k) + 18*hxx(i,j,k+dk)   &
                        - 6*hxx(i,j,k+2*dk) + hxx(i,j,k+3*dk)) / dz12
        ad1_hh(1,1) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * ( -3*hxy(i-di,j,k) - 10*hxy(i,j,k) + 18*hxy(i+di,j,k)   &
                        - 6*hxy(i+2*di,j,k) + hxy(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*hxy(i,j-dj,k) - 10*hxy(i,j,k) + 18*hxy(i,j+dj,k)   &
                        - 6*hxy(i,j+2*dj,k) + hxy(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*hxy(i,j,k-dk) - 10*hxy(i,j,k) + 18*hxy(i,j,k+dk)   &
                        - 6*hxy(i,j,k+2*dk) + hxy(i,j,k+3*dk)) / dz12
        ad1_hh(1,2) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * ( -3*hxz(i-di,j,k) - 10*hxz(i,j,k) + 18*hxz(i+di,j,k)   &
                        - 6*hxz(i+2*di,j,k) + hxz(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*hxz(i,j-dj,k) - 10*hxz(i,j,k) + 18*hxz(i,j+dj,k)   &
                        - 6*hxz(i,j+2*dj,k) + hxz(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*hxz(i,j,k-dk) - 10*hxz(i,j,k) + 18*hxz(i,j,k+dk)   &
                        - 6*hxz(i,j,k+2*dk) + hxz(i,j,k+3*dk)) / dz12
        ad1_hh(1,3) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * ( -3*hyy(i-di,j,k) - 10*hyy(i,j,k) + 18*hyy(i+di,j,k)   &
                        - 6*hyy(i+2*di,j,k) + hyy(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*hyy(i,j-dj,k) - 10*hyy(i,j,k) + 18*hyy(i,j+dj,k)   &
                        - 6*hyy(i,j+2*dj,k) + hyy(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*hyy(i,j,k-dk) - 10*hyy(i,j,k) + 18*hyy(i,j,k+dk)   &
                        - 6*hyy(i,j,k+2*dk) + hyy(i,j,k+3*dk)) / dz12
        ad1_hh(2,2) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * ( -3*hyz(i-di,j,k) - 10*hyz(i,j,k) + 18*hyz(i+di,j,k)   &
                        - 6*hyz(i+2*di,j,k) + hyz(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*hyz(i,j-dj,k) - 10*hyz(i,j,k) + 18*hyz(i,j+dj,k)   &
                        - 6*hyz(i,j+2*dj,k) + hyz(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*hyz(i,j,k-dk) - 10*hyz(i,j,k) + 18*hyz(i,j,k+dk)   &
                        - 6*hyz(i,j,k+2*dk) + hyz(i,j,k+3*dk)) / dz12
        ad1_hh(2,3) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * ( -3*hzz(i-di,j,k) - 10*hzz(i,j,k) + 18*hzz(i+di,j,k)   &
                        - 6*hzz(i+2*di,j,k) + hzz(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*hzz(i,j-dj,k) - 10*hzz(i,j,k) + 18*hzz(i,j+dj,k)   &
                        - 6*hzz(i,j+2*dj,k) + hzz(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*hzz(i,j,k-dk) - 10*hzz(i,j,k) + 18*hzz(i,j,k+dk)   &
                        - 6*hzz(i,j,k+2*dk) + hzz(i,j,k+3*dk)) / dz12
        ad1_hh(3,3) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        ad1_hh(2,1) = ad1_hh(1,2)
        ad1_hh(3,1) = ad1_hh(1,3)
        ad1_hh(3,2) = ad1_hh(2,3)

        ! ad1_trk
        d1_f(1) =di*(-3*tracek(i-di,j,k) - 10*tracek(i,j,k) + 18*tracek(i+di,j,k)&
                    - 6*tracek(i+2*di,j,k) + tracek(i+3*di,j,k)) / dx12
        d1_f(2) =dj*(-3*tracek(i,j-dj,k) - 10*tracek(i,j,k) + 18*tracek(i,j+dj,k)&
                    - 6*tracek(i,j+2*dj,k) + tracek(i,j+3*dj,k)) / dy12
        d1_f(3) =dk*(-3*tracek(i,j,k-dk) - 10*tracek(i,j,k) + 18*tracek(i,j,k+dk)&
                    - 6*tracek(i,j,k+2*dk) + tracek(i,j,k+3*dk)) / dz12
        ad1_trk = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        ! ad1_aa(3,3)
        d1_f(1) = di * ( -3*axx(i-di,j,k) - 10*axx(i,j,k) + 18*axx(i+di,j,k)   &
                        - 6*axx(i+2*di,j,k) + axx(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*axx(i,j-dj,k) - 10*axx(i,j,k) + 18*axx(i,j+dj,k)   &
                        - 6*axx(i,j+2*dj,k) + axx(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*axx(i,j,k-dk) - 10*axx(i,j,k) + 18*axx(i,j,k+dk)   &
                        - 6*axx(i,j,k+2*dk) + axx(i,j,k+3*dk)) / dz12
        ad1_aa(1,1) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * ( -3*axy(i-di,j,k) - 10*axy(i,j,k) + 18*axy(i+di,j,k)   &
                        - 6*axy(i+2*di,j,k) + axy(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*axy(i,j-dj,k) - 10*axy(i,j,k) + 18*axy(i,j+dj,k)   &
                        - 6*axy(i,j+2*dj,k) + axy(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*axy(i,j,k-dk) - 10*axy(i,j,k) + 18*axy(i,j,k+dk)   &
                        - 6*axy(i,j,k+2*dk) + axy(i,j,k+3*dk)) / dz12
        ad1_aa(1,2) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * ( -3*axz(i-di,j,k) - 10*axz(i,j,k) + 18*axz(i+di,j,k)   &
                        - 6*axz(i+2*di,j,k) + axz(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*axz(i,j-dj,k) - 10*axz(i,j,k) + 18*axz(i,j+dj,k)   &
                        - 6*axz(i,j+2*dj,k) + axz(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*axz(i,j,k-dk) - 10*axz(i,j,k) + 18*axz(i,j,k+dk)   &
                        - 6*axz(i,j,k+2*dk) + axz(i,j,k+3*dk)) / dz12
        ad1_aa(1,3) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * ( -3*ayy(i-di,j,k) - 10*ayy(i,j,k) + 18*ayy(i+di,j,k)   &
                        - 6*ayy(i+2*di,j,k) + ayy(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*ayy(i,j-dj,k) - 10*ayy(i,j,k) + 18*ayy(i,j+dj,k)   &
                        - 6*ayy(i,j+2*dj,k) + ayy(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*ayy(i,j,k-dk) - 10*ayy(i,j,k) + 18*ayy(i,j,k+dk)   &
                        - 6*ayy(i,j,k+2*dk) + ayy(i,j,k+3*dk)) / dz12
        ad1_aa(2,2) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * ( -3*ayz(i-di,j,k) - 10*ayz(i,j,k) + 18*ayz(i+di,j,k)   &
                        - 6*ayz(i+2*di,j,k) + ayz(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*ayz(i,j-dj,k) - 10*ayz(i,j,k) + 18*ayz(i,j+dj,k)   &
                        - 6*ayz(i,j+2*dj,k) + ayz(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*ayz(i,j,k-dk) - 10*ayz(i,j,k) + 18*ayz(i,j,k+dk)   &
                        - 6*ayz(i,j,k+2*dk) + ayz(i,j,k+3*dk)) / dz12
        ad1_aa(2,3) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) = di * ( -3*azz(i-di,j,k) - 10*azz(i,j,k) + 18*azz(i+di,j,k)   &
                        - 6*azz(i+2*di,j,k) + azz(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3*azz(i,j-dj,k) - 10*azz(i,j,k) + 18*azz(i,j+dj,k)   &
                        - 6*azz(i,j+2*dj,k) + azz(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3*azz(i,j,k-dk) - 10*azz(i,j,k) + 18*azz(i,j,k+dk)   &
                        - 6*azz(i,j,k+2*dk) + azz(i,j,k+3*dk)) / dz12
        ad1_aa(3,3) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        ad1_aa(2,1) = ad1_aa(1,2)
        ad1_aa(3,1) = ad1_aa(1,3)
        ad1_aa(3,2) = ad1_aa(2,3)

        ! ad1_gammat(3)
        d1_f(1) =di*(-3*gammatx(i-di,j,k) - 10*gammatx(i,j,k) + 18*gammatx(i+di,j,k)&
                    - 6*gammatx(i+2*di,j,k) + gammatx(i+3*di,j,k)) / dx12
        d1_f(2) =dj*(-3*gammatx(i,j-dj,k) - 10*gammatx(i,j,k) + 18*gammatx(i,j+dj,k)&
                    - 6*gammatx(i,j+2*dj,k) + gammatx(i,j+3*dj,k)) / dy12
        d1_f(3) =dk*(-3*gammatx(i,j,k-dk) - 10*gammatx(i,j,k) + 18*gammatx(i,j,k+dk)&
                    - 6*gammatx(i,j,k+2*dk) + gammatx(i,j,k+3*dk)) / dz12
        ad1_gammat(1) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) =di*(-3*gammaty(i-di,j,k) - 10*gammaty(i,j,k) + 18*gammaty(i+di,j,k)&
                    - 6*gammaty(i+2*di,j,k) + gammaty(i+3*di,j,k)) / dx12
        d1_f(2) =dj*(-3*gammaty(i,j-dj,k) - 10*gammaty(i,j,k) + 18*gammaty(i,j+dj,k)&
                    - 6*gammaty(i,j+2*dj,k) + gammaty(i,j+3*dj,k)) / dy12
        d1_f(3) =dk*(-3*gammaty(i,j,k-dk) - 10*gammaty(i,j,k) + 18*gammaty(i,j,k+dk)&
                    - 6*gammaty(i,j,k+2*dk) + gammaty(i,j,k+3*dk)) / dz12
        ad1_gammat(2) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        d1_f(1) =di*(-3*gammatz(i-di,j,k) - 10*gammatz(i,j,k) + 18*gammatz(i+di,j,k)&
                    - 6*gammatz(i+2*di,j,k) + gammatz(i+3*di,j,k)) / dx12
        d1_f(2) =dj*(-3*gammatz(i,j-dj,k) - 10*gammatz(i,j,k) + 18*gammatz(i,j+dj,k)&
                    - 6*gammatz(i,j+2*dj,k) + gammatz(i,j+3*dj,k)) / dy12
        d1_f(3) =dk*(-3*gammatz(i,j,k-dk) - 10*gammatz(i,j,k) + 18*gammatz(i,j,k+dk)&
                    - 6*gammatz(i,j,k+2*dk) + gammatz(i,j,k+3*dk)) / dz12
        ad1_gammat(3) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

        ! ad1_alph
        if (evolve_alp) then
           d1_f(1) = di * ( -3*alp(i-di,j,k) - 10*alp(i,j,k) + 18*alp(i+di,j,k)   &
                - 6*alp(i+2*di,j,k) + alp(i+3*di,j,k)) / dx12
           d1_f(2) = dj * ( -3*alp(i,j-dj,k) - 10*alp(i,j,k) + 18*alp(i,j+dj,k)   &
                - 6*alp(i,j+2*dj,k) + alp(i,j+3*dj,k)) / dy12
           d1_f(3) = dk * ( -3*alp(i,j,k-dk) - 10*alp(i,j,k) + 18*alp(i,j,k+dk)   &
                - 6*alp(i,j,k+2*dk) + alp(i,j,k+3*dk)) / dz12
           ad1_alph = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)
        end if

        ! ad1_beta(3)
        if (evolve_beta) then
           d1_f(1) = di * ( -3*betax(i-di,j,k) - 10*betax(i,j,k) + 18*betax(i+di,j,k)   &
                - 6*betax(i+2*di,j,k) + betax(i+3*di,j,k)) / dx12
           d1_f(2) = dj * ( -3*betax(i,j-dj,k) - 10*betax(i,j,k) + 18*betax(i,j+dj,k)   &
                - 6*betax(i,j+2*dj,k) + betax(i,j+3*dj,k)) / dy12
           d1_f(3) = dk * ( -3*betax(i,j,k-dk) - 10*betax(i,j,k) + 18*betax(i,j,k+dk)   &
                - 6*betax(i,j,k+2*dk) + betax(i,j,k+3*dk)) / dz12
           ad1_beta(1) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

           d1_f(1) = di * ( -3*betay(i-di,j,k) - 10*betay(i,j,k) + 18*betay(i+di,j,k)   &
                - 6*betay(i+2*di,j,k) + betay(i+3*di,j,k)) / dx12
           d1_f(2) = dj * ( -3*betay(i,j-dj,k) - 10*betay(i,j,k) + 18*betay(i,j+dj,k)   &
                - 6*betay(i,j+2*dj,k) + betay(i,j+3*dj,k)) / dy12
           d1_f(3) = dk * ( -3*betay(i,j,k-dk) - 10*betay(i,j,k) + 18*betay(i,j,k+dk)   &
                - 6*betay(i,j,k+2*dk) + betay(i,j,k+3*dk)) / dz12
           ad1_beta(2) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)

           d1_f(1) = di * ( -3*betaz(i-di,j,k) - 10*betaz(i,j,k) + 18*betaz(i+di,j,k)   &
                - 6*betaz(i+2*di,j,k) + betaz(i+3*di,j,k)) / dx12
           d1_f(2) = dj * ( -3*betaz(i,j-dj,k) - 10*betaz(i,j,k) + 18*betaz(i,j+dj,k)   &
                - 6*betaz(i,j+2*dj,k) + betaz(i,j+3*dj,k)) / dy12
           d1_f(3) = dk * ( -3*betaz(i,j,k-dk) - 10*betaz(i,j,k) + 18*betaz(i,j,k+dk)   &
                - 6*betaz(i,j,k+2*dk) + betaz(i,j,k+3*dk)) / dz12
           ad1_beta(3) = beta_l(1)*d1_f(1) + beta_l(2)*d1_f(2) + beta_l(3)*d1_f(3)
        end if
      else

        ! ad1_ww
        ad1_ww  = beta_l(1)*d1_ww(1) + beta_l(2)*d1_ww(2) + beta_l(3)*d1_ww(3)
        
        if (evolve_scalar) then
        ! ad1_phi
        ad1_lphi = beta_l(1)*d1_lphi(1) + beta_l(2)*d1_lphi(2) + beta_l(3)*d1_lphi(3)

        ! ad1_Kphi
        ad1_lKphi = beta_l(1)*d1_lKphi(1) + beta_l(2)*d1_lKphi(2) + beta_l(3)*d1_lKphi(3)
        end if

        ! ad1_hh(3,3)
        ad1_hh(1,1) = beta_l(1)*d1_hh11(1) + beta_l(2)*d1_hh11(2) + beta_l(3)*d1_hh11(3)
        ad1_hh(1,2) = beta_l(1)*d1_hh12(1) + beta_l(2)*d1_hh12(2) + beta_l(3)*d1_hh12(3)
        ad1_hh(1,3) = beta_l(1)*d1_hh13(1) + beta_l(2)*d1_hh13(2) + beta_l(3)*d1_hh13(3)
        ad1_hh(2,2) = beta_l(1)*d1_hh22(1) + beta_l(2)*d1_hh22(2) + beta_l(3)*d1_hh22(3)
        ad1_hh(2,3) = beta_l(1)*d1_hh23(1) + beta_l(2)*d1_hh23(2) + beta_l(3)*d1_hh23(3)
        ad1_hh(3,3) = beta_l(1)*d1_hh33(1) + beta_l(2)*d1_hh33(2) + beta_l(3)*d1_hh33(3)

        ad1_hh(2,1) = ad1_hh(1,2)
        ad1_hh(3,1) = ad1_hh(1,3)
        ad1_hh(3,2) = ad1_hh(2,3)

        ! ad1_trk
        ad1_trk = beta_l(1)*d1_trk(1) + beta_l(2)*d1_trk(2) + beta_l(3)*d1_trk(3)

        ! ad1_aa(3,3)
        ad1_aa(1,1) = beta_l(1)*d1_aa11(1) + beta_l(2)*d1_aa11(2) + beta_l(3)*d1_aa11(3)
        ad1_aa(1,2) = beta_l(1)*d1_aa12(1) + beta_l(2)*d1_aa12(2) + beta_l(3)*d1_aa12(3)
        ad1_aa(1,3) = beta_l(1)*d1_aa13(1) + beta_l(2)*d1_aa13(2) + beta_l(3)*d1_aa13(3)
        ad1_aa(2,2) = beta_l(1)*d1_aa22(1) + beta_l(2)*d1_aa22(2) + beta_l(3)*d1_aa22(3)
        ad1_aa(2,3) = beta_l(1)*d1_aa23(1) + beta_l(2)*d1_aa23(2) + beta_l(3)*d1_aa23(3)
        ad1_aa(3,3) = beta_l(1)*d1_aa33(1) + beta_l(2)*d1_aa33(2) + beta_l(3)*d1_aa33(3)

        ad1_aa(2,1) = ad1_aa(1,2)
        ad1_aa(3,1) = ad1_aa(1,3)
        ad1_aa(3,2) = ad1_aa(2,3)

        ! ad1_gammat(3)
        ad1_gammat(1) = beta_l(1)*d1_gammat1(1) + beta_l(2)*d1_gammat1(2) + beta_l(3)*d1_gammat1(3)
        ad1_gammat(2) = beta_l(1)*d1_gammat2(1) + beta_l(2)*d1_gammat2(2) + beta_l(3)*d1_gammat2(3)
        ad1_gammat(3) = beta_l(1)*d1_gammat3(1) + beta_l(2)*d1_gammat3(2) + beta_l(3)*d1_gammat3(3)

        ! ad1_alph
        if (evolve_alp) then
           ad1_alph = beta_l(1)*d1_alph(1) + beta_l(2)*d1_alph(2) + beta_l(3)*d1_alph(3)
        end if

        ! ad1_beta(3)
        if (evolve_beta) then
           ad1_beta(1) = beta_l(1)*d1_beta1(1) + beta_l(2)*d1_beta1(2) + beta_l(3)*d1_beta1(3)
           ad1_beta(2) = beta_l(1)*d1_beta2(1) + beta_l(2)*d1_beta2(2) + beta_l(3)*d1_beta2(3)
           ad1_beta(3) = beta_l(1)*d1_beta3(1) + beta_l(2)*d1_beta3(2) + beta_l(3)*d1_beta3(3)
        end if

        !if( abs(y(i,j,k)) < 1.0d-05 .and. abs(z(i,j,k)) < 1.0d-05 ) then
        !  write(*,*) 'i, j, k, x = ', i, j, k, x(i,j,k)
        !  write(*,*) 'ad1_ww   = ', ad1_ww
        !  write(*,*) 'ad1_hh   = ', ad1_hh
        !  write(*,*) 'ad1_trk  = ', ad1_trk
        !  write(*,*) 'ad1_aa   = ', ad1_aa
        !  write(*,*) 'ad1_alph = ', ad1_alph
        !  write(*,*) 'ad1_gam  = ', ad1_gammat
        !  call flush(6)
        !end if

      end if ! use_advection_stencils /= 0

    end if

    !-------------------------------------------
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

    call LeanBSSN_apply_jacobian2(d1_alph, d2_alph, jac, hes)
    call LeanBSSN_apply_jacobian2(d1_ww, d2_ww, jac, hes)
    call LeanBSSN_apply_jacobian2(d1_beta1, d2_beta1, jac, hes)
    call LeanBSSN_apply_jacobian2(d1_beta2, d2_beta2, jac, hes)
    call LeanBSSN_apply_jacobian2(d1_beta3, d2_beta3, jac, hes)
    call LeanBSSN_apply_jacobian2(d1_hh11, d2_hh11, jac, hes)
    call LeanBSSN_apply_jacobian2(d1_hh12, d2_hh12, jac, hes)
    call LeanBSSN_apply_jacobian2(d1_hh13, d2_hh13, jac, hes)
    call LeanBSSN_apply_jacobian2(d1_hh22, d2_hh22, jac, hes)
    call LeanBSSN_apply_jacobian2(d1_hh23, d2_hh23, jac, hes)
    call LeanBSSN_apply_jacobian2(d1_hh33, d2_hh33, jac, hes)
  end if

    d1_beta(1,:) = d1_beta1(:)
    d1_beta(2,:) = d1_beta2(:)
    d1_beta(3,:) = d1_beta3(:)

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

    d2_beta(1,:,:) = d2_beta1(:,:)
    d2_beta(2,:,:) = d2_beta2(:,:)
    d2_beta(3,:,:) = d2_beta3(:,:)

    d2_hh(1,1,:,:) = d2_hh11(:,:)
    d2_hh(1,2,:,:) = d2_hh12(:,:)
    d2_hh(1,3,:,:) = d2_hh13(:,:)
    d2_hh(2,2,:,:) = d2_hh22(:,:)
    d2_hh(2,3,:,:) = d2_hh23(:,:)
    d2_hh(3,3,:,:) = d2_hh33(:,:)
    d2_hh(2,1,:,:) = d2_hh(1,2,:,:)
    d2_hh(3,1,:,:) = d2_hh(1,3,:,:)
    d2_hh(3,2,:,:) = d2_hh(2,3,:,:)

    !-------------------------------------------


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
    cd2_ww   = d2_ww
    cd2_alph = d2_alph
    if (evolve_scalar) then 
       cd2_lphi = d2_lphi
    end if
    do a = 1, 3
      do b = a, 3
        do m = 1, 3
          cd2_ww(a,b)   = cd2_ww(a,b)   - cf2(m,a,b) * d1_ww(m)
          cd2_alph(a,b) = cd2_alph(a,b) - cf2(m,a,b) * d1_alph(m)
          if (evolve_scalar) then 
                 cd2_lphi(a,b) = cd2_lphi(a,b) - cf2(m,a,b) * d1_lphi(m)
          end if
        end do
      end do
    end do
    cd2_ww(2,1)   = cd2_ww(1,2)
    cd2_ww(3,1)   = cd2_ww(1,3)
    cd2_ww(3,2)   = cd2_ww(2,3)
    cd2_alph(2,1) = cd2_alph(1,2)
    cd2_alph(3,1) = cd2_alph(1,3)
    cd2_alph(3,2) = cd2_alph(2,3)
    if (evolve_scalar) then 
            cd2_lphi(2,1) = cd2_lphi(1,2)
            cd2_lphi(3,1) = cd2_lphi(1,3)
            cd2_lphi(3,2) = cd2_lphi(2,3)
    end if

    !---Calculate scalar field cov derivative W factor
    if (evolve_scalar) then 
            aux = 0
            do a=1,3
               do b=1,3
                  cW_dphi(a,b) = d1_lphi(a)*d1_ww(b) + d1_lphi(b)*d1_ww(a) 
                  aux = aux + hu(a,b)*d1_lphi(a)*d1_ww(b)
               end do
            end do 

            cW_dphi = (cW_dphi - aux*hh)/ww
     end if 

    !-------------------------------------------
    !------------ Scalar field terms ----------
    ! First we transform from phi -> varphi 
    ! For DEF
    ! Bphi      = exp(varphi*varphi/2)
    ! BKphi     = Bphi * varphi * Kphi
    ! d1_Bphi   = Bphi * varphi * d1_phi
    ! cd2_Bphi   = Bphi * (1 + varphi*varphi) * d1_varphi * d1_varphi + varphi * Bphi * covD_phi
    ! d1_BKphi  = Bphi*Kphi * (1 + varphi*varphi) * d1_varphi + varphi * Bphi * d1_Kphi
    ! omega     = 2/(B * varphi) -3/2 
    ! d1_omega  = -4/(B*Bphi * varphi^4)
    !
    ! For BD 
    ! Bphi      = exp(varphi)
    ! BKphi     = Bphi*Kphi
    ! d1_Bphi   = Bphi* d_phi
    ! cd2_Bphi  = Bphi * (d_phi * d_phi + covd_phi)
    ! d1_BKphi  = Bphi * ( Kphi * d_phi + d_kphi)
    ! omega     = 1/(2*k0**2) -3/2 
    ! d1_omega  = 0   ! omega is constant
    !

    if (evolve_scalar) then      ! Begin evolveScalar
       if (CCTK_EQUALS(theory,"BDold") .OR. CCTK_EQUALS(theory,"BD")  .OR.  CCTK_EQUALS(theory,"onlymetric") ) then                       !---define BD
            Bphi        = exp(lphi)
            BKphi       = Bphi * lKphi
            omega_Bphi  = 1.0d0/(2.0d0*k0BD*k0BD)-3.0d0/2.0d0
            domega_Bphi = 0.0d0 
            
       end if                                                   ! -- end  Define BD
       
       if (CCTK_EQUALS(theory,"DEF") .OR. CCTK_EQUALS(theory,"DEFold") .OR. CCTK_EQUALS(theory,"DEFdecoupling") .OR. CCTK_EQUALS(theory,"onlySF")) then                      !--- Define DEF
               Bphi        = exp(lphi2/2.0d0)
               BKphi       = Bphi*lphi*lKphi
               omega_Bphi  = 2.0d0/(B_DEF*lphi2)-3.0d0/2.0d0
               domega_Bphi = -4.0d0/(B_DEF*lphi2*lphi2*Bphi)               
       end if                                                   ! -- end define DEF

            !---- Calculate traces             

                tr_dalp_dphi = 0
                tr_cd2_phi   = 0
                tr_cd2_phi_new = 0 
                tr_dww_dphi  = 0
                tr_dphi_dphi = 0
                do a =1,3
                   do b = 1,3
                        tr_dalp_dphi = tr_dalp_dphi + hu(a,b) * d1_alph(a)*d1_lphi(b)
                        tr_cd2_phi   = tr_cd2_phi   + hu(a,b) * cd2_lphi(a,b)
                        tr_cd2_phi_new   = tr_cd2_phi_new   + hu(a,b) * (cd2_lphi(a,b)+cW_dphi(a,b) )*ww*ww
                        tr_dww_dphi  = tr_dww_dphi  + hu(a,b) * d1_ww(a) * d1_lphi(b)
                        tr_dphi_dphi = tr_dphi_dphi + hu(a,b) * d1_lphi(a)*d1_lphi(b)
                    end do
                end do

            if (r(i,j,k)==r_debug .AND. show_debug/=0) THEN
            write(*,*) "---------------------------" 
            write(*,*) "---------------------------"
            write(*,*) " r= ", r(i,j,k) 
            write(*,*) "Bphi ", Bphi 
            write(*,*) "BKphi ", BKphi
            write(*,*) "d1_lphi ", d1_lphi
            write(*,*) "d1_lKphi ", d1_lKphi
            write(*,*) "d2_lphi ", d2_lphi
            write(*,*) "phi ", lphi
            write(*,*) "kphi ", lKphi
            end if

    end if   !--- End if evolve_scalar

    !------------ Ricci Tensor -----------------
    ! Note: we implement W^2 R_{ij}
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
    ! Note: we implement W^2 R_{ij}
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
    ! Note: we implement W^2 R_{ij}
    c_ri_hh = ww*ww * (ri_1 + ri_2 + ri_3)
    c_ri    = c_ri_ww + c_ri_hh

    c_ri(2,1) = c_ri(1,2)
    c_ri(3,1) = c_ri(1,3)
    c_ri(3,2) = c_ri(2,3)
    !-------------------------------------------


    !------------ DD_alpha ---------------------
    ! Note: we implement W^2 D_{a}D_{b}\alpha
    c_ll  = 0
    aux = 0
    do a = 1, 3
      do b = 1, 3
        c_ll(a,b) = d1_alph(a) * d1_ww(b) + d1_alph(b) * d1_ww(a)
        aux       = aux + hu(a,b) * d1_alph(a) * d1_ww(b)
      end do
    end do
    c_ll = ww * ( c_ll - hh * aux ) + ww*ww * cd2_alph
    !-------------------------------------------

    !------------ Advection and Twist terms ----
    divbeta = 0
    do m = 1, 3
      divbeta = divbeta + d1_beta(m,m)
    end do

    ! rhs_ww
    rhs_ww = ad1_ww - ww * divbeta / 3

    ! rhs_hh
    hhdbeta = 0
    do a = 1, 3
      do b = 1, 3
        do m = 1, 3
          hhdbeta(a,b) = hhdbeta(a,b) + hh(a,m) * d1_beta(m,b)
        end do
      end do
    end do

    rhs_hh = ad1_hh
    do a = 1, 3
      do b = 1, 3
        rhs_hh(a,b) = rhs_hh(a,b) + hhdbeta(a,b) + hhdbeta(b,a)                &
                      - 2 * hh(a,b) * divbeta / 3
      end do
    end do

    ! rhs_trk
    rhs_trk = ad1_trk

    ! rhs_aa
    aadbeta = 0
    do a = 1, 3
      do b = 1, 3
        do m = 1, 3
          aadbeta(a,b) = aadbeta(a,b) + aa(a,m) * d1_beta(m,b)
        end do
      end do
    end do

    rhs_aa = ad1_aa
    do a = 1, 3
      do b = 1, 3
        rhs_aa(a,b) = rhs_aa(a,b) + aadbeta(a,b) + aadbeta(b,a)                &
                      - 2 * aa(a,b) * divbeta / 3
      end do
    end do

    ! rhs_gammat
    rhs_gammat = ad1_gammat + 2 * gammat * divbeta / 3

    gamcon = gammat
    do a = 1, 3
      do m = 1, 3
        do n = 1, 3
          gamcon(a) = gamcon(a) - hu(m,n) * cf2(a,m,n)
        end do
      end do
    end do

    do a = 1, 3
      do m = 1, 3
        rhs_gammat(a) = rhs_gammat(a) - gammat(m) * d1_beta(a,m)              &
                       - (chi_gamma + 2.0d0/3.0d0) * gamcon(a) * d1_beta(m,m)
      end do
    end do

    !-------------------------------------------


    !------------ Source terms -----------------
    ! rhs_ww
    rhs_ww = rhs_ww + alph * ww * trk / 3

    ! rhs_hh
    rhs_hh = rhs_hh - 2 * alph * aa

    ! rhs_trk
    tr_ll = 0
    sq_aa = 0
    a2    = 0
    do m = 1, 3
      do n = 1, 3
        tr_ll = tr_ll + hu(m,n) * c_ll(m,n)
        do p = 1, 3
          do q = 1, 3
            a2(m,n) = a2(m,n) + hu(p,q) * aa(m,p) * aa(n,q)
          end do
        end do
        sq_aa = sq_aa + hu(m,n) * a2(m,n)
      end do
    end do

    ! Note: tr_ll is now $D^i D_i \alpha$
    ! Note: we implemented W^2 D_{a}D_{b}\alpha
    rhs_trk = rhs_trk - tr_ll + alph * (sq_aa + trk*trk / 3)

    ! rhs_aa
    ! Note: we implemented W^2 R_{ij}
    trr = 0
    do l = 1, 3
      do m = 1, 3
        trr = trr + hu(l,m) * c_ri(l,m)
      end do
    end do

    ! tr_ll = already calculated for rhs_trk
    ! sq_aa = already calculated for rhs_trk
    ! a2    = already calculated at rhs_trk
    tf_c_ll = c_ll - hh * tr_ll / 3    ! Note: hh = ww2 * gg
    tf_c_ri = c_ri - hh * trr / 3      ! Again hh = ww2 * gg
    rhs_aa = rhs_aa + (-tf_c_ll + alph * tf_c_ri)               &
             + alph * (trk * aa - 2 * a2)

    ! rhs_gammat
    au = 0
    do a = 1, 3
      do b = a, 3
        do m = 1, 3
          do n = 1, 3
            au(a,b) = au(a,b) + hu(a,m) * hu(b,n) * aa(m,n)
          end do
        end do
      end do
    end do
    au(2,1) = au(1,2)
    au(3,1) = au(1,3)
    au(3,2) = au(2,3)

    do a = 1, 3
      do m = 1, 3
        rhs_gammat(a) = rhs_gammat(a) - 4 * alph * hu(a,m) * d1_trk(m) / 3   &
                       - 2 * au(a,m) * ( d1_alph(m) + 3 * alph * d1_ww(m) / ww )
        do n = 1, 3
          rhs_gammat(a) = rhs_gammat(a) + 2 * alph * cf2(a,m,n) * au(m,n)    &
                         + hu(a,m) * d2_beta(n,n,m) / 3                    &
                         + hu(m,n) * d2_beta(a,m,n)
        end do
      end do
    end do

    !-------------------------------------------


    !------------ Matter terms -----------------
    !
    ! n_mu = (-alph, 0, 0, 0)
    ! n^mu = (1, -betax, -betay, -betaz)/alph
    !
    ! E   = n^mu n^nu T_{mu nu}
    !     = (T_{00} - 2 beta^i T_{i0} + beta^i beta^j T_{ij})/alph^2
    !
    ! j_a = -h_a^mu n^nu T_{mu nu}
    !     = -(T_{a 0} - beta^j T_{a j})/alph
    !
    ! S_{a b} = h_{a mu} h_{b nu} T^{mu nu} = T_{a b}

    ! stress-energy tensor variables
    Tab = 0
    if (stress_energy_state /= 0) then   !---- Start matter
               
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

               srcji = 0
               do a = 1, 3
                  do m = 1, 3
                     srcji(a) = srcji(a) + hu(a,m) * srcjdi(m)
                  end do
               end do


               do a = 1, 3
                  do b = 1, 3
                     srcSij(a,b) = Tab(a,b)
                  end do
               end do

               srcS_ww2 = 0
               do m = 1, 3
                  do n = 1, 3
             ! Contracting with conformal metric hu
                     srcS_ww2 = srcS_ww2 + hu(m,n) * srcSij(m,n)
                  end do
               end do
               srcSijTF = srcSij - srcS_ww2 * hh / 3
       
               src_trT = srcS_ww2 * ww*ww - srcE
               
                    if (r(i,j,k)==r_debug ) THEN
                    write(*,*) "------ Matter terms----------------"
                    write(*,*) "r = ", r_debug
                    write(*,*) "Tab", Tab 
                    write(*,*) "srcE", srcE
                    write(*,*) "srcS_ww2 ", srcS_ww2 
                    write(*,*) "srcSijTF ", srcSijTF
                    write(*,*) "srcji ", srcji          
                    write(*,*) "src_trT", src_trT
                    end if 

               if(JordanFrame .OR. CCTK_EQUALS(theory,"onlymetric")) then  ! begin evolve Jordan
                       !!! Normalize matter sources with the SF
                       srcE     = srcE/Bphi
                       srcS_ww2 = srcS_ww2/Bphi
                       srcSijTF = srcSijTF/Bphi
                       srcji    = srcji/Bphi
                   
                    if(CCTK_EQUALS(theory,"BD") .OR. CCTK_EQUALS(theory,"BDold") .OR. CCTK_EQUALS(theory,"onlymetric")) then   ! Init BD
                        F_phi = k0BD*k0BD*(8*pi*src_trT)/Bphi                        

                        trk_omega = omega_Bphi*lKphi*lKphi 
                        trk_phi   = ww*ww*(tr_cd2_phi+tr_dphi_dphi) - ww*tr_dww_dphi - trk*lKphi -1.5d0*F_phi                        
                    
                        gammat_phi   = 0.0d0
                        gammat_omega = 0.0d0
                        do a=1,3
                          do b=1,3
                            aa_omega(a,b)   = d1_lphi(a)*d1_lphi(b) - tr_dphi_dphi * hh(a,b)/3.0d0 
                            aa_phi(a,b)     = ww*ww*(d1_lphi(a)*d1_lphi(b) + cd2_lphi(a,b) + cW_dphi(a,b)) - aa(a,b)*lKphi 
                            aa_phi(a,b)     = aa_phi(a,b) - ww*hh(a,b)*(ww*tr_dphi_dphi + ww*tr_cd2_phi - tr_dww_dphi)/3.0d0   ! take trace
                            gammat_omega(a) = gammat_omega(a) + omega_Bphi * lKphi * hu(a,b)*d1_lphi(b) 
                            gammat_phi(a)   = gammat_phi(a) + hu(a,b)*(d1_lKphi(b)+lKphi*d1_lphi(b) - trk * d1_lphi(b)/3.0d0)  - au(a,b)*d1_lphi(b)
                           end do
                         end do

                 end if ! END BD
                 

                    if (CCTK_EQUALS(theory,"DEFold") .OR. CCTK_EQUALS(theory,"DEF")) then   ! init DEF
                         F_phi = 2.0d0*pi*B_DEF*src_trT*lphi2/Bphi  +  ww*ww* tr_dphi_dphi - lKphi*lKphi
                         
                         trk_omega = omega_Bphi*lKphi*lKphi*lphi2 
                         trk_phi   = ww*ww*(lphi*tr_cd2_phi+(1.0d0+lphi2)*tr_dphi_dphi) - lphi*ww*tr_dww_dphi - trk*lphi*lKphi -1.5d0*F_phi                        


                         gammat_phi   = 0.0d0
                         gammat_omega = 0.0d0
                         do a=1,3
                           do b=1,3
                            aa_omega(a,b)   = lphi2*d1_lphi(a)*d1_lphi(b) - lphi2*tr_dphi_dphi * hh(a,b)/3.0d0 
                            aa_phi(a,b)     = ww*ww*( d1_lphi(a)*d1_lphi(b)*(lphi2+1.0d0) + lphi *(cd2_lphi(a,b) + cW_dphi(a,b)) ) - aa(a,b)*lKphi*lphi 
                            aa_phi(a,b)     = aa_phi(a,b) - ww*hh(a,b)*( ww*(lphi2+1.0d0)*tr_dphi_dphi  + lphi*(ww*tr_cd2_phi - tr_dww_dphi) )/3.0d0   ! take trace
                                

                            !Need to check this next step
                            !gammat_omega(a) = gammat_omega(a) + omega_Bphi * lKphi * hu(a,b)*d1_lphi(b)*lphi2                           
                            gammat_omega(a) = gammat_omega(a) + (1.0d0+2.0/B_DEF-lphi2/2.0d0) * lKphi * hu(a,b)*d1_lphi(b) 
                            gammat_phi(a)   = gammat_phi(a) + hu(a,b)*(lphi*d1_lKphi(b)+lKphi*(lphi2+1.0d0)*d1_lphi(b) - trk * d1_lphi(b)*lphi/3.0d0)  -  au(a,b)*d1_lphi(b)*lphi
                            end do
                          end do


                    end if  !end DEF

                end if                   !--- End JordanFranme 


                    if (r(i,j,k)==r_debug) THEN
                    write(*,*) "---------------------------" 
                    write(*,*) "r = ", r(i,j,k)
                    if(JordanFrame) then  ! begin evolve Jordan
                            write(*,*) "trk_omega ", trk_omega 
                            write(*,*) "trk_phi ", trk_phi
                            write(*,*) "F_phi", F_phi
                            write(*,*) "aa_omega ", aa_omega
                            write(*,*) "aa_phi ", aa_phi
                            write(*,*) "gamma_omega ", gammat_omega                           
                            write(*,*) "gamma_phi ", gammat_phi 
                            write(*,*) "------ Matter terms----------------" 
                            write(*,*) "srcE", srcE
                            write(*,*) "srcS_ww2 ", srcS_ww2 
                            write(*,*) "srcSijTF ", srcSijTF
                            write(*,*) "srcji ", srcji          
                            write(*,*) "src_trT", src_trT 

                        end if 
                    end if     
!                    !
       !------------ Correct source terms ---------
       rhs_trk    = rhs_trk    + pi4  * alph * (srcE + ww*ww * srcS_ww2) !+ alph*(trk_omega + trk_phi)
       rhs_aa     = rhs_aa     - pi8  * alph * ww*ww * srcSijTF
       rhs_gammat = rhs_gammat - pi16 * alph * srcji

      if(JordanFrame .OR. CCTK_EQUALS(theory,"onlymetric")) then
       rhs_trk    = rhs_trk + alph*(  trk_omega + trk_phi)
       rhs_aa     = rhs_aa - alph * (ww*ww*aa_omega*Omega_Bphi + aa_phi )
       rhs_gammat = rhs_gammat - 2.0d0*alph*(gammat_omega + gammat_phi)
       end if
    end if  ! ---- end matter         


            if (r(i,j,k)==r_debug .AND. show_debug/=0) THEN
                    write(*,*) "----------  RHS --------------"
                    write(*,*) "rhs aa ", rhs_aa 
                    write(*,*) "rhs trk ", rhs_trk
                    write(*,*) "rhs gamma ", rhs_gammat
            end if            
    !------------ Write to grid functions ------
    rhs_conf_fac(i,j,k) = rhs_ww


    rhs_hxx(i,j,k) = rhs_hh(1,1)
    rhs_hxy(i,j,k) = rhs_hh(1,2)
    rhs_hxz(i,j,k) = rhs_hh(1,3)
    rhs_hyy(i,j,k) = rhs_hh(2,2)
    rhs_hyz(i,j,k) = rhs_hh(2,3)
    rhs_hzz(i,j,k) = rhs_hh(3,3)

    rhs_tracek(i,j,k) = rhs_trk

    rhs_axx(i,j,k) = rhs_aa(1,1)
    rhs_axy(i,j,k) = rhs_aa(1,2)
    rhs_axz(i,j,k) = rhs_aa(1,3)
    rhs_ayy(i,j,k) = rhs_aa(2,2)
    rhs_ayz(i,j,k) = rhs_aa(2,3)
    rhs_azz(i,j,k) = rhs_aa(3,3)

    rhs_gammatx(i,j,k) = rhs_gammat(1)
    rhs_gammaty(i,j,k) = rhs_gammat(2)
    rhs_gammatz(i,j,k) = rhs_gammat(3)

    if(CCTK_EQUALS(theory,"onlymetric")) then        
    rhs_axx(i,j,k) = 0.0
    rhs_axy(i,j,k) = 0.0
    rhs_axz(i,j,k) = 0.0
    rhs_ayy(i,j,k) = 0.0
    rhs_ayz(i,j,k) = 0.0
    rhs_azz(i,j,k) = 0.0
    end if

    if (CCTK_EQUALS(theory,"onlySF") .OR. CCTK_EQUALS(theory,"DEFdecouplingEF")) then
    rhs_tracek(i,j,k) = 0.0

    rhs_axx(i,j,k) = 0.0
    rhs_axy(i,j,k) = 0.0
    rhs_axz(i,j,k) = 0.0
    rhs_ayy(i,j,k) = 0.0
    rhs_ayz(i,j,k) = 0.0
    rhs_azz(i,j,k) = 0.0

    rhs_gammatx(i,j,k) = 0.0
    rhs_gammaty(i,j,k) = 0.0
    rhs_gammatz(i,j,k) = 0.0


    end if
    !-------------------------------------------

    !------------ Now for the lapse -----------
    if (evolve_alp) then
       rhs_alp(i,j,k) = zeta_alpha * ad1_alph - 2.d0 * alph * trk
    end if
    !-------------------------------------------


    ! ----------- And shift -------------------
    if (evolve_beta) then
       ! rhs_beta
       rhs_beta = zeta_beta * ad1_beta

       myeta = eta_beta

       if( eta_transition /= 0 ) then
          r2 = x(i,j,k)**2 + y(i,j,k)**2 + z(i,j,k)**2
          if( r2 < eps_r ) r2 = eps_r
          !write(*,*) 'eta_beta = ', eta_beta
          !call flush(6)
          myeta = eta_beta * eta_transition_r**2 / (r2 + eta_transition_r**2)
          !write(*,*) 'ijk = ', i, j, k
          !write(*,*) r2, eta_transition_r**2, myeta
          !call flush(6)
       end if

       rhs_beta = rhs_beta + beta_Gamma * alph**beta_Alp * gammat - myeta * beta

       rhs_betax(i,j,k) = rhs_beta(1)
       rhs_betay(i,j,k) = rhs_beta(2)
       rhs_betaz(i,j,k) = rhs_beta(3)
    end if
   
   if (evolve_scalar) then    ! Start evolve scalar field


           ! Advective derivatives
           rhs_lphi  = ad1_lphi
           rhs_lKphi = ad1_lKphi

           ! Add source to scalar field
           ! --- Equations in Einstein frame
           if (CCTK_EQUALS(theory,"BDdecouplingEF") .OR. CCTK_EQUALS(theory,"DEFdecouplingEF") ) then
                   rhs_lphi  = rhs_lphi-2.0d0 * alph * lKphi
                   rhs_lKphi = rhs_lKphi -0.5d0 * alph * ww * ( ww*tr_cd2_phi - tr_dww_dphi)
                   rhs_lKphi = rhs_lKphi + alph*trk*lKphi -0.5d0*ww*ww*tr_dalp_dphi                   
                   rhs_phi(i,j,k) = rhs_lphi
                   if (CCTK_EQUALS(theory,"BDdecouplingEF")) then 
                       rhs_lKphi = rhs_lKphi - 2*pi*alph*src_trT*k0BD    
                       rhs_Kphi(i,j,k) = rhs_lKphi
                       if (r(i,j,k)==r_debug .AND. show_debug/=0) THEN
                            write(*,*) "rhs phi ", rhs_lphi 
                            write(*,*) "rhs Kphi ", rhs_lKphi                            
                       end if        
                   end if
                   if (CCTK_EQUALS(theory,"DEFdecouplingEF")) then
                           rhs_Kphi(i,j,k) = rhs_lKphi - 2*pi*alph*src_trT*lphi*betaDEF
                   end if
           end if
           ! ----- End eqs in Einstein frame

           !----- Equations in Jordan frame
           if (CCTK_EQUALS(theory,"decouplingJ")) then 
                   rhs_lphi  = rhs_lphi- alph * lKphi
                   rhs_lKphi = rhs_lKphi - alph * ww * ( ww*tr_cd2_phi - tr_dww_dphi)
                   rhs_lKphi = rhs_lKphi + alph*trk*lKphi -ww*ww*tr_dalp_dphi
        
                   rhs_phi(i,j,k) = rhs_lphi
                   rhs_Kphi(i,j,k) = rhs_lKphi + 8*pi*alph*src_trT*(k0BD*k0BD)
           end if

            if (JordanFrame .OR. CCTK_EQUALS(theory,"onlySF")) then  ! Begin Jordan Frame SF ev
                   rhs_lphi  = rhs_lphi- alph * lKphi
                                  
                    
                    if (CCTK_EQUALS(theory,"BD") .OR. CCTK_EQUALS(theory,"BDold")) then  ! add coupling term if BD
                   
                        rhs_lKphi = rhs_lKphi - alph*tr_cd2_phi_new &
                               + alph*trk*lKphi -ww*ww*tr_dalp_dphi &
                               + alph*(lKphi*lKphi - ww*ww*tr_dphi_dphi) &
                               + 8.0d0*pi*alph*src_trT*(k0BD*k0BD)/Bphi                
                    end if 
              
                    if (CCTK_EQUALS(theory,"DEF") .OR. CCTK_EQUALS(theory,"DEFold") .OR. CCTK_EQUALS(theory,"onlySF")) then  ! add coupling if DEF                                                
                        rhs_lKphi = rhs_lKphi -ww*ww*tr_dalp_dphi &
                               + alph*( -tr_cd2_phi_new + trk*lKphi &                        
                               + lphi*(lKphi*lKphi - ww*ww*tr_dphi_dphi) &
                               + 2.0d0*pi*src_trT*B_DEF*lphi/Bphi )
                    end if 

                        rhs_phi(i,j,k) = rhs_lphi

                        rhs_Kphi(i,j,k) = rhs_lKphi                                        
                
           end if                   !End Jordan Frame SF ev

            if (r(i,j,k)==r_debug .AND. show_debug/=0) THEN
                    write(*,*) "---------- SF RHS --------------"
                    write(*,*) "r =,", r_debug
                    write(*,*) "tr_cd2_phi ", tr_cd2_phi 
                    write(*,*) "tr_cd2_phi_new ", tr_cd2_phi_new
                    write(*,*) "trace term ", -alph*tr_cd2_phi_new
                    write(*,*) "coupling ", 2.0d0*alph*pi*src_trT*B_DEF*lphi/Bphi !8*pi*alph*src_trT*(k0BD*k0BD)/Bphi
                    write(*,*) "tr_dww_dphi", tr_dww_dphi
                    write(*,*) "beta" , beta
                    write(*,*) "beta locl" , beta_l
                    write(*,*) "alpha", alph
                    write(*,*) "ww", ww
                    write(*,*) "trk", trk                    
                    write(*,*) "tr_dalph_dphi", tr_dalp_dphi
                    write(*,*) "ad1_phi", ad1_lphi
                    write(*,*) "ad1_Kphi", ad1_lKphi
                    write(*,*) "tr_dphi_dphi", tr_dphi_dphi
                    write(*,*) "src_trT", src_trT                    
                    write(*,*) "rhs phi ", rhs_lphi 
                    write(*,*) "rhs Kphi ", rhs_lKphi
                    write(*,*) "---------------------"
                    write(*,*) "---------------------"
            end if            
             

    end if       !! end evolve scalar field

  end do
  end do
  end do
  !$OMP END PARALLEL DO

end subroutine LeanBSSN_calc_bssn_rhs
!
!===========================================================================
!
subroutine LeanBSSN_calc_bssn_rhs_bdry( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL, parameter :: one  = 1.0d0
  CCTK_REAL, parameter :: zero = 0.0d0
  CCTK_REAL phi0,kphi0
  CCTK_INT ierr

  ! NewRad_Apply calling syntax is as follows
  ! NewRad_Apply(cctkGH, var, rhs, var0, v0, radpower)
  !
  ! where
  !
  !   var  =  var0 + u(r-v0*t)/r^radpower

  ierr = NewRad_Apply(cctkGH, conf_fac, rhs_conf_fac, one, one, n_conf_fac)

  ierr = NewRad_Apply(cctkGH, hxx, rhs_hxx, one , one, n_hij)
  ierr = NewRad_Apply(cctkGH, hxy, rhs_hxy, zero, one, n_hij)
  ierr = NewRad_Apply(cctkGH, hxz, rhs_hxz, zero, one, n_hij)
  ierr = NewRad_Apply(cctkGH, hyy, rhs_hyy, one , one, n_hij)
  ierr = NewRad_Apply(cctkGH, hyz, rhs_hyz, zero, one, n_hij)
  ierr = NewRad_Apply(cctkGH, hzz, rhs_hzz, one , one, n_hij)

  ierr = NewRad_Apply(cctkGH, tracek, rhs_tracek, zero, one, n_trk)

  ierr = NewRad_Apply(cctkGH, axx, rhs_axx, zero, one, n_aij)
  ierr = NewRad_Apply(cctkGH, axy, rhs_axy, zero, one, n_aij)
  ierr = NewRad_Apply(cctkGH, axz, rhs_axz, zero, one, n_aij)
  ierr = NewRad_Apply(cctkGH, ayy, rhs_ayy, zero, one, n_aij)
  ierr = NewRad_Apply(cctkGH, ayz, rhs_ayz, zero, one, n_aij)
  ierr = NewRad_Apply(cctkGH, azz, rhs_azz, zero, one, n_aij)

  ierr = NewRad_Apply(cctkGH, gammatx, rhs_gammatx, zero, one, n_gammat)
  ierr = NewRad_Apply(cctkGH, gammaty, rhs_gammaty, zero, one, n_gammat)
  ierr = NewRad_Apply(cctkGH, gammatz, rhs_gammatz, zero, one, n_gammat)

  if (CCTK_EQUALS(lapse_evolution_method, "LeanBSSNMoL")) then
     ierr = NewRad_Apply(cctkGH, alp, rhs_alp, one, one, n_alpha)
  end if

  if (CCTK_EQUALS(shift_evolution_method, "LeanBSSNMoL")) then
     ierr = NewRad_Apply(cctkGH, betax, rhs_betax, zero, one, n_beta)
     ierr = NewRad_Apply(cctkGH, betay, rhs_betay, zero, one, n_beta)
     ierr = NewRad_Apply(cctkGH, betaz, rhs_betaz, zero, one, n_beta)
  end if
 
  if (CCTK_EQUALS(scalar_evolution_method, "LeanBSSNMoL")) then
     if (CCTK_EQUALS(theory, "decouplingEF") .OR. CCTK_EQUALS(theory,"DEFdecouplingEF") ) then
             phi0  = phi1_0
             kphi0 = Kphi1_0
     end if
     if (CCTK_EQUALS(theory, "decouplingJ")) then
             phi0  = 1.0
             kphi0 = Kphi1_0
     end if
     if (CCTK_EQUALS(theory, "onlySF") .OR. CCTK_EQUALS(theory,"BDold") .OR. CCTK_EQUALS(theory,"DEFold") .OR. CCTK_EQUALS(theory,"BD") .OR. CCTK_EQUALS(theory,"DEF")) then             
             phi0  = phi_at_inf
             kphi0 = Kphi1_0
     end if
     ierr = NewRad_Apply(cctkGH, phi1, rhs_phi, phi0, one, n_phi1)
     ierr = NewRad_Apply(cctkGH, Kphi1, rhs_Kphi, kphi0, one, n_Kphi1)
  end if
  
end subroutine LeanBSSN_calc_bssn_rhs_bdry
!
!===========================================================================
!
subroutine LeanBSSN_calc_bssn_rhs_bdry_sph( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT  i, j, k
  CCTK_REAL odr2

  CCTK_REAL rr
  CCTK_REAL alph, beta(3)
  CCTK_REAL ww, hh(3,3), trk, aa(3,3), gammat(3)

  CCTK_REAL dr_alph, dr_beta(3)
  CCTK_REAL dr_ww, dr_hh(3,3), dr_trk, dr_aa(3,3), dr_gammat(3)

  CCTK_INT  reflevel, map

  odr2 = 1.0d0 / (2.0d0*CCTK_DELTA_SPACE(3))

  reflevel = GetRefinementLevel(cctkGH)
  map      = MultiPatch_GetMap(cctkGH)

  ! Apply only on the coarsest level and in the spherical shell. Points marked
  ! with cctk_bbox(6) == 0 are inter-processor boundaries, so we also do not
  ! want those
  if (reflevel /= 0 .or. map == 0 .or. cctk_bbox(6) == 0) return
  ! write(*,*) 'map = ', map

  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE( i, j, k, rr, ww, hh, trk, aa, gammat, &
  !$OMP dr_ww, dr_hh, dr_trk, dr_aa, dr_gammat )
  do k = cctk_lsh(3)-cctk_nghostzones(3)+1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           rr        = r(i,j,k)

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

           dr_ww       = (conf_fac(i,j,k-2) - 4*conf_fac(i,j,k-1) + 3*conf_fac(i,j,k))*odr2

           dr_hh(1,1)  = (hxx(i,j,k-2) - 4*hxx(i,j,k-1) + 3*hxx(i,j,k))*odr2
           dr_hh(1,2)  = (hxy(i,j,k-2) - 4*hxy(i,j,k-1) + 3*hxy(i,j,k))*odr2
           dr_hh(1,3)  = (hxz(i,j,k-2) - 4*hxz(i,j,k-1) + 3*hxz(i,j,k))*odr2
           dr_hh(2,2)  = (hyy(i,j,k-2) - 4*hyy(i,j,k-1) + 3*hyy(i,j,k))*odr2
           dr_hh(2,3)  = (hyz(i,j,k-2) - 4*hyz(i,j,k-1) + 3*hyz(i,j,k))*odr2
           dr_hh(3,3)  = (hzz(i,j,k-2) - 4*hzz(i,j,k-1) + 3*hzz(i,j,k))*odr2

           dr_trk  = (tracek(i,j,k-2) - 4*tracek(i,j,k-1) + 3*tracek(i,j,k))*odr2

           dr_aa(1,1)  = (axx(i,j,k-2) - 4*axx(i,j,k-1) + 3*axx(i,j,k))*odr2
           dr_aa(1,2)  = (axy(i,j,k-2) - 4*axy(i,j,k-1) + 3*axy(i,j,k))*odr2
           dr_aa(1,3)  = (axz(i,j,k-2) - 4*axz(i,j,k-1) + 3*axz(i,j,k))*odr2
           dr_aa(2,2)  = (ayy(i,j,k-2) - 4*ayy(i,j,k-1) + 3*ayy(i,j,k))*odr2
           dr_aa(2,3)  = (ayz(i,j,k-2) - 4*ayz(i,j,k-1) + 3*ayz(i,j,k))*odr2
           dr_aa(3,3)  = (azz(i,j,k-2) - 4*azz(i,j,k-1) + 3*azz(i,j,k))*odr2

           dr_gammat(1)  = (gammatx(i,j,k-2) - 4*gammatx(i,j,k-1) + 3*gammatx(i,j,k))*odr2
           dr_gammat(2)  = (gammaty(i,j,k-2) - 4*gammaty(i,j,k-1) + 3*gammaty(i,j,k))*odr2
           dr_gammat(3)  = (gammatz(i,j,k-2) - 4*gammatz(i,j,k-1) + 3*gammatz(i,j,k))*odr2

           ! FIXME: change the wave speeds below

           rhs_conf_fac(i,j,k)  = -dr_ww - (ww - 1.0d0) / rr

           rhs_hxx(i,j,k)  = -dr_hh(1,1) - (hh(1,1) - 1.0d0) / rr
           rhs_hxy(i,j,k)  = -dr_hh(1,2) -  hh(1,2)        / rr
           rhs_hxz(i,j,k)  = -dr_hh(1,3) -  hh(1,3)        / rr
           rhs_hyy(i,j,k)  = -dr_hh(2,2) - (hh(2,2) - 1.0d0) / rr
           rhs_hyz(i,j,k)  = -dr_hh(2,3) -  hh(2,3)        / rr
           rhs_hzz(i,j,k)  = -dr_hh(3,3) - (hh(3,3) - 1.0d0) / rr

           rhs_tracek(i,j,k) = -dr_trk - trk / rr

           rhs_axx(i,j,k)  = -dr_aa(1,1) - aa(1,1) / rr
           rhs_axy(i,j,k)  = -dr_aa(1,2) - aa(1,2) / rr
           rhs_axz(i,j,k)  = -dr_aa(1,3) - aa(1,3) / rr
           rhs_ayy(i,j,k)  = -dr_aa(2,2) - aa(2,2) / rr
           rhs_ayz(i,j,k)  = -dr_aa(2,3) - aa(2,3) / rr
           rhs_azz(i,j,k)  = -dr_aa(3,3) - aa(3,3) / rr

           rhs_gammatx(i,j,k) = -dr_gammat(1) - gammat(1) / rr
           rhs_gammaty(i,j,k) = -dr_gammat(2) - gammat(2) / rr
           rhs_gammatz(i,j,k) = -dr_gammat(3) - gammat(3) / rr

        end do
     end do
  end do

  if (CCTK_EQUALS(lapse_evolution_method, "LeanBSSNMoL")) then
     !$OMP PARALLEL DO COLLAPSE(3) &
     !$OMP PRIVATE( i, j, k, rr, dr_alph )
     do k = cctk_lsh(3)-cctk_nghostzones(3)+1, cctk_lsh(3)
        do j = 1, cctk_lsh(2)
           do i = 1, cctk_lsh(1)
              rr        = r(i,j,k)
              dr_alph  = (alp(i,j,k-2) - 4*alp(i,j,k-1) + 3*alp(i,j,k))*odr2
              rhs_alp(i,j,k)  = -dr_alph - (alp(i,j,k) - 1.0) / rr
           end do
        end do
     end do
  end if

  if (CCTK_EQUALS(shift_evolution_method, "LeanBSSNMoL")) then
     !$OMP PARALLEL DO COLLAPSE(3) &
     !$OMP PRIVATE( i, j, k, rr, dr_beta )
     do k = cctk_lsh(3)-cctk_nghostzones(3)+1, cctk_lsh(3)
        do j = 1, cctk_lsh(2)
           do i = 1, cctk_lsh(1)
              rr        = r(i,j,k)

              dr_beta(1)  = (betax(i,j,k-2) - 4*betax(i,j,k-1) + 3*betax(i,j,k))*odr2
              dr_beta(2)  = (betay(i,j,k-2) - 4*betay(i,j,k-1) + 3*betay(i,j,k))*odr2
              dr_beta(3)  = (betaz(i,j,k-2) - 4*betaz(i,j,k-1) + 3*betaz(i,j,k))*odr2

              rhs_betax(i,j,k) = -dr_beta(1) - betax(i,j,k) / rr
              rhs_betay(i,j,k) = -dr_beta(2) - betay(i,j,k) / rr
              rhs_betaz(i,j,k) = -dr_beta(3) - betaz(i,j,k) / rr
           end do
        end do
     end do
  end if


end subroutine LeanBSSN_calc_bssn_rhs_bdry_sph
