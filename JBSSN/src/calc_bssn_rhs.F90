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


subroutine JBSSN_calc_bssn_rhs( CCTK_ARGUMENTS )

!  use finite_difference_mod      
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
  CCTK_REAL                lphi_delta, lphi_delta2, delta_p
  CCTK_REAL                d1_lphi(3), d1_lKphi(3)
  CCTK_REAL                omega_Bphi, domega_Bphi, F_phi
  CCTK_REAL                d2_lphi(3,3), cd2_lphi(3,3), cW_dphi(3,3) 
  CCTK_REAL                trk_phi,trk_omega,trk_mass
  CCTK_REAL                aa_omega(3,3), aa_phi(3,3), gammat_omega(3), gammat_phi(3)

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
  CCTK_REAL                c_dalp_dphi, c_cd2_phi, c_cd2_phi_new, c_dww_dphi, c_dphi_dphi, c_lKphi  

  ! Matter variables
  CCTK_REAL                srcE, srcjdi(3), srcji(3), srcSij(3,3),    &
                           srcSijTF(3,3), srcS_ww2, src_trT

  ! Right hand sides
  CCTK_REAL                rhs_ww, rhs_hh(3,3), rhs_trk, rhs_aa(3,3),      &
                           rhs_gammat(3), rhs_beta(3), rhs_lphi, rhs_lKphi

  ! Misc variables
  CCTK_REAL                dx, dy, dz                   
  CCTK_REAL                dx12, dy12, dz12, dxsq12, dysq12, dzsq12,         &
                           dxdy144, dxdz144, dydz144
  CCTK_INT                 i, j, k, nx, ny, nz
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
  logical                   evolve_Jordan, JordanFrame, cowling 

  evolve_alp    = CCTK_EQUALS(lapse_evolution_method, "JBSSN")
  evolve_beta   = CCTK_EQUALS(shift_evolution_method, "JBSSN")
  evolve_scalar = CCTK_EQUALS(scalar_evolution_method, "JBSSN")
  evolve_Jordan = CCTK_EQUALS(theory,"DEFold") .OR. CCTK_EQUALS(theory,"BDold")
  JordanFrame   = CCTK_EQUALS(theory,"DEF") .OR. CCTK_EQUALS(theory,"BD") .OR. evolve_Jordan .OR. CCTK_EQUALS(theory,"full")
  cowling       = CCTK_EQUALS(theory,"onlySF") .OR. CCTK_EQUALS(theory,"BDdecouplingEF") .OR. CCTK_EQUALS(theory,"DEFdecouplingEF")

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

  dx = CCTK_DELTA_SPACE(1)
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)
  


  dx12 = 12*CCTK_DELTA_SPACE(1)
  dy12 = 12*CCTK_DELTA_SPACE(2)
  dz12 = 12*CCTK_DELTA_SPACE(3)

  dxsq12 = 12*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(1)
  dysq12 = 12*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(2)
  dzsq12 = 12*CCTK_DELTA_SPACE(3)*CCTK_DELTA_SPACE(3)

  dxdy144 = 144*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2)
  dxdz144 = 144*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(3)
  dydz144 = 144*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(3)




  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE( i, j, k, di, dj, dk, nx,ny,nz,  &
  !$OMP ww, hh, trk, aa, gammat, alph, beta, Tab, dethh, hu, beta_l, &
  !$OMP lphi, lKphi, Bphi, BKphi, omega_Bphi, domega_Bphi, F_phi, lphi2, lphi_delta, lphi_delta2,&
  !$OMP trk_omega, trk_phi, trk_mass, aa_omega, aa_phi, gammat_omega, gammat_phi, delta_p,&
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
  !$OMP c_dalp_dphi, c_cd2_phi, c_cd2_phi_new, c_dww_dphi, c_dphi_dphi, c_lKphi,& 
  !$OMP srcE, srcjdi, srcji, srcSij, srcS_ww2, srcSijTF, src_trT, &
  !$OMP a, b, c, l, m, n, p, q, jac, hes)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

    !------------ Get local variables ----------
    nx = size(conf_fac,1)
    ny = size(conf_fac,2)
    nz = size(conf_fac,3)

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

     !  derivs_order == 4

      !------------ Centered 1st derivatives -----      

! d1_ww(3)     

      call calc_d1(conf_fac,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz,d1_ww)            

      if (evolve_scalar) then 
! d1_phi1(3)  
       call  calc_d1(phi1,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz,d1_lphi)     

! d1_Kphi1(3)  
        call calc_d1(Kphi1,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_lKphi)             
      end if 
      
! d1_hh(3,3,3) 
      call calc_d1(hxx,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz,d1_hh11)     
      call calc_d1(hxy,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_hh12)     
      call calc_d1(hxz,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_hh13)     
      call calc_d1(hyy,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_hh22)     
      call calc_d1(hyz,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_hh23)     
      call calc_d1(hzz,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_hh33)     
      
! d1_trk(3)  
      call calc_d1(tracek,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_trk)           

! d1_aa(3,3,3)  

      call calc_d1(axx,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_aa11)     
      call calc_d1(axy,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_aa12)     
      call calc_d1(axz,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_aa13)     
      call calc_d1(ayy,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_aa22)     
      call calc_d1(ayz,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_aa23)     
      call calc_d1(azz,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_aa33)           

! d1_gammat(3,3)  
      call calc_d1(gammatx,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_gammat1)     
      call calc_d1(gammaty,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_gammat2)     
      call calc_d1(gammatz,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_gammat3)           

! d1_alph(3)
      call calc_d1(alp,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_alph)           

! d1_beta(3,3) 
      call calc_d1(betax,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_beta1)     
      call calc_d1(betay,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_beta2)     
      call calc_d1(betaz,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz, d1_beta3)     

      !--------------------------------------------------


      !------------- Centered 2nd derivatives -----------

! d2_ww(3,3)
      call calc_d2(conf_fac,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz,k_sum/=0, d2_ww)

!d2_phi(3,3)  
      if (evolve_scalar) then              
      call calc_d2(phi1,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz,k_sum/=0, d2_lphi)
      end if
      
! d2_hh(3,3,3,3)
      call calc_d2(hxx,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz,k_sum/=0, d2_hh11)
      call calc_d2(hxy,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz,k_sum/=0, d2_hh12)
      call calc_d2(hxz,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz,k_sum/=0, d2_hh13)
      call calc_d2(hyy,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz,k_sum/=0, d2_hh22)
      call calc_d2(hyz,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz,k_sum/=0, d2_hh23)
      call calc_d2(hzz,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz,k_sum/=0, d2_hh33)

! d2_alph(3,3)  
      call calc_d2(alp,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz,k_sum/=0,d2_alph)      


! d2_beta(3,3,3) 
      call calc_d2(betax,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz,k_sum/=0,d2_beta1)      
      call calc_d2(betay,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz,k_sum/=0,d2_beta2)      
      call calc_d2(betaz,nx,ny,nz, i,j,k, derivs_order, dx,dy,dz,k_sum/=0,d2_beta3)      

      !--------------------------------------------------


      !------------ Advection derivatives --------
      if( use_advection_stencils /= 0 ) then

        di = int( sign( one, beta_l(1) ) )
        dj = int( sign( one, beta_l(2) ) )
        dk = int( sign( one, beta_l(3) ) )        

! ad1_ww 
        call calc_ad1(conf_fac,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_ww)                 

        if (evolve_scalar) then
! ad1_phi       
        call calc_ad1(phi1,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_lphi)                                  

! ad1_Kphi        
!       NEW
        call calc_ad1(Kphi1,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_lKphi)                        
        end if  

! ad1_hh(3,3)
!       NEW

        call  calc_ad1(hxx,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_hh(1,1))                
        call  calc_ad1(hxy,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_hh(1,2))                       
        call  calc_ad1(hxz,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_hh(1,3))        
        call  calc_ad1(hyy,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_hh(2,2))        
        call  calc_ad1(hyz,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_hh(2,3))        
        call  calc_ad1(hzz,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_hh(3,3))       

        ad1_hh(2,1) = ad1_hh(1,2)
        ad1_hh(3,1) = ad1_hh(1,3)
        ad1_hh(3,2) = ad1_hh(2,3)


! ad1_trk        

        call calc_ad1(tracek,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_trk)  


! ad1_aa(3,3) 
        call calc_ad1(axx,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_aa(1,1))                
        call calc_ad1(axy,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_aa(1,2))                
        call calc_ad1(axz,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_aa(1,3))                
        call calc_ad1(ayy,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_aa(2,2))                
        call calc_ad1(ayz,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_aa(2,3))                
        call calc_ad1(azz,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_aa(3,3))                       
!
        ad1_aa(2,1) = ad1_aa(1,2)
        ad1_aa(3,1) = ad1_aa(1,3)
        ad1_aa(3,2) = ad1_aa(2,3)
        

! ad1_gammat(3)
        call calc_ad1(gammatx,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk,ad1_gammat(1))  
        call calc_ad1(gammaty,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk,ad1_gammat(2))  
        call calc_ad1(gammatz,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk,ad1_gammat(3))          

!ad1_alph!       
        if (evolve_alp) then
        call calc_ad1(alp,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_alph)          
        end if
        

! ad1_beta(3)
       if (evolve_beta) then
        call calc_ad1(betax,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_beta(1))  
        call calc_ad1(betay,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_beta(2))  
        call calc_ad1(betaz,beta_l,nx,ny,nz,i,j,k,dx,dy,dz,di,dj,dk, ad1_beta(3))                 
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

      end if ! use_advection_stencils /= 0

!    end if

    !-------------------------------------------
  if (use_jacobian) then
    call JBSSN_apply_jacobian(d1_trk, jac)
    call JBSSN_apply_jacobian(d1_gammat1, jac)
    call JBSSN_apply_jacobian(d1_gammat2, jac)
    call JBSSN_apply_jacobian(d1_gammat3, jac)
    call JBSSN_apply_jacobian(d1_aa11, jac)
    call JBSSN_apply_jacobian(d1_aa12, jac)
    call JBSSN_apply_jacobian(d1_aa13, jac)
    call JBSSN_apply_jacobian(d1_aa22, jac)
    call JBSSN_apply_jacobian(d1_aa23, jac)
    call JBSSN_apply_jacobian(d1_aa33, jac)

    call JBSSN_apply_jacobian2(d1_alph, d2_alph, jac, hes)
    call JBSSN_apply_jacobian2(d1_ww, d2_ww, jac, hes)
    call JBSSN_apply_jacobian2(d1_beta1, d2_beta1, jac, hes)
    call JBSSN_apply_jacobian2(d1_beta2, d2_beta2, jac, hes)
    call JBSSN_apply_jacobian2(d1_beta3, d2_beta3, jac, hes)
    call JBSSN_apply_jacobian2(d1_hh11, d2_hh11, jac, hes)
    call JBSSN_apply_jacobian2(d1_hh12, d2_hh12, jac, hes)
    call JBSSN_apply_jacobian2(d1_hh13, d2_hh13, jac, hes)
    call JBSSN_apply_jacobian2(d1_hh22, d2_hh22, jac, hes)
    call JBSSN_apply_jacobian2(d1_hh23, d2_hh23, jac, hes)
    call JBSSN_apply_jacobian2(d1_hh33, d2_hh33, jac, hes)
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

    !---Calculate scalar field cov derivative W factor
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
    ! U(Bphi)   = 2/B m^2 Bphi^2 varphi^2
    !
    ! For BD 
    ! Bphi      = exp(varphi)
    ! BKphi     = Bphi*Kphi
    ! d1_Bphi   = Bphi* d_phi
    ! cd2_Bphi  = Bphi * (d_phi * d_phi + covd_phi)
    ! d1_BKphi  = Bphi * ( Kphi * d_phi + d_kphi)
    ! omega     = 1/(2*k0**2) -3/2 
    ! d1_omega  = 0   ! omega is constant
    ! U(Bphi)   = 1/2 m^2/k0^2 Bphi^2 varphi^2
    !
    ! For DEF+BD
    ! Bphi      = exp(varphi*varphi/2+delta*varphi)
    ! BKphi     = Bphi * (varphi+delta) * Kphi
    ! d1_Bphi   = Bphi * (varphi+delta) * d1_phi
    ! cd2_Bphi   = Bphi * (1 + (varphi+delta)^2 ) * d1_varphi * d1_varphi + (varphi+delta) * Bphi * covD_phi
    ! d1_BKphi  = Bphi*Kphi * (1 + (varphi+delta)^2) * d1_varphi + (varphi+delta) * Bphi * d1_Kphi
    ! omega     = 2/(B * (varphi+delta)**2) -3/2 
    ! d1_omega  = -4/(B*Bphi * (varphi+delta)^4)
    ! U(Bphi)   = 2/B m^2 Bphi^2 varphi^2
    !

    if (evolve_scalar) then      ! Begin evolveScalar
       if (CCTK_EQUALS(theory,"BDold") .OR. CCTK_EQUALS(theory,"BD")    ) then                       !---define BD
            Bphi        = exp(lphi)
            BKphi       = Bphi * lKphi
            omega_Bphi  = 1.0d0/(2.0d0*k0BD*k0BD)-3.0d0/2.0d0
            domega_Bphi = 0.0d0
            
       end if                                                   ! -- end  Define BD
       
       if (CCTK_EQUALS(theory,"DEF") .OR. CCTK_EQUALS(theory,"DEFold") .OR. CCTK_EQUALS(theory,"DEFdecoupling")  ) then                      !--- Define DEF
               Bphi        = exp(lphi2/2.0d0)
               BKphi       = Bphi*lphi*lKphi
               omega_Bphi  = 2.0d0/(B_DEF*lphi2)-3.0d0/2.0d0
               domega_Bphi = -4.0d0/(B_DEF*lphi2*lphi2*Bphi)               
       end if                                                   ! -- end define DEF

       if (CCTK_EQUALS(theory,"full") .OR. CCTK_EQUALS(theory,"onlySF")  .OR.  CCTK_EQUALS(theory,"onlymetric")   ) then                      !--- Define DEF+BD
               delta_p     = -2*k0BD/sqrt(B_DEF)               
               lphi_delta  = lphi+delta_p
               lphi_delta2 = lphi_delta*lphi_delta
               Bphi        = exp(lphi2/2.0d0+delta_p*lphi)
               BKphi       = Bphi*lphi_delta*lKphi               
               omega_Bphi  = 2.0d0/(B_DEF*lphi_delta2)-3.0d0/2.0d0
               domega_Bphi = -4.0d0/(B_DEF*lphi_delta2*lphi_delta2*Bphi)               
       end if                                                   ! -- end define DEF+BD

            !---- Calculate traces             

                tr_dalp_dphi = 0.0d0
                tr_cd2_phi   = 0.0d0
                tr_cd2_phi_new = 0.0d0
                tr_dww_dphi  = 0.0d0
                tr_dphi_dphi = 0.0d0

                c_dalp_dphi   = 0.0d0
                c_cd2_phi     = 0.0d0
                c_cd2_phi_new = 0.0d0
                c_dww_dphi    = 0.0d0
                c_dphi_dphi   = 0.0d0 
                do a =1,3
                   do b = 1,3                                                
                          if (k_sum/=0) then  ! with summation compensation
                                call kahan_sum(tr_dalp_dphi,c_dalp_dphi,hu(a,b)*d1_alph(a)*d1_lphi(b))
                                call kahan_sum(tr_cd2_phi,c_cd2_phi,hu(a,b)*cd2_lphi(a,b))
                                call kahan_sum(tr_cd2_phi_new,c_cd2_phi_new,hu(a,b)*(cd2_lphi(a,b)+cW_dphi(a,b))*ww*ww)
                                call kahan_sum(tr_dww_dphi,c_dww_dphi,hu(a,b)*d1_ww(a)*d1_lphi(b))
                                call kahan_sum(tr_dphi_dphi,c_dphi_dphi,hu(a,b)*d1_lphi(a)*d1_lphi(b))
                           else 
                                tr_dalp_dphi = tr_dalp_dphi + hu(a,b) * d1_alph(a)*d1_lphi(b)
                                tr_cd2_phi   = tr_cd2_phi   + hu(a,b) * cd2_lphi(a,b)
                                tr_cd2_phi_new   = tr_cd2_phi_new   + hu(a,b) * (cd2_lphi(a,b)+cW_dphi(a,b) )*ww*ww
                                tr_dww_dphi  = tr_dww_dphi  + hu(a,b) * d1_ww(a) * d1_lphi(b)
                                tr_dphi_dphi = tr_dphi_dphi + hu(a,b) * d1_lphi(a)*d1_lphi(b)
                           end if         
                    end do
                end do

           !--- Debug options
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
               
                    if (r(i,j,k)==r_debug .AND. show_debug/=0) THEN
                    write(*,*) "------ Matter terms----------------"
                    write(*,*) "r = ", r_debug
                    write(*,*) "Tab", Tab 
                    write(*,*) "srcE", srcE
                    write(*,*) "srcS_ww2 ", srcS_ww2 
                    write(*,*) "srcSijTF ", srcSijTF
                    write(*,*) "srcji ", srcji          
                    write(*,*) "src_trT", src_trT
                    end if 

                ! if Scalar field is present 
               if(JordanFrame .OR. CCTK_EQUALS(theory,"onlymetric")) then  ! begin evolve Jordan
                       !!! Normalize matter sources with the SF
                       srcE     = srcE/Bphi
                       srcS_ww2 = srcS_ww2/Bphi
                       srcSijTF = srcSijTF/Bphi
                       srcji    = srcji/Bphi
                   
                    if(CCTK_EQUALS(theory,"BD") )  then   ! Init BD
                        F_phi = k0BD*k0BD*(8*pi*src_trT)/Bphi + mass_phi*mass_phi*lphi*Bphi                       

                        trk_omega = omega_Bphi*lKphi*lKphi 
                        trk_phi   = ww*ww*(tr_cd2_phi+tr_dphi_dphi) - ww*tr_dww_dphi - trk*lKphi -1.5d0*F_phi                                                
                        trk_mass  = (1.0d0/(2.0d0*k0BD*k0BD))*mass_phi*mass_phi*lphi2*Bphi 

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
                 

                    if ( CCTK_EQUALS(theory,"DEF") ) then   ! init DEF
                         F_phi = 2.0d0*pi*B_DEF*src_trT*lphi2/Bphi  +  ww*ww* tr_dphi_dphi - lKphi*lKphi + mass_phi*mass_phi*lphi2*Bphi
                         
                         trk_omega = omega_Bphi*lKphi*lKphi*lphi2 
                         trk_phi   = ww*ww*(lphi*tr_cd2_phi+(1.0d0+lphi2)*tr_dphi_dphi) - lphi*ww*tr_dww_dphi - trk*lphi*lKphi -1.5d0*F_phi 
                         trk_mass  = (1.0d0/B_DEF)*mass_phi*mass_phi*lphi2*Bphi 

                         gammat_phi   = 0.0d0
                         gammat_omega = 0.0d0
                         do a=1,3
                           do b=1,3
                            aa_omega(a,b)   = lphi2*d1_lphi(a)*d1_lphi(b) - lphi2*tr_dphi_dphi * hh(a,b)/3.0d0 
                            aa_phi(a,b)     = ww*ww*( d1_lphi(a)*d1_lphi(b)*(lphi2+1.0d0) + lphi *(cd2_lphi(a,b) + cW_dphi(a,b)) ) - aa(a,b)*lKphi*lphi 
                            aa_phi(a,b)     = aa_phi(a,b) - ww*hh(a,b)*( ww*(lphi2+1.0d0)*tr_dphi_dphi  + lphi*(ww*tr_cd2_phi - tr_dww_dphi) )/3.0d0   ! take trace
                                
                            !Need to check this next step                            
                            gammat_omega(a) = gammat_omega(a) + (1.0d0+2.0d0/B_DEF-lphi2/2.0d0) * lKphi * hu(a,b)*d1_lphi(b) 
                            gammat_phi(a)   = gammat_phi(a) + hu(a,b)*(lphi*d1_lKphi(b) - trk * d1_lphi(b)*lphi/3.0d0)  -  au(a,b)*d1_lphi(b)*lphi
                            end do
                          end do
                    end if  !end DEF

                    if ( CCTK_EQUALS(theory,"full") .OR. CCTK_EQUALS(theory,"onlymetric") ) then   ! init DEF+BD
                         F_phi = 2.0d0*pi*B_DEF*src_trT*lphi_delta2/Bphi  +  ww*ww* tr_dphi_dphi - lKphi*lKphi + mass_phi*mass_phi* Bphi * lphi*lphi_delta
                         
                         trk_omega = omega_Bphi*lKphi*lKphi*lphi_delta2 
                         trk_phi   = ww*ww*(lphi_delta*tr_cd2_phi+(1.0d0+lphi_delta2)*tr_dphi_dphi) - lphi_delta*ww*tr_dww_dphi - trk*lphi_delta*lKphi -1.5d0*F_phi 
                         trk_mass  = (1.0d0/B_DEF)*mass_phi*mass_phi*lphi2*Bphi 

                         gammat_phi   = 0.0d0
                         gammat_omega = 0.0d0
                         do a=1,3
                           do b=1,3
                            aa_omega(a,b)   = lphi_delta2*d1_lphi(a)*d1_lphi(b) - lphi_delta2*tr_dphi_dphi * hh(a,b)/3.0d0 
                            aa_phi(a,b)     = ww*ww*( d1_lphi(a)*d1_lphi(b)*(lphi_delta2+1.0d0) + lphi_delta *(cd2_lphi(a,b) + cW_dphi(a,b)) ) - aa(a,b)*lKphi*lphi_delta 
                            aa_phi(a,b)     = aa_phi(a,b) - ww*hh(a,b)*( ww*(lphi_delta2+1.0d0)*tr_dphi_dphi  + lphi_delta*(ww*tr_cd2_phi - tr_dww_dphi) )/3.0d0   ! take trace
                                
                            !Need to check this next step                            
                            gammat_omega(a) = gammat_omega(a) + (1.0d0+2.0d0/B_DEF-lphi_delta2/2.0d0) * lKphi * hu(a,b)*d1_lphi(b) 
                            gammat_phi(a)   = gammat_phi(a) + hu(a,b)*(lphi_delta*d1_lKphi(b) - trk * d1_lphi(b)*lphi_delta/3.0d0) -  au(a,b)*d1_lphi(b)*lphi_delta
                            end do
                          end do
                    end if  !end full DEF+BD

                end if                   !--- End JordanFranme 


                ! DEBUG
                    if (r(i,j,k)==r_debug .AND. show_debug/=0) THEN
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
       rhs_trk    = rhs_trk    + pi4  * alph * (srcE + ww*ww * srcS_ww2) 
       rhs_aa     = rhs_aa     - pi8  * alph * ww*ww * srcSijTF
       rhs_gammat = rhs_gammat - pi16 * alph * srcji
        ! ---- add scalar field terms as source too
      if(JordanFrame .OR. CCTK_EQUALS(theory,"onlymetric")) then
       rhs_trk    = rhs_trk + alph*(  trk_omega + trk_phi - trk_mass)
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

!    if(CCTK_EQUALS(theory,"onlymetric")) then        
!           rhs_tracek(i,j,k) = 0.0
!           rhs_axy(i,j,k)    = 0.0
!           rhs_axz(i,j,k)    = 0.0
!           rhs_ayy(i,j,k)    = 0.0
!           rhs_ayz(i,j,k)    = 0.0
!           rhs_azz(i,j,k)    = 0.0
!    end if

    if (cowling) then
    rhs_tracek(i,j,k) = 0.0d0

    rhs_axx(i,j,k) = 0.0d0
    rhs_axy(i,j,k) = 0.0d0
    rhs_axz(i,j,k) = 0.0d0
    rhs_ayy(i,j,k) = 0.0d0
    rhs_ayz(i,j,k) = 0.0d0
    rhs_azz(i,j,k) = 0.0d0

    rhs_gammatx(i,j,k) = 0.0d0
    rhs_gammaty(i,j,k) = 0.0d0
    rhs_gammatz(i,j,k) = 0.0d0


    end if
    !-------------------------------------------

    !------------ Now for the lapse -----------
    if (evolve_alp) then
            if (cowling) then
               rhs_alp(i,j,k) = 0.0d0
            else
               rhs_alp(i,j,k) = zeta_alpha * ad1_alph - 2.d0 * alph * trk
            end if
    end if
    !-------------------------------------------


    ! ----------- And shift -------------------
    if (evolve_beta) then
       ! rhs_beta
            if (cowling) then
               rhs_betax(i,j,k) = 0.0d0
               rhs_betay(i,j,k) = 0.0d0
               rhs_betaz(i,j,k) = 0.0d0
           else 
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
    end if
  
!!!!!!!!!!!!!!!!!!!!!!!!!
!  Scalar field evolution
!!!!!!!!!!!!!!!!!!!!!!!!
   if (evolve_scalar) then    ! Start evolve scalar field


           ! Advective derivatives
           rhs_lphi  = ad1_lphi
           rhs_lKphi = ad1_lKphi
           

           ! Add source to scalar field
! --- Equations in Einstein frame
           if (CCTK_EQUALS(theory,"BDdecouplingEF") .OR. CCTK_EQUALS(theory,"DEFdecouplingEF") ) then
                   rhs_lphi  = rhs_lphi-2.0d0 * alph * lKphi
                   rhs_phi(i,j,k) = rhs_lphi
                   rhs_lKphi = rhs_lKphi - 0.5d0 * alph * ww * ( ww*tr_cd2_phi - tr_dww_dphi)
                   rhs_lKphi = rhs_lKphi + alph*trk*lKphi -0.5d0*ww*ww*tr_dalp_dphi                   
                   rhs_lKphi = rhs_lKphi - 2*pi*alph*src_trT*(k0BD - sqrt(-betaDEF)* tanh(sqrt(-betaDEF)/5.0 *lphi)  ) !+betaDEF*lphi)     
                   rhs_lKphi = rhs_lKphi + 0.5d0 * alph * lphi * mass_phi*mass_phi                
                   rhs_Kphi(i,j,k) = rhs_lKphi
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
                                  
                    
                    if (CCTK_EQUALS(theory,"BD") ) then  ! add coupling term if BD
                        rhs_lKphi = rhs_lKphi - alph*tr_cd2_phi_new &
                               + alph*trk*lKphi -ww*ww*tr_dalp_dphi &
                               + alph*(lKphi*lKphi - ww*ww*tr_dphi_dphi) &
                               + 8.0d0*pi*alph*src_trT*(k0BD*k0BD)/Bphi &
                               + alph*lphi*mass_phi*mass_phi*Bphi                
                    end if 
              
                    if (CCTK_EQUALS(theory,"DEF") ) then  ! add coupling if DEF
                            if (k_sum/=0) then  ! with summation compensation
                                c_lKphi = 0.0d0
                                call kahan_sum(rhs_lKphi,c_lKphi,-ww*ww*tr_dalp_dphi)
                                call kahan_sum(rhs_lKphi,c_lKphi,-alph*tr_cd2_phi_new)
                                call kahan_sum(rhs_lKphi,c_lKphi,alph*trk*lKphi)
                                call kahan_sum(rhs_lKphi,c_lKphi,alph*lphi*lKphi*lKphi)
                                call kahan_sum(rhs_lKphi,c_lKphi,-alph*lphi*ww*ww*tr_dphi_dphi)
                                call kahan_sum(rhs_lKphi,c_lKphi,alph*lphi*2.0d0*pi*src_trT*B_DEF/Bphi)
                                call kahan_sum(rhs_lKphi,c_lKphi,alph*lphi*mass_phi*mass_phi*Bphi)
                             else
                                rhs_lKphi = rhs_lKphi -ww*ww*tr_dalp_dphi &
                                       + alph*( -tr_cd2_phi_new + trk*lKphi &                        
                                       + lphi*(lKphi*lKphi - ww*ww*tr_dphi_dphi &
                                       + 2.0d0*pi*src_trT*B_DEF/Bphi)) + alph*lphi*mass_phi*mass_phi*Bphi
                             end if 
                    end if              ! end DEF 
                    
                    if (CCTK_EQUALS(theory,"full") .OR. CCTK_EQUALS(theory,"onlySF") ) then  ! add coupling term if full                   
                        rhs_lKphi = rhs_lKphi - alph*tr_cd2_phi_new &
                               + alph*trk*lKphi -ww*ww*tr_dalp_dphi &
                               + alph*lphi_delta*(lKphi*lKphi - ww*ww*tr_dphi_dphi) &
                               + 2.0d0*pi*alph*src_trT*B_DEF*lphi_delta/Bphi &
                               + alph * mass_phi*mass_phi * Bphi * lphi               
                    end if      ! end full

                        


                    rhs_phi(i,j,k) = rhs_lphi
                    rhs_Kphi(i,j,k) = rhs_lKphi                                        
                
           end if               
!-------- End Jordan Frame SF ev

            if (r(i,j,k)==r_debug .AND. show_debug/=0) THEN
                    write(*,*) "---------- SF RHS --------------"
                    write(*,*) "r =,", r_debug                    
                    write(*,*) "tr_cd2_phi_new ", tr_cd2_phi_new
                    write(*,*) "trace term ", -alph*tr_cd2_phi_new
                    write(*,*) "coupling ", 2.0d0*pi*alph*src_trT*B_DEF*lphi_delta/Bphi  !2.0d0*alph*pi*src_trT*B_DEF*lphi/Bphi !8*pi*alph*src_trT*(k0BD*k0BD)/Bphi
                    write(*,*) "coup+trace ", -alph*tr_cd2_phi_new + 2.0d0*pi*alph*src_trT*B_DEF*lphi_delta/Bphi                     
                    write(*,*) "mass term ", alph*lphi*mass_phi*mass_phi*Bphi
                    write(*,*) "coup+trace+mass",  -alph*tr_cd2_phi_new + 2.0d0*pi*alph*src_trT*B_DEF*lphi_delta/Bphi + alph*lphi*mass_phi*mass_phi*Bphi
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
             
           if( CCTK_EQUALS(theory,"onlymetric")  ) then                   
                    rhs_phi(i,j,k) = 0.0d0
                    rhs_Kphi(i,j,k) = 0.0d0                                                            
           end if

    end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! end evolve scalar field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end do
  end do
  end do
  !$OMP END PARALLEL DO

end subroutine JBSSN_calc_bssn_rhs
!
!===========================================================================
!
subroutine JBSSN_calc_bssn_rhs_bdry( CCTK_ARGUMENTS )

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

  if (CCTK_EQUALS(lapse_evolution_method, "JBSSN")) then
     ierr = NewRad_Apply(cctkGH, alp, rhs_alp, one, one, n_alpha)
  end if

  if (CCTK_EQUALS(shift_evolution_method, "JBSSN")) then
     ierr = NewRad_Apply(cctkGH, betax, rhs_betax, zero, one, n_beta)
     ierr = NewRad_Apply(cctkGH, betay, rhs_betay, zero, one, n_beta)
     ierr = NewRad_Apply(cctkGH, betaz, rhs_betaz, zero, one, n_beta)
  end if
 
  if (CCTK_EQUALS(scalar_evolution_method, "JBSSN")) then
     if (CCTK_EQUALS(theory, "decouplingEF") .OR. CCTK_EQUALS(theory,"DEFdecouplingEF") ) then
             phi0  = phi1_0
             kphi0 = Kphi1_0
     end if
     if (CCTK_EQUALS(theory, "decouplingJ")) then
             phi0  = 1.0
             kphi0 = Kphi1_0
     end if
     if (CCTK_EQUALS(theory, "onlySF") .OR. CCTK_EQUALS(theory,"full") .OR. CCTK_EQUALS(theory,"DEFold") .OR. CCTK_EQUALS(theory,"BD") .OR. CCTK_EQUALS(theory,"DEF")) then             
             phi0  = phi_at_inf
             kphi0 = Kphi1_0
     end if
     ierr = NewRad_Apply(cctkGH, phi1, rhs_phi, phi0, one, n_phi1)
     ierr = NewRad_Apply(cctkGH, Kphi1, rhs_Kphi, kphi0, one, n_Kphi1)
  end if
  
end subroutine JBSSN_calc_bssn_rhs_bdry
!
!===========================================================================
!
subroutine JBSSN_calc_bssn_rhs_bdry_sph( CCTK_ARGUMENTS )

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

  if (CCTK_EQUALS(lapse_evolution_method, "JBSSN")) then
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

  if (CCTK_EQUALS(shift_evolution_method, "JBSSN")) then
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


end subroutine JBSSN_calc_bssn_rhs_bdry_sph
