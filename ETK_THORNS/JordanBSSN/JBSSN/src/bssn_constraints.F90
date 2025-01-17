! bssn_constraints.F90
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine JBSSN_bssn_constraints( CCTK_ARGUMENTS )


  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  ! Metric variables
  CCTK_REAL                alph, beta(3)
  CCTK_REAL                ww, hh(3,3), hu(3,3), trk, aa(3,3), gammat(3),   &
                           dethh, Tab(4,4)

  ! Scalar field 
  CCTK_REAL                lphi,lphi2, lKphi, Bphi, BKphi, delta_p, lphi_delta, lphi_delta2
  CCTK_REAL                tr_dalp_dphi, tr_cd2_phi, tr_cd2_phi_new, tr_dww_dphi, tr_dphi_dphi
  CCTK_REAL                ham_phi, mom_phi(3)

  ! First derivatives
  CCTK_REAL                d1_hh11(3), d1_hh12(3), d1_hh13(3), d1_hh22(3), d1_hh23(3), d1_hh33(3)
  CCTK_REAL                d1_aa11(3), d1_aa12(3), d1_aa13(3), d1_aa22(3), d1_aa23(3), d1_aa33(3)
  CCTK_REAL                d1_gammat1(3), d1_gammat2(3), d1_gammat3(3)
  CCTK_REAL                d1_lphi(3), d1_lKphi(3)
  CCTK_REAL                d1_ww(3), d1_hh(3,3,3), d1_trk(3), d1_aa(3,3,3),&
                           d1_gammat(3,3)

  ! Second derivatives
  CCTK_REAL                d2_hh11(3,3), d2_hh12(3,3), d2_hh13(3,3), d2_hh22(3,3), d2_hh23(3,3), d2_hh33(3,3)
  CCTK_REAL                d2_ww(3,3), d2_hh(3,3,3,3)
  CCTK_REAL                d2_lphi(3,3)

  ! Covariant derivatives
  CCTK_REAL                cd2_ww(3,3), cd1_aa(3,3,3), cd2_lphi(3,3), cW_dphi(3,3)

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
  CCTK_REAL                dx, dy, dz, aux          
  CCTK_REAL, parameter ::  one  = 1
  CCTK_REAL, parameter ::  pi   = acos(-one)
  CCTK_REAL, parameter ::  pi8  = 8*pi
  CCTK_REAL, parameter ::  pi16 = 16*pi
  CCTK_INT                 i, j, k, nx, ny, nz
  CCTK_INT                 a, b, c, l, m, n, p, q

  ! Jacobian
  CCTK_REAL                jac(3,3), hes(3,3,3)

  integer                  istat
  logical                  use_jacobian
  logical                  evolve_scalar
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

  evolve_scalar = CCTK_EQUALS(scalar_evolution_method, "JBSSN")
  
  dx = CCTK_DELTA_SPACE(1)
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)  

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
  !$OMP lphi, lphi2, lKphi, d1_lphi, d1_lKphi, d2_lphi,         &
  !$OMP tr_dalp_dphi, tr_cd2_phi, tr_cd2_phi_new, tr_dww_dphi,  &
  !$OMP delta_p, lphi_delta, lphi_delta2,                       &
  !$OMP tr_dphi_dphi, cd2_lphi, cW_dphi, ham_phi, mom_phi,      &
  !$OMP i, j, k, a, b, c, l, m, n, p, q, jac, hes,nx,ny,nz)
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

      !------------ Centered 1st derivatives -----
      ! d1_ww(3)

      ! d1_ww(3)           
      call calc_d1(conf_fac, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,d1_ww)

      if (evolve_scalar) then 
    !   d1_phi1(3)
        call calc_d1(phi1, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,d1_lphi)     

      ! d1_Kphi1(3)
        call calc_d1(Kphi1, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,d1_lKphi)     
      end if 
      
      ! d1_hh(3,3,3)

      call calc_d1(hxx, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz, d1_hh11)     
      call calc_d1(hxy, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,d1_hh12)     
      call calc_d1(hxz, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,d1_hh13)     
      call calc_d1(hyy, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,d1_hh22) 
      call calc_d1(hyz, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,d1_hh23)     
      call calc_d1(hzz, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,d1_hh33)    
      
      ! d1_trk(3)
      call calc_d1(tracek, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,d1_trk)     
      

      ! d1_aa(3,3,3)

      call calc_d1(axx, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,d1_aa11)     
      call calc_d1(axy, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz, d1_aa12)     
      call calc_d1(axz, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz, d1_aa13)     
      call calc_d1(ayy, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz, d1_aa22)     
      call calc_d1(ayz, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz, d1_aa23)     
      call calc_d1(azz, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz, d1_aa33)     



      ! d1_gammat(3,3)
      call calc_d1(gammatx, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz, d1_gammat1)     
      call calc_d1(gammaty, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz, d1_gammat2)     
      call calc_d1(gammatz, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz, d1_gammat3)     




      !------------- Centered 2nd derivatives -----------

      ! d2_ww(3,3)

      call calc_d2(conf_fac, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,k_sum/=0, d2_ww)

      if (evolve_scalar) then
      !d2_phi(3,3)
      call calc_d2(phi1, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,k_sum/=0, d2_lphi)
      end if


      ! d2_hh(3,3,3,3)

      call calc_d2(hxx, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,k_sum/=0, d2_hh11)
      call calc_d2(hxy, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,k_sum/=0, d2_hh12)
      call calc_d2(hxz, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,k_sum/=0, d2_hh13)
      call calc_d2(hyy, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,k_sum/=0, d2_hh22)
      call calc_d2(hyz, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,k_sum/=0, d2_hh23)
      call calc_d2(hzz, nx,ny,nz,i,j,k, derivs_order, dx,dy,dz,k_sum/=0,d2_hh33)      
      
    !--------------------------------------------------

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

      call JBSSN_apply_jacobian2(d1_ww, d2_ww, jac, hes)
      call JBSSN_apply_jacobian2(d1_hh11, d2_hh11, jac, hes)
      call JBSSN_apply_jacobian2(d1_hh12, d2_hh12, jac, hes)
      call JBSSN_apply_jacobian2(d1_hh13, d2_hh13, jac, hes)
      call JBSSN_apply_jacobian2(d1_hh22, d2_hh22, jac, hes)
      call JBSSN_apply_jacobian2(d1_hh23, d2_hh23, jac, hes)
      call JBSSN_apply_jacobian2(d1_hh33, d2_hh33, jac, hes)
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

    if (evolve_scalar) then 
       cd2_lphi = d2_lphi
    end if

    do a = 1, 3
      do b = a, 3
        do l = 1, 3
          cd2_ww(a,b) = cd2_ww(a,b) - cf2(l,a,b) * d1_ww(l)

          if (evolve_scalar) then 
                 cd2_lphi(a,b) = cd2_lphi(a,b) - cf2(l,a,b) * d1_lphi(l)
          end if

           do m = 1, 3
            cd1_aa(a,b,l) = cd1_aa(a,b,l) - cf2(m,a,l) * aa(b,m) - cf2(m,b,l) * aa(a,m)
           end do
          end do   ! l
        end do     ! b
      end do       ! a 
    cd2_ww(2,1)   = cd2_ww(1,2)
    cd2_ww(3,1)   = cd2_ww(1,3)
    cd2_ww(3,2)   = cd2_ww(2,3)
    cd1_aa(2,1,:) = cd1_aa(1,2,:)
    cd1_aa(3,1,:) = cd1_aa(1,3,:)
    cd1_aa(3,2,:) = cd1_aa(2,3,:)
    if (evolve_scalar) then 
            cd2_lphi(2,1) = cd2_lphi(1,2)
            cd2_lphi(3,1) = cd2_lphi(1,3)
            cd2_lphi(3,2) = cd2_lphi(2,3)
    end if
    !-------------------------------------------

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

    !-------------------------------------------
    !------------ Scalar field terms ----------
    
    if (evolve_scalar) then      ! Begin evolveScalar

            !---- Calculate traces                                     
        tr_cd2_phi_new = 0         
        tr_dphi_dphi = 0
        
        do a =1,3
           do b = 1,3                      
              tr_cd2_phi_new   = tr_cd2_phi_new   + hu(a,b) * (cd2_lphi(a,b)+cW_dphi(a,b) )*ww*ww              
              tr_dphi_dphi = tr_dphi_dphi + hu(a,b) * d1_lphi(a)*d1_lphi(b)
           end do           
        end do
        tr_cd2_phi_new = tr_cd2_phi_new * ww*ww         
        tr_dphi_dphi = tr_dphi_dphi*ww*ww
        
        !!! Choose theory

       if (  CCTK_EQUALS(theory,"BD") ) then                       !---define BD
            Bphi        = exp(lphi)
            BKphi       = Bphi * lKphi

            ham_phi = 0.0d0 
            mom_phi(1) = 0.d0
            mom_phi(2) = 0.d0
            mom_phi(3) = 0.d0

            ham_phi =  0.5d0 * (1.0d0 /(k0BD*k0BD) - 3.0d0) * lKphi*lKphi
            ham_phi = ham_phi +  0.5d0 * ( 1.0d0/(k0BD*k0BD)+1.0d0) *tr_dphi_dphi
            ham_phi = ham_phi + 2.0d0*(-trk * lKphi + tr_cd2_phi_new)             
            

            do a=1,3
              do b=1,3
                do c=1,3                   
                   mom_phi(a) = mom_phi(a) +  hu(b,c)*aa(a,b)*d1_lphi(c)
                end do
              end do
            end do
            ! write like this to avoid error accumulation for multiple substractions
            mom_phi = ( 0.5d0 * (1.0d0/(k0BD*k0BD) -1.0d0)*lKphi -trk/3.0d0 )*d1_lphi + d1_lKphi - mom_phi

       end if                                                   ! -- end  Define BD
       
       if (CCTK_EQUALS(theory,"DEF") .OR. CCTK_EQUALS(theory,"DEFdecoupling")) then                      !--- Define DEF
            ham_phi = 0.0d0 
            mom_phi(1) = 0.d0
            mom_phi(2) = 0.d0
            mom_phi(3) = 0.d0
             Bphi        = exp(lphi2/2.0d0)
             BKphi       = Bphi*lphi*lKphi

            ham_phi =  (2.0d0 /B_DEF - 1.5d0*lphi2) * lKphi*lKphi
            !ham_phi = ham_phi +  0.5d0 * (1.0d0 + 2.0d0/B_DEF - lphi2/2.0d0) *tr_dphi_dphi
            ham_phi = ham_phi +  (2.0d0 + 2.0d0/B_DEF + lphi2/2.0d0) *tr_dphi_dphi
            ham_phi = ham_phi + 2.0d0*lphi*(-trk * lKphi + tr_cd2_phi_new) 
            ham_phi = ham_phi + 2.0d0 * mass_phi*mass_phi * lphi2 *Bphi/B_DEF           
                       
            do a=1,3
              do b=1,3
                do c=1,3                   
                   mom_phi(a) = mom_phi(a) + lphi*hu(b,c)*aa(a,b)*d1_lphi(c)
                end do
              end do
            end do
            ! write like this to avoid error accumulation for multiple substractions
            mom_phi =  ( (2.0d0/B_DEF +1.0d0-lphi2/2.0d0)*lKphi - trk/3.0d0 )*d1_lphi + lphi*d1_lKphi - mom_phi

       end if                                                   ! -- end define DEF

       if (CCTK_EQUALS(theory,"full")  .OR. CCTK_EQUALS(theory,"onlySF")) then               !--- Define DEF+BD
            ham_phi = 0.0d0 
            mom_phi(1) = 0.d0
            mom_phi(2) = 0.d0
            mom_phi(3) = 0.d0
            delta_p     = -2*k0BD/sqrt(B_DEF)               
            lphi_delta  = lphi+delta_p
            lphi_delta2 = lphi_delta*lphi_delta
            Bphi        = exp(lphi2/2.0d0+delta_p*lphi)
            BKphi       = Bphi*lphi_delta*lKphi               


            ham_phi =  (2.0d0 /B_DEF - 1.5d0*lphi_delta2) * lKphi*lKphi
            !ham_phi = ham_phi +  0.5d0 * (1.0d0 + 2.0d0/B_DEF - lphi2/2.0d0) *tr_dphi_dphi
            ham_phi = ham_phi +  (2.0d0 + 2.0d0/B_DEF + lphi_delta2/2.0d0) *tr_dphi_dphi
            ham_phi = ham_phi + 2.0d0*lphi_delta2*(-trk * lKphi + tr_cd2_phi_new) 
            ham_phi = ham_phi + 2.0d0 * mass_phi*mass_phi * lphi_delta2 *Bphi/B_DEF           
                       
            do a=1,3
              do b=1,3
                do c=1,3                   
                   mom_phi(a) = mom_phi(a) + lphi_delta*hu(b,c)*aa(a,b)*d1_lphi(c)
                end do
              end do
            end do
            ! write like this to avoid error accumulation for multiple substractions
            mom_phi =  ( (2.0d0/B_DEF +1.0d0-lphi_delta2/2.0d0)*lKphi - trk/3.0d0 )*d1_lphi + lphi_delta*d1_lKphi - mom_phi
       end if                                                   ! -- end define DEF+BD


        
        ham = ham - ham_phi
        mom = mom - mom_phi

else 
        Bphi = 1.0d0       
   
end if  !! End evolveScalar


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




!! Matter sources
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
       ham = ham - pi16 * srcE/Bphi
       mom = mom - pi8  * srcjdi/Bphi

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

end subroutine JBSSN_bssn_constraints
