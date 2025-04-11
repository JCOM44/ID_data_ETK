! boundaries.F90
!
!  JBSSN_Boundaries
!  JBSSN_Constraints_Boundaries
!  JBSSN_calc_bssn_rhs_bdry
!  JBSSN_calc_bssn_rhs_bdry_sph
!  JBSSN_calc_sf_rhs_bdry_sph
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine JBSSN_Boundaries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr

  CCTK_INT, parameter :: one = 1


  ! The outgoing (radiative) boundary conditions are being handled from the rhs
  ! routine through calls to the NewRad infrastructure. Here we register all
  ! BCs as 'none', which enforces all the symmetry BCs.

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "JBSSN::conf_fac", "none")
  if (ierr < 0)                                                            &
       call CCTK_ERROR("Failed to register BC for JBSSN::conf_fac!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "JBSSN::hmetric", "none")
  if (ierr < 0)                                                            &
       call CCTK_ERROR("Failed to register BC for JBSSN::hmetric!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "JBSSN::hcurv", "none")
  if (ierr < 0)                                                            &
       call CCTK_ERROR("Failed to register BC for JBSSN::hcurv!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "JBSSN::trk", "none")
  if (ierr < 0)                                                            &
       call CCTK_ERROR("Failed to register BC for JBSSN::trk!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "JBSSN::gammat", "none")
  if (ierr < 0)                                                            &
       call CCTK_ERROR("Failed to register BC for JBSSN::gammat!")

  if (CCTK_EQUALS(lapse_evolution_method, "JBSSN")) then
     ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,   &
          "ADMBase::lapse", "none")
     if (ierr < 0)                                                         &
          call CCTK_ERROR("Failed to register BC for ADMBase::lapse!")
  end if

  if (CCTK_EQUALS(shift_evolution_method, "JBSSN")) then
     ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,   &
          "ADMBase::shift", "none")
     if (ierr < 0)                                                         &
          call CCTK_ERROR("Failed to register BC for ADMBase::shift!")
  end if

  if (CCTK_EQUALS(scalar_evolution_method, "JBSSN")) then
     ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,   &
          "ScalarBase::phi", "none")
     if (ierr < 0)                                                         &
          call CCTK_ERROR("Failed to register BC for ScalarBase::phi!")
     ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,   &
          "ScalarBase::Kphi", "none")
     if (ierr < 0)                                                         &
          call CCTK_ERROR("Failed to register BC for ScalarBase::Kphi!")
  end if

end subroutine JBSSN_Boundaries
!
!=============================================================================
!
subroutine JBSSN_Constraints_Boundaries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr
  CCTK_INT, parameter :: one = 1
  CCTK_INT bndsize

  if (calculate_constraints_every .le. 0) then
     return
  end if

  if (MOD(cctk_iteration, calculate_constraints_every) .ne. 0 ) then
     return
  endif

  if (derivs_order == 6) then
     bndsize = 5
  else if (derivs_order == 4) then
     bndsize = 3
  else
     call CCTK_ERROR("derivs_order not yet implemented.")
  end if

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, bndsize, -one, &
       "JBSSN::ham", "flat")
  if (ierr < 0)                                                           &
       call CCTK_ERROR("Failed to register BC for JBSSN::ham!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, bndsize, -one, &
       "JBSSN::mom", "flat")
  if (ierr < 0)                                                           &
       call CCTK_ERROR("Failed to register BC for JBSSN::mom!")

end subroutine JBSSN_Constraints_Boundaries

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
 
  if (CCTK_EQUALS(scalar_evolution_method, "JBSSN") .AND. mass_phi==0) then
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
!
!===========================================================================
!
subroutine JBSSN_calc_sf_rhs_bdry_sph( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT  i, j, k
  CCTK_REAL odr2

  CCTK_REAL rr, expmass
  CCTK_REAL lphi, lKphi

  CCTK_REAL dr_lphi, dr_lKphi
  CCTK_REAL dr_rlphi
  CCTK_REAL k_wave, vel_g 

  CCTK_INT  reflevel

  odr2 = 1.0d0 / (2.0d0*CCTK_DELTA_SPACE(3))

  reflevel = GetRefinementLevel(cctkGH)


  ! Apply only on the coarsest level and in the spherical shell. Points marked
  ! with cctk_bbox(6) == 0 are inter-processor boundaries, so we also do not
  ! want those
  if (reflevel /= 0  .or. cctk_bbox(6) == 0) return


  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE( i, j, k, rr, lphi, lKphi, &
  !$OMP dr_lphi, dr_lKphi, expmass, dr_rlphi, k_wave, vel_g )
  do k = cctk_lsh(3)-cctk_nghostzones(3)+1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           rr        = r(i,j,k)
           lphi      = phi1(i,j,k)           
           lKphi     = Kphi1(i,j,k)

           expmass   = exp(-mass_phi*rr)

           dr_lphi   = (phi1(i,j,k-2) - 4*phi1(i,j,k-1) + 3*phi1(i,j,k))*odr2
           dr_rlphi   = ( r(i,j,k-2)*phi1(i,j,k-2) - 4*r(i,j,k-1)*phi1(i,j,k-1) + 3*r(i,j,k)*phi1(i,j,k) )*odr2

           dr_lKphi  = (Kphi1(i,j,k-2) - 4*Kphi1(i,j,k-1) + 3*Kphi1(i,j,k))*odr2
           
           ! In massive ST, we have a dispersion relation w^2 = k^2 + m^2, which results in a 
           ! group velocity smaller than 1 in some cases. To calculate the velocity we need to
           ! calculate it for every k.  Considering a plane wave Yukawa like solution 
           ! where phi ~ 1/r e^i(kr + omega t) the value of 
           ! k can be approximated by k ~ d/dr (r phi) / (r phi) 
           
           k_wave = abs(dr_rlphi/(rr*lphi)) 
           vel_g  = k_wave/sqrt(k_wave*k_wave + mass_phi*mass_phi)


           ! FIX: check the wave speeds below

           rhs_phi(i,j,k)  = - vel_g *  (dr_lphi + (lphi - phi1_0) / rr + mass_phi * lphi) 
           ! Not sure if Kphi has the same velocity
           rhs_Kphi(i,j,k) = - vel_g *  (dr_lKphi + (lKphi - Kphi1_0)  / rr + mass_phi * lKphi) 

        end do
     end do
  end do

end subroutine JBSSN_calc_sf_rhs_bdry_sph
