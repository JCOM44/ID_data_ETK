! Basegrid.F90 : Register symmetries of the grid functions
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine JBSSN_symmetries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  if (CCTK_EQUALS(lapse_evolution_method, "JBSSN")) then
     call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::rhs_alp" )
  end if

  if (CCTK_EQUALS(scalar_evolution_method, "JBSSN")) then
     call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::rhs_phi" )
     call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::rhs_Kphi" )
  end if

  if (CCTK_EQUALS(shift_evolution_method, "JBSSN")) then
     call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "JBSSN::rhs_betax" )
     call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "JBSSN::rhs_betay" )
     call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "JBSSN::rhs_betaz" )
  end if

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::conf_fac" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::rhs_conf_fac" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::hxx" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1, 1/), "JBSSN::hxy" )
  call SetCartSymVN( ierr, cctkGH, (/-1, 1,-1/), "JBSSN::hxz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::hyy" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1,-1/), "JBSSN::hyz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::hzz" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::rhs_hxx" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1, 1/), "JBSSN::rhs_hxy" )
  call SetCartSymVN( ierr, cctkGH, (/-1, 1,-1/), "JBSSN::rhs_hxz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::rhs_hyy" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1,-1/), "JBSSN::rhs_hyz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::rhs_hzz" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::axx" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1, 1/), "JBSSN::axy" )
  call SetCartSymVN( ierr, cctkGH, (/-1, 1,-1/), "JBSSN::axz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::ayy" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1,-1/), "JBSSN::ayz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::azz" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::rhs_axx" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1, 1/), "JBSSN::rhs_axy" )
  call SetCartSymVN( ierr, cctkGH, (/-1, 1,-1/), "JBSSN::rhs_axz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::rhs_ayy" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1,-1/), "JBSSN::rhs_ayz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::rhs_azz" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::tracek" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::rhs_tracek" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "JBSSN::gammatx" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "JBSSN::gammaty" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "JBSSN::gammatz" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "JBSSN::rhs_gammatx" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "JBSSN::rhs_gammaty" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "JBSSN::rhs_gammatz" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JBSSN::hc" )
  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "JBSSN::mcx" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "JBSSN::mcy" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "JBSSN::mcz" )

end subroutine JBSSN_symmetries
