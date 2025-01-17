! Basegrid.F90 : Register symmetries of the grid functions
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine LeanBSSN_symmetries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  if (CCTK_EQUALS(lapse_evolution_method, "LeanBSSNMoL")) then
     call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::rhs_alp" )
  end if

  if (CCTK_EQUALS(scalar_evolution_method, "LeanBSSNMoL")) then
     call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::rhs_phi" )
     call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::rhs_Kphi" )
  end if

  if (CCTK_EQUALS(shift_evolution_method, "LeanBSSNMoL")) then
     call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "LeanBSSNMoL::rhs_betax" )
     call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "LeanBSSNMoL::rhs_betay" )
     call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "LeanBSSNMoL::rhs_betaz" )
  end if

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::conf_fac" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::rhs_conf_fac" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::hxx" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1, 1/), "LeanBSSNMoL::hxy" )
  call SetCartSymVN( ierr, cctkGH, (/-1, 1,-1/), "LeanBSSNMoL::hxz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::hyy" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1,-1/), "LeanBSSNMoL::hyz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::hzz" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::rhs_hxx" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1, 1/), "LeanBSSNMoL::rhs_hxy" )
  call SetCartSymVN( ierr, cctkGH, (/-1, 1,-1/), "LeanBSSNMoL::rhs_hxz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::rhs_hyy" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1,-1/), "LeanBSSNMoL::rhs_hyz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::rhs_hzz" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::axx" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1, 1/), "LeanBSSNMoL::axy" )
  call SetCartSymVN( ierr, cctkGH, (/-1, 1,-1/), "LeanBSSNMoL::axz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::ayy" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1,-1/), "LeanBSSNMoL::ayz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::azz" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::rhs_axx" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1, 1/), "LeanBSSNMoL::rhs_axy" )
  call SetCartSymVN( ierr, cctkGH, (/-1, 1,-1/), "LeanBSSNMoL::rhs_axz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::rhs_ayy" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1,-1/), "LeanBSSNMoL::rhs_ayz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::rhs_azz" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::tracek" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::rhs_tracek" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "LeanBSSNMoL::gammatx" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "LeanBSSNMoL::gammaty" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "LeanBSSNMoL::gammatz" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "LeanBSSNMoL::rhs_gammatx" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "LeanBSSNMoL::rhs_gammaty" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "LeanBSSNMoL::rhs_gammatz" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "LeanBSSNMoL::hc" )
  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "LeanBSSNMoL::mcx" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "LeanBSSNMoL::mcy" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "LeanBSSNMoL::mcz" )

end subroutine LeanBSSN_symmetries
