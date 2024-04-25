! Basegrid.F90 : Register symmetries of the grid functions
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine JordanFBSSN_symmetries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  if (CCTK_EQUALS(lapse_evolution_method, "JordanFBSSN")) then
     call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::rhs_alp" )
  end if

  if (CCTK_EQUALS(scalar_evolution_method, "JordanFBSSN")) then
     call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::rhs_phi" )
     call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::rhs_kphi" )
  end if

  if (CCTK_EQUALS(shift_evolution_method, "JordanFBSSN")) then
     call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "JordanFBSSN::rhs_betax" )
     call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "JordanFBSSN::rhs_betay" )
     call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "JordanFBSSN::rhs_betaz" )
  end if

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::conf_fac" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::rhs_conf_fac" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::hxx" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1, 1/), "JordanFBSSN::hxy" )
  call SetCartSymVN( ierr, cctkGH, (/-1, 1,-1/), "JordanFBSSN::hxz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::hyy" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1,-1/), "JordanFBSSN::hyz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::hzz" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::rhs_hxx" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1, 1/), "JordanFBSSN::rhs_hxy" )
  call SetCartSymVN( ierr, cctkGH, (/-1, 1,-1/), "JordanFBSSN::rhs_hxz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::rhs_hyy" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1,-1/), "JordanFBSSN::rhs_hyz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::rhs_hzz" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::axx" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1, 1/), "JordanFBSSN::axy" )
  call SetCartSymVN( ierr, cctkGH, (/-1, 1,-1/), "JordanFBSSN::axz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::ayy" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1,-1/), "JordanFBSSN::ayz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::azz" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::rhs_axx" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1, 1/), "JordanFBSSN::rhs_axy" )
  call SetCartSymVN( ierr, cctkGH, (/-1, 1,-1/), "JordanFBSSN::rhs_axz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::rhs_ayy" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1,-1/), "JordanFBSSN::rhs_ayz" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::rhs_azz" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::tracek" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::rhs_tracek" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "JordanFBSSN::gammatx" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "JordanFBSSN::gammaty" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "JordanFBSSN::gammatz" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "JordanFBSSN::rhs_gammatx" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "JordanFBSSN::rhs_gammaty" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "JordanFBSSN::rhs_gammatz" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "JordanFBSSN::hc" )
  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "JordanFBSSN::mcx" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "JordanFBSSN::mcy" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "JordanFBSSN::mcz" )



end subroutine JordanFBSSN_symmetries
