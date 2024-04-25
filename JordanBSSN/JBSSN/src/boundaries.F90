! boundaries.F90
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine JordanFBSSN_Boundaries( CCTK_ARGUMENTS )

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
       "JordanFBSSN::conf_fac", "none")
  if (ierr < 0)                                                            &
       call CCTK_ERROR("Failed to register BC for JordanFBSSN::conf_fac!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "JordanFBSSN::hmetric", "none")
  if (ierr < 0)                                                            &
       call CCTK_ERROR("Failed to register BC for JordanFBSSN::hmetric!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "JordanFBSSN::hcurv", "none")
  if (ierr < 0)                                                            &
       call CCTK_ERROR("Failed to register BC for JordanFBSSN::hcurv!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "JordanFBSSN::trk", "none")
  if (ierr < 0)                                                            &
       call CCTK_ERROR("Failed to register BC for JordanFBSSN::trk!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "JordanFBSSN::gammat", "none")
  if (ierr < 0)                                                            &
       call CCTK_ERROR("Failed to register BC for JordanFBSSN::gammat!")

  if (CCTK_EQUALS(lapse_evolution_method, "JordanFBSSN")) then
     ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,   &
          "ADMBase::lapse", "none")
     if (ierr < 0)                                                         &
          call CCTK_ERROR("Failed to register BC for ADMBase::lapse!")
  end if

  if (CCTK_EQUALS(shift_evolution_method, "JordanFBSSN")) then
     ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,   &
          "ADMBase::shift", "none")
     if (ierr < 0)                                                         &
          call CCTK_ERROR("Failed to register BC for ADMBase::shift!")
  end if

  if (CCTK_EQUALS(scalar_evolution_method, "JordanFBSSN")) then
     ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,   &
          "ScalarBase::phi", "none")
     if (ierr < 0)                                                         &
          call CCTK_ERROR("Failed to register BC for ScalarBase::phi1!")
     ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,   &
          "ScalarBase::Kphi", "none")
     if (ierr < 0)                                                         &
          call CCTK_ERROR("Failed to register BC for ScalarBase::Kphi1!")
  end if

end subroutine JordanFBSSN_Boundaries
!
!=============================================================================
!
subroutine JordanFBSSN_Constraints_Boundaries( CCTK_ARGUMENTS )

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
       "JordanFBSSN::ham", "flat")
  if (ierr < 0)                                                           &
       call CCTK_ERROR("Failed to register BC for JordanFBSSN::ham!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, bndsize, -one, &
       "JordanFBSSN::mom", "flat")
  if (ierr < 0)                                                           &
       call CCTK_ERROR("Failed to register BC for JordanFBSSN::mom!")

end subroutine JordanFBSSN_Constraints_Boundaries
