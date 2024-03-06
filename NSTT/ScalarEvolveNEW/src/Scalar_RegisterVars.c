
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void Scalar_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT use_jacobian;
  CCTK_INT ierr = 0, group, rhs;

  // register evolution and rhs gridfunction groups with MoL

  /* ADM metric and extrinsic curvature */
  group = CCTK_GroupIndex("ADMBase::lapse");
  ierr += MoLRegisterSaveAndRestoreGroup(group);
  group = CCTK_GroupIndex("ADMBase::shift");
  ierr += MoLRegisterSaveAndRestoreGroup(group);
  group = CCTK_GroupIndex("ADMBase::metric");
  ierr += MoLRegisterSaveAndRestoreGroup(group);
  group = CCTK_GroupIndex("ADMBase::curv");
  ierr += MoLRegisterSaveAndRestoreGroup(group);

  /* phi and rhs_phi */
  group = CCTK_GroupIndex("ScalarBase::phi");
  rhs   = CCTK_GroupIndex("ScalarEvolveNEW::rhs_phi");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Kphi and rhs_Kphi */
  group = CCTK_GroupIndex("ScalarBase::Kphi");
  rhs   = CCTK_GroupIndex("ScalarEvolveNEW::rhs_Kphi");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

/* if ((check_rho == 1) && (use_jacobian == 0)){
#  group = CCTK_GroupIndex("ScalarEvolveNEW::matter_density");
#  ierr += MoLRegisterSaveAndRestoreGroup(group);
# } */

  if (ierr) CCTK_ERROR("Problems registering with MoL");

}
