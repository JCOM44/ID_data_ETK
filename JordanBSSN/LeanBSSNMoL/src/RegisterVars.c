
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

void LeanBSSN_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, var, rhs;

  // register evolution and rhs gridfunction groups with MoL

  group = CCTK_GroupIndex("ADMBase::metric");
  ierr += MoLRegisterConstrainedGroup(group);
  group = CCTK_GroupIndex("ADMBase::curv");
  ierr += MoLRegisterConstrainedGroup(group);


  if (CCTK_EQUALS(lapse_evolution_method, "LeanBSSNMoL")){
    var   = CCTK_VarIndex("ADMBase::alp");
    rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_alp");
    ierr += MoLRegisterEvolved(var, rhs);
  }

  if (CCTK_EQUALS(scalar_evolution_method, "LeanBSSNMoL")){
    var   = CCTK_VarIndex("ScalarBase::phi1");
    rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_phi");
    ierr += MoLRegisterEvolved(var, rhs);
    
    var   = CCTK_VarIndex("ScalarBase::Kphi1");
    rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_Kphi");
    ierr += MoLRegisterEvolved(var, rhs);
  }

  if (CCTK_EQUALS(shift_evolution_method, "LeanBSSNMoL")){
      var   = CCTK_VarIndex("ADMBase::betax");
      rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_betax");
      ierr += MoLRegisterEvolved(var, rhs);

      var   = CCTK_VarIndex("ADMBase::betay");
      rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_betay");
      ierr += MoLRegisterEvolved(var, rhs);

      var   = CCTK_VarIndex("ADMBase::betaz");
      rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_betaz");
      ierr += MoLRegisterEvolved(var, rhs);
  }

  var   = CCTK_VarIndex("LeanBSSNMoL::conf_fac");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_conf_fac");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("LeanBSSNMoL::hxx");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_hxx");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("LeanBSSNMoL::hxy");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_hxy");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("LeanBSSNMoL::hxz");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_hxz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("LeanBSSNMoL::hyy");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_hyy");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("LeanBSSNMoL::hyz");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_hyz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("LeanBSSNMoL::hzz");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_hzz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("LeanBSSNMoL::axx");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_axx");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("LeanBSSNMoL::axy");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_axy");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("LeanBSSNMoL::axz");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_axz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("LeanBSSNMoL::ayy");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_ayy");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("LeanBSSNMoL::ayz");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_ayz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("LeanBSSNMoL::azz");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_azz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("LeanBSSNMoL::tracek");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_tracek");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("LeanBSSNMoL::gammatx");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_gammatx");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("LeanBSSNMoL::gammaty");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_gammaty");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("LeanBSSNMoL::gammatz");
  rhs   = CCTK_VarIndex("LeanBSSNMoL::rhs_gammatz");
  ierr += MoLRegisterEvolved(var, rhs);


  if (ierr) CCTK_ERROR("Problems registering with MoL");

}
