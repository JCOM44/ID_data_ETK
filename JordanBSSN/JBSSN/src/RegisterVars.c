
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

void JBSSN_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, var, rhs;

  // register evolution and rhs gridfunction groups with MoL

  group = CCTK_GroupIndex("ADMBase::metric");
  ierr += MoLRegisterConstrainedGroup(group);
  group = CCTK_GroupIndex("ADMBase::curv");
  ierr += MoLRegisterConstrainedGroup(group);


  if (CCTK_EQUALS(lapse_evolution_method, "JBSSN")){
    var   = CCTK_VarIndex("ADMBase::alp");
    rhs   = CCTK_VarIndex("JBSSN::rhs_alp");
    ierr += MoLRegisterEvolved(var, rhs);
  }

  if (CCTK_EQUALS(scalar_evolution_method, "JBSSN")){
    var   = CCTK_VarIndex("ScalarBase::phi1");
    rhs   = CCTK_VarIndex("JBSSN::rhs_phi");
    ierr += MoLRegisterEvolved(var, rhs);
    
    var   = CCTK_VarIndex("ScalarBase::Kphi1");
    rhs   = CCTK_VarIndex("JBSSN::rhs_Kphi");
    ierr += MoLRegisterEvolved(var, rhs);
  }

  if (CCTK_EQUALS(shift_evolution_method, "JBSSN")){
      var   = CCTK_VarIndex("ADMBase::betax");
      rhs   = CCTK_VarIndex("JBSSN::rhs_betax");
      ierr += MoLRegisterEvolved(var, rhs);

      var   = CCTK_VarIndex("ADMBase::betay");
      rhs   = CCTK_VarIndex("JBSSN::rhs_betay");
      ierr += MoLRegisterEvolved(var, rhs);

      var   = CCTK_VarIndex("ADMBase::betaz");
      rhs   = CCTK_VarIndex("JBSSN::rhs_betaz");
      ierr += MoLRegisterEvolved(var, rhs);
  }

  var   = CCTK_VarIndex("JBSSN::conf_fac");
  rhs   = CCTK_VarIndex("JBSSN::rhs_conf_fac");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JBSSN::hxx");
  rhs   = CCTK_VarIndex("JBSSN::rhs_hxx");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JBSSN::hxy");
  rhs   = CCTK_VarIndex("JBSSN::rhs_hxy");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JBSSN::hxz");
  rhs   = CCTK_VarIndex("JBSSN::rhs_hxz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JBSSN::hyy");
  rhs   = CCTK_VarIndex("JBSSN::rhs_hyy");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JBSSN::hyz");
  rhs   = CCTK_VarIndex("JBSSN::rhs_hyz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JBSSN::hzz");
  rhs   = CCTK_VarIndex("JBSSN::rhs_hzz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JBSSN::axx");
  rhs   = CCTK_VarIndex("JBSSN::rhs_axx");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JBSSN::axy");
  rhs   = CCTK_VarIndex("JBSSN::rhs_axy");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JBSSN::axz");
  rhs   = CCTK_VarIndex("JBSSN::rhs_axz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JBSSN::ayy");
  rhs   = CCTK_VarIndex("JBSSN::rhs_ayy");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JBSSN::ayz");
  rhs   = CCTK_VarIndex("JBSSN::rhs_ayz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JBSSN::azz");
  rhs   = CCTK_VarIndex("JBSSN::rhs_azz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JBSSN::tracek");
  rhs   = CCTK_VarIndex("JBSSN::rhs_tracek");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JBSSN::gammatx");
  rhs   = CCTK_VarIndex("JBSSN::rhs_gammatx");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JBSSN::gammaty");
  rhs   = CCTK_VarIndex("JBSSN::rhs_gammaty");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JBSSN::gammatz");
  rhs   = CCTK_VarIndex("JBSSN::rhs_gammatz");
  ierr += MoLRegisterEvolved(var, rhs);


  if (ierr) CCTK_ERROR("Problems registering with MoL");

}
