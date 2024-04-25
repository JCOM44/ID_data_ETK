
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

void JordanFBSSN_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, var, rhs;

  // register evolution and rhs gridfunction groups with MoL

  group = CCTK_GroupIndex("ADMBase::metric");
  ierr += MoLRegisterConstrainedGroup(group);
  group = CCTK_GroupIndex("ADMBase::curv");
  ierr += MoLRegisterConstrainedGroup(group);


  if (CCTK_EQUALS(lapse_evolution_method, "JordanFBSSN")){
    var   = CCTK_VarIndex("ADMBase::alp");
    rhs   = CCTK_VarIndex("JordanFBSSN::rhs_alp");
    ierr += MoLRegisterEvolved(var, rhs);
  }

  if (CCTK_EQUALS(scalar_evolution_method, "JordanFBSSN")){
    var   = CCTK_VarIndex("ScalarBase::phi1");
    rhs   = CCTK_VarIndex("JordanFBSSN::rhs_phi");
    ierr += MoLRegisterEvolved(var, rhs);

    var   = CCTK_VarIndex("ScalarBase::Kphi1");
    rhs   = CCTK_VarIndex("JordanFBSSN::rhs_kphi");
    ierr += MoLRegisterEvolved(var, rhs);
   }

  if (CCTK_EQUALS(shift_evolution_method, "JordanFBSSN")){
      var   = CCTK_VarIndex("ADMBase::betax");
      rhs   = CCTK_VarIndex("JordanFBSSN::rhs_betax");
      ierr += MoLRegisterEvolved(var, rhs);

      var   = CCTK_VarIndex("ADMBase::betay");
      rhs   = CCTK_VarIndex("JordanFBSSN::rhs_betay");
      ierr += MoLRegisterEvolved(var, rhs);

      var   = CCTK_VarIndex("ADMBase::betaz");
      rhs   = CCTK_VarIndex("JordanFBSSN::rhs_betaz");
      ierr += MoLRegisterEvolved(var, rhs);
  }

  var   = CCTK_VarIndex("JordanFBSSN::conf_fac");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_conf_fac");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JordanFBSSN::hxx");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_hxx");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JordanFBSSN::hxy");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_hxy");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JordanFBSSN::hxz");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_hxz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JordanFBSSN::hyy");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_hyy");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JordanFBSSN::hyz");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_hyz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JordanFBSSN::hzz");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_hzz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JordanFBSSN::axx");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_axx");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JordanFBSSN::axy");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_axy");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JordanFBSSN::axz");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_axz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JordanFBSSN::ayy");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_ayy");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JordanFBSSN::ayz");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_ayz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JordanFBSSN::azz");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_azz");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JordanFBSSN::tracek");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_tracek");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JordanFBSSN::gammatx");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_gammatx");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JordanFBSSN::gammaty");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_gammaty");
  ierr += MoLRegisterEvolved(var, rhs);

  var   = CCTK_VarIndex("JordanFBSSN::gammatz");
  rhs   = CCTK_VarIndex("JordanFBSSN::rhs_gammatz");
  ierr += MoLRegisterEvolved(var, rhs);



  if (ierr) CCTK_ERROR("Problems registering with MoL");

}
