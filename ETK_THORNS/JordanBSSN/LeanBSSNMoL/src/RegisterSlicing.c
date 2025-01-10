#include <cctk.h>
#include <Slicing.h>

int LeanBSSN_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("LeanBSSNMoL");
  return 0;
}
