#include <cctk.h>
#include <Slicing.h>

int JBSSN_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("JBSSN");
  return 0;
}
