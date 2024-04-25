#include <cctk.h>
#include <Slicing.h>

int JordanFBSSN_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("JordanFBSSN");
  return 0;
}
