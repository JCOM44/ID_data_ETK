#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <stddef.h>
#include <stdlib.h>


#define IMAX(a,b) ( a>b ? a : b ) 
#define IMIN(a,b) ( a<b ? a : b ) 
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SQ(x) ((x)*(x))              /* square macro */
#define DBL_EPSILON 1e-15


void huntSF(double xx[], 
          int n, 
          double x, 
          int *jlo);


double interp_4SF(double xp[5], 
                double yp[5], 
                int    np ,
                double xb);


void grid_interpSF(double **old, 
                 double *s_gp, 
                 double *mu, 
                 double r_e, 
                 int nx,
                 int ny,
                 int nz,  
                 double *x_grid,
                 double *y_grid,
                 double *z_grid, 
                 int i,
                 int j,
                 int k, 
                 double *new,
                 int sign);
 
