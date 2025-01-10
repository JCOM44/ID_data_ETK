#include "interpolations.h"
#include "memory.h"

/********************************/
/*  Driver to find the interval of our */
/* grid where a point is contained, returns the lower index that*/
/* delimits the domain where the point is contained  */
/*********************************/
void huntSF(double xx[], int n, double x, int *jlo)
{ 
	int jm,jhi,inc,ascnd;

	ascnd=(xx[n] > xx[1]);
	if (*jlo <= 0 || *jlo > n) {
		*jlo=0;
		jhi=n+1;
	} else {
		inc=1;
		if (x >= xx[*jlo] == ascnd) {
			if (*jlo == n) return;
			jhi=(*jlo)+1;
			while (x >= xx[jhi] == ascnd) {
				*jlo=jhi;
				inc += inc;
				jhi=(*jlo)+inc;
				if (jhi > n) {
					jhi=n+1;
					break;
				}
			}
		} else {
			if (*jlo == 1) {
				*jlo=0;
				return;
			}
			jhi=(*jlo);
			*jlo -= 1;
			while (x < xx[*jlo] == ascnd) {
				jhi=(*jlo);
				inc += inc;
				*jlo=jhi-inc;
				if (*jlo < 1) {
					*jlo=0;
					break;
				}
			}
		}
	}
	while (jhi-(*jlo) != 1) {
		jm=(jhi+(*jlo)) >> 1;
		if (x > xx[jm] == ascnd)
			*jlo=jm;
		else
			jhi=jm;
	}
}


/*************************************************************************/
/* Driver for the interpolation routine. First we find the tab. point    */
/* nearest to xb, then we interpolate using four points around xb.       */  
/*************************************************************************/
double interp_4SF(double xp[5], 
                double yp[5], 
                int    np ,
                double xb)
{ 
 int k=1;      /* index of 1st point */ 
 double y;     /* intermediate value */


 if( xb==xp[k] ||  xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3]) 
    xb += DBL_EPSILON;

 y= (xb-xp[k+1])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k]/
        ((xp[k]-xp[k+1])*(xp[k]-xp[k+2])*(xp[k]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k+1]/
       ((xp[k+1]-xp[k])*(xp[k+1]-xp[k+2])*(xp[k+1]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+3])*yp[k+2]/
       ((xp[k+2]-xp[k])*(xp[k+2]-xp[k+1])*(xp[k+2]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+2])*yp[k+3]/
       ((xp[k+3]-xp[k])*(xp[k+3]-xp[k+1])*(xp[k+3]-xp[k+2]));

 return (y);
 
}


/*************************************************************************
* Interpolation between two different grids
*************************************************************************/
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
                 int sign) {

int s,
    m,
    s_nearest,
    m_nearest,                /* nearest points in interpolation */
    k_s,                      /* first s point in interpolation */
    k_m;                      /* first s point in interpolation */


double r_c,                   /* r of cartesian x,y,z point */
       s_c,                   /* s of cartesian x,y,z point */
       mu_c,                  /* mu of cartesian x,y,z point */
       s_4[5],                /* s of the 4 nearest points */
       mu_4[5],               /* mu of the 4 nearest points */
       old_s[5],              /* old at 4 nearest constant s points */
       old_m[5];              /* old at 4 nearest constant mu points */  


	    
      r_c = sqrt(   SQ(x_grid[i-1+nx*(j-1+ny*(k-1))])   
                       + SQ(y_grid[i-1+nx*(j-1+ny*(k-1))]) 
                       + SQ(z_grid[i-1+nx*(j-1+ny*(k-1))]) ); 
      s_c = r_c/(r_e+r_c);
      


      if(r_c==0.0 ) 
        mu_c = 0.0;
      else
        mu_c = fabs(z_grid[i-1+nx*(j-1+ny*(k-1))])/r_c;
 
      s_nearest = 0; m_nearest = 0;

      huntSF(s_gp, SDIV, s_c, &s_nearest);
      huntSF(mu, MDIV, mu_c, &m_nearest); 
  
      k_s = IMIN(IMAX((s_nearest)-(4-1)/2,1),SDIV+1-4);
      k_m = IMIN(IMAX((m_nearest)-(4-1)/2,1),MDIV+1-4);
      	       


      for(s=1;s<=4;s++) 
            s_4[s] = s_gp[k_s-1+s];


      for(m=1;m<=4;m++) 
            mu_4[m] = mu[k_m-1+m]; 
			

      for(s=1;s<=4;s++) {
          for(m=1;m<=4;m++) {
                  old_s[m] = old[k_s-1+s][k_m-1+m];
          }
          old_m[s] = interp_4SF(mu_4, old_s, 4, mu_c);  
      }

/*     if (z_grid[i-1+nx*(j-1+ny*(k-1))]<=0.0)                               
              (*new) = (1.0*sign)*interp_4SF(s_4, old_m, 4, s_c); 
     else */
              (*new) = interp_4SF(s_4, old_m, 4, s_c); 

	
}


