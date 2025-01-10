/* intial data thorn: ID_SF_Read2D */
/*======================================================*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "ID_SF_utils.h"



#include "memory.h"
#include "readhdf5.h"
#include "interpolations.h"



  /*Reading related parameters */

/* -------------------------------------------------------------------*/
void ID_SF_Read2D(CCTK_ARGUMENTS);
void
ID_SF_Read2D (CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO("=== EBegin ID_SF_Read initial data ===");

  /*=== define parameters ====================================*/
  /*==========================================================*/


  /*Reading related parameters */
int m,s;
/*char filename[] = "/home/jolivera/python/phi.h5";*/
int sdiv,mdiv;
double  betaphi,axes_ratio,r_e;
double *s_gp,*mu_gp;
double **phi_td;




 // Allocate memory for s_gp and mu
      /* SET UP GRID */
      s_gp=malloc((SDIV+1)*sizeof(double));
      mu_gp=malloc((MDIV+1)*sizeof(double));

      /* ALLLOCATE MEMORY */

      phi_td = array_allocateSF(1,SDIV,1,MDIV);


      /* INITIALIZE VARIABLES WITH ZERO */

      for(s=1;s<=SDIV;s++)
        for(m=1;m<=MDIV;m++) {
	phi_td[s][m] = 0.0e0;
		}

/*=============Read from the HDF5 file and store in local variables====*/
hdf5_read_varSF(&sdiv,&mdiv,SFmodel_file,&axes_ratio,&r_e,&betaphi,s_gp,mu_gp,phi_td);


int i,j,k;
int ind;
int nx=cctkGH->cctk_lsh[0]; int ny=cctkGH->cctk_lsh[1]; int nz=cctkGH->cctk_lsh[2];

 CCTK_REAL *x_coord=0, *y_coord=0, *z_coord=0;  
 double phi_ijk;
CCTK_REAL *x_grid,*y_grid,*z_grid;
double x_i,y_j,z_k;


/*============ set up initial cartesian grid data==========*/

  x_coord = x;
  y_coord = y;
  z_coord = z;


for(i=1;i<=nx;i++)
        for(j=1;j<=ny;j++)
                for(k=1;k<=nz;k++) {
			ind = i-1+nx*(j-1+ny*(k-1));

			grid_interpSF(phi_td,s_gp,mu_gp,r_e,nx,ny,nz,x_coord,y_coord,z_coord,i,j,k,&phi_ijk,1);
		/*	if(phi_ijk==0.0){
				phi_ijk = phi_at_inf; 	
			}*/

			phi1[ind] = phi_ijk;
			/*phi2[ind] = 0.0;*/
			Kphi1[ind] = 0.0;
			/*Kphi2[ind] = 0.0;*/

					}

  CCTK_INFO("=== End ID_SF_Read2D initial data ===");
 /*=========== derivatives and advective derivatives ====*/
  CCTK_REAL     dphi_x,dphi_y,dphi_z;

  CCTK_REAL dx12,dy12,dz12;

  CCTK_INFO("=== Setting differences... ===");
  /*========= Define delta_xi =========*/
  dx12 = 12.0*CCTK_DELTA_SPACE(0);
  dy12 = 12.0*CCTK_DELTA_SPACE(1);
  dz12 = 12.0*CCTK_DELTA_SPACE(2);

printf("\n%f",dx12);
printf("\n%f",dy12);
printf("\n%f",dz12);
  CCTK_INFO("=== Calculating bounds ... ===");
  
  CCTK_INT imin[3], imax[3];
  for (int d = 0; d < 3; ++ d)
  {
    imin[d] = 0;
    imax[d] = cctk_lsh[d];
  }
  for (int k = imin[2]; k < imax[2]; ++k)
  {
   for (int j = imin[1]; j < imax[1]; ++j)
   {
    for (int i = imin[0]; i < imax[0]; ++i)
    {

     const int ind = CCTK_GFINDEX3D (cctkGH, i, j, k);



  /*====== calculate derivatives =======*/



   dphi_x = (   8*phi1[CCTK_GFINDEX3D (cctkGH, i+1, j, k)]
                - phi1[CCTK_GFINDEX3D (cctkGH, i+2, j, k)]
                - 8*phi1[CCTK_GFINDEX3D (cctkGH, i-1, j, k)]
                + phi1[CCTK_GFINDEX3D (cctkGH, i-2, j, k)]   ) / dx12;

   dphi_y = (   8*phi1[CCTK_GFINDEX3D (cctkGH, i, j+1, k)]
                - phi1[CCTK_GFINDEX3D (cctkGH, i, j+2, k)]
                - 8*phi1[CCTK_GFINDEX3D (cctkGH, i, j-1, k)]
                + phi1[CCTK_GFINDEX3D (cctkGH, i, j-2, k)]   ) / dy12;

   dphi_z = (   8*phi1[CCTK_GFINDEX3D (cctkGH, i, j, k+1)]
                - phi1[CCTK_GFINDEX3D (cctkGH, i, j, k+2)]
                - 8*phi1[CCTK_GFINDEX3D (cctkGH, i, j, k-1)]
                + phi1[CCTK_GFINDEX3D (cctkGH, i, j, k-2)]   ) / dz12;


	if (CCTK_EQUALS (theory, "BDdecouplingEF")){
    Kphi1[ind] = 1/(2*alp[ind]) * (betax[ind] * dphi_x 
                         + betay[ind] * dphi_y
                        + betaz[ind] * dphi_z );
						}
	else    {

    Kphi1[ind] = 1/(alp[ind]) * (betax[ind] * dphi_x 
                         + betay[ind] * dphi_y
                        + betaz[ind] * dphi_z );
		}
        }
     }
  }

}
/* -------------------------------------------------------------------*/
