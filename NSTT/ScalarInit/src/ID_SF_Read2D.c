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
#include <hdf5.h>


#include "memory.h"
#include "readhdf5.h"
#include "interpolations.h"



  /*Reading related parameters */

/* -------------------------------------------------------------------*/
void ID_SF_Read2D(CCTK_ARGUMENTS);
void
ID_SF_Read2D (CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INFO("=== Begin ID_SF_Read2D initial data ===");

  /*=== define parameters ====================================*/
  /*==========================================================*/


  /*====================Reading related parameters ===========*/
  int m,s;
  int sdiv,mdiv;
  double  betaphi,axes_ratio,r_e;
  double *s_gp,*mu;
  double **phi_td;



 /*======== Allocate memory for s_gp and mu=====*/
      /* SET UP GRID */
      s_gp=malloc((SDIV+1)*sizeof(double));
      mu=malloc((MDIV+1)*sizeof(double));

      /* ALLLOCATE MEMORY */

      phi_td = array_allocate(1,SDIV,1,MDIV);


      /* INITIALIZE VARIABLES WITH ZERO */

      for(s=1;s<=SDIV;s++)
        for(m=1;m<=MDIV;m++) {
	phi_td[s][m] = 0.0e0;
		}

/*=============Read from the HDF5 file and store in local variables====*/
hdf5_read_var(&sdiv,&mdiv,SFmodel_file,&axes_ratio,&r_e,&betaphi,s_gp,mu,phi_td);


int i,j,k;
int ind;
int nx=cctkGH->cctk_lsh[0]; int ny=cctkGH->cctk_lsh[1]; int nz=cctkGH->cctk_lsh[2];

 CCTK_REAL *x_coord=0, *y_coord=0, *z_coord=0;  
 double phi_ijk;

/*============ set up initial cartesian grid data==========*/

  x_coord = x;
  y_coord = y;
  z_coord = z;


for(i=1;i<=nx;i++)
        for(j=1;j<=ny;j++)
                for(k=1;k<=nz;k++) {
			ind = i-1+nx*(j-1+ny*(k-1));
			x_i = x_grid[ind]; 
			y_j = y_grid[ind];
			z_k = z_grid[ind]; 


			grid_interp(phi_td,s_gp,mu,r_e,nx,ny,nz,x_coord,y_coord,z_coord,i,j,k,&phi_ijk,-1);
			if(z_k==0)
				printf("\n phi( %.2f, %.2f ) = %.2f ",x_i,y_j,phi_ijk);
			
			

			phi1[ind] = phi_ijk;
			phi2[ind] = 0.0;
			Kphi1[ind] = 1.0;
			Kphi2[ind] = 0.0;

					}


  CCTK_INFO("=== Ending ID_SF_Read2D initial data ===");

  /*=======================================*/
}

	

				
/* -------------------------------------------------------------------*/
