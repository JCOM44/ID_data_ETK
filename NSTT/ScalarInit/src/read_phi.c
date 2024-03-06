#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <hdf5.h>


#include "memory.h"
#include "readhdf5.h"
#include "interpolations.h"

int main() {

  /*Reading related parameters */
int m,s;
char filename[] = "/home/jolivera/python/phi.h5";
int sdiv,mdiv;
double  betaphi,axes_ratio,r_e;
double *s_gp,*mu;
double **phi;

 // Allocate memory for s_gp and mu
      /* SET UP GRID */
      s_gp=malloc((SDIV+1)*sizeof(double));
      mu=malloc((MDIV+1)*sizeof(double));

      /* ALLLOCATE MEMORY */
         
      phi = array_allocate(1,SDIV,1,MDIV);        


      /* INITIALIZE VARIABLES WITH ZERO */

      for(s=1;s<=SDIV;s++)
        for(m=1;m<=MDIV;m++) {
          phi[s][m] = 0.0e0;
        }

/*Read from the HDF5 file and store in local variables*/
hdf5_read_var(&sdiv,&mdiv,filename,&axes_ratio,&r_e,&betaphi,s_gp,mu,phi);


printf("\nThe scalard  value at s=%f and mu=%f is %f \n\n",s_gp[2],mu[2],phi[1][1]);
printf("Define new values to interpolate");

/*Define new grid in which to interpolate*/

int i,j,k;
int nx,ny,nz;  /*Number of grid points, this is going to be given by the ETK variables*/
nx=10;
ny=10;
nz=10;
double Lz,Ly,Lx; /*Total size of the domains*/
Lz = 5.0;
Ly = 5.0;
Lx = 5.0;
double *x_coord,*y_coord,*z_coord;
double phi_ijk;

printf("Values defined");

x_coord = malloc(nx * ny * nz * sizeof(double));
y_coord = malloc(nx * ny * nz * sizeof(double));
z_coord = malloc(nx * ny * nz * sizeof(double));

if (x_coord == NULL || y_coord == NULL || z_coord == NULL) {
    printf("Memory allocation failed\n");
    exit(1);
}

printf("\n Grid contains the following values\n");

for(i=1;i<=nx;i++)
        for(j=1;j<=ny;j++)
                for(k=1;k<=nz;k++) {
          x_coord[i-1+nx*(j-1+ny*(k-1))] = Lx*i/nx;
          y_coord[i-1+nx*(j-1+ny*(k-1))] = Ly*j/ny;
          z_coord[i-1+nx*(j-1+ny*(k-1))] = Lz*k/nz;
        }

grid_interp(phi,s_gp,mu,r_e,nx,ny,nz,x_coord,y_coord,z_coord,0,1,1,&phi_ijk,-1);

printf("\n phi=%f\n",phi_ijk);

return 0;
}