
#define READHDF5_H


/*#define H5Acreate_vers 2*/
/*#define H5Dcreate_vers 2*/
#define H5Dopen_vers 1 


void hdf5_read_varSF(int *sdiv, int *mdiv,
		   char fname[],
		   /* ------------------- */
                   double *r_ratio,
                   double *r_e,	
		   /* ------------------- */
                   double *betaphi,
       /* ------------------- */
                   double *s_qp,  /*SDIV+1*/
                   double *mu,    /*MDIV+1*/
		   /* ------------------- */
                   double **phi);


