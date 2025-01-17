#include<stdio.h>
#include"memory.h"
#include<stdlib.h>
void AllocationErrorSF(char error_text[])
{
	fprintf(stderr,"RNS RUNTIME ERROR: %s\n Exiting the system ",error_text);
	exit(1);
}

double **array_allocateSF(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+ARRAY_END)*sizeof(double*)));
	if (!m) AllocationError("allocation failure 1 in matrix()");
	m += ARRAY_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+ARRAY_END)*sizeof(double)));
	if (!m[nrl]) AllocationError("allocation failure 2 in matrix()");
	m[nrl] += ARRAY_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}
 
 

void array_freeSF(double **m, long nrl, long ncl)
/* free a double matrix allocated by array_allocate() */
{
	free((FREE_ARG) (m[nrl]+ncl-ARRAY_END));
	free((FREE_ARG) (m+nrl-ARRAY_END));
}
 
void check_status_error(int status)
/* Control variable to monitor opening and closing of HDF5 files */
{
	if (!status) printf("ERROR WHEN ACCESSING...");
}
