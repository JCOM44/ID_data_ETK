#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


#define H5_USE_16_API 1
#include <hdf5.h>

#include <string.h>
#include <stdlib.h>

#include "readhdf5.h"


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
                   double **phi)
{

  hid_t       file_id;   /* file identifier */
  herr_t      status;

  hid_t       attribute_id; /* identifiers for attributes*/
  hid_t       attributeH5type;
  hid_t       datasetH5type;
  hid_t       dataset_id;  /* identifiers for dsets*/

  int         i,j;
  double      *dset_data;
  double      **var;
  int         varIndex;

  /* =========================================== */
  /* Create a new file using default properties. */
  /* =========================================== */
  /* file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT); */
  printf("Trying to open...\n");
  file_id = H5Fopen(fname, H5F_ACC_RDONLY,H5F_ACC_RDONLY);

  /* ============================================= */
  /*   Create and save the data attributes         */ 
  /* ============================================= */
  int attrblistLEN=5;
  struct {char  name[128];hid_t H5type; void * data;} attrblist[] ={
    {"axis_ratio"  , H5T_NATIVE_DOUBLE, r_ratio},
    {"re"          , H5T_NATIVE_DOUBLE, r_e},
    {"SDIV"        , H5T_NATIVE_INT, sdiv},
    {"MDIV"         ,H5T_NATIVE_INT, mdiv},
    {"alphaphi"     ,H5T_NATIVE_DOUBLE, betaphi}
  };
  /* ============================================= */
  /*  Read the data set for variables and save it  */ 
  /* ============================================= */
  int varlistLEN = 1;
  struct {char  name[128];double ** data;char description[128];} varlist[] ={
    {"/phi"     ,phi     ,"Values for RNSID variable phi     "}
  };

  for(varIndex = 0; varIndex< attrblistLEN;varIndex++) {
    if ( attrblist[varIndex].H5type == H5T_C_S1 ) {
        attributeH5type = H5Tcopy(H5T_C_S1);
  	H5Tset_size(attributeH5type,strlen(attrblist[varIndex].data)+1);
    } else {
      attributeH5type=attrblist[varIndex].H5type;
      
    }

    attribute_id  = H5Aopen(file_id,attrblist[varIndex].name,H5P_DEFAULT);
    attributeH5type = H5Aget_type(attribute_id);
    attributeH5type = H5Tget_native_type(attributeH5type,H5T_DIR_DEFAULT);

    /* -------------------------------------------------------------
    if ( attrblist[varIndex].H5type == H5T_C_S1 ) {
    } else {
      if ( attributeH5type != attrblist[varIndex].H5type) {
	printf("Attribute type missmatch %d %d %d !\n",attributeH5type,attrblist[varIndex].H5type,H5T_IEEE_F64LE); 
	attributeH5type = attrblist[varIndex].H5type;
      }
    } 
    -------------------------------------------------------------- */
    //    attribute_id = H5Acreate (file_id, 
    //	  	       attrblist[varIndex].name, attributeH5type,
    //                  dataspace_id,  H5P_DEFAULT, H5P_DEFAULT);
    status = H5Aread(attribute_id,attributeH5type, attrblist[varIndex].data);
    status = H5Aclose(attribute_id);


    if (attrblist[varIndex].H5type == H5T_NATIVE_DOUBLE)
      printf("Read attribute %s value is %5.4e\n",attrblist[varIndex].name,*((double *)attrblist[varIndex].data));
    else if (attrblist[varIndex].H5type == H5T_NATIVE_INT) 
      printf("Read attribute %s value is %d\n",attrblist[varIndex].name,*((int *)attrblist[varIndex].data));
    else if (attrblist[varIndex].H5type == H5T_C_S1) 
      printf("Read attribute %s value is %s\n",attrblist[varIndex].name,(char*)attrblist[varIndex].data);
    else 
      printf("Read attribute %s value is not knowns\n",attrblist[varIndex].name);

    /// status = H5Sclose(dataspace_id);
    if ( attrblist[varIndex].H5type == H5T_C_S1 ) {
      status= H5Tclose(attributeH5type);
      
    }
  } 

  /* ============================================= */
  /*  Read the data set for GRID points            */
  /* ============================================= */
  dataset_id =H5Dopen(file_id,"/s_qp");

  datasetH5type =H5Dget_type(dataset_id);
  datasetH5type =H5Tget_native_type(datasetH5type,H5T_DIR_DEFAULT); 
  status = H5Dread(dataset_id,datasetH5type, H5S_ALL, H5S_ALL, H5P_DEFAULT,s_qp);
  status = H5Dclose(dataset_id);

  dataset_id =H5Dopen(file_id,"/mu");
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,mu);
  status = H5Dclose(dataset_id);


  /* ============================================= */
  /*  Read the data set for variables              */
  /* ============================================= */
  for(varIndex = 0; varIndex< varlistLEN;varIndex++) {
    dset_data = malloc( *sdiv * *mdiv * sizeof(double) );
    dataset_id =H5Dopen(file_id,varlist[varIndex].name);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
    //status = H5Sclose(dataspace_id);
    var = varlist[varIndex].data;
    double minval,maxval;
    minval = 100.0;
    maxval =0.0;
    for (i = 0; i < *sdiv; i++) {
      for (j = 0; j < *mdiv; j++) {
	var[i+1][j+1] = dset_data[i*(*mdiv)+j] ;
 	if (minval > dset_data[i*(*mdiv)+j])
	  minval = dset_data[i*(*mdiv)+j];
 	if (maxval < dset_data[i*(*mdiv)+j])
	  maxval = dset_data[i*(*mdiv)+j];
      }
    }
    status = H5Dclose(dataset_id);
    printf("Read variable %s %g %g \n",varlist[varIndex].name,minval,maxval);
    free(dset_data);

  }

  printf("All done, variables and attributes read correctly!\n");
  /* Terminate access to the file. */
  status = H5Fclose(file_id); 
}
