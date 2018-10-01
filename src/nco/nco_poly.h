#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef NCO_POLY_H
#define NCO_POLY_H 


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define CORNER_MAX 100
  

  typedef struct{
  
    double *dp_x;    /* x  vertices */
    double *dp_y;    /* y vertices */
    double *dp_xyz;  /* maybe useful for 3D stuff */ 
    double area;  
  
    int  stat;     
    int crn_nbr; /* number of vertices */
    int mem_flg; /* [flg]    */ 

  } poly_sct;   


  poly_sct *
  nco_poly_free
  (poly_sct *pl);

  poly_sct *   
  nco_poly_init
  (void);

  poly_sct *
  nco_poly_init_lst
  (int arr_nbr,
   double *dp_x_in,
   double *dp_y_in);

  poly_sct *
  nco_poly_init_crn
  (int crn_nbr_in);
  
  poly_sct*
  nco_poly_dpl
  (poly_sct *pl);

  void
  nco_poly_prn
  (int style,
   poly_sct *pl);


   poly_sct**             /* [O] [nbr] Array of poly_sct */   
   nco_poly_lst_mk(
   double *area, /* I [sr] Area of source grid */
   int *msk, /* I [flg] Mask on source grid */
   double *lat_ctr, /* I [dgr] Latitude  centers of source grid */
   double *lon_ctr, /* I [dgr] Longitude centers of source grid */
   double *lat_crn, /* I [dgr] Latitude  corners of source grid */
   double *lon_crn, /* I [dgr] Longitude corners of source grid */
   size_t grd_sz, /* I [nbr] Number of elements in single layer of source grid */
   long grd_crn_nbr, /* I [nbr] Maximum number of corners in source gridcell */
   int *pl_nbr);    /* O [nbr] size  poly_sct */  
		   
   poly_sct **
   nco_poly_lst_free(
   poly_sct **pl_lst,
   int arr_nbr);


  
#ifdef __cplusplus
} /* end extern "C" */
#endif /* __cplusplus */

#endif /* NCO_POLY_H  */
