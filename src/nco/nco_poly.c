#include "nco_poly.h"


poly_sct *
nco_poly_free
(poly_sct *pl)
{

  /* mem flag set -so pointers from external array */
  if( pl->mem_flg ==1 )
  {
    pl->dp_x=(double*)NULL_CEWI;
    pl->dp_y=(double*)NULL_CEWI;

  }
  else
  {  
    pl->dp_x=(double*)nco_free(pl->dp_x);
    pl->dp_y=(double*)nco_free(pl->dp_y);
  }

  if(pl->dp_xyz)
    pl->dp_xyz=(double*)nco_free(pl->dp_xyz);

  
    
}  


poly_sct *   
nco_poly_init
(void)
{  
  poly_sct *pl;

  pl=(poly_sct*)nco_malloc( sizeof(poly_sct));

  pl->dp_x=(double*)NULL_CEWI;
  pl->dp_y=(double*)NULL_CEWI;
  pl->dp_xyz=(double*)NULL_CEWI;

  pl->dp_x_minmax[0]=0.0;
  pl->dp_x_minmax[1]=0.0;

  pl->dp_y_minmax[0]=0.0;
  pl->dp_y_minmax[1]=0.0;

  
  pl->stat=0;
  pl->area=0.0;
  pl->crn_nbr=0;
  pl->mem_flg=0;

  return pl;
}

poly_sct*
nco_poly_dpl
(poly_sct *pl)
{

  poly_sct *pl_cpy;
  int crn_nbr_in;
  
  pl_cpy=nco_poly_init();

  crn_nbr_in=pl->crn_nbr;

  pl_cpy->stat=pl->stat;
  pl_cpy->area=pl->area;
  pl_cpy->crn_nbr=crn_nbr_in;

  /* mem flag is ALWAYS 0 for a copy  */
  pl_cpy->mem_flg=0;

  pl_cpy->dp_x=(double*)nco_malloc((size_t)crn_nbr_in* sizeof(double));
  pl_cpy->dp_y=(double*)nco_malloc((size_t)crn_nbr_in* sizeof(double));

  memcpy(pl_cpy->dp_x, pl->dp_x, (size_t)crn_nbr_in* sizeof(double));
  memcpy(pl_cpy->dp_y, pl->dp_y, (size_t)crn_nbr_in* sizeof(double));  

  pl->dp_x_minmax[0];
  pl->dp_x_minmax[1]=0.0;

  pl->dp_y_minmax[0]=0.0;
  pl->dp_y_minmax[1]=0.0;


  pl_cpy->dp_x_minmax[0] = pl->dp_x_minmax[0];
  pl_cpy->dp_x_minmax[1] = pl->dp_x_minmax[1];

  pl_cpy->dp_y_minmax[0] = pl->dp_y_minmax[0];
  pl_cpy->dp_y_minmax[1] = pl->dp_y_minmax[1];


  


  
  return pl_cpy;
  
} 
  
poly_sct *
nco_poly_init_crn
(int crn_nbr_in)
{
  poly_sct *pl;
  pl=nco_poly_init();

  pl->crn_nbr=crn_nbr_in;

  pl->dp_x=(double*)nco_calloc((size_t)crn_nbr_in, sizeof(double));
  pl->dp_y=(double*)nco_calloc((size_t)crn_nbr_in, sizeof(double));

  pl->mem_flg=0;
  
  return pl;
}
  

poly_sct *
nco_poly_init_lst
(int arr_nbr,
 double *dp_x_in,
 double *dp_y_in)
{
 int idx;
 int sz;

 poly_sct *pl;


 /* less than a triangle */
 if (arr_nbr <3 )
   return (poly_sct*)NULL_CEWI;   


 /* check repeated points at end of arrray - nb must be an exact match */
 for(idx=1; idx<arr_nbr; idx++ )
   if( dp_x_in[idx] == dp_x_in[idx-1] && dp_y_in[idx] == dp_y_in[idx-1] )
     break;

 if(idx < 3 )
     return (poly_sct*)NULL_CEWI;   

 /* we have at least a triangle */ 
 pl=nco_poly_init();
 
 /* dont free  pointers */
 pl->mem_flg=1;
 pl->crn_nbr=idx;
 
 pl->dp_x=dp_x_in;
 pl->dp_y=dp_y_in;
 
 
 
 return pl;
 

}  

void nco_poly_add_minmax
(poly_sct *pl)
{  
  
  int idx;
  int sz;


  sz=pl->crn_nbr; 
  
  pl->dp_x_minmax[0]=DBL_MAX;
  pl->dp_x_minmax[1]=-DBL_MAX;

  pl->dp_y_minmax[0]=DBL_MAX;
  pl->dp_y_minmax[1]=-DBL_MAX;


  
  for(idx=0; idx<sz;idx++)
  {
    /* min */
    if( pl->dp_x[idx] < pl->dp_x_minmax[0] )
      pl->dp_x_minmax[0] = pl->dp_x[idx]; 

    /* max */
    if( pl->dp_x[idx] > pl->dp_x_minmax[1] )
          pl->dp_x_minmax[1] = pl->dp_x[idx];  

    /* min */
    if( pl->dp_y[idx] < pl->dp_y_minmax[0] )
      pl->dp_y_minmax[0] = pl->dp_y[idx]; 

    /* max */
    if( pl->dp_y[idx] > pl->dp_y_minmax[1] )
          pl->dp_y_minmax[1] = pl->dp_y[idx];  

    
    
  }

  return; 
  

}  




void
nco_poly_prn
(int style,
 poly_sct *pl)
{
  int idx;


  switch(style){ 

    case 0:
      (void)fprintf(stdout,"\n%s: crn_nbr=%d stat=%d mem_flg=%d area=%f\n", nco_prg_nm_get(), pl->crn_nbr, pl->stat, pl->mem_flg, pl->area);      
      (void)fprintf(stdout,"dp_x ");
      for(idx=0; idx<pl->crn_nbr; idx++)
	(void)fprintf(stdout,"%20.14f, ",pl->dp_x[idx]);
      (void)fprintf(stdout,"\n");		  

      (void)fprintf(stdout,"dp_y ");
      for(idx=0; idx<pl->crn_nbr; idx++)
	(void)fprintf(stdout,"%20.14f, ",pl->dp_y[idx]);
      (void)fprintf(stdout,"\n");

      (void)fprintf(stdout,"min/max x( %g, %g) y(%g %g)\n", pl->dp_x_minmax[0], pl->dp_x_minmax[1], pl->dp_y_minmax[0], pl->dp_y_minmax[1]);       
      
      break;

   case 1:  
   default:
     (void)fprintf(stdout,"%s: crn_nbr=%d\n", nco_prg_nm_get(), pl->crn_nbr);
     
     for(idx=0; idx<pl->crn_nbr; idx++)
        (void)fprintf(stdout,"{ %20.14f, %20.14f }\n",pl->dp_x[idx], pl->dp_y[idx]);

     break;
  }

  return;
     
}


poly_sct*
nco_poly_do_vrl(
poly_sct *pl_in,
poly_sct *pl_out){

  
 poly_sct *pl_vrl;
  

  /* for now just copy pl_in so  we can test other functions */
 pl_vrl=nco_poly_dpl( pl_in);

 return pl_vrl;
  
}  





/************************ functions that manipulate lists of polygons ****************************************************/

poly_sct **             /* [O] [nbr]  size of array */   
nco_poly_lst_mk(
double *area, /* I [sr] Area of source grid */
int *msk, /* I [flg] Mask on source grid */
double *lat_ctr, /* I [dgr] Latitude  centers of source grid */
double *lon_ctr, /* I [dgr] Longitude centers of source grid */
double *lat_crn, /* I [dgr] Latitude  corners of source grid */
double *lon_crn, /* I [dgr] Longitude corners of source grid */
size_t grd_sz, /* I [nbr] Number of elements in single layer of source grid */
long grd_crn_nbr, /* I [nbr] Maximum number of corners in source gridcell */
int *pl_nbr)
{

    int idx=0;
    int idx_cnt=0;

    double *lat_ptr=lat_crn;
    double *lon_ptr=lon_crn;
    poly_sct *pl;
    poly_sct **pl_lst;


    pl_lst=(poly_sct**)nco_malloc( (size_t)grd_sz * sizeof (poly_sct*) );
    
    printf("About to print poly sct\n");
    for(idx=0;idx<grd_sz; idx++) 
    {
      /* check mask and area */
      if( msk[idx]==0 || area[idx] == 0.0d)
	continue;

      
      pl=nco_poly_init_lst( grd_crn_nbr, lon_ptr, lat_ptr);
      lon_ptr+=(size_t)grd_crn_nbr;
      lat_ptr+=(size_t)grd_crn_nbr;

      /* if poly is less  than a triangle then  null is returned*/
      if(!pl)
	continue;

      /* add min max */
      nco_poly_add_minmax(pl);
      
      pl_lst[idx_cnt++]=pl;     
      
    }
    
    /* realloc if ncessary */
    if(idx_cnt< grd_sz)
      pl_lst=(poly_sct**)nco_realloc( pl_lst, (size_t)idx_cnt * sizeof (poly_sct*) );
    
    *pl_nbr=idx_cnt;
     
    return pl_lst;
 
}  

poly_sct **
nco_poly_lst_free(
poly_sct **pl_lst,
int arr_nbr)
{
  int idx;

   for(idx=0; idx<arr_nbr; idx++)
     pl_lst[idx]=nco_poly_free(pl_lst[idx]);

   pl_lst=(poly_sct**)nco_free(pl_lst);

   return pl_lst;

}  


void
nco_poly_set_priority(
int nbr_lst,		      
KDPriority *list){		      

int idx;

 for(idx=0;idx<nbr_lst;idx++){

   list[idx].dist = 1.1;
   list[idx].elem = (KDElem*)NULL;
 }  

 return ; 

}
  


poly_sct **
nco_poly_mk_vrl_lst(   /* create overlap mesh */
 poly_sct ** pl_lst_in,
 int pl_cnt_in,
 poly_sct ** pl_lst_out,
 int pl_cnt_out,
 int *pl_cnt_vrl_ret){

/* just duplicate output list to overlap */

 int idx;
 int jdx;
 int sz;
 int max_nbr_vrl=200; 
 int pl_cnt_vrl=0;
 
 char *chr_ptr;
 char fnc_nm[]="nco_poly_mk_vrl()";  

 kd_box size;

 poly_sct ** pl_lst_vrl=NULL_CEWI;
 
 KDElem *my_elem;
 KDTree *rtree;

 KDPriority *list;

  list = (KDPriority *)nco_calloc(sizeof(KDPriority),(size_t)max_nbr_vrl); 
 
  printf("INFO - entered function nco_poly_mk_vrl\n"); 
 
  /* create kd_tree from output polygons */
  rtree=kd_create();


   /* populate kd_tree */
  for(idx=0 ; idx<pl_cnt_out;idx++){
    
       
    my_elem=(KDElem*)nco_calloc((size_t)1,sizeof (KDElem) );
   
    size[KD_LEFT]  =  pl_lst_out[idx]->dp_x_minmax[0];
    size[KD_RIGHT] =  pl_lst_out[idx]->dp_x_minmax[1];

    size[KD_BOTTOM] = pl_lst_out[idx]->dp_y_minmax[0];
    size[KD_TOP]    = pl_lst_out[idx]->dp_y_minmax[1];    

    //chr_ptr=(char*)pl_lst_out[idx];

    kd_insert(rtree, (kd_generic)pl_lst_out[idx], size, (char*)my_elem);

  }

  /* rebuild tree for faster access */
  kd_rebuild(rtree);

  //printf("about to output kd_tree\n");
  // kd_print(rtree);


  
/* start main loop over input polygons */ 
 for(idx=0 ; idx<pl_cnt_in ;idx++ )
 { 
   int cnt_vrl=0;

   (void)nco_poly_set_priority(max_nbr_vrl,list); 
   /* get bounds of polygon in */   
    size[KD_LEFT]  =  pl_lst_in[idx]->dp_x_minmax[0];
    size[KD_RIGHT] =  pl_lst_in[idx]->dp_x_minmax[1];

    size[KD_BOTTOM] = pl_lst_in[idx]->dp_y_minmax[0];
    size[KD_TOP]    = pl_lst_in[idx]->dp_y_minmax[1];    

    /* find overlapping polygons */
    cnt_vrl=kd_nearest_intersect(rtree, size, max_nbr_vrl,list );

    nco_poly_prn(0, pl_lst_in[idx] );
    fprintf(stdout,"%s: number of overlaps=%d -overlapping polygons to follow\n/**************************/\n"  , fnc_nm,  cnt_vrl);

   
    /* for testing purposes just use first overlap polygon */
    cnt_vrl= ( cnt_vrl  ? 1: 0);
    
    for(jdx=0; jdx <cnt_vrl ;jdx++){

      poly_sct *pl_vrl=(poly_sct*)NULL_CEWI;	 
      poly_sct *pl_out=(poly_sct*)list[jdx].elem->item;           ;

      nco_poly_prn(0, pl_out);           
     
  
      pl_vrl=nco_poly_do_vrl(pl_lst_in[idx], pl_out);

      if(pl_vrl){
	pl_lst_vrl=(poly_sct**)nco_realloc(pl_lst_vrl, sizeof(poly_sct*) * (pl_cnt_vrl+1));
	pl_lst_vrl[pl_cnt_vrl]=pl_vrl;
	pl_cnt_vrl++;
      }

    }


    
 }   


 kd_destroy(rtree,NULL);

 list = (KDPriority *)nco_free(list);

 /* return size of list */
 *pl_cnt_vrl_ret=pl_cnt_vrl;

 
 return pl_lst_vrl;

}  


