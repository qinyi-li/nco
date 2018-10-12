#ifndef NCO_VRL_H /* Contents have not yet been inserted in current source file */
#define NCO_VRL_H


#include        <stdlib.h>
#include        <stdio.h>
#include        <math.h>

/* Personal headers */
#include "nco.h" /* netCDF Operator (NCO) definitions */
#include "nco_mmr.h" /* Memory management */
#include "nco_omp.h" /* OpenMP utilities */
#include "nco_rgr.h" /* Regridding */
#include "nco_sld.h" /* Swath-Like Data */
#include "nco_sng_utl.h" /* String utilities */



/* Dimension of points */
#define DIM 2

#define DSIGMA 1.0e-14d

/* define minimium area in AreaSign (cross-product) */
#define DAREA  1.0e-28d  

#define VP_MAX    1000            /* Max # of pts in polygon */



#ifdef __cplusplus
/* Use C-bindings so C++-compiled and C-compiled libraries are compatible */
extern "C" {
#endif /* !__cplusplus */



  

  

typedef enum { Pin, Qin, Unknown } tInFlag;
typedef int     tPointi[DIM];   /* type integer point */
typedef double  tPointd[DIM];   /* type double point */


typedef tPointi tPolygoni[VP_MAX]; /* type integer polygon */
typedef tPointd tPolygond[VP_MAX]; /* type integer polygon */





/*---------------------------------------------------------------------
Function prototypes.
---------------------------------------------------------------------*/

int    ConvexIntersect( tPolygond P, tPolygond Q, tPolygond R, int n, int m, int *r );
char    SegSegInt( tPointd a, tPointd b, tPointd c, tPointd d, tPointd p, tPointd q );
char    ParallelInt( tPointd a, tPointd b, tPointd c, tPointd d, tPointd p, tPointd q );
int	AreaSign( tPointd a, tPointd b, tPointd c );
nco_bool  Between( tPointd a, tPointd b, tPointd c );

double  Dot( tPointd a, tPointd b );
void    SubVec( tPointd a, tPointd b, tPointd c );
void    Adi( tPointd p, tPointd a );
void    AddPoint( tPolygond R, int *r, tPointd p); 
nco_bool  Collinear( tPointd a, tPointd b, tPointd c );

nco_bool  LeftOn( tPointd a, tPointd b, tPointd c );
nco_bool  Left( tPointd a, tPointd b, tPointd c );
tInFlag InOut( tPointd p, tInFlag inflag, int aHB, int bHA );



void    ClosePostscript( void );
void	PrintSharedSeg( tPointd p, tPointd q );
void    PrintPoly( tPolygond P, int n );
//void    PrintPolyd( int r, tPolygond R );


int     Advance( int a, int *aa, int n, int inside, tPointi v );
void	OutputPolygons( tPolygond P, tPolygond Q, int n, int m );
int     ReadPoly( tPolygond P );

/*-------------------------------------------------------------------*/

#ifdef __cplusplus
} /* end extern "C" */
#endif /* !__cplusplus */

#endif /* NCO_VRL_H */
