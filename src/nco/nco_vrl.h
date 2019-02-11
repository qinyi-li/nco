#ifndef NCO_VRL_H /* Contents have not yet been inserted in current source file */
#define NCO_VRL_H


#include        <stdlib.h>
#include        <stdio.h>
#include        <math.h>

/* Personal headers */
#include "nco.h" /* netCDF Operator (NCO) definitions */
#include "nco_mmr.h"     /* Memory management */
#include "nco_omp.h"     /* OpenMP utilities */
#include "nco_rgr.h"     /* Regridding */
#include "nco_sld.h"     /* Swath-Like Data */
#include "nco_sng_utl.h" /* String utilities */
#include "nco_poly.h"    /* poly sct stuff */


/* Dimension of points */
#define DIM 2

#define DSIGMA 1.0e-14d

/* define minimium area in AreaSign (cross-product) */
#define DAREA  1.0e-28d  

#define VP_MAX    1000            /* Max # of pts in polygon */

#define ARC_MIN_LENGTH (1.0e-20d)

/* this is 1.0e-20 * PI / 180.0 */
#define ARC_MIN_LENGTH_RAD (1.0e-15d)

/* if true then longitude 0-360 */
/* we need this to convert 3D back to 2D */
#define IS_LON_360 (1)

#define DEBUG_VRL (1)


#ifdef __cplusplus
/* Use C-bindings so C++-compiled and C-compiled libraries are compatible */
extern "C" {
#endif /* !__cplusplus */



typedef enum { Pin, Qin, Unknown } tInFlag;
typedef int     tPointi[DIM];   /* type integer point */
typedef double  tPointd[DIM];   /* type double point */
typedef double  tPointds[5];    /* type spherical double point */

typedef tPointi tPolygoni[VP_MAX]; /* type integer polygon */
typedef tPointd tPolygond[VP_MAX]; /* type integer polygon */
typedef tPointds tPolygonds[VP_MAX]; /* 3D sperical coords */


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



// void    ClosePostscript( void );
// void	PrintSharedSeg( tPointd p, tPointd q );
void    PrintPoly( tPolygond P, int n );
const char * prnInFlag(tInFlag in);

int     Advance( int a, int *aa, int n, int inside, tPointi v );
// void	OutputPolygons( tPolygond P, tPolygond Q, int n, int m );
// int     ReadPoly( tPolygond P );

/*-------------------------------------------------------------------*/
/* spherical methods */

int  sConvexIntersect( tPolygonds P, tPolygonds Q, tPolygonds R, int n, int m, int *r );

int  snewConvexIntersect( poly_sct *sP, poly_sct * sQ, poly_sct *sR, int *r);


char    sSegSegInt( tPointds a, tPointds b, tPointds c, tPointds d, tPointds p, tPointds q );

int sLHS(tPointds Pi, tPointds Qi );
nco_bool sFace( int iLHS, int iRHS, int jRHS  );

double  sDot( tPointds a, tPointds b );
double  sCross(tPointds a, tPointds b, tPointds c);
double sRadius(tPointds a);

double  sxCross( tPointds a, tPointds b, tPointds c );
void    sAdi(tPointds a, tPointds b );
void    sph2crt(tPointds a,  double *lon, double *lat, nco_bool bDeg);
void    crt2sph(tPointd a, tPointds b);
void    sphAddcrt(tPointds ds);
void    sAddPoint( tPolygonds R , int *r, tPointds P);

void sphAddPoint(double **sph, int *r, double *P );

double latCorrect( double lat1, double lon1, double lon2  );

void getLatCorrect_old(tPointds a, tPointds b, double *dp_min, double *dp_max );

void getLatCorrect(double lon1, double lat1, double lon2, double lat2, double *dp_min, double *dp_max, nco_bool bDeg);

nco_bool iBetween(double a, double b, double x  );
nco_bool sLatLonBetween(tPointds a, tPointds b, tPointds x);
char    sParallelDouble( tPointds a, tPointds b, tPointds c, tPointds d, tPointds p, tPointds q );
void prnPoint(const char *sMsg, tPointds p, int style, nco_bool bRet );

nco_bool sConvex(tPolygonds sP, int np);
void sPrintPoly(tPolygonds sR, int r, int istyle);

nco_bool sPointInPolygon( tPolygonds sP, int n, tPointds pControl, tPointds pVertex);

void setStaticGlobals(double lon_min_rad, double lon_max_rad, double lat_min_rad, double lat_max_rad   );

/*-------------------------------------------------------------------*/



#ifdef __cplusplus
} /* end extern "C" */
#endif /* !__cplusplus */

#endif /* NCO_VRL_H */
