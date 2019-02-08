/*
This code is described in "Computational Geometry in C" (Second Edition),
Chapter 7.  It is not written to be comprehensible without the
explanation in that book.

Written by Joseph O'Rourke.
Last modified: December 1997
Questions to orourke@cs.smith.edu.
--------------------------------------------------------------------
This code is Copyright 1997 by Joseph O'Rourke.  It may be freely
redistributed in its entirety provided that this copyright notice is
not removed.
--------------------------------------------------------------------
*/

#include "nco_vrl.h"

/* global variables for latitude, longitude in RADIANS
   these may be set in nco_poly.c or
   should be safe with OPenMP  ? */

static double LAT_MIN_RAD;
static double LAT_MAX_RAD;

static double LON_MIN_RAD;
static double LON_MAX_RAD;



/*---------------------------------------------------------------------
---------------------------------------------------------------------*/
int ConvexIntersect( tPolygond P, tPolygond Q, tPolygond R, int n, int m, int *r )
{
   int lcl_dbg=0; 
   nco_bool FirstPoint=True;    /*s this the first point? (used to initialize).*/   
   int     a=0, b=0;            /* indices on P and Q (resp.) */
   int     a1, b1;              /* a-1, b-1 (resp.) */
   int     aa=0, ba=0;          /* # advances on a & b indices (after 1st inter.) */
   int     cross;               /* sign of z-component of A x B */
   int     bHA, aHB;            /* b in H(A); a in H(b). */
   int     code;                /* SegSegInt return code. */
   
   tPointd A, B;                /* directed edges on P and Q (resp.) */
   tPointd Origin = {0.0,0.0};  /* (0,0) */
   tPointd p0;                  /* The first point. */
   tPointd p;                   /* double point of intersection */
   tPointd q;                   /* second point of intersection */

   
   tInFlag inflag = Unknown; /* {Pin, Qin, Unknown}: which inside */

   do {

     
      /* Computations of key variables. */
      a1 = (a + n - 1) % n;
      b1 = (b + m - 1) % m;

      SubVec( P[a], P[a1], A );
      SubVec( Q[b], Q[b1], B );
      cross = AreaSign( Origin, A, B );
      aHB   = AreaSign( Q[b1], Q[b], P[a] );
      bHA   = AreaSign( P[a1], P[a], Q[b] );

      /* If A & B intersect, update inflag. */
      code = SegSegInt( P[a1], P[a], Q[b1], Q[b], p, q );

      if(DEBUG_VRL)
        (void)fprintf(stdout, "%s: cross=%d, aHB=%d, bHA=%d code = %c\n", nco_prg_nm_get(),cross, aHB, bHA, code );

      if ( code == '1' || code == 'v' ) {
         if ( inflag == Unknown && FirstPoint ) {
	    aa = 0;
	    ba = 0;
            FirstPoint = False ;
	    Adi(p0,p);
	    AddPoint(R,r, p0);  
         }

         inflag = ( aHB >0 ? Pin : bHA >0 ? Qin : inflag );
	 
	 AddPoint(R,r, p);

         if(DEBUG_VRL)
	      (void)fprintf(stdout, "%s: InOut sets inflag=%d\n", nco_prg_nm_get(),  inflag);
      }

      /*-----Advance rules-----*/

      /* Special case: A & B overlap and oppositely oriented. */
      if ( code == 'e' && Dot( A, B ) < 0  )
      {	
	   AddPoint(R,r,p );
      	   AddPoint(R,r,q );
           return EXIT_FAILURE;  
      }
	   
      /* Special case: A & B parallel and separated. */
      if ( (cross == 0) && ( aHB < 0) && ( bHA < 0 ) )
	{

          if(DEBUG_VRL)
              (void)fprintf(stdout, "%s: P and Q are disjoint\n", nco_prg_nm_get());
	  
	  return EXIT_FAILURE;
      }
      /* Special case: A & B collinear. */
      else if ( (cross == 0) && ( aHB == 0) && ( bHA == 0 ) )
      {
            /* Advance but do not output point. */
            if ( inflag == Pin )
	    {
		// b = Advance( b, &ba, m, inflag == Qin, Q[b] );
	        b++; ba++;
	    } 
            else
	    {
		//a = Advance( a, &aa, n, inflag == Pin, P[a] );
	        a++; aa++;
	    } 
      }
      /* Generic cases. */
      else if ( cross >= 0 )
      {
         if ( bHA > 0)
	 {   //a = Advance( a, &aa, n, inflag == Pin, P[a] );
	   if( inflag == Pin ) AddPoint(R,r, P[a]);

	   a++; aa++;

	 }  
         else
	 {  
	   // b = Advance( b, &ba, m, inflag == Qin, Q[b] );
           if( inflag == Qin) AddPoint(R,r, Q[b]);

	   b++; ba++;
	 }   
      }
      
      else /* if ( cross < 0 ) */
      {
         if ( aHB > 0)
	 {  
	   //b = Advance( b, &ba, m, inflag == Qin, Q[b] );
           if(inflag == Qin ) AddPoint(R,r, Q[b]);
	   
           b++; ba++;
	 }   
         else
	 {  
	   //a = Advance( a, &aa, n, inflag == Pin, P[a] );
	   if( inflag == Pin ) AddPoint(R,r, P[a]);
	   
	   a++; aa++;
	 }    
      }

      /* normalize counters */
      a%=n;
      b%=m;

     if(DEBUG_VRL)
       (void)fprintf(stdout, "%s: Before Advances:a=%d, b=%d; aa=%d, ba=%d; inflag=%d\n", nco_prg_nm_get(),   a, b, aa, ba, inflag);


   /* Quit when both adv. indices have cycled, or one has cycled twice. */
   } while ( ((aa < n) || (ba < m)) && (aa < 2*n) && (ba < 2*m) );

   if ( !FirstPoint ) 
   {
      if(DEBUG_VRL)
         (void)fprintf(stdout, "%s: no points output\n", nco_prg_nm_get());
      
      return EXIT_FAILURE;

   }
   
   /* Deal with special cases: not implemented. */
   if ( inflag == Unknown)
   {

      if(DEBUG_VRL)
         (void)fprintf(stdout, "The boundaries of P and Q do not cross.\n", nco_prg_nm_get());
      
      return EXIT_FAILURE;



   }
   
   return EXIT_SUCCESS;
   
}

/*---------------------------------------------------------------------
SegSegInt: Finds the point of intersection p between two closed
segments ab and cd.  Returns p and a char with the following meaning:
   'e': The segments collinearly overlap, sharing a point.
   'v': An endpoint (vertex) of one segment is on the other segment,
        but 'e' doesn't hold.
   '1': The segments intersect properly (i.e., they share a point and
        neither 'v' nor 'e' holds).
   '0': The segments do not intersect (i.e., they share no points).
Note that two collinear segments that share just one point, an endpoint
of each, returns 'e' rather than 'v' as one might expect.
---------------------------------------------------------------------*/
char SegSegInt( tPointd a, tPointd b, tPointd c, tPointd d, tPointd p, tPointd q )
{
   double  s, t;       /* The two parameters of the parametric eqns. */
   double num, denom;  /* Numerator and denoninator of equations. */
   char code = '?';    /* Return char characterizing intersection. */

   /*printf("%%SegSegInt: a,b,c,d: (%d,%d), (%d,%d), (%d,%d), (%d,%d)\n",
	a[0],a[1], b[0],b[1], c[0],c[1], d[0],d[1]);*/

   denom = a[0] * ( d[1] - c[1] ) +
           b[0] * ( c[1] - d[1] ) +
           d[0] * ( b[1] - a[1] ) +
           c[0] * ( a[1] - b[1] );

   /* If denom is zero, then segments are parallel: handle separately. */
   if (denom == 0.0)
      return  ParallelInt(a, b, c, d, p, q);

   num =    a[0] * ( d[1] - c[1] ) +
            c[0] * ( a[1] - d[1] ) +
            d[0] * ( c[1] - a[1] );
   
   if ( num == 0.0 || num == denom )
     code = 'v';
   
   s = num / denom;
   /*printf("num=%lf, denom=%lf, s=%lf\n", num, denom, s);*/

   num = -( a[0] * ( c[1] - b[1] ) +
            b[0] * ( a[1] - c[1] ) +
            c[0] * ( b[1] - a[1] ) );
   
   if ( num == 0.0 || num == denom )
     code = 'v';
   
   t = num / denom;
   /*printf("num=%lf, denom=%lf, t=%lf\n", num, denom, t);*/

   if(  s >0.0 && s < 1.0  &&  t >0.0 && t < 1.0  )
     code = '1';
   else  if(  s <0.0 || s > 1.0 || t <0.0 || t > 1.0  )
     code = '0';

   p[0] = a[0] + s * ( b[0] - a[0] );
   p[1] = a[1] + s * ( b[1] - a[1] );

   return code;
}
char   ParallelInt( tPointd a, tPointd b, tPointd c, tPointd d, tPointd p, tPointd q )
{
/*   
   printf("ParallelInt: a,b,c,d: (%d,%d), (%d,%d), (%d,%d), (%d,%d)\n",
	a[0],a[1], b[0],b[1], c[0],c[1], d[0],d[1]);
*/
  /* Check if collinear */
   if ( AreaSign( a, b, c) == 0  )
      return '0';

   if ( Between( a, b, c ) && Between( a, b, d ) ) {
      Adi( p, c );
      Adi( q, d );
      return 'e';
   }
   if ( Between( c, d, a ) && Between( c, d, b ) ) {
      Adi( p, a );
      Adi( q, b );
      return 'e';
   }
   if ( Between( a, b, c ) && Between( c, d, b ) ) {
      Adi( p, c );
      Adi( q, b );
      return 'e';
   }
   if ( Between( a, b, c ) && Between( c, d, a ) ) {
      Adi( p, c );
      Adi( q, a );
      return 'e';
   }
   if ( Between( a, b, d ) && Between( c, d, b ) ) {
      Adi( p, d );
      Adi( q, b );
      return 'e';
   }
   if ( Between( a, b, d ) && Between( c, d, a ) ) {
      Adi( p, d );
      Adi( q, a );
      return 'e';
   }
   return '0';
}

/*---------------------------------------------------------------------
Returns the dot product of the two input vectors.
---------------------------------------------------------------------*/
double  Dot( tPointd a, tPointd b )
{
    int i;
    double sum = 0.0;

    for( i = 0; i < DIM; i++ )
       sum += a[i] * b[i];

    return  sum;
}

/*---------------------------------------------------------------------
a - b ==> c.
---------------------------------------------------------------------*/
void SubVec( tPointd a, tPointd b, tPointd c )
{
   int i;

   for( i = 0; i < DIM; i++ )
      c[i] = a[i] - b[i];
}


void  Adi( tPointd p, tPointd a )
{
  p[0]=a[0];
  p[1]=a[1];
  /*
   int i;
   for ( i = 0; i < DIM; i++ )
      p[i] = a[i];
  */
}

/*---------------------------------------------------------------------
Prints out the double point of intersection, and toggles in/out flag.
---------------------------------------------------------------------*/
tInFlag InOut( tPointd p, tInFlag inflag, int aHB, int bHA )
{
  //printf("%8.2lf %8.2lf lineto\n", p[0], p[1] );

   /* Update inflag. */
   if      ( aHB > 0)
      return Pin;
   else if ( bHA > 0)
      return Qin;
   else    /* Keep status quo. */
      return inflag;
}

/*---------------------------------------------------------------------
   Advances and prints out an inside vertex if appropriate.
---------------------------------------------------------------------*/
int     Advance( int a, int *aa, int n, nco_bool inside, tPointi v )
{
   if ( inside )
      printf("%5d    %5d    lineto\n", v[0], v[1] );
   (*aa)++;
   return  (a+1) % n;
}


int AreaSign( tPointd a, tPointd b, tPointd c )
{
    double area2;

    area2 = ( b[0] - a[0] ) * ( c[1] - a[1] ) -
            ( c[0] - a[0] ) * ( b[1] - a[1] );

    /* The area should be an integer. */
    /*
    if      ( area2 >  0.5 ) return  1;
    else if ( area2 < -0.5 ) return -1;
    else                     return  0;
    */
    if      ( area2 >  DAREA ) return  1;
    else if ( area2 < -DAREA ) return -1;
    else                       return  0;

    
}



/*
   Returns true iff c is strictly to the left of the directed
   line through a to b.
*/
nco_bool Left( tPointd a, tPointd b, tPointd c )
{
        return  AreaSign( a, b, c ) > 0;
}

nco_bool LeftOn( tPointd a, tPointd b, tPointd c )
{
        return  AreaSign( a, b, c ) >= 0;
}

nco_bool Collinear( tPointd a, tPointd b, tPointd c )
{
        return  AreaSign( a, b, c ) == 0;
}


/*---------------------------------------------------------------------
Returns TRUE iff point c lies on the closed segement ab.
Assumes it is already known that abc are collinear.
---------------------------------------------------------------------*/
nco_bool Between( tPointd a, tPointd b, tPointd c )
{
   tPointd      ba, ca;

   /* If ab not vertical, check betweenness on x; else on y. */
   if ( a[0] != b[0] )
     return (a[0] <= c[0] && c[0] <= b[0])  || (a[0] >= c[0] && c[0] >= b[0] ) ;
   else
     return (a[1] <= c[1] && c[1] <= b[1]) || (a[1] >= c[1] && c[1] >= b[1] ) ;

}




void AddPoint( tPolygond R, int *r, tPointd P)
{

  
  /* only add  point if its distance from  from previous point is more than DSIGMA */ 
  if ( *r == 0  ||  (  pow( (R[*r-1][0] - P[0]),2.0 )  + pow( (R[*r-1][1] - P[1]),2.0) > DAREA  ) )
  {  

    R[*r][0] = P[0];
    R[*r][1] = P[1];
    (*r)++;     
    
  }
    
}

const char * prnInFlag(tInFlag in)
{
  if(in == Pin)
    return "Pin";
  else if(in == Qin)
    return "Qin";
  else if(in == Unknown)
    return "Unknown";
}


/*---------------------------------------------------------------------
Polygon I/O functions  
---------------------------------------------------------------------*/



void PrintPoly(tPolygond R, int r)
{
  int idx;



   printf("Polygon R:\n");
   
   for( idx = 0; idx < r; idx++ )
      printf("%20.14f %20.14f\n", R[idx][0], R[idx][1]);
   
   printf("End Polygon\n");


}

void sPrintPoly(tPolygonds sR, int r, int istyle)
{
  int idx;


  printf("\nSpherical Polygon\n");

  for( idx = 0; idx < r; idx++ )
    prnPoint(">",sR[idx],istyle,True );
    //printf("%20.14f %20.14f\n", sR[idx][0], sR[idx][1]);

  printf("End Polygon\n");


}




/* spherical functions */
int sConvexIntersect( tPolygonds P, tPolygonds Q, tPolygonds R, int n, int m, int *r ) {

   nco_bool flg_dbg=True;

   nco_bool qpFace = False;
   nco_bool pqFace = False;
   nco_bool isGeared = False;

   int numIntersect=0;

   int a = 0, a1 = 0, aa=0;
   int b = 0, b1 = 0, bb=0;


   int ipqLHS = 0;
   int ip1qLHS = 0 ;
   int iqpLHS = 0;
   int iq1pLHS = 0 ;

   double nx1;
   double nx2;
   double nx3;

   char code='0';

   tPointds p;
   tPointds q;

   tInFlag inflag= Unknown;

   if(DEBUG_VRL)
     fprintf(stdout, "%s: just entered sConvexIntersect()\n", nco_prg_nm_get() );


   do{


      a1 = (a + n - 1) % n;
      b1 = (b + m - 1) % m;


      tPointds Pcross;
      tPointds Qcross;
      tPointds Xcross;

      nx1=sCross(P[a1], P[a], Pcross);
      nx2=sCross(Q[b1], Q[b], Qcross);

      nx3=sCross(Pcross,Qcross,Xcross);


      ipqLHS = sLHS(P[a], Qcross);
      ip1qLHS = sLHS(P[a1], Qcross);


      /* imply rules facing if 0 */

      if(ipqLHS==0 && ip1qLHS!=0)
         ipqLHS=ip1qLHS*-1;
      else if( ipqLHS != 0 && ip1qLHS == 0 )
         ip1qLHS=ipqLHS*-1;


      iqpLHS = sLHS(Q[b], Pcross);
      iq1pLHS = sLHS(Q[b1], Pcross);

      /* imply rules facing if 0 */

      if(iqpLHS == 0 && iq1pLHS != 0)
         iqpLHS=iq1pLHS*-1;
      else if(iqpLHS != 0 && iq1pLHS == 0)
         iq1pLHS=iqpLHS*-1;


      /* now calculate face rules */
      qpFace = sFace(ip1qLHS, ipqLHS, iqpLHS);
      pqFace = sFace(iq1pLHS, iqpLHS, ipqLHS);

      /* Xcross product near zero !! so make it zero*/
      if(nx3< 1.0e-10)
      {
         ip1qLHS=0;
         ipqLHS=0;
         iq1pLHS=0;
         iqpLHS=0;
         qpFace=0;
         pqFace=0;
      }



      if( isGeared == False)
      {
         if(  (ipqLHS == 1 && iqpLHS == 1) ||  ( qpFace && pqFace )     )
         {
            aa++;a++;
         }
         else
         {
            isGeared = True;
         }
      }

      if(isGeared) {
         code = sSegSegInt(P[a1], P[a], Q[b1], Q[b], p, q);


         if (code == '1' || code == 'e') {
            if(DEBUG_VRL)
               prnPoint("(): intersect", p,3,True  );

            sAddPoint(R, r, p);

            /*
            if(code=='e')
              sAddPoint(R, r, q);
            */

            if (numIntersect++ == 0) {
               /* reset counters */
               aa = 0;
               bb = 0;
            }



            inflag = ( ipqLHS ==1 ? Pin : iqpLHS ==1 ? Qin : inflag );


            if(DEBUG_VRL)
              printf("%%InOut sets inflag=%s\n", prnInFlag(inflag));

         }

         if(DEBUG_VRL)
            printf("numIntersect=%d code=%c (ipqLHS=%d, ip1qLHS=%d), (iqpLHS=%d, iq1pLHS=%d), (qpFace=%d pqFace=%d)\n",numIntersect, code, ipqLHS, ip1qLHS,  iqpLHS,iq1pLHS, qpFace,pqFace);



         if (qpFace && pqFace)  {

            /* Advance either P or Q which has previously arrived ? */
            if(inflag == Pin) sAddPoint(R,r, P[a]);

            aa++;a++;


         } else if (qpFace) {
            if(inflag == Qin) sAddPoint(R,r, Q[b]);

            bb++;b++;


            /* advance q */
         } else if (pqFace) {
            /* advance p */
            if(inflag == Pin) sAddPoint(R,r,P[a]);

            aa++;a++;

         } else if (iqpLHS == -1) {
            /* advance q */
            //if(inflag== Qin) sAddPoint(R,r,Q[b]);
            bb++;b++;

            /* cross product zero  */
         } else if( ipqLHS==0 && ip1qLHS==0 && iq1pLHS ==0 && iqpLHS ==0   ){
            if(inflag==Pin)
            {bb++;b++;}
            else
            {aa++;a++;}

         }



         else {
            /* catch all */
            if(inflag==Pin) sAddPoint(R,r,P[a]);
            aa++;a++;

         }

      }

      a%=n;
      b%=m;

      if(DEBUG_VRL)
         fprintf(stdout, "\ndebug isGeared=%d a=%d aa=%d b=%d bb=%d \n",isGeared, a, aa, b, bb);

      /* quick exit if current point is same a First point  - nb an exact match ?*/
      if( *r >3 &&  R[0][3]==R[*r-1][3] && R[0][4]==R[*r-1][4] )
      {
         --*r;
         break;
      }


   } while ( ((aa < n) || (bb < m)) && (aa < 2*n) && (bb < 2*m) );

   return EXIT_SUCCESS;

}

char  sSegSegInt( tPointds a, tPointds b, tPointds c, tPointds d, tPointds p, tPointds q )
{
   int flg_dbg=1;
   int flg_sx=0;

   double nx1;
   double nx2;
   double nx3;

   double darc;

   tPointds Pcross;
   tPointds Qcross;
   tPointds Icross;



   if(flg_sx) {
      nx1=sxCross(a, b, Pcross);
      nx2=sxCross(c, d, Qcross);

      sphAddcrt(Pcross);
      sphAddcrt(Qcross);

      nx3=sCross(Pcross, Qcross, Icross);
      sphAddcrt(Icross);
   }
   else
   {
      nx1=sCross(a, b, Pcross);
      nx2=sCross(c, d, Qcross);

      nx3=sCross(Pcross, Qcross, Icross);
      sphAddcrt(Icross);
   }

   darc=atan(nx3);

   if(DEBUG_VRL) {
      prnPoint("sSegSegInt(): intersection", Icross, 3, True);
      printf("sSegSegInt(): ||Pcross||=%.20g ||Qcross||=%.20g ||Icross||=%.20g arc=%.20g\n", nx1, nx2, nx3, darc);
   }

   /* Icross is zero, should really have a range rather than an explicit zero */
   if( nx3 < 1.0e-15)
      return sParallelDouble(a,b,c,d,p,q);


   if( sLatLonBetween(a,b, Icross ) && sLatLonBetween(c,d, Icross) )
   {
      memcpy(p,Icross, sizeof(tPointds));
      return '1';
   }

   /* try antipodal point */
   Icross[0]*= -1.0;
   Icross[1]*= -1.0;
   Icross[2]*= -1.0;

   sphAddcrt(Icross);

   if( sLatLonBetween(a,b, Icross ) && sLatLonBetween(c,d, Icross) )
   {

      memcpy(p,Icross, sizeof(tPointds));
      return '1';
   }

   return '0';





}


/* takes a point and a cross product representing the normal to the arc plane */
/* returns 1 if point on LHS of arc plane */
/* returns -1 if point on RHS of arc plane */
/* return 0 if point on the arc - (given suitable tolerances ) */
int sLHS(tPointds Pi, tPointds Qcross )
{
   double ds;

   ds=sDot(Pi,Qcross);

   if(ds  > 0.0)
      return 1;
   else if(ds <0.0)
      return -1;
   else
      return 0;


   /*
   ds=acos( sDot(Pi,Qcross) );

   if( ds < M_PI_2 - ARC_MIN_LENGTH )
     return 1;
   else if ( ds > M_PI_2 + ARC_MIN_LENGTH )
     return -1;
   else
     return 0;

   */
}

/* implement face rules */
nco_bool sFace( int iLHS, int iRHS, int jRHS  )
{
   if( iLHS == 1 && iRHS == -1 && jRHS == -1 )
      return True;

   if( iLHS == -1 && iRHS == 1 && jRHS == 1  )
      return True;

   return False;


}



double  sDot( tPointds a, tPointds b )
{
   int idx;
   double sum=0.0;

   for(idx=0; idx<3; idx++)
      sum+=a[idx]*b[idx];

   return sum;


}

double  sCross(tPointds a, tPointds b, tPointds c)
{
   //
   int flg_dbg=0;
   double n1;

   c[0]=a[1]*b[2]-a[2]*b[1];
   c[1]=a[2]*b[0]-a[0]*b[2];
   c[2]=a[0]*b[1]-a[1]*b[0];

   // normalize vector
   n1=sqrt( c[0]*c[0]+c[1]*c[1] + c[2]*c[2] );

   if( n1 >  0.0 && n1 != 1.0  )
   {
      c[0] /= n1;
      c[1] /= n1;
      c[2] /= n1;
   }

   if(DEBUG_VRL)
      printf("sCross(): n1=%f (%f, %f %f)\n", n1, c[0],c[1], c[2]);

   return n1;

}

double sRadius(tPointds a){
  double n1;

  n1=sqrt( a[0]*a[0]+a[1]*a[1] + a[2]*a[2] );

  return n1;
}


/* new method for calculating cross product */
double sxCross( tPointds a, tPointds b, tPointds c  )
{
   int flg_dbg=0;

   double n1;
   double lon1;
   double lon2;

   double lat1;
   double lat2;

   lon1=a[3] * M_PI /180.0;
   lat1=a[4] * M_PI /180.0;

   lon2=b[3] * M_PI /180.0;
   lat2=b[4] * M_PI /180.0;



   c[0] =   sin(lat1+lat2) * cos( (lon1+lon2) / 2.0) * sin( (lon1-lon2)/2.0)
            - sin(lat1-lat2) * sin ((lon1+lon2) / 2.0) * cos( (lon1-lon2)/2.0);

   c[1] =   sin(lat1+lat2) * sin( (lon1+lon2) / 2.0) * sin( (lon1-lon2)/2.0)
            + sin(lat1-lat2) * cos ((lon1+lon2) / 2.0) * cos( (lon1-lon2)/2.0);



   c[2]=cos(lat1) * cos(lat2) * sin(lon2-lon1);


   // normalize vector
   n1=sqrt( c[0]*c[0]+c[1]*c[1] + c[2]*c[2] );

   if( n1 != 0.0 && n1 !=1.0  )
   {
      c[0] /= n1;
      c[1] /= n1;
      c[2] /= n1;
   }

   if(DEBUG_VRL)
      printf("sxCross(): n1=%f (%f, %f %f)\n", n1, c[0],c[1], c[2]);

   return n1;

}


void  sAdi(tPointds a, tPointds b )
{
   (void)memcpy(a,b, sizeof(tPointds));
}


void  sph2crt(tPointds a,  double *lon, double *lat, nco_bool bDeg)
{

   /* nb this returns range (-180, 180) */
   *lon = atan2(a[1],a[0]) ;
   if( *lon < 0.0 && IS_LON_360)
      *lon+= (M_PI*2);

   // b[1]= asin(a[2]) * 180.0 /M_PI;
   *lat=atan2( a[2], sqrt( a[0]*a[0]+a[1]*a[1] ) ) ;

   /* convert to degrees if required */
   if(bDeg)
   {
      *lon*=(180.0 / M_PI );
      *lat*=(180.0 / M_PI );

   }

   return;
}

void crt2sph(tPointd a, tPointds b)
{
   double lon;
   double lat;

   lon=a[0] * M_PI / 180.0;
   lat=a[1] * M_PI / 180.0;


   b[0] = cos(lat) * cos(lon);
   b[1] = cos(lat) * sin(lon);
   b[2] = sin(lat);

   /* lat lon - we need this for bounding box */
   b[3] = lon;
   b[4] = lat;

}


void sphAddcrt(tPointds ds)
{
   sph2crt(ds, &ds[3], &ds[4],0 );
}

void sAddPoint( tPolygonds R, int *r, tPointds P)
{
   int flg_dbg=1;

   double delta;


   delta = ( *r==0 ? 0.0 :   2.0 *asin(    sqrt( pow( R[*r-1][0] - P[0],2 ) + pow( R[*r-1][1] - P[1],2 ) + pow( R[*r-1][2] - P[2],2 )  ) /2.0) );

   if(DEBUG_VRL)
      prnPoint("aAddPoint():", P,3,True );



   /* only add  point if its distinct from previous point */
   if ( *r==0 ||  delta > ARC_MIN_LENGTH_RAD )
   {

      memcpy(R[*r], P, sizeof(tPointds));
      (*r)++;
   }


}

nco_bool iBetween(double a, double b, double x  )
{

   nco_bool sdiff=False;
   int flg_dbg=0;

   if(DEBUG_VRL)
      printf("iBetween(): a=%.20f, b=%.20f, x=%.20f\n", a, b, x);

   if(fabs(b-a) < DSIGMA  )
      sdiff=True;
   else
      sdiff=False;

   if(sdiff) {
      if (fabs(x - a) < DSIGMA || fabs(x - b) < DSIGMA)
         return True;
      else
         return False;
   }

   if(  b >a &&  x>= a && x<=b  || b<a && x>=b && x<=a    )
      return True;
   else
      return False;

}

/* assume latitude -90,90 */
double latCorrect( double lat1, double lon1, double lon2  )
{

   double dp;

   if( lon1 == lon2  || lat1==0.0 || lat1 == M_PI /2.0   || lat1 == -M_PI/2.0  )
      return lat1;

   //lat1=lat1*M_PI / 180.0;

   dp= tan(lat1) / cos ( lon2-lon1 ) ;

   dp=atan(dp);


   return dp;


}

void getLatCorrect(double lon1, double lat1, double lon2, double lat2, double *dp_min, double *dp_max, nco_bool bDeg)
{


  if( lat2 >lat1 )
  {
    double dswp;

    dswp=lat1;
    lat1=lat2;
    lat2=dswp;

  }

  if(bDeg)
  {
    lat1 *= M_PI / 180.0;
    lat2 *= M_PI / 180.0;
    lon1 *= M_PI / 180.0;
    lon2 *= M_PI / 180.0;
  }




  /* lat1 & lat2 >0.0 */
  if( lat1>0.0 && lat2 >=0.0)
  {
    *dp_max = latCorrect(lat1, lon1, lon2);
    *dp_min = lat2;
  }
  else if( lat1 <= 0.0 && lat2<0.0 )
  {
    *dp_max = lat1;
    *dp_min = latCorrect(lat2, lon1, lon2);
  }

  else if( lat1 >0.0 && lat2  < 0.0)
  {
    *dp_max=latCorrect(lat1, lon1, lon2);
    *dp_min=latCorrect(lat2, lon1, lon2);

  }
  else
  {
    *dp_max=0.0;
    *dp_min=0.0;

  }

  /* convert back to degrees */
  if(bDeg)
  {
    *dp_max *= 180.0 / M_PI;
    *dp_min *= 180.0 / M_PI;
  }

  return;



}


void getLatCorrect_old(tPointds a, tPointds b, double *dp_min, double *dp_max )
{

   double lat1;
   double lat2;

   if( a[4] >= b[4] )
   {
      lat1 = a[4];
      lat2 = b[4];
   }
   else
   {
      lat1 = b[4];
      lat2 = a[4];
   }


   /* lat1 & lat2 >0.0 */
   if( lat1>0.0 && lat2 >=0.0)
   {
      *dp_max = latCorrect(lat1, a[3], b[3]);
      *dp_min = lat2;
   }
   else if( lat1 <= 0.0 && lat2<0.0 )
   {
      *dp_max = lat1;
      *dp_min = latCorrect(lat2, a[3], b[3]);
   }

   else if( lat1 >0.0 && lat2  < 0.0)
   {
      *dp_max=latCorrect(lat1, a[3], b[3]);
      *dp_min=latCorrect(lat2, a[3], b[3]);

   }
   else
   {
      *dp_max=0.0;
      *dp_min=0.0;

   }


   return;

}



/* use crt coords to check bounds */
nco_bool sLatLonBetween(tPointds a, tPointds b, tPointds x)
{

   /* working in radians here */
   nco_bool bDeg=False;
   int flg_dbg=1;

   double lat_min;
   double lat_max;

   if(  iBetween(a[3], b[3], x[3]) == False )
      return False;

   /* special lat check */
   //getLatCorrect(a,b, &lat_min,&lat_max);
   getLatCorrect(a[3],a[4],b[3],b[4], &lat_min,&lat_max, bDeg);


   if(DEBUG_VRL)
      printf("sBetween(): lat_min=%f lat_max=%f lat=%f\n", lat_min, lat_max, x[4]);

   if( x[4]>=lat_min && x[4]<=lat_max )
      return True;
   else
      return False;

   return False;


}

nco_bool sxBetween(tPointds a, tPointds b, tPointds c)
{

   if ( a[3] != b[3] )
      return (  c[3] >= a[3] && c[3] <=b[3] ||  c[3] <= a[3] && c[3] >= b[3] ) ;
   else
      return (  c[4] >= a[4] && c[4] <=b[4] ||  c[4] <= a[4] && c[4] >= b[4] ) ;


   /*
   if ( a[3] != b[3] )
     return (  a[3] <= c[3] && b[3] >= c[3] ||  a[3] >= c[3] && b[3] <= c[3] ) ;
   else
     return (  a[4] <= c[4] && b[4] >=c[4]   ||  a[4] >= c[4] && b[4] <= c[4] ) ;
   */

}


char sParallelDouble( tPointds a, tPointds b, tPointds c, tPointds d, tPointds p, tPointds q )
{

   int flg_dbg=1;

   char code='0';
   char *ptype="none";

   if( sxBetween( a, b, c ) && sxBetween( a, b, d ) ) {
      sAdi( p, c );
      sAdi( q, d );
      ptype="abc-abd";
      code= 'e';
   }
   else if( sxBetween( c, d, a ) && sxBetween( c, d, b ) ) {
      sAdi( p, a );
      sAdi( q, b );
      ptype="cda-cdb";
      code= 'e';
   }
   else if( sxBetween( a, b, c ) && sxBetween( c, d, b ) ) {
      sAdi( p, c );
      sAdi( q, b );
      ptype="abc-cdb";
      code= 'e';
   }
   else if( sxBetween( a, b, c ) && sxBetween( c, d, a ) ) {
      sAdi( p, c );
      sAdi( q, a );
      ptype="abc-cda";
      code= 'e';
   }
   else if( sxBetween( a, b, d ) && sxBetween( c, d, b ) ) {
      sAdi( p, d );
      sAdi( q, b );
      ptype="abd-cdb";
      code= 'e';
   }
   else if( sxBetween( a, b, d ) && sxBetween( c, d, a ) ) {
      sAdi( p, d );
      sAdi( q, a );
      ptype="abd-cda";
      code= 'e';
   }

   if(DEBUG_VRL)
      printf("sParallelDouble(): code=%c type=%s\n", code, ptype);

   return code;
}



void prnPoint(const char *sMsg, tPointds p, int style, nco_bool bRet )
{

   printf("%s ", sMsg);

   switch(style)
   {
      case 0:
      default:
         printf( "(dx=%.20f, dy=%.20f, dz=%.20f), (lon=%.20f,lat=%.20f)",p[0], p[1], p[2], p[3], p[4] );
       break;

      case 1:
         printf( "(dx=%.20f, dy=%.20f, dz=%.20f)",p[0], p[1], p[2] );
       break;

      case 2:
         printf( "(lon=%.20f,lat=%.20f)",p[3], p[4] );
       break;

      case 3:
         printf( "(lon=%.20f,lat=%.20f)",p[3] *180.0/M_PI,  p[4]*180/M_PI );
       break;

      case 4:
         printf( "(dx=%.20f, dy=%.20f, dz=%.20f), (lon=%.20f,lat=%.20f)",p[0], p[1], p[2], p[3] *180.0/M_PI,  p[4]*180/M_PI);
       break;

      case 5:
         printf( "(dx=%f, dy=%f, dz=%f), (lon=%f,lat=%f)",p[0], p[1], p[2], p[3] *180.0/M_PI,  p[4]*180/M_PI);
       break;



   }

   if(bRet)
      printf("\n");
   else
      printf(" * ");

}

nco_bool sConvex(tPolygonds sP, int nbr_vrt)
{

int flg_dbg=1;
int idx;
int idx_pre;
int idx_nex;


double n1;
double n2;

double dp;
double theta;
double rad1=1.0;
double rad=1.0;

tPointds aCross;
tPointds bCross;

for(idx=0; idx<nbr_vrt;idx++)
{
  idx_pre=(idx + nbr_vrt -1)% nbr_vrt;
  idx_nex=(idx + nbr_vrt +1)% nbr_vrt;

  n1=sxCross(sP[idx], sP[idx_pre],aCross);
  n2=sxCross(sP[idx], sP[idx_nex], bCross);

  //rad1 = sRadius(aCross);
  //rad  = sRadius(bCross);
  dp=sDot( aCross, bCross);


  // dp=sDot(sP[idx1], sP[idx]) / rad1 /rad;
  theta=acos(dp);

  if(DEBUG_VRL)
    printf("sConvex():, %d angle=%f n1=%.15g n2=%.15g\n", idx, theta*180.0/M_PI, n1, n2);

  if(theta > 2.0*M_PI   )
     return False;

}

return True;


}


/* works by counting the number of intersections of the
   line (pControl, pVertex) and each edge in sP
   pControl is chosen so that it is OUTSIDE sP
 */
nco_bool sPointInPolygon( tPolygonds sP, int n, tPointds pControl, tPointds pVertex)
{

  char code;
  int idx;
  int idx1=0;
  int numIntersect=0;

  tPointds p;
  tPointds q;


  /* count number of intersections */
  for(idx=0; idx< n ; idx++)
  {
    idx1=(idx+n -1) % n ;

    code=sSegSegInt( sP[idx1], sP[idx], pControl, pVertex, p, q );

    if(code=='1' || code=='v' || code == 'e')
      numIntersect++;


  }

  /* for any polygon (convex or concave)
    an odd  number of crossings means that the point is inside
    while an even number means that it is outside */

  return (numIntersect % 2  );





}

void sMakeControlPoint(tPolygonds P, int n, tPointds pControl)
{



}

/* set static globals */
void setStaticGlobals(double lon_min_rad, double lon_max_rad, double lat_min_rad, double lat_max_rad   )
{

  LON_MIN_RAD=lon_min_rad;
  LON_MAX_RAD=lon_max_rad;

  LAT_MIN_RAD=lat_min_rad;
  LAT_MAX_RAD=lat_max_rad;

  return;

}







