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

// #include "nco_vrl.h"


/*
int     	n, m;
tPolygoni	P, Q;
*/

/*
int main(int argc, char **argv)
{


  int     	n, m;
  tPolygond	P, Q;

  int r=0;
  tPolygond     R;
  
   n = ReadPoly( P );
   m = ReadPoly( Q );
   OutputPolygons(P, Q, n, m);
   ConvexIntersect( P, Q, R, n, m, &r);

   if(r >1) 
     PrintPoly(r,R); 
   
   ClosePostscript();
}

*/

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

      if(lcl_dbg)
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

         if(lcl_dbg)
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

          if(lcl_dbg)
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

     if(lcl_dbg)
       (void)fprintf(stdout, "%s: Before Advances:a=%d, b=%d; aa=%d, ba=%d; inflag=%d\n", nco_prg_nm_get(),   a, b, aa, ba, inflag);


   /* Quit when both adv. indices have cycled, or one has cycled twice. */
   } while ( ((aa < n) || (ba < m)) && (aa < 2*n) && (ba < 2*m) );

   if ( !FirstPoint ) 
   {
      if(lcl_dbg)
         (void)fprintf(stdout, "%s: no points output\n", nco_prg_nm_get());
      
      return EXIT_FAILURE;

   }
   
   /* Deal with special cases: not implemented. */
   if ( inflag == Unknown)
   {

      if(lcl_dbg)
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

  
  /* only add  point if its distinct from previous point */ 
  if ( *r == 0  ||    (fabs(R[*r-1][0] - P[0]) >DSIGMA || fabs(R[*r-1][1] - P[1])>DSIGMA) )
  {  

    R[*r][0] = P[0];
    R[*r][1] = P[1];
    (*r)++;     
    
  }
    
}



/*---------------------------------------------------------------------
Polygon I/O functions  
---------------------------------------------------------------------*/


/*
   Reads in the coordinates of the vertices of a polygon from stdin,
   puts them into P, and returns n, the number of vertices.
   Formatting conventions: etc.
*/
int   ReadPoly( tPolygond P )
{
   int   n = 0;
   int   nin;

   scanf("%d", &nin);
   /*printf("%%Polygon:\n");
   printf("%%  i   x   y\n");*/
   while ( (n < nin) && (scanf("%lf %lf",&P[n][0],&P[n][1]) != EOF) ) {
      /*printf("%%%3d%4d%4d\n", n, P[n][0], P[n][1]);*/
      ++n;
   }
/*
   if (n < PMAX)
      printf("%%n = %3d vertices read\n",n);
   else   printf("Error in read_poly:  too many points; max is %d\n", PMAX);
   putchar('\n');
*/

   return   n;
}

void PrintPoly(tPolygond R, int r)
{
  int idx;
  
   printf("Polygon R:\n");
   
   for( idx = 0; idx < r; idx++ )
      printf("%20.14f %20.14f\n", R[idx][0], R[idx][1]);
   
   printf("End Polygon\n");


}  


void   OutputPolygons(tPolygond P, tPolygond Q, int n, int m )
{
   int i;
   double xmin, ymin, xmax, ymax;

   /* Compute Bounding Box for Postscript header. */
   xmin = xmax = P[0][0];
   ymin = ymax = P[0][1];
   for (i = 1; i < n; i++) {
      if      ( P[i][0] > xmax ) xmax = P[i][0];
      else if ( P[i][0] < xmin ) xmin = P[i][0];
      if      ( P[i][1] > ymax ) ymax = P[i][1];
      else if ( P[i][1] < ymin ) ymin = P[i][1];
   }
   for (i = 0; i < m; i++) {
      if      ( Q[i][0] > xmax ) xmax = Q[i][0];
      else if ( Q[i][0] < xmin ) xmin = Q[i][0];
      if      ( Q[i][1] > ymax ) ymax = Q[i][1];
      else if ( Q[i][1] < ymin ) ymin = Q[i][1];
   }


   /* PostScript header */
   printf("%%!PS\n");
   printf("%%%%Creator: convconv.c (Joseph O'Rourke)\n");
   printf("%%%%BoundingBox: %f %f %f %f\n", xmin, ymin, xmax, ymax);
   printf("%%%%EndComments\n");
   printf(".00 .00 setlinewidth\n");
   printf("%d %d translate\n", -xmin+100, -ymin+100 );
   /* The +100 shifts the figure from the lower left corner. */

   printf("\n%%Polygon P:\n");
   printf("newpath\n");
   printf("%f\t%f\tmoveto\n", P[0][0], P[0][1]);
   for( i = 1; i <= n; i++ )
      printf("%f\t%f\tlineto\n", P[i][0], P[i][1]);
   
   printf("closepath stroke\n");

   printf("\n%%Polygon Q:\n");
   printf("newpath\n");
   printf("%f\t%f\tmoveto\n", Q[0][0], Q[0][1]);
   for( i = 1; i <= m; i++ )
      printf("%f\t%f\tlineto\n", Q[i][0], Q[i][1]);
   
   printf("closepath stroke\n");

   printf("2 2 setlinewidth\n");
   printf("newpath\n");
 }

void	PrintSharedSeg( tPointd p, tPointd q )
{
   printf("%%A int B:\n");
   printf("%8.2lf %8.2lf moveto\n", p[0], p[1] );
   printf("%8.2lf %8.2lf lineto\n", q[0], q[1] );
   ClosePostscript();
}

void   ClosePostscript( void )
{
   printf("closepath stroke\n");
   printf("showpage\n%%%%EOF\n");
}

  
