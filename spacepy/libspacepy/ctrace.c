/*************************************************************************/
/*Fast, accurate tracing routines for tracing a streamline through a 2D
vector field.  The subroutines included here are meant to be called from
python routines.

Test compilation/execution:
gcc -DCTRACE_TEST ctrace.c -lm
./a.out

Copyright 2010 - 2014  Los Alamos National Security, LLC. */

/*************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Bilinear interpolation for x1,y1=0 and x2,y2=1     */
/* Q's are surrounding points such that Q00 = F[0,0], */
/* Q10 = F[1,0], etc.                                 */
static double bilin_reg(double x, double y, double Q00, 
			double Q01, double Q10, double Q11)
{
  double fout;
  //printf("Input = %.3f, %.3f\n",x, y);
  fout = 
    Q00*(1.0-x)*(1.0-y) +
    Q10* x *    (1.0-y) +
    Q01* y *    (1.0-x) +
    Q11* x * y;

//  printf("points = %.5f, %.5f, %.5f, %.5f; ans = %.5f\n",
//	 Q00, Q10, Q01, Q11, fout);

  return fout;
}

/* Trilinear interpolation for x1,y1,z1=0 and x2,y2,z2=1 */
/* Q's are surrounding points such that Q000 = F[0,0,0], */
/* Q100 = F[1,0,0], etc.                                 */
static double trilin_reg(double x, double y, double z,
			 double Q000, double Q001, double Q010, double Q011,
			 double Q100, double Q101, double Q110, double Q111)
{
  double fout;

  fout = 
    Q000*(1.0-x)*(1.0-y)*(1.0-z) +
    Q001*(1.0-x)*(1.0-y)*     z  +
    Q010*(1.0-x)*     y *(1.0-z) +
    Q011*(1.0-x)*     y *     z  +
    Q100*     x *(1.0-y)*(1.0-z) +
    Q101*     x *(1.0-y)*     z  +
    Q110*     x *     y *(1.0-z) +
    Q111*     x *     y *     z  ;

  return fout;
}

/* Check to see if we should break out of an integration */
/* As set now, code will extrapolate up to 1 point outside of current block. */
int DoBreak(int iloc, int jloc, int iSize, int jSize)
{
  int ibreak = 0;
  //printf("At %i, %i with limits %i, %i\n", iloc, jloc, iSize, jSize);
  if (iloc >= iSize || jloc >= jSize)
    ibreak = 1;
  if (iloc < -1 || jloc < -1)
    ibreak = 1;

  return ibreak;
}

/* Check to see if we should break out of an integration */
/* Same as above but for 3D fields */
int DoBreak3d(int iloc, int jloc, int kloc, int iSize, int jSize, int kSize)
{
  int ibreak = 0;
  if (iloc >= iSize || jloc >= jSize || kloc >= kSize)
    ibreak = 1;
  if (iloc < -1 || jloc < -1 || kloc < -1)
    ibreak = 1;

  return ibreak;
}

/* Create unit vectors of 2D field */
static void make_unit(int iSize, int jSize, double *ux, double *uy)
{
  int i;
  double magnitude;

  for (i=0; i< iSize * jSize; i++) {
	magnitude = sqrt(pow(ux[i],2) + pow(uy[i],2));
	ux[i] /= magnitude;
	uy[i] /= magnitude;
  }

  return;
}

/* Create unit vectors of 3D field */
static void make_unit3d(int iSize, int jSize, int kSize,
			double *ux, double *uy, double *uz)
{
  int i;
  double magnitude;

  for (i=0; i< iSize * jSize * kSize; i++) {
    magnitude = sqrt( pow(ux[i],2) + pow(uy[i],2) + pow(uz[i],2) );
    ux[i] /= magnitude;
    uy[i] /= magnitude;
    uz[i] /= magnitude;
  }

  return;
}

/* Min and MAX: */
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

/*Interpolate a value from a 2D grid
 *double x, double y: X and Y of location to interp
 *double *field: Field to interpolate on
 *int xloc, int yloc: position of grid point below x,y
 *int xsize, int ysize: Size of field in X and y
 */
#define grid_interp(x, y, field, xloc, yloc, xsize, ysize) \
  bilin_reg(x-xloc, y-yloc, *(field+yloc*xsize+xloc), \
	    *(field+(yloc+1)*xsize+xloc), *(field+yloc*xsize+xloc+1), \
	    *(field+(yloc+1)*xsize+xloc+1))

/*Interpolate a value from a 3D grid
 *double x, double y, double z: X, Y, Z of location to interp
 *double *field: Field to interpolate on
 *int xloc, int yloc, int zloc: position of grid point below x,y,z
 *int xsize, int ysize, int zsize: Size of field in X, Y, and Z
 */
#define grid_interp3d(x, y, z, field, xloc, yloc, zloc, xsize, ysize, zsize) \
  trilin_reg(x-xloc, y-yloc, z-zloc,           			  \
	     *(field+ xloc   *ysize*zsize+ yloc   *zsize+zloc),	  \
	     *(field+ xloc   *ysize*zsize+ yloc   *zsize+zloc+1), \
	     *(field+ xloc   *ysize*zsize+(yloc+1)*zsize+zloc),	  \
	     *(field+ xloc   *ysize*zsize+(yloc+1)*zsize+zloc+1), \
	     *(field+(xloc+1)*ysize*zsize+ yloc   *zsize+zloc),   \
	     *(field+(xloc+1)*ysize*zsize+ yloc   *zsize+zloc+1), \
	     *(field+(xloc+1)*ysize*zsize+(yloc+1)*zsize+zloc),   \
	     *(field+(xloc+1)*ysize*zsize+(yloc+1)*zsize+zloc+1))

//	     *(field+ zloc   *ysize*xsize+ yloc   *xsize+xloc),	  \
//	     *(field+(zloc+1)*ysize*xsize+ yloc   *xsize+xloc),	  \
//	     *(field+ zloc   *ysize*xsize+(yloc+1)*xsize+xloc),	  \
//	     *(field+(zloc+1)*ysize*xsize+(yloc+1)*xsize+xloc),	  \
//	     *(field+ zloc   *ysize*xsize+ yloc   *xsize+xloc+1), \
//	     *(field+(zloc+1)*ysize*xsize+ yloc   *xsize+xloc+1), \
//	     *(field+ zloc   *ysize*xsize+(yloc+1)*xsize+xloc+1), \
//	     *(field+(zloc+1)*ysize*xsize+(yloc+1)*xsize+xloc+1))

/* Simple tracing using Euler's method. */
/* Super fast but not super accurate.   */
int cEuler(int iSize, int jSize,           /* Grid size and max steps */
	   int maxstep, double ds,         /* Maxsteps and step size  */
	   double xstart, double ystart,   /* Starting locations      */
	   double xGrid[], double yGrid[], /* Actual coord system     */
	   double *ux, double *uy,         /* Field to trace through  */
	   double x[], double y[])         /* x, y of result stream   */
{
  /* Variable declarations */
  int n=0, i, xloc, yloc;
  double dx, dy, fx, fy;

  /* Get starting points in normalized/array coordinates */
  dx = xGrid[1] - xGrid[0];
  dy = yGrid[1] - yGrid[0];
  x[0] = (xstart-xGrid[0]) / dx;
  y[0] = (ystart-yGrid[0]) / dy;
  
  /* Create unit vectors from full vector field */
  make_unit(iSize, jSize, ux, uy);

  /* Perform tracing using Euler's method */
  for(n=0; n<maxstep-1; n++)
    {
      /* Find surrounding points */
      xloc = floor(x[n]);
      yloc = floor(y[n]);

      /* Break if we leave the domain */
      if (DoBreak(xloc, yloc, iSize, jSize))
      break;
      if(xloc>iSize-2) xloc=iSize-2;
      if(xloc<0)       xloc=0;
      if(yloc>iSize-2) yloc=iSize-2;
      if(yloc<0)       yloc=0;
      /* Interpolate unit vectors to current location */
      fx = grid_interp(x[n], y[n], ux, xloc, yloc, iSize, jSize);
      fy = grid_interp(x[n], y[n], uy, xloc, yloc, iSize, jSize);
      
      /* Detect NaNs in function values */
      //if (isnan(fx) || isnan(fy) || isinf(fx) || isinf(fy))
      //break;

      /* Perform single step */
      x[n+1] = x[n] + ds * fx;
      y[n+1] = y[n] + ds * fy;
    }

  /* Return traced points to original coordinate system. */
  for (i=0; i<=n; i++)
    {
      x[i] = x[i]*dx + xGrid[0];
      y[i] = y[i]*dy + yGrid[0];
    }
  //printf("Used %i of %i points\n",n,maxstep);
  //printf("Breaking at %.3f, %.3f\n", x[n], y[n]);
  return n;
}

/* Fast and reasonably accurate tracing with */
/* Runge-Kutta 4 method (constant step size) */
int cRk4(int iSize, int jSize,             /* Grid size and max steps */
	 int maxstep, double ds,           /* Maxsteps and step size  */
	 double xstart, double ystart,     /* Starting locations      */
	 double xGrid[], double yGrid[],   /* Actual coord system     */
	 double *ux, double *uy,           /* Field to trace through  */
	 double x[], double y[])           /* x, y of result stream   */
{
  /* Variable declarations */
  int n=0, i, xloc, yloc;
  double dx, dy, xpos, ypos,
    f1x, f2x, f3x, f4x, f1y, f2y, f3y, f4y;

  /* Get starting points in normalized/array coordinates */
  dx = xGrid[1] - xGrid[0];
  dy = yGrid[1] - yGrid[0];
  x[0] = (xstart-xGrid[0]) / dx;
  y[0] = (ystart-yGrid[0]) / dy;
  
  /* Create unit vectors from full vector field */
  make_unit(iSize, jSize, ux, uy);

  /* Perform tracing using RK4 */
  for(n=0; n<maxstep-1; n++){
    /* See Euler's method for more descriptive comments. */
    /* SUBSTEP #1 */
    xloc = floor(x[n]);
    yloc = floor(y[n]);
    if (DoBreak(xloc, yloc, iSize, jSize))
      break;
    if(xloc>iSize-2) xloc=iSize-2;
    if(xloc<1)       xloc=1;
    if(yloc>iSize-2) yloc=iSize-2;
    if(yloc<1)       yloc=1;
    //xloc=MAX(0, MIN(iSize-2, xloc));
    //yloc=MAX(0, MIN(jSize-2, yloc));
    f1x = grid_interp(x[n], y[n], ux, xloc, yloc, iSize, jSize);
    f1y = grid_interp(x[n], y[n], uy, xloc, yloc, iSize, jSize);
    if (isnan(f1x) || isnan(f1y) || isinf(f1x) || isinf(f1y))
      break;

    /* SUBSTEP #2 */
    xpos = x[n]+f1x*ds/2.0;
    ypos = y[n]+f1y*ds/2.0;
    xloc = floor(xpos);
    yloc = floor(ypos);
    if (DoBreak(xloc, yloc, iSize, jSize))
      break;
    if(xloc>iSize-2) xloc=iSize-2;
    if(xloc<1)       xloc=1;
    if(yloc>iSize-2) yloc=iSize-2;
    if(yloc<1)       yloc=1;
    //xloc=MAX(0, MIN(iSize-2, xloc));
    //yloc=MAX(0, MIN(jSize-2, yloc));
    f2x = grid_interp(xpos, ypos, ux, xloc, yloc, iSize, jSize);
    f2y = grid_interp(xpos, ypos, uy, xloc, yloc, iSize, jSize);

    if (isnan(f2x) || isnan(f2y) || isinf(f2x) || isinf(f2y))
      break;

    /* SUBSTEP #3 */
    xpos = x[n]+f2x*ds/2.0;
    ypos = y[n]+f2y*ds/2.0;
    xloc = floor(xpos);
    yloc = floor(ypos);
    if (DoBreak(xloc, yloc, iSize, jSize))
      break;
    if(xloc>iSize-2) xloc=iSize-2;
    if(xloc<1)       xloc=1;
    if(yloc>iSize-2) yloc=iSize-2;
    if(yloc<1)       yloc=1;
    //xloc=MAX(0, MIN(iSize-2, xloc));
    //yloc=MAX(0, MIN(jSize-2, yloc));
    f3x = grid_interp(xpos, ypos, ux, xloc, yloc, iSize, jSize);
    f3y = grid_interp(xpos, ypos, uy, xloc, yloc, iSize, jSize);
    if (isnan(f3x) || isnan(f3y) || isinf(f3x) || isinf(f3y))
      break;

    /* SUBSTEP #4 */
    xpos = x[n]+f3x*ds;
    ypos = y[n]+f3y*ds;
    xloc = floor(xpos);
    yloc = floor(ypos);
    if (DoBreak(xloc, yloc, iSize, jSize))
      break;
    if(xloc>iSize-2) xloc=iSize-2;
    if(xloc<1)       xloc=1;
    if(yloc>iSize-2) yloc=iSize-2;
    if(yloc<1)       yloc=1;
    //xloc=MAX(0, MIN(iSize-2, xloc));
    //yloc=MAX(0, MIN(jSize-2, yloc));
    f4x = grid_interp(xpos, ypos, ux, xloc, yloc, iSize, jSize);
    f4y = grid_interp(xpos, ypos, uy, xloc, yloc, iSize, jSize);
    if (isnan(f4x) || isnan(f4y) || isinf(f4x) || isinf(f4y))
      break;

    /* Peform the full step using all substeps */
    x[n+1] = (x[n] + ds/6.0 * (f1x + f2x*2.0 + f3x*2.0 + f4x));
    y[n+1] = (y[n] + ds/6.0 * (f1y + f2y*2.0 + f3y*2.0 + f4y));

  }

  /* Return traced points to original coordinate system. */
  for (i=0; i<=n; i++)
    {
      x[i] = x[i]*dx + xGrid[0];
      y[i] = y[i]*dy + yGrid[0];
    }
  return n;

}

/* Fast and reasonably accurate tracing with */
/* Runge-Kutta 4 method (constant step size) */
int cRk4_3d(int iSize, int jSize, int kSize,    /* Grid size and max steps */
	    int maxstep, double ds,             /* Maxsteps and step size  */
	    double xstart, double ystart,       /* Starting locations      */
	    double zstart,                      /* Starting locations      */
	    double xGrid[], double yGrid[],     /* Actual coord system     */
	    double zGrid[],                     /* Actual coord system     */
	    double *ux, double *uy, double *uz, /* Field to trace through  */
	    double x[], double y[], double z[]) /* x, y of result stream   */
{
  /* Variable declarations */
  int n=0, i, xloc, yloc, zloc;
  double dx, dy, dz, xpos, ypos, zpos,
    f1x, f2x, f3x, f4x, f1y, f2y, f3y, f4y, f1z, f2z, f3z, f4z;

  printf("STARTING XYZ = %.4f, %.4f, %.4f\n", xstart, ystart, zstart);
  
  /* Get starting points in normalized/array coordinates */
  dx = xGrid[1] - xGrid[0];
  dy = yGrid[1] - yGrid[0];
  dz = zGrid[1] - zGrid[0];
  x[0] = (xstart-xGrid[0]) / dx;
  y[0] = (ystart-yGrid[0]) / dy;
  z[0] = (zstart-zGrid[0]) / dz;

  printf("dX, dY, dZ = %.4f, %.4f, %.4f\n", dx, dy, dz);
  printf("Grid starts (XYZ) = %.4f, %.4f, %.4f\n",xGrid[0], yGrid[0], zGrid[0]);
  printf("Normalized starting XYZ = %.4f, %.4f, %.4f\n",x[0],y[0],z[0]);
  
  /* Create unit vectors from full vector field */
  printf("Full starting U = %.4f, %.4f, %.4f\n",ux[0],uy[0],uz[0]);
  make_unit3d(iSize, jSize, kSize, ux, uy, uz);
  printf("Normalized starting U = %.4f, %.4f, %.4f\n",ux[0],uy[0],uz[0]);
  
  /* Perform tracing using RK4 */
  for(n=0; n<maxstep-1; n++){
    
    /* SUBSTEP #1 */
    /* Find surrounding points */
    xloc = floor(x[n]);
    yloc = floor(y[n]);
    zloc = floor(z[n]);

    printf("XYZ Locs: %d, %d, %d\n", xloc, yloc, zloc);
    
    /* Break if we leave the domain */
    if (DoBreak3d(xloc, yloc, zloc, iSize, jSize, kSize))
      break;
    if(xloc>iSize-2) xloc=iSize-2;
    if(xloc<1)       xloc=1;
    if(yloc>jSize-2) yloc=jSize-2;
    if(yloc<1)       yloc=1;
    if(zloc>kSize-2) zloc=kSize-2;
    if(zloc<1)       zloc=1;

    /* Interpolate unit vectors to current location */
    f1x = grid_interp3d(x[n],y[n],z[n], ux, xloc,yloc,zloc, iSize,jSize,kSize);
    f1y = grid_interp3d(x[n],y[n],z[n], uy, xloc,yloc,zloc, iSize,jSize,kSize);
    f1z = grid_interp3d(x[n],y[n],z[n], uz, xloc,yloc,zloc, iSize,jSize,kSize);
    if ( isnan(f1x) || isnan(f1y) || isnan(f1z) ||
	 isinf(f1x) || isinf(f1y) || isinf(f1z) )
      break;

    printf("xny[n] = %.4f, %.4f, %.4f\n", x[n],y[n],z[n]);
    printf("F1xyz = %.4f, %.4f, %.4f\n", f1x, f1y, f1z);
    
    /* SUBSTEP #2 */
    xpos = x[n]+f1x*ds/2.0;
    ypos = y[n]+f1y*ds/2.0;
    zpos = z[n]+f1z*ds/2.0;
    xloc = floor(xpos);
    yloc = floor(ypos);
    zloc = floor(zpos);
    
    if (DoBreak3d(xloc, yloc, zloc, iSize, jSize, kSize))
      break;
    if(xloc>iSize-2) xloc=iSize-2;
    if(xloc<1)       xloc=1;
    if(yloc>jSize-2) yloc=jSize-2;
    if(yloc<1)       yloc=1;
    if(zloc>kSize-2) zloc=kSize-2;
    if(zloc<1)       zloc=1;

    f2x = grid_interp3d(xpos,ypos,zpos, ux, xloc,yloc,zloc, iSize,jSize,kSize);
    f2y = grid_interp3d(xpos,ypos,zpos, uy, xloc,yloc,zloc, iSize,jSize,kSize);
    f2z = grid_interp3d(xpos,ypos,zpos, uz, xloc,yloc,zloc, iSize,jSize,kSize);
    if ( isnan(f2x) || isnan(f2y) || isnan(f2z) ||
	 isinf(f2x) || isinf(f2y) || isinf(f2z) )
      break;

    /* SUBSTEP #3 */
    xpos = x[n]+f2x*ds/2.0;
    ypos = y[n]+f2y*ds/2.0;
    zpos = z[n]+f2z*ds/2.0;
    xloc = floor(xpos);
    yloc = floor(ypos);
    zloc = floor(zpos);

    if (DoBreak3d(xloc, yloc, zloc, iSize, jSize, kSize))
      break;
    if(xloc>iSize-2) xloc=iSize-2;
    if(xloc<1)       xloc=1;
    if(yloc>jSize-2) yloc=jSize-2;
    if(yloc<1)       yloc=1;
    if(zloc>kSize-2) zloc=kSize-2;
    if(zloc<1)       zloc=1;
    
    f3x = grid_interp3d(xpos,ypos,zpos, ux, xloc,yloc,zloc, iSize,jSize,kSize);
    f3y = grid_interp3d(xpos,ypos,zpos, uy, xloc,yloc,zloc, iSize,jSize,kSize);
    f3z = grid_interp3d(xpos,ypos,zpos, uz, xloc,yloc,zloc, iSize,jSize,kSize);
    if ( isnan(f3x) || isnan(f3y) || isnan(f3z) ||
	 isinf(f3x) || isinf(f3y) || isinf(f3z) )
      break;
    
    /* SUBSTEP #4 */
    xpos = x[n]+f3x*ds;
    ypos = y[n]+f3y*ds;
    zpos = z[n]+f3z*ds;
    xloc = floor(xpos);
    yloc = floor(ypos);
    zloc = floor(zpos);

    if (DoBreak3d(xloc, yloc, zloc, iSize, jSize, kSize))
      break;
    if(xloc>iSize-2) xloc=iSize-2;
    if(xloc<1)       xloc=1;
    if(yloc>jSize-2) yloc=jSize-2;
    if(yloc<1)       yloc=1;
    if(zloc>kSize-2) zloc=kSize-2;
    if(zloc<1)       zloc=1;

    f4x = grid_interp3d(xpos,ypos,zpos, ux, xloc,yloc,zloc, iSize,jSize,kSize);
    f4y = grid_interp3d(xpos,ypos,zpos, uy, xloc,yloc,zloc, iSize,jSize,kSize);
    f4z = grid_interp3d(xpos,ypos,zpos, uz, xloc,yloc,zloc, iSize,jSize,kSize);
    if ( isnan(f4x) || isnan(f4y) || isnan(f4z) ||
	 isinf(f4x) || isinf(f4y) || isinf(f4z) )
      break;

    /* Peform the full step using all substeps */
    x[n+1] = (x[n] + ds/6.0 * (f1x + f2x*2.0 + f3x*2.0 + f4x));
    y[n+1] = (y[n] + ds/6.0 * (f1y + f2y*2.0 + f3y*2.0 + f4y));
    z[n+1] = (z[n] + ds/6.0 * (f1z + f2z*2.0 + f3z*2.0 + f4z));

    printf("Advanced to:\n\tX=%.4f\n\tY=%.4f\n\tZ=%.4f\n",
	   x[n+1]*dx + xGrid[0], y[n+1]*dy + yGrid[0], z[n+1]*dz + zGrid[0]);
    //if(n>3) break;
  }

  /* Return traced points to original coordinate system. */
  for (i=0; i<=n; i++)
    {
      x[i] = x[i]*dx + xGrid[0];
      y[i] = y[i]*dy + yGrid[0];
      z[i] = z[i]*dz + zGrid[0];
    }
  return n;

}

#ifdef CTRACE_TEST
static int test_arrays2d(int iSize, int jSize, double xGrid[], double yGrid[],
			 double *ux, double *uy){
  printf("Seems to work.\n");
  return 1;
}

static int test_arrays3d(int iSize, int jSize, int kSize,
			 double xGrid[], double yGrid[], double zGrid[],
			 double *ux, double *uy){
  printf("Seems to work.\n");
  return 1;
}


/* A utility for printing test results. */
void test_check(double result, double answer)
{
  double diff, thresh;

  thresh = 0.00001;
  diff = 100.0*abs(result - answer)/answer;

  if (diff < thresh)
    printf("TEST PASSED! (Result=%.4f, Answer=%.4f)\n", result, answer);
  else{
    printf("TEST FAILED!\n");
    printf("Result %.20f differs from answer %.20f\n", result, answer);
    printf("Difference of %.3f%% is over threshold of %.5f%%\n", diff, thresh);
  }
}

/* Main simply tests the functionality of all funcs in this file.*/
int main()
{
  /* Declarations */
  double out=0, x=0.1, y=0.2, z=0.3, Q00=3.0, Q01=5.0, Q10=40.0, Q11=60.0,
    Q000=100., Q001=10000., Q010=300., Q011=30000.0,
    Q100=2000., Q101=200000.0, Q110=4000., Q111=400000.0;
  double xall[27] = {1, 1, 1, 1, 1, 1, 1, 1, 1,
		     2, 2, 2, 2, 2, 2, 2, 2, 2,
		     3, 3, 3, 3, 3, 3, 3, 3, 3};
  double yall[27] = {4, 4, 4, 5, 5, 5, 6, 6, 6,
		     4, 4, 4, 5, 5, 5, 6, 6, 6,
		     4, 4, 4, 5, 5, 5, 6, 6, 6};
  double zall[27] = {7, 8, 9, 7, 8, 9, 7, 8, 9,
		     7, 8, 9, 7, 8, 9, 7, 8, 9,
		     7, 8, 9, 7, 8, 9, 7, 8, 9};
  double fall[27] = {30, 33, 36, 32, 35, 38, 34, 37, 40,
		     31, 34, 37, 33, 36, 39, 35, 38, 41,
		     32, 35, 38, 34, 37, 40, 36, 39, 42};
  double sol1 = 7.460;
  double sol2 = 11236.2;
  double sol3 = 39.0;
  int l[1];

  /* Test bilin_reg */
  printf("Testing bilin_reg\n");
  out = bilin_reg(x, y, Q00, Q01, Q10, Q11);
  printf("TEST 1: ");
  test_check(out, sol1);

  /* Test Trilin_reg */
  printf("Testing Triilin_reg\n");
  out = trilin_reg(x,y,z, Q000,Q001,Q010,Q011,Q100,Q101,Q110,Q111);
  printf("TEST 2: ");
  test_check(out, sol2);

  /* Test grid_interp3d */
  printf("Testing Grid_Interp3d\n");
  out = grid_interp3d(2.5-xall[0], 5.5-yall[0], 8.5-zall[0],
		      fall, 1,1,1, 3,3,3);
  printf("TEST 3: ");
  test_check(out, sol3);
  
  /* Test cEuler 1 */
  int i, j, nx=841, ny=121, maxstep=10000, npoints;
  double xgrid[nx], ygrid[nx], *ux, *uy, xt[maxstep], yt[maxstep], ds=1.0;
  ux = malloc(nx * ny * sizeof(double));
  uy = malloc(nx * ny * sizeof(double));
  
  for (i=0; i<nx; i++)
    {
      xgrid[i] = -10.0+0.25*i;
      ygrid[i] = xgrid[i];
    }

  for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
      {
	*(ux+i*ny+j) = xgrid[i];
	*(uy+i*ny+j) = -1.0 * ygrid[j];
      }

  npoints = cEuler(nx, ny, maxstep, ds, 1.0, 10.0, xgrid, ygrid,
			 ux, uy, xt, yt);
  printf("Npoints = %i\n", npoints);
  printf("Grid goes from %.2f to %.2f\n", xgrid[0], xgrid[nx-1]);
  printf("Our trace starts at %.2f, %.2f\n", xt[0], yt[0]);
  printf("...and ends at      %.2f, %.2f\n", xt[npoints], yt[npoints]);

  npoints = cRk4(nx, ny, maxstep, ds, 1.0, 10.0, xgrid, ygrid,
		 ux, uy, xt, yt);
  printf("Npoints = %i\n", npoints);
  printf("Grid goes from %.2f to %.2f\n", xgrid[0], xgrid[nx-1]);
  printf("Our trace starts at %.2f, %.2f\n", xt[0], yt[0]);
  printf("...and ends at      %.2f, %.2f\n", xt[npoints], yt[npoints]);  

  /*for (i=0; i<npoints; i++)
    printf("%.3f, %.3f\n", xt[i], yt[i]);*/
  return 0;
}
#endif /*CTRACE_TEST*/
