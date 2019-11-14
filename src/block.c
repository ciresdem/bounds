/*------------------------------------------------------------
 * block.c
 *
 * This file is part of BOUNDS
 *
 * Copyright (c) 2016, 2018, 2019 Matthew Love <matthew.love@colorado.edu>
 * BOUNDS is liscensed under the GPL v.2 or later and 
 * is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details. 
 * <http://www.gnu.org/licenses/> 
 *--------------------------------------------------------------*/

#include "bounds.h"

/* Get the minimum and maximum values from a set of points
 */
void
minmax (point_t* points, int npoints, region_t *xyzi) 
{
  int p;

  xyzi->xmin = points[0].x, xyzi->ymin = points[0].y;
  xyzi->xmax = points[0].x, xyzi->ymax = points[0].y;

  for (p = 1; p < npoints; p++) 
    {
      if (points[p].x <= xyzi->xmin) 
	xyzi->xmin = points[p].x;
      if (points[p].x >= xyzi->xmax) 
	xyzi->xmax = points[p].x;
      if (points[p].y <= xyzi->ymin) 
	xyzi->ymin = points[p].y;
      if (points[p].y >= xyzi->ymax) 
	xyzi->ymax = points[p].y;
    }
}

/* Return 1 if p1 and p2 are equal
 */
int
pnts_equal_p (point_t p1, point_t p2)
{
  if (fabs (p1.x - p2.x) > FLT_EPSILON) 
    return 0;
  if (fabs (p1.y - p2.y) > FLT_EPSILON) 
    return 0;
  return 1;
}

int
int_in_range (int val, int min, int max) 
{
  if (val < min) 
    return 0;
  if (val > max) 
    return 0;
  return 1;
}

/* "Bounding Block"
 * Generates a grid at `inc` cell-size and polygonizes it.
 */
int
bbe_block (point_t* points, int npoints, double inc, region_t region, int vflag) 
{
  int i, j, xpos, ypos, edge, lxi, lyi;
  int bb = 0, done = 0, bcount = 0, fcount = 0, bp = 1;
  point_t bb1, bb2;
  region_t xyi;

  /* Gather region info 
   */
  if (region_valid_p (&region)) 
    {
      xyi = region;
      if (vflag > 0) 
	fprintf (stderr, "bounds: Using user supplied region: %f/%f/%f/%f\n", 
		 region.xmin, region.xmax, region.ymin, region.ymax);
    } else minmax (points, npoints, &xyi);

  /* Set the rows and columns of the internal grid 
   */
  int ysize = fabs ((xyi.ymax - xyi.ymin) / inc) + 1;
  int xsize = fabs ((xyi.xmax - xyi.xmin) / inc) + 1;

  /* `xys` is the size of the bbarray array, this has to be large enough
   * to hold all 4 sides of every other cell in the grid. 
   */
  ssize_t xys = (xsize * ysize) * 2;
  
  if (vflag > 0) 
    fprintf (stderr, "bounds: Size of internal grid: %d/%d\n", 
	     ysize, xsize);

  /* Allocate memory for arrays 
   */  
  int** blockarray;
  blockarray = (int**) malloc (ysize * sizeof (int*));  
  for (i = 0; i < ysize; i++) 
    blockarray[i] = (int*) malloc (xsize * sizeof (int));  
  
  if (!blockarray)
    {
      if (vflag > 0) 
	fprintf (stderr, "bounds: Failed to allocate needed memory, try increasing the distance value (%f)\n", inc);
      exit (EXIT_FAILURE);
    }

  edge_t *bbarray;
  bbarray = (edge_t *) malloc (sizeof (edge_t) * xys);

  if (!bbarray) 
    {
      if (vflag > 0) 
	fprintf (stderr, "bounds: Failed to allocate needed memory, try increasing the distance value (%f)\n", inc);
      exit (EXIT_FAILURE);
    }

  /* Loop through the point records and grid them 
   */
  for (i = 0; i < npoints; i++) 
    {
      xpos = (points[i].x - xyi.xmin) / inc;
      ypos = (points[i].y - xyi.ymin) / inc;
      if (xpos >= 0 && xpos < xsize)
	if (ypos >= 0 && ypos < ysize)
	  blockarray[ypos][xpos] = 1;
    }

  /* Loop through the grid and find and record the edge lines 
   */
  for (i = 0; i < ysize; i++)
    for (j = 0; j < xsize; j++)
      if (blockarray[i][j] == 1) 
	{
	  // bottom edge (records left to right)
	  if (i == 0 || blockarray[i-1][j] != 1) 
	    {
	      bb1.x = ((j * inc) + xyi.xmin), bb1.y = ((i * inc) + xyi.ymin);
	      bb2.x = ((j * inc) + xyi.xmin) + inc, bb2.y = ((i * inc) + xyi.ymin);
	      bbarray[bb].p1 = bb1, bbarray[bb].p2 = bb2;
	      bbarray[bb].xi = j, bbarray[bb].yi = i;
	      bb++;
	    }
	  // left edge (records top to bottom)
	  if (j == 0 || blockarray[i][j-1] != 1) 
	    {
	      bb1.x = ((j * inc) + xyi.xmin), bb1.y = ((i * inc) + xyi.ymin) + inc;
	      bb2.x = ((j * inc) + xyi.xmin), bb2.y = ((i * inc) + xyi.ymin);
	      bbarray[bb].p1 = bb1, bbarray[bb].p2 = bb2;
	      bbarray[bb].xi = j, bbarray[bb].yi = i;
	      bb++;
	    }
	  // top edge (records right to left)
	  if (i == ysize-1 || blockarray[i+1][j] != 1) 
	    {
	      bb1.x = ((j * inc) + xyi.xmin) + inc, bb1.y = ((i * inc) + xyi.ymin) + inc;
	      bb2.x = ((j * inc) + xyi.xmin), bb2.y = ((i * inc) + xyi.ymin) + inc;
	      bbarray[bb].p1 = bb1, bbarray[bb].p2 = bb2;
	      bbarray[bb].xi = j, bbarray[bb].yi = i;
	      bb++;
	    }
	  // right edge (records bottom to top)
	  if (j == xsize-1 || blockarray[i][j+1] != 1) 
	    {
	      bb1.x = ((j * inc) + xyi.xmin) + inc, bb1.y = ((i * inc) + xyi.ymin) + inc;
	      bb2.x = ((j * inc) + xyi.xmin) + inc, bb2.y = ((i * inc) + xyi.ymin);
	      bbarray[bb].p1 = bb1, bbarray[bb].p2 = bb2;
	      bbarray[bb].xi = j, bbarray[bb].yi = i;
	      bb++;
	    }
	}

  for (i = 0; i < ysize; i++) 
    {
      free (blockarray[i]);
      blockarray[i] = NULL;
    }
  free (blockarray);
  blockarray = NULL;

  /* Allocate memory for the edge points. 
   * There may be more boundary points than input points, 
   * so a new array is used. 
   * Add one to the point count to repeat the first point.
  */
  point_t* bnds;  
  bnds = (point_t*) malloc (sizeof (point_t) * (bb + 1));

  if (vflag > 0) 
    {
      fprintf (stderr, "bounds: Sorting %d edge lines into a boundary\n", bb);
      fprintf (stderr, "bounds: %-10s\t %-10s\t %-10s\n", "Boundaries", "Points Sorted", "Remaining Edges");
      fprintf (stderr, "\rbounds: %-10d\t %-10d\t %-10d", 0, fcount, bb);
    }

  /* Sort the edge lines into polygons 
   * Do for each polygon...finishes (done == 1) when finds a matching point.
   */
  while (bb >= 4) 
    {
      if (bcount == 0) 
	{
	  //bnds[0] = bbarray[bb-1].p1, bnds[1] = bbarray[bb-1].p2;
	  bnds[0] = bbarray[0].p1, bnds[1] = bbarray[0].p2;
	  lxi = bbarray[0].xi, lyi = bbarray[0].yi;
	  bbarray[0] = bbarray[bb-1];
	  bcount=2, bb--;
	}

      while (done == 0) 
	{
	  /* Loop through the edges in bbarray and record the next point. 
	   */
	  for (edge = 0; edge < bb; edge++) 
	    if (int_in_range (bbarray[edge].xi, lxi-1, lxi+1))
	      if (int_in_range (bbarray[edge].yi, lyi-1, lyi+1))
		{
		  if (pnts_equal_p (bnds[bcount-1],  bbarray[edge].p1)) 
		    {
		      if (pnts_equal_p (bnds[0], bbarray[edge].p2)) 
			done = 1;
		      bnds[bcount] = bbarray[edge].p2;
		      lxi = bbarray[edge].xi, lyi = bbarray[edge].yi;
		      bbarray[edge] = bbarray[bb-1];
		      bb--, bcount++;
		      break;
		    }
		  if (pnts_equal_p (bnds[bcount-1], bbarray[edge].p2)) 
		    {
		      if (pnts_equal_p (bnds[0], bbarray[edge].p1)) 
			done = 1;
		      bnds[bcount] = bbarray[edge].p1;
		      lxi = bbarray[edge].xi, lyi = bbarray[edge].yi;
		      bbarray[edge] = bbarray[bb-1];
		      bb--, bcount++;
		      break;
		    }
		}
	}
      
      /* Print out the polygon 
       */
      for (edge = 0; edge < bcount; edge++) 
	  printf ("%f %f\n", bnds[edge].x, bnds[edge].y);
      
      /* Reset some values 
       */
      fcount += bcount;
      bcount = 0, done = 0;
      if (bb >= 4) 
	{
	  bp++;
	  printf( ">\n" );
	}

      if (vflag > 0) 
	{
	  fprintf (stderr, "\rbounds: %-10d\t %-10d\t %-10d", bp, fcount, bb);
	  fflush (stderr);
	}
    }
  
  if (vflag > 0) 
    fprintf (stderr, "\nbounds: Found %d boundaries using %d points\n", bp, fcount);

  /* Free remaining allocations
   */
  free (bbarray);
  bbarray = NULL;
  free (bnds);
  bnds = NULL;
  
  return (1);  
}
