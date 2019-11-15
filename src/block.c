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

/* Return 1 if `val` is is between `min` and `max`
 */
int
int_in_range (int val, int min, int max) 
{
  if (val < min) 
    return 0;
  if (val > max) 
    return 0;
  return 1;
}

/* Return 0 if `val` is less than zero, else return `val`
 */
int
int_or_zero (int val) 
{
  if (val < 0) return 0;
  else return val;
}

/* Return `max` if `val` is greater than `max`, else return `val`
 */
int
int_or_max (int val, int vmax) 
{
  if (val > vmax) return vmax;
  else return val;
}

/* Return a point based on cell location, increment and region (of a grid array)
 */
point_t
pixel_to_point (int xi, int yi, double inc, region_t xyi)
{
  point_t p1;
  p1.x = ((xi * inc) + xyi.xmin);
  p1.y = ((yi * inc) + xyi.ymin);
  return p1;
}

/* Return 1 if `ge` contains an edge
 */
int
g_edges_p (g_edges_t ge)
{
  if (ge.t == 1)
    return 1;
  else if (ge.b == 1)
    return 1;
  else if (ge.l == 1)
    return 1;
  else if (ge.r == 1)
    return 1;
  else
    return 0;
}

/* Return 1 if point `p1` equals point `p2` and return 2 if `p1` equals point `p3`, else return 0 (no match).
 */
int
pl_match (point_t p1, point_t p2, point_t p3)
{
  if (pnts_equal_p (p1,  p2))
    return 1;
  else if (pnts_equal_p (p1, p3))
    return 2;
  else
    return 0;
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
  int ysize = fabs ((xyi.ymax - xyi.ymin) / inc);
  int xsize = fabs ((xyi.xmax - xyi.xmin) / inc);

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

/* "Bounding Block"
 * Generates a grid at `inc` cell-size and polygonizes it into a boundary.
 */
int
bbp_block(point_t* points, int npoints, double inc, region_t region, int vflag) {
  int i, j, xpos, ypos, edge, l, lxi, lyi;
  int bp = 0, done = 0, pdone = 0, bcount = 0, fcount = 0, fyi = 0;
  point_t bb1, bb2, bb3;
  region_t xyi;

  /* Gather region info 
   */
  if (region_valid_p(&region)) 
    {
      xyi = region;
      if (vflag > 0) fprintf (stderr, "bounds: Using user supplied region: %f/%f/%f/%f\n", 
			      xyi.xmin, xyi.xmax, xyi.ymin, xyi.ymax);
    } else minmax(points, npoints, &xyi);

  /* Set the rows and columns of the internal grid 
   */
  int ysize = fabs((xyi.ymax - xyi.ymin) / inc);
  int xsize = fabs((xyi.xmax - xyi.xmin) / inc);
  ssize_t xys = (xsize * ysize) * 2;

  if (vflag > 0) fprintf(stderr,"bounds: Size of internal grid: %d/%d\n", 
			 ysize, xsize);

  /* Allocate memory for arrays 
   * `blockarray` will hold the point data location information
   * `edgearray` wll hold the edge information for each edge cell
   */
  int** blockarray;
  blockarray = (int**) malloc (ysize * sizeof(int*));
  for (i = 0; i < ysize; i++) blockarray[i] = (int*) malloc (xsize * sizeof (int));
  
  if (!blockarray)
    {
      if (vflag > 0) 
	fprintf (stderr,"bounds: Failed to allocate needed memory, try increasing the distance value (%f)\n", inc);
      exit (EXIT_FAILURE);
    }

  g_edges_t** edgearray;
  edgearray = (g_edges_t**) malloc (ysize * sizeof (g_edges_t*));
  for (i = 0; i < ysize; i++) edgearray[i] = (g_edges_t*) malloc (xsize * sizeof (g_edges_t));

  if (!edgearray) 
    {
      if (vflag > 0) 
	fprintf (stderr,"bounds: Failed to allocate needed memory, try increasing the distance value (%f)\n", inc);
      exit (EXIT_FAILURE);
    }

  if (vflag > 0) fprintf(stderr,"bounds: Gridding points\n");

  /* Assign all blockarray values to zero
   */
  for (i = 0; i < ysize; i++) 
    for (j = 0; j < xsize; j++) 
      blockarray[i][j] = 0;

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

  if (vflag > 0) 
    fprintf (stderr,"bounds: points gridded\nbounds: recording edges from grid\n");

  /* Allocate memory for the edge points.
   * There may be more boundary points than input points,
   * so a new array is used
   */
  point_t* bnds;
  bnds = (point_t*) malloc (sizeof (point_t) * xys);

  /* Loop through the grid and record the cell edges into edgearray
   */
  for (i = 0; i < ysize; i++) 
    for (j = 0; j < xsize; j++) 
      if (blockarray[i][j] == 1) 
	{
	  if (i == 0 || blockarray[i-1][j] == 0) 
	    edgearray[i][j].b = 1;

	  if (j == 0 || blockarray[i][j-1] == 0) 
	    edgearray[i][j].l = 1;

	  if (i == ysize-1 || blockarray[i+1][j] == 0) 
	    edgearray[i][j].t = 1;

	  if (j == xsize-1 || blockarray[i][j+1] == 0) 
	    edgearray[i][j].r = 1;

	  /* Set the cell to 2 if there are no edges recorded.
	   */
	  if (g_edges_p (edgearray[i][j]) == 0)
	    blockarray[i][j] = 2;
	}

  if (vflag > 0) 
    {
      fprintf (stderr, "bounds: %-10s\t %-10s\n", "Boundaries", "Points Sorted");
      fprintf (stderr, "\rbounds: %-10d\t %-10d", 0, fcount);
    }

  /* Scan the edgearray and arrange and output the edges into polygon(s).
   * pdone is 1 when we can't find any more edge cells.
   */
  while (pdone == 0)
    {

      /* Find the first edge in the polygon; add it to bnds and record it's position.
       */
      for (i = int_or_zero (fyi - 2); i < ysize; i++) 
	for (j = 0; j < xsize; j++) 
	  if (blockarray[i][j] == 1)
	    {
	      if (edgearray[i][j].b == 1)
		{
		  bb3 = pixel_to_point (j, i, inc, xyi);
		  bnds[0].x = bb3.x, bnds[0].y = bb3.y;
		  bnds[1].x = bb3.x + inc, bnds[1].y = bb3.y;
		  edgearray[i][j].b = 0, bcount = 2;
		  fyi = i, lxi = j, lyi = i, i = ysize, j = xsize;
		}
	      else if (edgearray[i][j].l == 1)
		{
		  bb3 = pixel_to_point (j, i, inc, xyi);
		  bnds[0].x = bb3.x, bnds[0].y = bb3.y + inc;
		  bnds[1].x = bb3.x, bnds[1].y = bb3.y;
		  edgearray[i][j].l = 0, bcount = 2;
		  fyi = i, lxi = j, lyi = i, i = ysize, j = xsize;
		}
	      else if (edgearray[i][j].t == 1)
		{
		  bb3 = pixel_to_point (j, i, inc, xyi);
		  bnds[0].x = bb3.x + inc, bnds[0].y = bb3.y + inc;
		  bnds[1].x = bb3.x, bnds[1].y = bb3.y + inc;
		  edgearray[i][j].t = 0, bcount = 2;
		  fyi = i, lxi = j, lyi = i, i = ysize, j = xsize;
		}
	      else if (edgearray[i][j].r == 1)
		{
		  bb3 = pixel_to_point (j, i, inc, xyi);
		  bnds[0].x = bb3.x + inc, bnds[0].y = bb3.y + inc;
		  bnds[1].x = bb3.x + inc, bnds[1].y = bb3.y;
		  edgearray[i][j].r = 0, bcount = 2;
		  fyi = i, lxi = j, lyi = i, i = ysize, j = xsize;
		}
	    }
    
      if (bcount != 2) 
	{
	  done = 1;
	  pdone = 1;
	}
      
      /* Scan the nearby cells in the edgearray and build polygons.
       * done is 1 when we match the first point found above.
       */
      while (done == 0)
	{
	  for (i = int_or_zero (lyi-1); i < int_or_max (lyi + 2, ysize); i++)
	    for (j = int_or_zero (lxi-1); j < int_or_max (lxi + 2, xsize); j++)
	      if (blockarray[i][j] == 1)
		{
		  if (edgearray[i][j].b == 1)
		    {
		      bb3 = pixel_to_point (j, i, inc, xyi);
		      bb1.x = bb3.x, bb1.y = bb3.y;
		      bb2.x = bb3.x + inc, bb2.y = bb3.y;
		      l = pl_match (bnds[bcount - 1], bb1, bb2);
		      if (l)
			{
			  if (pl_match(bnds[0], bb1, bb2)) done = 1;
			  if (l == 1) bnds[bcount].x = bb2.x, bnds[bcount].y = bb2.y;
			  else bnds[bcount].x = bb1.x, bnds[bcount].y = bb1.y;

			  edgearray[i][j].b = 0;
			  lxi = j, lyi = i;
			  bcount++;

			  if (g_edges_p (edgearray[i][j]) == 0) blockarray[i][j] = 2;
			}
		    }
		  if (edgearray[i][j].l == 1)
		    {
		      bb3 = pixel_to_point (j, i, inc, xyi);
		      bb1.x = bb3.x, bb1.y = bb3.y + inc;
		      bb2.x = bb3.x, bb2.y = bb3.y;
		      l = pl_match (bnds[bcount - 1], bb1, bb2);
		      if (l)
			{
			  if (pl_match(bnds[0], bb1, bb2)) done = 1;
			  if (l == 1) bnds[bcount].x = bb2.x, bnds[bcount].y = bb2.y;
			  else bnds[bcount].x = bb1.x, bnds[bcount].y = bb1.y;

			  edgearray[i][j].l = 0;
			  lxi = j, lyi = i;
			  bcount++;
			  
			  if (g_edges_p (edgearray[i][j]) == 0) blockarray[i][j] = 2;
			}
		    }
		  if (edgearray[i][j].t == 1)
		    {
		      bb3 = pixel_to_point (j, i, inc, xyi);
		      bb1.x = bb3.x + inc, bb1.y = bb3.y + inc;
		      bb2.x = bb3.x, bb2.y = bb3.y + inc;
		      l = pl_match (bnds[bcount - 1], bb1, bb2);
		      if (l)
			{
			  if (pl_match(bnds[0], bb1, bb2)) done = 1;
			  if (l == 1) bnds[bcount].x = bb2.x, bnds[bcount].y = bb2.y;
			  else bnds[bcount].x = bb1.x, bnds[bcount].y = bb1.y;

			  edgearray[i][j].t = 0;
			  lxi = j, lyi = i;
			  bcount++;

			  if (g_edges_p (edgearray[i][j]) == 0) blockarray[i][j] = 2;
			}
		    }
		  if (edgearray[i][j].r == 1)
		    {
		      bb3 = pixel_to_point (j, i, inc, xyi);
		      bb1.y = bb3.y + inc, bb1.x = bb3.x + inc;
		      bb2.x = bb3.x + inc, bb2.y = bb3.y;
		      l = pl_match (bnds[bcount - 1], bb1, bb2);
		      if (l)
			{
			  if (pl_match(bnds[0], bb1, bb2)) done = 1;
			  if (l == 1) bnds[bcount].x = bb2.x, bnds[bcount].y = bb2.y;
			  else bnds[bcount].x = bb1.x, bnds[bcount].y = bb1.y;

			  edgearray[i][j].r = 0;
			  lxi = j, lyi = i;
			  bcount++;

			  if (g_edges_p (edgearray[i][j]) == 0) blockarray[i][j] = 2;
			}
		    }
		}
	}

      for (edge = 0; edge < bcount; edge++) 
	printf("%f %f\n", bnds[edge].x, bnds[edge].y);
  
      /* Reset some values 
       */
      fcount = fcount + bcount;
      bcount = 0, done = 0;
      
      if (pdone == 0) 
	{
	  bp++;
	  printf( ">\n" );
	  if (vflag > 0) 
	    {
	      fprintf (stderr, "\rbounds: %-10d\t %-10d", bp, fcount);
	      fflush (stderr);
	    }
	}
    }

  /* Cleanup up and return.
   */
  for (i = 0; i < ysize; i++) 
    {
      free (blockarray[i]);
      free (edgearray[i]);
      blockarray[i] = NULL;
      edgearray[i] = NULL;
    }

  free (blockarray);
  blockarray = NULL;
  free (edgearray);
  edgearray = NULL;
  
  if (vflag > 0) 
    fprintf (stderr,"\nbounds: Found %d total boundary points\n", fcount);

  free (bnds);
  bnds = NULL;
  
  return (1);
}
