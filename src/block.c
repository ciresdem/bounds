/*------------------------------------------------------------
 * block.c
 *
 * This file is part of BOUNDS
 *
 * Copyright (c) 2016 - 2021 Matthew Love <matthew.love@colorado.edu>
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

int
int_or_min (int val, int min) 
{
  if (val < min) return min;
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
 * Generates a grid at `inc` cell-size and polygonizes it into a boundary.
 */
int
bbs_block(FILE *infile, double inc, region_t region, int vflag, int jflag) {
  int i, j, xpos, ypos, edge, l, lxi, lyi;
  int bp = 0, done = 0, pdone = 0, bcount = 0, fcount = 0, fyi = 0, dflag = 0;
  point_t rpnt, bb1, bb2, bb3;
  region_t xyi;
  ssize_t npr = 0;
  char* delim;
  char* ptrec = "xy";
  int forp_flag = 0;
  point_t* pnts;
  pnts = (point_t*) malloc (sizeof (point_t));
  
  /* Gather region info 
   */
  if (region_valid_p(&region)) 
    {
      xyi = region;
      if (vflag > 0) fprintf (stderr, "bounds: using user supplied region: %f/%f/%f/%f\n", 
			      xyi.xmin, xyi.xmax, xyi.ymin, xyi.ymax);
    }
  else
    {
      fprintf (stderr, "bounds: scanning xy data for region\n");
      load_pnts (infile, &pnts, &npr, ptrec, vflag);
      minmax(pnts, npr, &xyi);
      forp_flag = 1;
    }
  
  /* Set the rows and columns of the internal grid 
   */
  int ysize = fabs((xyi.ymax - xyi.ymin) / inc);// + 1;
  int xsize = fabs((xyi.xmax - xyi.xmin) / inc);// + 1;
  ssize_t xys = ((xsize * ysize) * 3) + 1;

  if (vflag > 0) {
    fprintf(stderr, "bounds: region is %f/%f/%f/%f\n", xyi.xmin, xyi.xmax, xyi.ymin, xyi.ymax);
    fprintf(stderr,"bounds: size of internal grid: %d/%d\n", 
	    ysize, xsize);
  }

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
	fprintf (stderr,"bounds: failed to allocate needed memory, try increasing the distance value (%f)\n", inc);
      exit (EXIT_FAILURE);
    }

  g_edges_t** edgearray;
  edgearray = (g_edges_t**) malloc (ysize * sizeof (g_edges_t*));
  for (i = 0; i < ysize; i++) edgearray[i] = (g_edges_t*) malloc (xsize * sizeof (g_edges_t));

  if (!edgearray) 
    {
      if (vflag > 0) 
	fprintf (stderr,"bounds: failed to allocate needed memory, try increasing the distance value (%f)\n", inc);
      exit (EXIT_FAILURE);
    }

  /* Assign all blockarray values to zero
   */
  for (i = 0; i < ysize; i++) 
    for (j = 0; j < xsize; j++) 
      blockarray[i][j] = 0;

  if (vflag > 0) fprintf(stderr,"bounds: gridding points\n");

  if (forp_flag == 0)
    {
      while (read_point(infile, &rpnt, &delim, ptrec, dflag, vflag) == 0) 
	{
	  xpos = (rpnt.x - xyi.xmin) / inc;
	  ypos = (rpnt.y - xyi.ymin) / inc;
	  if (xpos >= 0 && xpos < xsize)
	    if (ypos >= 0 && ypos < ysize)
	      blockarray[ypos][xpos] = 1;
	  npr++;
	  if (npr>0)
	    dflag++;
	}
    }
  else
    {
      /* Loop through the point records and grid them 
       */
      for (i = 0; i < npr; i++) 
	{
	  xpos = (pnts[i].x - xyi.xmin) / inc;
	  ypos = (pnts[i].y - xyi.ymin) / inc;
	  if (xpos >= 0 && xpos < xsize)
	    if (ypos >= 0 && ypos < ysize)
	      blockarray[ypos][xpos] = 1;
	}
    }
  
  if (vflag > 0) 
    fprintf (stderr,"bounds: %d points gridded\nbounds: recording edges from grid\n", npr);

  /* Allocate memory for the edge points.
   * There may be more boundary points than input points,
   * so a new array is used
   */
  point_t* bnds;
  bnds = (point_t*) malloc (sizeof (point_t) * xys);

  if (!bnds) 
    {
      if (vflag > 0) 
	fprintf (stderr,"bounds: failed to allocate needed memory, try increasing the distance value (%f) or shrinking the region\n", inc);
      exit (EXIT_FAILURE);
    }
  
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
	  done = 1, pdone = 1;
      else
	if (bp > 0)
	  if (jflag > 0)
	    printf ( "]],[[" );
	  else
	    printf ( ">\n" );
      
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
			  lxi = j, lyi = i, bcount++;

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
			  lxi = j, lyi = i, bcount++;
			  
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
			  lxi = j, lyi = i, bcount++;

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
			  lxi = j, lyi = i, bcount++;

			  if (g_edges_p (edgearray[i][j]) == 0) blockarray[i][j] = 2;
			}
		    }
		}
	}

      if (jflag > 0)
	{
	  for (edge = 0; edge < bcount-1; edge++) 
	    printf("[%.10f, %.10f],", bnds[edge].x, bnds[edge].y);
	  //fprintf (stderr, "%d", bcount-1);
	  if ( bcount > 1 )
	    {
	      printf("[%.10f, %.10f]", bnds[bcount-1].x, bnds[bcount-1].y);
	    }
	}
      else
	{
	  for (edge = 0; edge < bcount; edge++) 
	    printf("%.10f %.10f\n", bnds[edge].x, bnds[edge].y);
	}

      //if (pdone == 0)
      //printf( ">\n" );
      
      /* Reset some values 
       */
      fcount = fcount + bcount;
      bcount = 0, done = 0;
      
      if (pdone == 0) 
	bp++;
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
    fprintf (stderr,"bounds: found %d total boundary points\n", fcount);

  free (bnds);
  bnds = NULL;
  
  return (0);
}
