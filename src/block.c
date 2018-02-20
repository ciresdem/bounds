/*------------------------------------------------------------
 * block.c
 *
 * This file is part of BOUNDS
 *
 * Copyright (c) 2016 Matthew Love <matthew.love@colorado.edu>
 * BOUNDS is liscensed under the GPL v.2 or later and 
 * is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details. 
 * <http://www.gnu.org/licenses/> 
 *--------------------------------------------------------------*/

#include "bounds.h"

/* 
 * "Bounding Block"
 * Generates a grid at `inc` cell-size and polygonizes it.
 */
int
block_pts(point_t* points, int npoints, double inc) {

  int i, j, xpos, ypos, edge;
  int bb = 0, done = 0, bcount = 0, l = 0, fcount = 0;
  double x, y;
  point_t bb1, bb2;

  /* Gather xyzinfo */
  xyz_info xyzi;
  minmax(points, npoints, &xyzi);

  /* Set the rows and columns of the internal grid */
  int ysize = fabs((xyzi.ymax - xyzi.ymin) / inc) + 1;
  int xsize = fabs((xyzi.xmax - xyzi.xmin) / inc) + 1;
  /* `xys` is the size of the bbarray array, this has to be large enough
   to hold all 4 sides of every point in the grid. */
  int xys = (xsize*ysize)*4;
  
  fprintf(stderr,"bounds: Size of internal grid: %d/%d\n",ysize,xsize);

  /* Allocate memory for arrays */  
  int** blockarray;
  blockarray = (int**) malloc(ysize * sizeof(int*));  
  for (i = 0; i < ysize; i++) blockarray[i] = (int*) malloc(xsize * sizeof(int));  
  
  if (!blockarray){
    fprintf(stderr,"bounds: Failed to allocate needed memory, try increasing the distance value (%f)\n",inc);
    exit(EXIT_FAILURE);
  }

  //fprintf(stderr,"bounds: Initializing bbarray with a value of %d\n",xys);
  line_t *bbarray;
  bbarray = (line_t *)malloc(sizeof(line_t)*xys);

  if (!bbarray) {
    fprintf(stderr,"bounds: Failed to allocate needed memory, try increasing the distance value (%f)\n",inc);
    exit(EXIT_FAILURE);
  }
  //fprintf(stderr,"bounds: bbarray Initialized\n");
  
  /* Loop through the point records and grid them */
  for (i = 0; i < npoints; i++) {
    xpos = fabs((points[i].x - xyzi.xmin) / inc);
    ypos = fabs((points[i].y - xyzi.ymin) / inc);
    blockarray[ypos][xpos] = 1;
  }
  
  /* Loop through the grid and find and record the edge lines */
  for (i = 0; i < ysize; i++) {
    for (j = 0; j < xsize; j++) {
      if (blockarray[i][j] == 1) {
	if (i == 0 || blockarray[i-1][j] != 1) {
	  bb1.x = ((j * inc) + xyzi.xmin), bb1.y = ((i * inc) + xyzi.ymin);
	  bb2.x = ((j * inc) + xyzi.xmin) + inc, bb2.y = ((i * inc) + xyzi.ymin);
	  bbarray[bb].p1 = bb1, bbarray[bb].p2 = bb2;
	  bb++;
	}
	if (j == 0 || blockarray[i][j-1] != 1) {
	  bb1.x = ((j * inc) + xyzi.xmin), bb1.y = ((i * inc) + xyzi.ymin);
	  bb2.x = ((j * inc) + xyzi.xmin), bb2.y = ((i * inc) + xyzi.ymin) + inc;
	  bbarray[bb].p1 = bb1, bbarray[bb].p2 = bb2;
	  bb++;
	}
	if (i == ysize-1 || blockarray[i+1][j] != 1) {
	  bb1.x = ((j * inc) + xyzi.xmin) + inc, bb1.y = ((i * inc) + xyzi.ymin) + inc;
	  bb2.x = ((j * inc) + xyzi.xmin), bb2.y = ((i * inc) + xyzi.ymin) + inc;
	  bbarray[bb].p1 = bb1, bbarray[bb].p2 = bb2;
	  bb++;
	}
	if (j == xsize-1 || blockarray[i][j+1] != 1) {
	  bb1.x = ((j * inc) + xyzi.xmin) + inc, bb1.y = ((i * inc) + xyzi.ymin);
	  bb2.x = ((j * inc) + xyzi.xmin) + inc, bb2.y = ((i * inc) + xyzi.ymin) + inc;
	  bbarray[bb].p1 = bb1, bbarray[bb].p2 = bb2;
	  bb++;
	}
      }
    }
  }

  /* Sort the edge lines into polygons */
  while (bb >= 4) {
    /* Do for each polygon...finishes (done == 1) when finds a matching point. */
    while (done == 0) {
      if (bcount == 0) {
	points[0] = bbarray[bb-1].p1, points[1] = bbarray[bb-1].p2;
	bcount=2, bb--;
      }
      for (edge = 0; edge < bb; edge++) {
	if (fabs( points[bcount-1].x - bbarray[edge].p1.x ) < FLT_EPSILON && fabs( points[bcount-1].y - bbarray[edge].p1.y) < FLT_EPSILON ) {
	  if (fabs( points[0].x - bbarray[edge].p2.x) < FLT_EPSILON && fabs (points[0].y - bbarray[edge].p2.y) < FLT_EPSILON) l = 1;
	  points[bcount] = bbarray[edge].p2;
	  bbarray[edge] = bbarray[bb-1];
	  bb--, bcount++, done = l;
	  break;
	}
	if (fabs (points[bcount-1].x - bbarray[edge].p2.x) < FLT_EPSILON && fabs (points[bcount-1].y - bbarray[edge].p2.y) < FLT_EPSILON) {
	  if (fabs (points[0].x - bbarray[edge].p1.x) < FLT_EPSILON && fabs (points[0].y - bbarray[edge].p1.y) < FLT_EPSILON) l = 1;
	  points[bcount] = bbarray[edge].p1;
	  bbarray[edge] = bbarray[bb-1];  
	  bb--, bcount++, done = l;
	  break;
	}
      }
    }

    /* Print out the polygon */
    for (edge = 0; edge < bcount; edge++) {
      printf("%f %f\n", points[edge].x, points[edge].y);
    }

    /* Reset some values */
    fcount = fcount + (bcount-1);
    bcount = 0, done = 0, l = 0;
    if (bb >= 4) printf(">\n");
  }
  
  fprintf(stderr,"bounds: Found %d total boundary points\n", fcount);

  /* Free allocated grid 'blockarray' from memory */
  for (i = 0; i < ysize; i++) free(blockarray[i]);
  free(blockarray);
  
  /* Free allocated grid bbarray from memory */
  free(bbarray);
  
  return(1);  
}
