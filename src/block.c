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

int
pnts_equal_p(point_t p1, point_t p2){
  return fabs( p1.x - p2.x ) < FLT_EPSILON && fabs( p1.y - p2.y ) < FLT_EPSILON;
}

int
xyz_info_valid_p(xyz_info *region) {
  if (region->xmin >= region->xmax) return 0;
  if (region->ymin >= region->ymax) return 0;
  return 1;
}

/* 
 * "Bounding Block"
 * Generates a grid at `inc` cell-size and polygonizes it.
 */
int
block_pts(point_t* points, int npoints, double inc, xyz_info region, int vflag) {

  int i, j, xpos, ypos, edge;
  int bb = 0, done = 0, bcount = 0, l = 0, fcount = 0;
  double x, y;
  point_t bb1, bb2;

  /* Gather xyzinfo */
  xyz_info xyzi;

  if (xyz_info_valid_p(&region)) {
    xyzi = region;
    if (vflag > 0) fprintf(stderr,"bounds: Using user supplied region: %f/%f/%f/%f\n",region.xmin,region.xmax,region.ymin,region.ymax);
  } else {
    minmax(points, npoints, &xyzi);
  }

  /* Set the rows and columns of the internal grid */
  int ysize = fabs((xyzi.ymax - xyzi.ymin) / inc);
  int xsize = fabs((xyzi.xmax - xyzi.xmin) / inc);
  /* `xys` is the size of the bbarray array, this has to be large enough
   to hold all 4 sides of every point in the grid. */
  ssize_t xys = (xsize*ysize)*4;
  
  if (vflag > 0) fprintf(stderr,"bounds: Size of internal grid: %d/%d\n",ysize,xsize);

  /* Allocate memory for arrays */  
  int** blockarray;
  blockarray = (int**) malloc(ysize * sizeof(int*));  
  for (i = 0; i < ysize; i++) blockarray[i] = (int*) malloc(xsize * sizeof(int));  
  
  if (!blockarray){
    if (vflag > 0) fprintf(stderr,"bounds: Failed to allocate needed memory, try increasing the distance value (%f)\n",inc);
    exit(EXIT_FAILURE);
  }

  if (vflag > 0) fprintf(stderr,"bounds: Initializing bbarray with a value of %d\n",xys);

  line_t *bbarray;
  bbarray = (line_t *) malloc(sizeof(line_t)*xys);

  if (!bbarray) {
    if (vflag > 0) fprintf(stderr,"bounds: Failed to allocate needed memory, try increasing the distance value (%f)\n",inc);
    exit(EXIT_FAILURE);
  }

  if (vflag > 0) fprintf(stderr,"bounds: Gridding points\n");

  /* Loop through the point records and grid them */
  for (i = 0; i < npoints; i++) {
    xpos = fabs((points[i].x - xyzi.xmin) / inc);
    ypos = fabs((points[i].y - xyzi.ymin) / inc);
    if (xpos < xsize && ypos < ysize) blockarray[ypos][xpos] = 1;
  }

  if (vflag > 0) fprintf(stderr,"bounds: points gridded\n");
  
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

  for (i = 0; i < ysize; i++) {
    free(blockarray[i]);
    blockarray[i]=NULL;
  }
  free(blockarray);
  blockarray=NULL;
  
  if (vflag > 0) fprintf(stderr,"bounds: grid edges recorded and blockarray freed\n");

  /* Allocate memory for the edge points. 
     There may be more boundary points than input points, 
     so a new array is used. 
     Add one to the point count to repeat the first point.
  */
  point_t* bnds;  
  bnds = (point_t*) malloc(sizeof(point_t) * (bb+1));

  if (vflag > 0) fprintf(stderr,"bounds: sorting %d points into boundary\n", bb);

  /* Sort the edge lines into polygons */
  while (bb >= 4) {
    /* Do for each polygon...finishes (done == 1) when finds a matching point. */
    while (done == 0) {
      if (bcount == 0) {
	bnds[0] = bbarray[bb-1].p1, bnds[1] = bbarray[bb-1].p2;
	bcount=2, bb--;
      }
      for (edge = 0; edge < bb; edge++) {
	if (pnts_equal_p( bnds[bcount-1],  bbarray[edge].p1 )) {
	  if (pnts_equal_p( bnds[0], bbarray[edge].p2 )) l = 1;
	  bnds[bcount] = bbarray[edge].p2;
	  bbarray[edge] = bbarray[bb-1];
	  bb--, bcount++, done = l;
	  break;
	}
	else if (pnts_equal_p( bnds[bcount-1], bbarray[edge].p2 )) {
	  if (pnts_equal_p( bnds[0], bbarray[edge].p1 )) l = 1;
	  bnds[bcount] = bbarray[edge].p1;
	  bbarray[edge] = bbarray[bb-1];  
	  bb--, bcount++, done = l;
	  break;
	}
      }
    }

    /* Print out the polygon */
    for (edge = 0; edge < bcount; edge++) {
      printf("%f %f\n", bnds[edge].x, bnds[edge].y);
    }

    /* Reset some values */
    fcount = fcount + bcount;
    bcount = 0, done = 0, l = 0;
    if (bb >= 4) printf(">\n");
  }
  
  if (vflag > 0) fprintf(stderr,"bounds: Found %d total boundary points\n", fcount);

  /* Free allocated grid bbarray from memory */
  free(bbarray);
  bbarray=NULL;
  free(bnds);
  bnds=NULL;
  
  return(1);  
}
