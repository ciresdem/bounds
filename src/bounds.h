/*------------------------------------------------------------
 * bounds : version 0.4.3
 *
 * This file is part of BOUNDS
 *
 * Copyright (c) 2011, 2012, 2016, 2018 Matthew Love <matthew.love@colorado.edu>
 * BOUNDS is liscensed under the GPL v.2 or later and 
 * is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details. 
 * <http://www.gnu.org/licenses/> 
 *--------------------------------------------------------------*/

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#define BOUNDS_VERSION "0.4.3"

#define MAX_RECORD_LENGTH 1024
#define MAX_HULLS 100000

#ifndef INFINITY
#define INFINITY (1.0 / 0.0)
#endif 

#ifndef NAN
#define NAN (0.0 / 0.0)
#endif 

typedef struct {
  double x;
  double y;
} point_t;
 
typedef point_t* point_ptr_t;

typedef struct {
  point_t p1, p2;
} line_t;

typedef line_t* line_ptr_t;

typedef struct {
  double xmax;
  double ymax;
  double xmin;
  double ymin;
} xyz_info;

/* Line-Count 
 */
ssize_t
linecnt(FILE *in);

/* Scan an xyz line and attempt to determine the delimiter and record length 
 */
int
scanline(char* infile);

/* Read a point record from the given infile
 */
int
read_point(FILE *infile, point_t *rpnt, char* delimiter, char* pnt_recr);

void
minmax(point_t* points, int npoints, xyz_info *xyzi);

/* Get the euclidean distance between p1 and p2 in meters
 */
static float
dist_euclid(point_t* p1, point_t* p2);

/* Calculate the angle between two points, from the previous point
 * p1 = current point, p2 = next point, p0 = previous point
 */
static float
atheta(point_t* p1, point_t* p2, point_t* p0);

static float
atheta_degrees(point_t p1, point_t p2, point_t p0);

/* Return the angle between the two given points 
 */
static float
theta(point_t* p1, point_t* p2);

static float
theta_degrees(point_t* p1, point_t* p2);

/* Return 1 if l1 and l2 intersect, -1 if they don't
 */
int
intersect(line_t l1, line_t l2, double thresh);

/* Return 1 if p1 is inside poly
 */
int
inside(point_t* p1, point_t* poly, ssize_t hullsize, double dist);

/* A 'package-wrap' concavehull 
 * -- Retruns the number of points in the boundary;
 * The hull makes up the begining of the points array.
 * Generate the concave hull using the given distance threshold
 */
ssize_t
dpw_concave(point_t* points, int npoints, double d);

/* Returns a list of points on the convex hull in counter-clockwise order.
 * Note: the last point in the returned list is the same as the first one. 
 */
void
mc_convex(point_t* points, ssize_t npoints, point_ptr_t** out_hull, ssize_t* out_hullsize);

/* A 'package-wrap' convexhull 
 * -- Retruns the number of points in the hull;
 * The hull makes up the begining of the points array.
 */
ssize_t
pw_convex(point_t* points, ssize_t npoints);

/* Generate a boundary with blocks
 * `inc` is the blocksize in the same units
 * as the input points.
 */
int
block_pts(point_t* points, int npoints, double inc);

// End
