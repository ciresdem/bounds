/*------------------------------------------------------------
 * hull.c
 *
 * Copyright (c) 2011-2018 Matthew Love <matthew.love@colorado.edu>
 * This file is liscensed under the GPL v.2 or later and 
 * is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details. 
 * <http://www.gnu.org/licenses/> 
 *--------------------------------------------------------------*/

#include "bounds.h"

/* Get the euclidean distance between p1 and p2 in meters
 */
static float
dist_euclid(point_t* p1, point_t* p2) {
  return sqrt( ((p2->y - p1->y) * (p2->y - p1->y)) + ((p2->x - p1->x) * (p2->x - p1->x)) );
}

/* Calculate the angle between two points, from the previous point
 * p1 = current point, p2 = next point, p0 = previous point
 */
static float
atheta(point_t* p1, point_t* p2, point_t* p0) {
  float t, t2, t3, d12x, d12y, d13x, d13y;

  d12x = p1->x - p0->x, d12y = p1->y - p0->y; 
  d13x = p1->x - p2->x, d13y = p1->y - p2->y; 

  t = atan2(d12y, d12x), t2 = atan2(d13y, d13x);

  t3 = (t2 - t);

  if (t3 < 0) t3 = t3 + (2 * M_PI);
  return t3;
}

static float
atheta_degrees(point_t p1, point_t p2, point_t p0) {
  float t, t2, t3, d12x, d12y, d13x, d13y;

  d12x = p1.x - p0.x, d12y = p1.y - p0.y; 
  d13x = p1.x - p2.x, d13y = p1.y - p2.y; 

  t = atan2(d12y, d12x), t2 = atan2(d13y, d13x);

  t3 = (t2 - t) * 180 / M_PI;

  if (t3 < 0) t3 = t3 + 360;
  return t3;
}

/* Return the angle between the two given points 
 */
static float
theta(point_t* p1, point_t* p2) {
  float t = atan2(p2->y - p1->y, p2->x - p1->x);;
  if (t < 0) t = t + (2 * M_PI);
  return t;
}

static float
theta_degrees(point_t* p1, point_t* p2) {
  float t = atan2(p2->y - p1->y, p2->x - p1->x) * 180 / M_PI;
  if (t < 0) t = t + 360;
  return t;
}
 
/* Three points are a counter-clockwise turn if ccw > 0, clockwise if
 * ccw < 0, and collinear if ccw = 0 because ccw is a determinant that
 * gives the signed area of the triangle formed by p1, p2 and p3.
 */
static double
ccw(point_t* p1, point_t* p2, point_t* p3) {
  return (p2->x - p1->x) * (p3->y - p1->y) - (p2->y - p1->y) * (p3->x - p1->x);
}

int
intersect(line_t l1, line_t l2, double thresh) {
  double a, b;
  a = ccw(&l1.p1, &l1.p2, &l2.p1) *ccw(&l1.p1, &l1.p2, &l2.p2);
  b = ccw(&l2.p1, &l2.p2, &l1.p1) *ccw(&l2.p1, &l2.p2, &l1.p2);

  if (a < 0-thresh && b < 0-thresh) return 1;
  else if (a >= 0-thresh && a <= thresh) return 0;
  else if (b >= 0-thresh && b <= thresh) return 0;
  else return -1;
}

/* Depreciated intersect function. 
 */
static int
intersectp(point_t* p1, point_t* p2, point_t* p3, point_t* p4) {
  float den = (p2->x - p1->x) * (p4->y - p3->y) - (p2->y - p1->y) * (p4->x - p3->x);

  //if (fabs(den) < FLT_EPSILON) return 0; // parallel lines
  if (!den) return 0; // parallel lines

  float a = ((p3->x - p1->x) * (p4->y - p3->y) - (p3->y - p1->y) * (p4->x - p3->x)) / den;
  if (a <= FLT_EPSILON  || a >= 1.0) return 0; // intersection point not between p1 and p2

  a = ((p3->x - p1->x) * (p2->y - p1->y) - (p3->y - p1->y) * (p2->x - p1->x)) / den;
  if (a <= FLT_EPSILON || a >= 1.0) return 0; // intersection point not between p3 and p4

  return 1;
}

int
inside(point_t* p1, point_t* poly, ssize_t hullsize, double dist) {

  int k = 0, l, l1;
  ssize_t i;
  point_t p2;
  line_t lt, lp;

  p2.y = p1->y, p2.x = INT_MAX;
  lt.p1 = *p1, lt.p2 = p2;
  
  for (i = 0; i < hullsize; i++) {
    lp.p1 = poly[i], lp.p2 = poly[i+1];

    l = intersect(lt,lp,FLT_EPSILON*dist);

    if (l == 0) return 1;
    if (l == 1) k++;
  }
  return k & 1;
}

/* A 'package-wrap' convexhull 
 * -- Retruns the number of points in the boundary;
 * The hull makes up the begining of the points array.
 * Generate the concave hull using the given distance threshold
 */
ssize_t
dpw_concave(point_t* points, int npoints, double d) {
  int i, min, M, k, j, fp=0;
  float th, cth;
  point_t t, t0;
  line_t l1, l2;

  /* Find the minimum y value in the dataset */
  for (min = 0, i = 1; i < npoints-1; i++)
    if (points[i].y < points[min].y) min = i;

  /* Loop through all the points, starting with the point found above. */
  for (M = 0; M < npoints; M++) {
    t = points[M], points[M] = points[min], points[min] = t;
    min = -1, k = 0, th = 2*M_PI;

    /* Set the previous point for angle calculations */
    if (M > 0) t0 = points[M - 1];

    /* 
     *  !! Need a better way to re-insert first point !! 
     */
    /* Re-insert the first point into the dataset */
    if (M > 10 || M > npoints*.1) points[npoints-1] = points[0];
 
    /* Loop through the remaining points and find the next concave point; 
       less than distance threshold -> less than working theta -> does not 
       intersect existing boundary */
    for (i = M + 1; i < npoints; i++)
      if (dist_euclid(&points[M], &points[i]) <= d) {
	if (M == 0) cth = theta(&points[M], &points[i]);
	else cth = atheta(&points[M], &points[i], &t0);
	if (cth > 0  && cth < th) {
	  /* Make sure the selected point doesn't create a line segment which
	     intersects with the existing hull. */
	  if (M > 0)
	    for (k = 0, j = 0; j < M; j++) {
	      l1.p1 = points[M], l1.p2 = points[i];
	      l2.p1 = points[j], l2.p2 = points[j+1];
	      if (intersect(l1,l2,DBL_EPSILON) > 0) k++;
	    }
	  /* If all criteria are met,, add this point to the hull */
	  if (k == 0) min = i, th = cth; 
	}
      }

    /* No point was found, try again with a larger distance threshhold. */
    if (min == -1) return min;

    /* The first point was found again, a successful hull was found! end here. */
    if (min == npoints-1) {
      points[M + 1] = points[min];
      return M+1;
    }
  }
}

/* 
 * A Monotone-Chain Convex Hull
 * -- Returns a list of points on the convex hull in counter-clockwise order.
 * Note: the last point in the returned list is the same as the first one. 
 */
void
mc_convex(point_t* points, ssize_t npoints, point_ptr_t** out_hull, ssize_t* out_hullsize) {
  point_ptr_t* hull;
  ssize_t i, t, k = 0;

  hull = *out_hull;

  /* lower hull */
  for (i = 0; i < npoints; ++i) {
    while (k >= 2 && ccw(hull[k - 2], hull[k - 1], &points[i]) <= 0) --k;
    hull[k++] = &points[i];
  }
 
  /* upper hull */
  for (i = npoints - 2, t = k + 1; i >= 0; --i) {
    while (k >= t && ccw(hull[k - 2], hull[k - 1], &points[i]) <= 0) --k;
    hull[k++] = &points[i];
  }
  
  *out_hull = hull;
  *out_hullsize = k;
}

/* A 'package-wrap' Convex Hull 
 * -- Retruns the number of points in the hull;
 * The hull makes up the begining of the points array.
 */
ssize_t
pw_convex(point_t* points, ssize_t npoints) {
  ssize_t i, min, M;
  float th = 0.0;
  float cth, v;
  point_t t;

  for (min = 0, i = 1; i < npoints; i++)
    if (points[i].y < points[min].y) min = i;

  for (M = 0; M < npoints; M++) {
    t = points[M], points[M] = points[min]; 
    points[min] = t, min = npoints; 
    v = th, th = 2 * M_PI;

    for (i = M+1; i <= npoints; i++) {
      cth = theta(&points[M], &points[i]);
      if (cth > v)
	if (cth < th)
	  min = i, th = cth; 
    }
    if (min == npoints) {
      points[M + 1] = points[0];
      return M + 2;
    }
  }
}
