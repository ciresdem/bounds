/*------------------------------------------------------------
 * hull.c
 *
 * This file is part of BOUNDS
 *
 * Copyright (c) 2011, 2012, 2018, 2019 Matthew Love <matthew.love@colorado.edu>
 * BOUNDS is liscensed under the GPL v.2 or later and 
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
dist_euclid (point_t* p1, point_t* p2) 
{
  return sqrt (((p2->y - p1->y) * (p2->y - p1->y)) + ((p2->x - p1->x) * (p2->x - p1->x)));
}

/* Calculate the angle between two points, from the previous point
 * p1 = current point, p2 = next point, p0 = previous point
 */
static float
atheta (point_t* p1, point_t* p2, point_t* p0) 
{
  float t, t2, t3, d12x, d12y, d13x, d13y;

  d12x = p1->x - p0->x, d12y = p1->y - p0->y; 
  d13x = p1->x - p2->x, d13y = p1->y - p2->y; 

  t = atan2 (d12y, d12x), t2 = atan2 (d13y, d13x);
  t3 = (t2 - t);

  if (t3 < 0) 
    t3 = t3 + (2 * M_PI);
  return t3;
}

static float
atheta_degrees (point_t p1, point_t p2, point_t p0) 
{
  float t, t2, t3, d12x, d12y, d13x, d13y;

  d12x = p1.x - p0.x, d12y = p1.y - p0.y; 
  d13x = p1.x - p2.x, d13y = p1.y - p2.y; 

  t = atan2 (d12y, d12x), t2 = atan2 (d13y, d13x);

  t3 = (t2 - t) * 180 / M_PI;

  if (t3 < 0) 
    t3 = t3 + 360;
  return t3;
}

/* Return the angle between the two given points 
 */
static float
theta (point_t* p1, point_t* p2) 
{
  float t = atan2 (p2->y - p1->y, p2->x - p1->x);;
  if (t < 0) 
    t = t + (2 * M_PI);
  return t;
}

static float
theta_degrees (point_t* p1, point_t* p2) 
{
  float t = atan2 (p2->y - p1->y, p2->x - p1->x) * 180 / M_PI;
  if (t < 0) 
    t = t + 360;
  return t;
}
 
/* Three points are a counter-clockwise turn if ccw > 0, clockwise if
 * ccw < 0, and collinear if ccw = 0 because ccw is a determinant that
 * gives the signed area of the triangle formed by p1, p2 and p3.
 */
static double
ccw2 (point_t* p1, point_t* p2, point_t* p3) 
{
  return (p2->x - p1->x) * (p3->y - p2->y) - (p2->y - p1->y) * (p3->x - p2->x);
}

int
ccw (point_t* p1, point_t* p2, point_t* p3) 
{
  float ccwv = (p2->y - p1->y) * (p3->x - p2->x) - (p2->x - p1->x) * (p3->y - p2->y);
  if (ccwv == 0) return 0;
  return (ccwv > 0)? 1: -1;
}

int
on_line_p (point_t* p1, line_t* l1)
{
  if (p1->x <= max (l1->p1.x, l1->p2.x) && p1->x >= min (l1->p1.x, l1->p2.x) &&
      p1->y <= max (l1->p1.y, l1->p2.y) && p1->y >= min (l1->p1.y, l1->p2.y))
    return 1;

  return 0;
}

/* int */
/* on_line_p (point_t* p1, point_t* p2, point_t* p3) */
/* { */
/*   if (p2->x <= max (p1->x, p3->x) && p2->x >= min (p1->x, p3->x) && */
/*       p2->y <= max (p1->y, p3->y) && p2->y >= min (p1->y, p3->y)) */
/*     return 1; */
/*   return 0; */
/* } */

int
intersect_p (line_t l1, line_t l2)
{
  int a, b, c, d;

  a = ccw (&l1.p1, &l1.p2, &l2.p1), b = ccw (&l1.p1, &l1.p2, &l2.p2);
  c = ccw (&l2.p1, &l2.p2, &l1.p1), d = ccw (&l2.p1, &l2.p2, &l1.p2);

  if (a != b && c != d) return 1;

  if (a == 0 && on_line_p (&l2.p1, &l1)) return 1;
  if (b == 0 && on_line_p (&l2.p2, &l1)) return 1;
  if (c == 0 && on_line_p (&l1.p1, &l2)) return 1;
  if (d == 0 && on_line_p (&l1.p2, &l2)) return 1;

  return 0;

}

/* Return 1 if point p1 is inside polygon poly, otherwise return 0
 */
int
inside_p (point_t* p1, point_t* poly, ssize_t hullsize) 
{
  int k = 0, l, l1;
  ssize_t i;
  point_t p2;
  line_t lt, lp;

  p2.y = p1->y, p2.x = FLT_MAX;
  lt.p1 = *p1, lt.p2 = p2;
  l = 0;

  for (i = 0; i < hullsize; i++) 
    {
      lp.p1 = poly[i], lp.p2 = poly[i+1];
    
      if (lt.p1.y == lp.p2.y || lt.p1.y == lp.p1.y) l++;

      if (intersect_p (lt, lp)) k++;
    }

  if (l == 2) return 1;

  return k & 1;
}

/* A 'package-wrap' concavehull 
 * -- Retruns the number of points in the boundary;
 * The hull makes up the begining of the points array.
 * Generate the concave hull using the given distance threshold
 * Note: The input points must be sorted first.
 */
ssize_t
dpw_concave (point_t* points, int npoints, double d) 
{
  int i, min, M, k, j, np;
  float th, cth, fth;
  point_t t, t0;
  line_t l1, l2;
  int rein = 0;

  /* Find the minimum y value in the dataset */
  for (min = 0, i = 1; i < npoints; i++)
    if (points[i].y < points[min].y) 
      min = i;

  points[npoints].x = -9999, points[npoints].y = -9999;
  np = npoints;

  /* Loop through all the points, starting with the point found above. */
  for (M = 0; M < np; M++) 
    {
      t = points[M], points[M] = points[min], points[min] = t;
      min = -1, k = 0, th = 2*M_PI;
      
      /* Set the previous point for angle calculations */
      if (M > 0) t0 = points[M - 1];
      
      /* Re-insert the first point into the dataset */
      if (rein == 0)
	{
	  //fth = atheta (&points[M], &points[0], &t0);
	  //if (fth > 0 && fth < th)
	  if (M > 10)
	    {
	      points[npoints] = points[0];
	      rein = 1, np++;
	    }
	}
      
      /* Loop through the remaining points and find the next concave point; 
	 less than distance threshold -> less than working theta -> does not 
	 intersect existing boundary */
      for (i = M + 1; i < np; i++) 
	if (dist_euclid (&points[M], &points[i]) <= d) 
	  {
	    if (M == 0) 
	      cth = theta (&points[M], &points[i]);
	    else 
	      cth = atheta (&points[M], &points[i], &t0);
	    if (cth > FLT_EPSILON && cth <= th) 
	      {
		/* Make sure the selected point doesn't create a line segment which
		   intersects with the existing hull. */
		//if (M > 0)
		l1.p1 = points[M], l1.p2 = points[i];

		for (k = 0, j = 1; j < M - 1; j++) 
		  {
		    l2.p1 = points[j], l2.p2 = points[j + 1];
		    if (intersect_p (l1, l2)) k++;

		  }

		/* If all criteria are met,, add this point to the hull */
		if (k == 0) 
		    min = i, th = cth; 
	      }
	  }

      /* No point was found, try again with a larger distance threshhold. */
      if (min == -1) 
	return min;
      
      /* The first point was found again, a successful hull was found! end here. */
      if (min == npoints) 
	{
	  /* matched the null point... */
	  if (points[min].x == -9999 || points[min].y == -9999)
	    return -1;
	  points[M + 1] = points[min];
	  return M + 1;
	}
    }
}

/* 
 * A Monotone-Chain Convex Hull
 * -- Returns a list of points on the convex hull in counter-clockwise order.
 * Note: The last point in the returned list is the same as the first one. 
 * Note: The input points must be sorted first.
 */
void
mc_convex (point_t* points, ssize_t npoints, point_ptr_t** out_hull, ssize_t* out_hullsize) 
{
  point_ptr_t* hull;
  ssize_t i, t, k = 0;

  hull = *out_hull;

  /* lower hull */
  for (i = 0; i < npoints; ++i) 
    {
      while (k >= 2 && ccw (hull[k - 2], hull[k - 1], &points[i]) <= 0) 
	--k;
      hull[k++] = &points[i];
    }
 
  /* upper hull */
  for (i = npoints - 2, t = k + 1; i >= 0; --i) 
    {
      while (k >= t && ccw (hull[k - 2], hull[k - 1], &points[i]) <= 0) 
	--k;
      hull[k++] = &points[i];
    }
  
  *out_hull = hull;
  *out_hullsize = k;
}

/* A 'package-wrap' Convex Hull 
 * -- Returns the number of points in the hull;
 * The hull makes up the begining of the points array.
 */
ssize_t
pw_convex (point_t* points, ssize_t npoints) 
{
  ssize_t i, min, M;
  float th = 0.0;
  float cth, v;
  point_t t;

  for (min = 0, i = 1; i < npoints; i++)
    if (points[i].y < points[min].y) min = i;

  for (M = 0; M < npoints; M++) 
    {
      t = points[M], points[M] = points[min]; 
      points[min] = t, min = npoints; 
      v = th, th = 2 * M_PI;
      
      for (i = M + 1; i <= npoints; i++) 
	{
	  cth = theta(&points[M], &points[i]);
	  if (cth > v)
	    if (cth < th)
	      min = i, th = cth; 
	}
      if (min == npoints) 
	{
	  points[M + 1] = points[0];
	  return M + 2;
	}
    }
}
