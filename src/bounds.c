/*------------------------------------------------------------
 * bounds.c
 * 
 * This file is part of BOUNDS
 *
 * Copyright (c) 2011, 2012, 2016 - 2023 Matthew Love <matthew.love@colorado.edu>
 * BOUNDS is liscensed under the GPL v.2 or later and 
 * is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details. 
 * <http://www.gnu.org/licenses/> 
 *--------------------------------------------------------------*/

#include <getopt.h>
#include "bounds.h"

/* Flags set by `--version' and `--help' */
static int version_flag;
static int help_flag;
static int verbose_flag;

static void
print_version(const char* command_name, const char* command_version) 
{
  fprintf (stderr, "%s %s \n\
Copyright © 2011, 2012, 2016 - 2023 Matthew Love <matthew.love@colorado.edu> \n\
%s is liscensed under the GPL v.2 or later and is \n\
distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;\n\
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A\n\
PARTICULAR PURPOSE.  See the GNU General Public License for more details.\n\
<http://www.gnu.org/licenses/>\n", command_name, command_version, command_name);
  exit (0);
}

static void 
usage () 
{
  fprintf (stderr, "\
Usage: bounds [OPTION]... [FILE]\n\
Generate a boundary of the set of xy points from FILE, or standard input, to standard output.\n\
\n\
  ---- xy i/o ----\n\n\
  -d, --delimiter\tThe input xy file record delimiter.\n\
                 \tIf omitted, the delimiter will be guessed from the first line read.\n\
  -g, --gmt\t\tFormat output as GMT vector multipolygon; Use twice to supress\n\
           \t\tthe initial header (e.g. -gg).\n\
  -j, --json\t\tFormat output as GeoJSON vector multipolygon; Use twice to supress\n\
           \t\tthe initial header (e.g. -jj).\n\
  -n, --name\t\tThe output layer name (only used with -g or -j).\n\
  -r, --record\t\tThe input record order, 'xy' should represent the locations\n\
              \t\tof the x and y records, respectively (e.g. --record zdyx).\n\
  -s, --skip\t\tThe number of lines to skip from the input.\n\n\
  ---- bounds ----\n\n\
  -b, --box\t\t'Bounding Box' boundary. \n\
  -k, --block\t\t'Bounding Block' boundary. Specify the blocking increment\n\
             \t\tin input units (e.g. --block 0.001). Specify a blocking region\n\
             \t\tafter the increment if desired (e.g. --block 0.001/west/east/south/north).\n\
  -v, --concave\t\t'Concave Hull' boundary using a distance weighted package wrap algorithm.\n\
               \t\tSpecify distance value or - to estimate appropriate distance.\n\
  -x, --convex\t\t'Convex Hull' boundary using a monotone chain algorithm. [default]\n\
              \t\tUse twice to use a package wrap algorithm (e.g. -xx).\n\n\
  ---- et cetra ----\n\n\
      --verbose\t\tincrease the verbosity.\n\
      --help\t\tprint this help menu and exit.\n\
      --version\t\tprint version information and exit.\n\n\
\n\
With no FILE, or when FILE is --, read standard input.\n\
All OPTION values must be in the same units as the input xy data.\n\n\
Examples:\n\
  bounds \t\toutput a convex hull from standard input.\n\
  bounds -g -k0.0001\toutput a GMT formatted 'block' boundary from standard input.\n\
  bounds -v10 -d,\toutput a concave hull from comma-delimited standard input.\n\
  bounds -v- in.xyz\toutput a concave hull from file in.xyz\n\
\n\
Report bugs to <matthew.love@colorado.edu>\n\
CIRES DEM home page: <http://ciresgroups.colorado.edu/coastalDEM>\n\
");
  exit (0);
}

/* Compare functions for use in qsort
 */
int
compare (const void* a, const void* b) 
{
  const point_t *elem1 = a;    
  const point_t *elem2 = b;

  if (elem1->x < elem2->x) 
    return -1;
  else if (elem1->x > elem2->x) 
    return 1;
  else 
    return 0;
}

int
sort_min_y (const void* a, const void* b) 
{
  const point_t *elem1 = a;    
  const point_t *elem2 = b;

  if (elem1->y < elem2->y) 
    return -1;
  else if (elem1->y > elem2->y) 
    return 1;
  else 
    return 0;
}

int
sort_max_y (const void* a, const void* b) 
{
  const point_t *elem1 = a;    
  const point_t *elem2 = b;

  if (elem1->y > elem2->y) 
    return -1;
  else if (elem1->y < elem2->y) 
    return 1;
  else 
    return 0;
}

int
sort_max_x (const void* a, const void* b) 
{
  const point_t *elem1 = a;    
  const point_t *elem2 = b;

  if (elem1->x > elem2->x) 
    return -1;
  else if (elem1->x < elem2->x) 
    return 1;
  else 
    return 0;
}

int
sort_min_x (const void* a, const void* b) 
{
  const point_t *elem1 = a;    
  const point_t *elem2 = b;

  if (elem1->x < elem2->x) 
    return -1;
  else if (elem1->x > elem2->x) 
    return 1;
  else 
    return 0;
}

int
region_valid_p (region_t *region) 
{
  if (region->xmin >= region->xmax) 
    return 0;

  if (region->ymin >= region->ymax) 
    return 0;

  return 1;
}

/* Quick and Dirty Density
 * TODO: use a convex hull instaed of bb.
 */
double
qadd (point_t* points, int npoints)
{
  int i, xmin, xmax, ymin, ymax;
  line_t l1, l2, l3, l4;
  double w, l;
  
  for (ymin = 0, ymax = 0, xmin = 0, xmax = 0, i = 1; i < npoints; i++) 
    {
      if (points[i].y < points[ymin].y) ymin = i;
      if (points[i].x < points[xmin].x) xmin = i;
      if (points[i].y > points[ymax].y)	ymax = i;
      if (points[i].x > points[xmax].x)	xmax = i;
    }

  w = points[xmax].x - points[xmin].x;
  l = points[ymax].y - points[ymin].y;

  return ((w * l) / npoints) ;
}

int
main (int argc, char **argv) 
{
  char* fn;
  FILE* fp;

  int c, i, status, min, j;
  int inflag = 0, vflag = 0, sflag = 0, dflag = 0, pc = 0, sl = 0;
  int cflag = 0, kflag = 0, bflag = 0, gmtflag = 0, jsonflag = 0, nflag = 0;
  double dist;

  point_t rpnt, pnt;
  point_ptr_t hull0[MAX_HULLS];
  point_ptr_t* hull = hull0;
  ssize_t hullsize;
  ssize_t npr = 0;

  //char* delim = " \t";
  char* delim;
  char* ptrec = "xy";
  char* kreg = "";
  char* lname = "bounds";
  
  while (1) 
    {
      static struct option long_options[] =
	{
	  /* These options set a flag. */
	  {"version", no_argument, &version_flag, 1},
	  {"help", no_argument, &help_flag, 1},
	  {"verbose", no_argument, &verbose_flag, 1},
	  /* These options don't set a flag.
	     We distinguish them by their indices. */
	  {"delimiter", required_argument, 0, 'd'},
	  {"skip", required_argument, 0, 's'},
	  {"name", required_argument, 0, 'n'},
	  {"gmt", no_argument, 0, 'g'},
	  {"json", no_argument, 0, 'j'},
	  {"record", required_argument, 0, 'r'},
	  {"box", no_argument, 0, 'b'},
	  {"block", required_argument, 0, 'k'},
	  {"convex", no_argument, 0, 'x'},
	  {"concave", required_argument, 0, 'v'},
	  {0, 0, 0, 0}
	};
      /* getopt_long stores the option index here. */
      int option_index = 0;
      
      c = getopt_long (argc, argv, "gjd:n:r:s:bk:xv:",
		       long_options, &option_index);
    
      /* Detect the end of the options. */
      if (c == -1) break;
      
      switch (c) {
      case 0:
	/* If this option set a flag, do nothing else now. */
	if (long_options[option_index].flag != 0)
	  break;
	printf ("option %s", long_options[option_index].name);
	if (optarg)
	  printf (" with arg %s", optarg);
	printf ("\n");
	break;
      case 'd':
	dflag++;
	delim = optarg;
	break;
      case 'r':
	ptrec = optarg;
	break;
      case 's':
	sflag++;
	sl = atoi(optarg);
	break;
      case 'n':
	nflag++;
	lname = optarg;
	break;
      case 'g':
	gmtflag++;
	break;
      case 'j':
	jsonflag++;
	break;
      case 'b':
	bflag++;
	break;
      case 'k':
	kflag++;
	kreg = optarg;
	//dist = atof(optarg);
	break;
      case 'x':
	cflag++;
	break;
      case 'v':
	vflag++;
	dist = atof(optarg);
	break;
	
      case '?':
	/* getopt_long already printed an error message. */
	fprintf (stderr, "Try 'bounds --help' for more information.\n");
	exit (1);
	break;
	
      default:
	abort ();
      }
    }
  
  if (version_flag)
    {
      //print_version("bounds", BOUNDS_VERSION);
      printf("%s\n", BOUNDS_VERSION);
      exit (0);
    }

  if (help_flag) 
    usage();

  fn = argv[optind];
  if (fn) 
    inflag++;

  if (inflag > 0) 
    {
      fp = fopen (fn, "r");
      if (!fp) 
	{
	  fprintf (stderr,"bounds: failed to open file: %s\n", fn);
	  exit (1);
	}
    } 
  else 
    {
      fp = stdin;
      fn = "stdin";
    }

  if (verbose_flag > 0) 
    fprintf (stderr, "bounds: working on file: %s\n", fn);
  
  /* Allocate memory for the `pnts` and `pnts2` array using the total number of points.
     `pnts2` is only used by concave and gets allocated there.
  */
  point_t* pnts;
  pnts = (point_t*) malloc (sizeof (point_t));

  /* This is for the GMT compatibility. More can be done here.
   */
  if (gmtflag == 1)
    {
      printf ("# @VGMT1.0 @GMULTIPOLYGON\n# @NName\n# @Tstring\n# FEATURE_DATA\n");
      printf (">\n# @D%s\n# @P\n", lname);
    }
  else if (gmtflag == 2)
    printf (">\n# @D%s\n# @P\n", lname);
  else if (gmtflag == 3)
    {
      printf ("# @VGMT1.0 @GMULTIPOLYGON\n# @NName\n# @Tstring\n# FEATURE_DATA\n");
      exit(0);
    }
  else if (jsonflag == 1)
    {
      printf ("{ \"type\": \"FeatureCollection\",\n\"features\": [\n");
      printf ("{ \"type\": \"Feature\", \"properties\": { \"Name\": \"%s\" },", lname);
      printf (" \"geometry\": { \"type\": \"MultiPolygon\",\n \"coordinates\": [[[");
    }
  else if (jsonflag == 2)
    {
      printf ("{ \"type\": \"Feature\", \"properties\": { \"Name\": \"%s\" },", lname);
      printf (" \"geometry\": { \"type\": \"MultiPolygon\",\n \"coordinates\": [[[");
    }
  else if (jsonflag == 3)
    {
      printf ("");
    }
  else
    printf (">\n");
    
  /* The default is a convex hull -- `cflag` */
  if (cflag == 0 && vflag == 0 && bflag == 0 && kflag == 0) 
    cflag++;
  
  /* Monotone Chain Convex Hull Algorithm - -*Default*-
   */
  if (cflag == 1) 
    {
      load_pnts (fp, &pnts, &npr, ptrec, verbose_flag);
      qsort (pnts, npr, sizeof (point_t), compare);
      mc_convex (pnts, npr, &hull, &hullsize);

      if (jsonflag > 0)
	{
	  for (i = 0; i < hullsize-1; i++)
	    printf ("[%f, %f], ", hull[i]->x, hull[i]->y);
	  printf("[%f, %f]", hull[hullsize-1]->x, hull[hullsize-1]->y);
	}
      else
	{
	  for (i = 0; i < hullsize; i++)
	    printf ("%f %f\n", hull[i]->x, hull[i]->y);
	}

      if (verbose_flag > 0) 
	fprintf (stderr, "bounds: found %d convex boundary points.\n", hullsize);
    }
  
  /* 'Package Wrap' Convex Hull Algorithm 
   */
  else if (cflag == 2)
    {
      load_pnts (fp, &pnts, &npr, ptrec, verbose_flag);
      hullsize = pw_convex (pnts, npr);

      if (jsonflag > 0)
	{
	  for (i = 0; i < hullsize; i++)
	    printf ("[%f, %f], ", pnts[i].x, pnts[i].y);
	  printf ("[%f, %f]", pnts[hullsize].x, pnts[hullsize].y);
	}
      else
	{
	  for (i = 0; i <= hullsize; i++)
	    printf ("%f %f\n", pnts[i].x, pnts[i].y);
	}
      
      if (verbose_flag > 0) 
	fprintf (stderr, "bounds: found %d convex boundary points.\n", hullsize);
    }
  
  /* Concave Hull - distance weighted pacakage wrap algorithm 
   * Note: will return a single polygon containing all the given
   * points, using whatever distance value is needed to
   * accomplish that.
   */
  else if (vflag == 1) 
    {
      load_pnts (fp, &pnts, &npr, ptrec, verbose_flag);

      /* The distance parameter can't be less than zero */
      if (!dist)
	dist = qadd (pnts, npr);
      if (dist > 0) hullsize = -1;
      else hullsize = 0;

      /* if (verbose_flag > 0) */
      /* 	fprintf (stderr,"bounds: %-10s %-10s\n\rbounds: %-10d %-10.12f", "iteration", "distance", pc, dist); */

      /* Keep a copy of the original point-set in `pts2` in-case
	 * we need to re-run with a higher `dist` value. */
      point_t* pnts2;  
      pnts2 = (point_t*) malloc ((npr + 1) * sizeof (point_t));
      memcpy (pnts2, pnts, sizeof (point_t) * (npr + 1));
      
      /* Find a boundary, inrease the distance variable until a boundary is found. */
      while (hullsize == -1) 
	{
	  hullsize = dpw_concave (pnts, npr, dist);
	  
	  /* If a hull was found, check that it gathered all the points.
	   * We want to have a 'hull' in that the output contains all points
	   * within the boundary polygon. Otherwise instead of increasing
	   * the boundary we could run on outside points to output a 
	   * multipolygon.
	   * If there are points still outside the boundary, increase the
	   * distance and retry.
	   */
	  if (hullsize >= 0)
	    for (i = hullsize+1; i < npr; i++)
	      if (!inside_p (&pnts[i], pnts, hullsize)) 
		{
		  hullsize = -1;
		  break;
		}

	  /* If a hull wasn't found, increase the `dist` and try again. */
	  if (hullsize == -1) 
	    {
	      dist += dist, pc++;
	      memcpy (pnts, pnts2, sizeof (point_t) * (npr + 1));
	      /* if (verbose_flag > 0)  */
	      /* 	{ */
	      /* 	  fprintf (stderr,"\rbounds: %-10d %-10.12f", pc, dist); */
	      /* 	  fflush (stderr); */
	      /* 	} */
	    }
	  
	  /* In case something funky happens; just make a convex hull */
	  if (dist == INFINITY || dist == NAN || dist < 0)
	    { 
	      memcpy (pnts, pnts2, sizeof (point_t) * (npr + 1));
	      hullsize = pw_convex (pnts, npr);
	    }
	}
      
      /* Print out the hull */
      if (jsonflag > 0)
	{
	  for (i = 0; i < hullsize; i++)
	    printf ("[%f, %f], ", pnts[i].x, pnts[i].y);
	  printf ("[%f, %f]", pnts[hullsize].x, pnts[hullsize].y);
	}
      else
	{
	  for (i = 0; i <= hullsize; i++)
	    printf ("%f %f\n", pnts[i].x, pnts[i].y);
	}
      
      if (verbose_flag > 0) 
	  fprintf (stderr, "bounds: found %d total boundary points\n", hullsize);

      free (pnts2);
      pnts2 = NULL;
  }
  
  /* Bounding Box - Generate a box around the given points.
   */
  else if (bflag == 1) 
    {
      double ymin, ymax, xmin, xmax;
      /* Read through the point records and find the min/max bounding box.
       */
      i = 0;
      while (read_point(fp, &rpnt, &delim, ptrec, dflag, verbose_flag) == 0) 
	{
	  if (sl > 0) 
	    sl--;
	  else 
	    {
	      npr++;
	      if (i==0)
		{
		  ymin = rpnt.y;
		  ymax = rpnt.y;
		  xmin = rpnt.x;
		  xmax = rpnt.x;
		}
	      else
		{
		  if (rpnt.y < ymin) 
		    ymin = rpnt.y;
		  if (rpnt.x < xmin) 
		    xmin = rpnt.x;
		  if (rpnt.y > ymax) 
		    ymax = rpnt.y;
		  if (rpnt.x > xmax) 
		    xmax = rpnt.x;
		}
	      i++;
	      if (i>1)
		dflag++;
	    }
	}

      if (jsonflag > 0)
	{
	  printf ("[%f,%f],,", xmin, ymin);
	  printf ("[%f,%f],", xmin, ymax);
	  printf ("[%f,%f],", xmax, ymax);
	  printf ("[%f,%f],", xmax, ymin);
	  printf ("[%f,%f]", xmin, ymin);
	}
      else
	{
	  printf ("%f %f\n", xmin, ymin);
	  printf ("%f %f\n", xmin, ymax);
	  printf ("%f %f\n", xmax, ymax);
	  printf ("%f %f\n", xmax, ymin);
	  printf ("%f %f\n", xmin, ymin);
	}
      
      if (verbose_flag > 0) fprintf (stderr, "bounds: processed %d points.\n", npr);
    }
  
  /* Bounding Block - polygonize a grid of the points using `dist` cell-size
   */
  else if (kflag > 0) 
    {
      int kr_length;
      region_t rgn;
      
      kr_length = strlen (kreg);
      char* p = strtok (kreg, "/");
      for (j = 0; j < kr_length; j++) 
	{ 
	  if (p != NULL) 
	    {
	      if (j == 0) 
		dist = atof (p);
	      if (j == 1) 
		rgn.xmin = atof (p);
	      if (j == 2) 
		rgn.xmax = atof (p);
	      if (j == 3) 
		rgn.ymin = atof (p);
	      if (j == 4) 
		rgn.ymax = atof (p);
	    }
	  p = strtok (NULL, "/");
	}

      /* The distance parameter can't be less than zero */
      if (dist > 0) bbs_block (fp, dist, rgn, verbose_flag, jsonflag);
    }

  free (pnts);
  pnts = NULL;
  fclose (fp);

  if (jsonflag == 1)
    {
      printf ("]]]}}]}\n");
    }
  else if (jsonflag == 2)
    {
      printf("]]]}}");
    }
  
  if (verbose_flag > 0)
    {
      fprintf (stderr, "bounds: done\n");
    }
  
  exit (0);
}
