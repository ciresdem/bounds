/*------------------------------------------------------------
 * bounds.c
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

#include <getopt.h>
#include "bounds.h"

/* Flags set by `--version' and `--help' */
static int version_flag;
static int help_flag;

static void
print_version(const char* command_name, const char* command_version) {
  fprintf(stderr, "%s version %s \n\
Copyright © 2011, 2012, 2016, 2018 Matthew Love <matthew.love@colorado.edu> \n\
%s is liscensed under the GPL v.2 or later and is \n\
distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;\n\
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A\n\
PARTICULAR PURPOSE.  See the GNU General Public License for more details.\n\
<http://www.gnu.org/licenses/>\n", command_name, command_version, command_name);
  exit (1);
}

static void 
usage() {
  fprintf(stderr, "\
bounds infile [Options]\n\n\
  Generate a boundary of the given xy point file(s) and print results \n\
  to stdout. Optionally, choose between different algorithms to use to \n\
  generate the boundary. The default behavior will generate a convex hull \n\
  boundary. \n\
\n\
  infile\t\tThe input point file. If not given, will read from stdin.\n\
\n\
 Options:\n\
  -d, --delimiter\tThe input xy file record delimiter\n\
     Note: The default will assume a space (\" \") is the delimiter.\n\
  -r, --record\t\tThe input record order, 'xy' should represent the locations\n\
              \t\tof the x and y records, respectively.\n\
     Note: The default will assume the record begins with 'xy'. If the x and \n\
           y records are located elsewhere, fill in records with any letter. i.e.\n\
           `-r zyx` for data where the record is `elevation,y-value,x-value`.\n\
  -s, --skip\t\tSkip header lines in the input xy file; specify the number of lines\n\
            \t\tto skip here.\n\n\
  -b, --box\t\tGenerate a 'bounding box' boundary. \n\
  -k, --block\t\tGenerate a 'block' boundary. \n\
             \t\tSpecify the blocking increment here. (e.g. --block 0.001)\n\
     Note: The blocking increment should be in the same units as the input xy data.\n\
  -x, --convex\t\tGenerate a 'convex hull' boundary using a monotone chain algorithm. (Default)\n\
  -v, --concave\t\tGenerate a 'concave hull' boundary using a distance weighted \n\
               \t\tpackage wrap algorithm. \n\
               \t\tSpecify the distance threshold used in generating the concave \n\
               \t\thull as an argument here. (e.g. --concave 0.01)\n\
     Note: The distance threshold should be in the same units as the input xy data.\n\
\n\
  --help\t\tPrint this help menu and exit.\n\
  --version\t\tPrint version information and exit.\n\n\
");
  exit(1);
}

/* Compare functions for use in qsort
 */
int
compare(const void* a, const void* b) {
  const point_t *elem1 = a;    
  const point_t *elem2 = b;

  if (elem1->x < elem2->x) return -1;
  else if (elem1->x > elem2->x) return 1;
  else return 0;
}

int
sort_min_y(const void* a, const void* b) {
  const point_t *elem1 = a;    
  const point_t *elem2 = b;

  if (elem1->y < elem2->y) return -1;
  else if (elem1->y > elem2->y) return 1;
  else return 0;
}

int
sort_max_y(const void* a, const void* b) {
  const point_t *elem1 = a;    
  const point_t *elem2 = b;

  if (elem1->y > elem2->y) return -1;
  else if (elem1->y < elem2->y) return 1;
  else return 0;
}

int
sort_max_x(const void* a, const void* b) {
  const point_t *elem1 = a;    
  const point_t *elem2 = b;

  if (elem1->x > elem2->x) return -1;
  else if (elem1->x < elem2->x) return 1;
  else return 0;
}

int
sort_min_x(const void* a, const void* b) {
  const point_t *elem1 = a;    
  const point_t *elem2 = b;

  if (elem1->x < elem2->x) return -1;
  else if (elem1->x > elem2->x) return 1;
  else return 0;
}

int
main (int argc, char **argv) {
  char* fn;
  FILE* fp;

  int c, i, status, min, j;
  int inflag = 0, vflag = 0, sflag = 0, pc = 0, sl = 0;
  int cflag = 0, kflag = 0, calg = 0, bflag = 0;
  double dist = 1;

  point_t rpnt, pnt;
  point_ptr_t hull0[MAX_HULLS];
  point_ptr_t* hull = hull0;
  ssize_t hullsize;

  char* delim = " ";
  char* ptrec = "xy";
  
  while (1) {
    static struct option long_options[] =
      {
	/* These options set a flag. */
	{"version", no_argument, &version_flag, 1},
	{"help", no_argument, &help_flag, 1},
	/* These options don't set a flag.
	   We distinguish them by their indices. */
	{"delimiter", required_argument, 0, 'd'},
	{"skip", required_argument, 0, 's'},
	{"record", required_argument, 0, 'r'},
	{"box", no_argument, 0, 'b'},
	{"block", required_argument, 0, 'k'},
	{"convex", no_argument, 0, 'x'},
	{"concave", required_argument, 0, 'v'},
	{0, 0, 0, 0}
      };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    c = getopt_long (argc, argv, "d:r:s:bk:xv:",
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
      delim = optarg;
      break;
    case 'r':
      ptrec = optarg;
      break;
    case 's':
      sflag++;
      sl = atoi(optarg);
      break;
    case 'b':
      bflag++;
      break;
    case 'k':
      kflag++;
      dist = atof(optarg);
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
      break;
      
    default:
      abort ();
    }
  }

  if (version_flag) print_version("bounds", BOUNDS_VERSION);
  if (help_flag) usage();
  //if (inflag == 0 && optind >= argc) usage ();

  fn = argv[optind];
  if (fn) inflag++;
  //if (inflag == 0) {
  //fp = fopen(fn, "r");
  if (inflag > 0) {
    fp = fopen(fn, "r");
  } else {
    fp = stdin;
    fn = "stdin";
  }
    
  if (!fp) {
    fprintf(stderr,"bounds: Failed to open file: %s\n", fn);
    exit(0);
  } else fprintf(stderr,"bounds: Working on file: %s ...\n", fn);

  ssize_t npr = 0;
  
  /* Allocate memory for the `pnts` and `pnts2` array using the total number of points.
   * `pnts2` is only used by concave.
   */
  point_t* pnts;  
  pnts = (point_t*) malloc(sizeof(point_t));
    
  /* Read through the point records and record them in `pnts`.
   */
  i = 0;
  while (read_point(fp, &rpnt, delim, ptrec) == 0) {
    if (sl > 0) sl--;
    else {
      npr++;
      point_t* temp = realloc(pnts, npr * sizeof(point_t));
      pnts = temp;
      pnts[i].x = rpnt.x, pnts[i].y = rpnt.y;
      i++;
    }
  }
  
  /* This is for the GMT compatibility. More can be done here.
   */
  printf("# @VGMT1.0 @GMULTIPOLYGON\n# @NName\n# @Tstring\n# FEATURE_DATA\n");
  printf(">\n# @D%s\n# @P\n", fn);
    
  /* The default is a convex hull -- `cflag` */
  if (cflag == 0 && vflag == 0 && bflag == 0 && kflag == 0) cflag++;
  
  /* 
   * Monotone Chain Convex Hull Algorithm - -*Default*-
   */
  if (cflag == 1) {
    qsort(pnts, npr, sizeof(point_t), compare);
    mc_convex(pnts, npr, &hull, &hullsize);
    
    for (i = 0; i < hullsize; i++)
      printf("%f %f\n", hull[i]->x, hull[i]->y);
    
    fprintf(stderr, "bounds: Found %d convex boundary points.\n", hullsize);
  }
  
  /* 
   * 'Package Wrap' Convex Hull Algorithm 
   */
  else if (cflag == 1 && calg == 2){
    hullsize = pw_convex(pnts, npr);
    
    for (i = 0; i < hullsize; i++)
      printf("%f %f\n", pnts[i].x, pnts[i].y);
    
    fprintf(stderr, "bounds: Found %d convex boundary points.\n", hullsize);
  }
  
  /* 
   * Concave Hull - distance weighted pacakage wrap algorithm 
   * Note: will return a single polygon containing all the given
   * points, using whatever distance value is needed to
   * accomplish that.
   */
  else if (vflag == 1) {
    /* The distance parameter can't be less than zero */
    if (dist > 0) hullsize = -1;
    else hullsize = 0;
    
    /* Keep a copy of the original point-set in `pts2` in-case
     * we need to re-run with a higher `dist` value. */

    point_t* pnts2;  
    pnts2 = (point_t*) malloc(npr*sizeof(point_t));
    memcpy(pnts2, pnts, sizeof(point_t)*npr);
    
    /* Find a boundary, inrease the distance variable until a boundary is found. */
    while (hullsize == -1) {
      qsort(pnts, npr, sizeof(point_t), sort_max_x);
      hullsize = dpw_concave(pnts, npr, dist);
      
      /* If a hull was found, check that it gathered all the points.
       * If there are points still outside the boundary, increase the
       * distance and retry.
       */
      if (hullsize >= 0)
	for (i = hullsize; i < npr; i++)
	  if (!inside(&pnts[i], pnts, hullsize, dist)) {
	    hullsize = -1;
	    break;
	  }
      
      /* If a hull wasn't found, increase the `dist` and try again. */
      if (hullsize == -1) {
	dist = dist + (dist * 0.1), pc++;
	fprintf(stderr,"\rbounds: Increasing distance ( %f )",dist);
	  fflush(stderr);
	  memcpy(pnts, pnts2, sizeof(point_t)*npr);
      }
      
      /* In case something funky happens. */
      if (dist == INFINITY || dist == NAN || dist <= FLT_EPSILON) hullsize = 0;
    }

    /* Print out the hull */
    for (i = 0; i <= hullsize; i++)
      printf("%.7f %.7f\n", pnts[i].x, pnts[i].y);
    
    if (pc > 0) fprintf(stderr,"\n");
    fprintf(stderr, "bounds: Found %d concave boundary points at a %.6f distance threshhold.\n", hullsize, dist);
  }
    
  /* 
   * Bounding Box - Generate a box around the given points.
   */
  else if (bflag == 1) {
    int xmin, xmax, ymin, ymax;
    
    /* Find the extent values in the dataset */
    for (ymin = 0, ymax = 0, xmin = 0, xmax = 0, i = 1; i < npr; i++) {
      if (pnts[i].y < pnts[ymin].y) ymin = i;
      if (pnts[i].x < pnts[xmin].x) xmin = i;
      if (pnts[i].y > pnts[ymax].y) ymax = i;
      if (pnts[i].x > pnts[xmax].x) xmax = i;
    }
    printf("%f %f\n", pnts[xmin].x, pnts[ymin].y);
    printf("%f %f\n", pnts[xmin].x, pnts[ymax].y);
    printf("%f %f\n", pnts[xmax].x, pnts[ymax].y);
    printf("%f %f\n", pnts[xmax].x, pnts[ymin].y);
    printf("%f %f\n", pnts[xmin].x, pnts[ymin].y);
  }
  
  /* 
   * Bounding Block - polygonize a grid of the points using `dist` cell-size
   */
  else if (kflag == 1) {
    /* The distance parameter can't be less than zero */
    if (dist > 0) block_pts(pnts, npr, dist);
  }
  free(pnts), fclose(fp);
  exit(1);
}
