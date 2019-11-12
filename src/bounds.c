/*------------------------------------------------------------
 * bounds.c
 * 
 * This file is part of BOUNDS
 *
 * Copyright (c) 2011, 2012, 2016, 2018, 2019 Matthew Love <matthew.love@colorado.edu>
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
print_version(const char* command_name, const char* command_version) {
  fprintf(stderr, "%s %s \n\
Copyright © 2011, 2012, 2016, 2018, 2019 Matthew Love <matthew.love@colorado.edu> \n\
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
bounds [OPTION]... [FILE]\n\
Generate a boundary of the set of xy points from FILE, or standard input, to standard output.\n\
\n\
  -d, --delimiter\tthe input xy file record delimiter\n\
  -r, --record\t\tthe input record order, 'xy' should represent the locations\n\
              \t\tof the x and y records, respectively. (e.g. --record zdyx)\n\
  -s, --skip\t\tspecify the number of lines to skip here.\n\
  -n, --name\t\tthe output layer name\n\
  -a, --append\t\tdon't print gmt-header\n\n\
  ---- bounds ----\n\n\
  -b, --box\t\t'bounding box' boundary. \n\
  -k, --block\t\t'block' boundary. specify the blocking increment here. (e.g. --block 0.001)\n\
             \t\tappend a region if desired (e.g. --block 0.001/west/east/south/north\n\
  -x, --convex\t\t'convex hull' boundary using a monotone chain algorithm. (default)\n\
  -v, --concave\t\t'concave hull' boundary using a distance weighted package wrap algorithm. \n\
               \t\tspecify the distance threshold here. (e.g. --concave 0.01)\n\
      --verbose\t\tincrease the verbosity.\n\
      --help\t\tprint this help menu and exit.\n\
      --version\t\tprint version information and exit.\n\n\
\n\
With no FILE, or when FILE is --, read standard input.\n\
All OPTION values must be in the same units as the input xy data.\n\n\
Examples:\n\
  bounds \t\toutput a convex hull from standard input.\n\
  bounds -k 0.0001\toutput a 'block' boundary from standard input.\n\
  bounds -v 10 -d \",\"\toutput a concave hull from comma-delimited standard input.\n\
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
  int cflag = 0, kflag = 0, calg = 0, bflag = 0, gmtflag = 0, nflag = 0;
  double dist = 1;

  point_t rpnt, pnt;
  point_ptr_t hull0[MAX_HULLS];
  point_ptr_t* hull = hull0;
  ssize_t hullsize;
  ssize_t npr = 0;

  char* delim = " \t";
  char* ptrec = "xy";
  char* kreg = "";
  char* lname = "bounds";
  
  while (1) {
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
	{"append", no_argument, 0, 'a'},
	{"record", required_argument, 0, 'r'},
	{"box", no_argument, 0, 'b'},
	{"block", required_argument, 0, 'k'},
	{"convex", no_argument, 0, 'x'},
	{"concave", required_argument, 0, 'v'},
	{0, 0, 0, 0}
      };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    c = getopt_long (argc, argv, "ad:n:r:s:bk:xv:",
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
    case 'n':
      nflag++;
      lname = optarg;
      break;
    case 'a':
      gmtflag++;
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
      fprintf(stderr,"Try 'bounds --help' for more information.\n");
      exit(0);
      break;
      
    default:
      abort ();
    }
  }

  if (version_flag) print_version("bounds", BOUNDS_VERSION);
  if (help_flag) usage();

  fn = argv[optind];
  if (fn) inflag++;

  if (inflag > 0) {
    fp = fopen(fn, "r");
    if (!fp) {
      if (verbose_flag > 0) fprintf(stderr,"bounds: Failed to open file: %s\n", fn);
      exit(EXIT_FAILURE);
    }
  } else {
    fp = stdin;
    fn = "stdin";
  }

  if (verbose_flag > 0) fprintf(stderr, "bounds: Working on file: %s\n", fn);
  
  /* Allocate memory for the `pnts` and `pnts2` array using the total number of points.
   * `pnts2` is only used by concave and gets allocated there.
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

  if (verbose_flag > 0) fprintf(stderr,"bounds: pts: %d\n",i);
 
  /* This is for the GMT compatibility. More can be done here.
   */
  if (gmtflag == 0) {
    printf("# @VGMT1.0 @GMULTIPOLYGON\n# @NName\n# @Tstring\n# FEATURE_DATA\n");
  }
  printf(">\n# @D%s\n# @P\n", lname);
    
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
    
    if (verbose_flag > 0) fprintf(stderr, "bounds: Found %d convex boundary points.\n", hullsize);
  }
  
  /* 
   * 'Package Wrap' Convex Hull Algorithm 
   */
  else if (cflag == 1 && calg == 2){
    hullsize = pw_convex(pnts, npr);
    
    for (i = 0; i < hullsize; i++)
      printf("%f %f\n", pnts[i].x, pnts[i].y);
    
    if (verbose_flag > 0) fprintf(stderr, "bounds: Found %d convex boundary points.\n", hullsize);
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
       * We want to have a 'hull' in that the output contains all points
       * within the boundary polygon. Otherwise instead of increasing
       * the boundary we could run on outside points to output a 
       * multipolygon.
       * If there are points still outside the boundary, increase the
       * distance and retry.
       */
      if (hullsize >= 0)
      	for (i = hullsize; i < npr; i++)
      	  if (!inside(&pnts[i], pnts, hullsize, dist)) {
      	    if (verbose_flag > 0) fprintf(stderr,"bounds: Failed to gather all points");
      	    hullsize = -1;
      	    break;
      	  }

      /* If a hull wasn't found, increase the `dist` and try again. */
      if (hullsize == -1) {
	dist = dist + (dist * 0.25), pc++;
	memcpy(pnts, pnts2, sizeof(point_t)*npr);
	if (verbose_flag > 0) {
	  fprintf(stderr,"\rbounds: Increasing distance ( %f )",dist);
	  fflush(stderr);
	}
      }
      
      /* In case something funky happens. */
      if (dist == INFINITY || dist == NAN || dist <= FLT_EPSILON) hullsize = 0;
    }

    /* Print out the hull */
    for (i = 0; i <= hullsize; i++)
      printf("%f %f\n", pnts[i].x, pnts[i].y);
    
    if (verbose_flag > 0) {
      if (pc > 0) fprintf(stderr,"\n");
      fprintf(stderr, "bounds: Found %d concave boundary points at a %.6f distance threshhold.\n", hullsize, dist);
    }
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
    int kr_length;
    xyz_info rgn;
    kr_length = strlen(kreg);

    char* p = strtok(kreg, "/");
    for (j = 0; j < kr_length; j++) { 
      if (p != NULL) {
	if (j == 0) dist = atof(p);
	if (j == 1) rgn.xmin = atof(p);
	if (j == 2) rgn.xmax = atof(p);
	if (j == 3) rgn.ymin = atof(p);
	if (j == 4) rgn.ymax = atof(p);
      }
      p = strtok(NULL, "/");
    }
    
    /* The distance parameter can't be less than zero */
    if (dist > 0) block_pts(pnts, npr, dist, rgn, verbose_flag);
  }

  free(pnts);
  pnts=NULL;
  fclose(fp);

  if (verbose_flag > 0) fprintf(stderr, "bounds: done\n");

  exit(1);
}
