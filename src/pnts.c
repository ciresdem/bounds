/*------------------------------------------------------------
 * pnts.c
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

/* Line-Count 
 */
ssize_t
linecnt(FILE *in) {
  char tmp[MAX_RECORD_LENGTH] = {0x0};
  int recordcnt = 0;

  while(fgets(tmp, sizeof(tmp), in) != 0) recordcnt++;
  rewind(in);

  return recordcnt; 
}

/* Scan an xyz line and attempt to determine the delimiter and record length 
 * TODO: guess the delimiter
 */
int
scanline(char* infile) {
  FILE *in = fopen(infile, "r");  
  int c, count;

  count = 0;
  for (;;) {
    c = fgetc(in);
    if (c == EOF || c == '\n') break;
    ++count;
  }
  return count;
}

/* Read a point record from the given infile
 */
int
read_point(FILE *infile, point_t *rpnt, char* delimiter, char* pnt_recr) {
  char tmp[MAX_RECORD_LENGTH] = {0x0};
  char pntp, *strs;

  int pf_length, j, status;

  /* Read in the file */
  if(infile == NULL) {
    perror("bounds: File open error");
    exit(EXIT_FAILURE);
  }

  /* read a record */
  if (fgets(tmp, sizeof(tmp), infile) !=0) {
    strs = malloc(strlen(tmp) + 1);
    strcpy(strs, tmp);
    status = 0;
  }
  else {
    status = -1;
    return status;
  }

  pf_length = strlen(pnt_recr);

  char* p = strtok(strs, delimiter);
  for (j = 0; j < pf_length; j++) {
    pntp = pnt_recr[j];
    if (p != NULL) {
      if (pntp == 'x') rpnt->x = atof(p);
      else if (pntp == 'y') rpnt->y = atof(p);
    }
    p = strtok(NULL, delimiter);
  }
  free(strs);  
  return status;
}

/*
 * Get the minimum and maximum values from a set of points
 */
void
minmax(point_t* points, int npoints, xyz_info *xyzi) {
  int i;

  for (i = 0; i < npoints; i++) {
    if (i == 0) {
      xyzi->xmin = points[i].x, xyzi->ymin = points[i].y;
      xyzi->xmax = points[i].x, xyzi->ymax = points[i].y;
    }
    else {
      if (points[i].x <= xyzi->xmin) xyzi->xmin = points[i].x;
      if (points[i].x >= xyzi->xmax) xyzi->xmax = points[i].x;
      if (points[i].y <= xyzi->ymin) xyzi->ymin = points[i].y;
      if (points[i].y >= xyzi->ymax) xyzi->ymax = points[i].y;
    }
  }
}
