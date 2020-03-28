/*------------------------------------------------------------
 * pnts.c
 *
 * This file is part of BOUNDS
 *
 * Copyright (c) 2011, 2012, 2013, 2014, 2015, 2018, 2019 Matthew Love <matthew.love@colorado.edu>
 * BOUNDS is liscensed under the GPL v.2 or later and 
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
linecnt (FILE *in) 
{
  char tmp[MAX_RECORD_LENGTH] = {0x0};
  int recordcnt = 0;
  
  while (fgets (tmp, sizeof (tmp), in) != 0) 
    recordcnt++;
  rewind (in);
  return recordcnt; 
}

/* Scan an xyz line and attempt to determine the record length 
 */
int
scanline (char* infile) 
{
  FILE *in = fopen (infile, "r");  
  int c, count;

  count = 0;
  for (;;) 
    {
      c = fgetc (in);
      if (c == EOF || c == '\n') 
	break;
      ++count;
    }
  return count;
}

/* Guess the delimiter from line.
 * The result will be recorded in `delimiter`
 */
int
auto_delim_l (char* inl, char** delimiter)
{
  char *strs;
  int cnt = 0, j;
  
  /* Known delimiters
   */
  int nkds = 3;
  char* kds[3] = {" \t", ",", "|"};
  
  strs = malloc (strlen (inl) + 1);
  strcpy (strs, inl);
  
  for (j = 0; j < nkds; j++)
    {
      cnt = 0;
      char* p = strtok (strs, kds[j]);
      while (p)
	{
	  cnt++;
	  p = strtok (NULL, kds[j]);
	}
      if (cnt > 1)
	{
	  *delimiter = kds[j];
	  j = nkds;
	}
    }
  
  free (strs);  
  return 1;
}

/* Read a point record from the given infile
 * If dflag < 1 then guess the delimiter and store it's value in `delimiter`.
 */
int
read_point (FILE *infile, point_t *rpnt, char** delimiter, char* pnt_recr, int dflag) 
{
  char tmp[MAX_RECORD_LENGTH] = {0x0};
  char pntp, *strs;
  int pf_length, j, status;

  /* Read in the file */
  if (infile == NULL) 
    {
      perror("bounds: File open error");
      exit(EXIT_FAILURE);
    }

  /* read a record */
  if (fgets (tmp, sizeof (tmp), infile) !=0) 
    {
      strs = malloc(strlen (tmp) + 1);
      strcpy (strs, tmp);
      status = 0;
    }
  else 
    {
      status = -1;
      return status;
    }

  pf_length = strlen (pnt_recr);

  if (!dflag)
    {
      dflag = auto_delim_l (strs, delimiter);
      fprintf(stderr,"bounds: delimiter is '%s'\n", *delimiter);
    }
  
  char* p = strtok (strs, *delimiter);
  for (j = 0; j < pf_length; j++) 
    {
      pntp = pnt_recr[j];
      if (p != NULL) 
	{
	  if (pntp == 'x') 
	    rpnt->x = atof (p);
	  else if (pntp == 'y') 
	    rpnt->y = atof (p);
	}
      p = strtok (NULL, *delimiter);
    }
  free (strs);  
  return status;
}

/* Load points
 */
int
load_pnts(FILE *infile, point_t **pnts, ssize_t *npr, char* pnt_recr)
{
  point_t rpnt;
  int i = 0, dflag = 0, sl = 0;
  char* delim;
    
  /* Read through the point records and record them in `pnts`.
   */
  *npr = 0;
  while (read_point(infile, &rpnt, &delim, pnt_recr, dflag) == 0)
    {
      if (sl > 0)
  	sl--;
      else
  	{
  	  *npr = *npr + 1;
  	  *pnts = realloc (*pnts, (*npr+1) * sizeof (point_t));
  	  (*pnts)[i].x = rpnt.x, (*pnts)[i].y = rpnt.y;
  	  i++;
  	  if (i>0)
  	    dflag++;
  	}
    }
  fprintf (stderr,"bounds: Processing %d points\n", *npr);
}
