## [ B O U N D S ]

BOUNDS is a command-line program for the GNU operating system that generates a boundary of a set of xy points.

With Bounds, generate an ogr-gmt compatible boundary vector of a given point set. There are a number of algorithms to choose from in generating the boundary, including convex hull, concave hull, bounding block and bounding box.

```
bounds infile [Options]

  Generate a boundary of the given xy point file(s) and print results 
  to stdout. Optionally, choose between different algorithms to use to 
  generate the boundary. The default behavior will generate a convex hull 
  boundary. 

  infile		The input point file.

 Options:
  -d, --delimiter	The input xy file record delimiter
     Note: The default will assume a space (" ") is the delimiter.
  -r, --record		The input record order, 'xy' should represent the locations
              		of the x and y records, respectively.
     Note: The default will assume the record begins with 'xy'. If the x and 
           y records are located elsewhere, fill in records with any letter. i.e.
           `-r zyx` for data where the record is `elevation,y-value,x-value`.
  -s, --skip		Skip header lines in the input xy file; specify the number of lines
            		to skip here.

  -b, --box		Generate a 'bounding box' boundary. 
  -k, --block		Generate a 'block' boundary. 
             		Specify the blocking increment here. (e.g. --block 0.001)
     Note: The blocking increment should be in the same units as the input xy data.
  -x, --convex		Generate a 'convex hull' boundary using a monotone chain algorithm. (Default)
  -v, --concave		Generate a 'concave hull' boundary using a distance weighted 
               		package wrap algorithm. 
               		Specify the distance threshold used in generating the concave 
               		hull as an argument here. (e.g. --concave 0.01)
     Note: The distance threshold should be in the same units as the input xy data.

  --help		Print this help menu and exit.
  --version		Print version information and exit.
```