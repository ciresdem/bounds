
#     [ B O U N D S ]

---------------------------

BOUNDS is a command-line program for the GNU operating system that generates a boundary of a set of xy points.

With Bounds, generate an ogr-gmt compatible boundary vector of a given point set. There are a number of algorithms to choose from in generating the boundary, including convex hull, concave hull, bounding block and bounding box.

---------------------------

```
bounds [OPTION]... [FILE]
Generate a boundary of the set of xy points from FILE, or standard input, to standard output.

  ---- xy i/o ----

  -d, --delimiter       The input xy file record delimiter.
                        If omitted, the delimiter will be guessed from the first line read.
  -g, --gmt             Format output as GMT vector multipolygon; Use twice to supress
                        the initial header (e.g. -gg).
  -n, --name            The output layer name (only used with -g).
  -r, --record          The input record order, 'xy' should represent the locations
                        of the x and y records, respectively (e.g. --record zdyx).
  -s, --skip            The number of lines to skip from the input.

  ---- bounds ----

  -b, --box             'Bounding Box' boundary. 
  -k, --block           'Bounding Block' boundary. Specify the blocking increment
                        in input units (e.g. --block 0.001). Specify a blocking region
                        after the increment if desired (e.g. --block 0.001/west/east/south/north).
  -x, --convex          'Convex Hull' boundary using a monotone chain algorithm. [default]
                        Use twice to use a package wrap algorithm (e.g. -xx).
  -v, --concave         'Concave Hull' boundary using a distance weighted package wrap algorithm.
                        Specify distance value or - to estimate appropriate distance.

  ---- et cetra ----

      --verbose         increase the verbosity.
      --help            print this help menu and exit.
      --version         print version information and exit.


With no FILE, or when FILE is --, read standard input.
All OPTION values must be in the same units as the input xy data.

Examples:
  bounds                output a convex hull from standard input.
  bounds -g -k0.0001    output a GMT formatted 'block' boundary from standard input.
  bounds -v10 -d,       output a concave hull from comma-delimited standard input.
  bounds -v- in.xyz     output a concave hull from file in.xyz
```

![](./media/bounds_box.jpg)

**Figure 1.** NOS survey H08891 bounded with 'box'

![](./media/bounds_block.jpg)

**Figure 2.** NOS survey H08891 bounded with 'block'

![](./media/bounds_convex.jpg)

**Figure 3.** NOS survey H08891 bounded with 'convex'

![](./media/bounds_concave.jpg)

**Figure 4.** NOS survey H08891 bounded with 'concave'