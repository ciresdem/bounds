## Headers
include_HEADERS=bounds.h

## Libraries
lib_LTLIBRARIES= libbounds.la
libbounds_la_SOURCES = hull.c pnts.c block.c bounds.h

## C Programs
bin_PROGRAMS = bounds
bounds_SOURCES = bounds.c
bounds_LDADD = libbounds.la
