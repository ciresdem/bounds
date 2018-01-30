#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([bounds], [0.4.1], [matthew.love@colorado.edu])
AC_CONFIG_AUX_DIR([./config])
AM_INIT_AUTOMAKE([-Wall -Werror])
#AC_CONFIG_SRCDIR([.])
AM_PROG_AR

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL

dnl AC_ARG_WITH(
dnl 	[gdal],
dnl 	AC_HELP_STRING([--with-gdal=path-to-gdal], [Specify the location of gdal]), 
dnl 	[ncloc=$withval
dnl 	CFLAGS="$CFLAGS -I$ncloc/include"
dnl 	LDFLAGS="$LDFLAGS -L$ncloc/lib"])


# Checks for libraries.

AC_CHECK_LIB([m], fabs, [], [
		  echo "Math library is required..."
		  exit -1])

dnl AC_CHECK_LIB([gdal], GDALAllRegister, [], [
dnl 		       echo "GDAL is required..."
dnl 		       exit -1])
# Checks for header files.
AC_CHECK_HEADERS([stdlib.h math.h float.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile src/Makefile doc/Makefile])
AC_OUTPUT
