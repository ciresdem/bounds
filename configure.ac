#                                               -*- Autoconf -*-

dnl Copyright (C) 2018 by Matthew Love

dnl bounds : version 0.4.2 : malo
dnl
dnl This file is part of BOUNDS
dnl
dnl Copyright (c) 2011, 2012, 2016, 2018 Matthew Love <matthew.love@colorado.edu>
dnl BOUNDS is liscensed under the GPL v.2 or later and 
dnl is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details. 
dnl <http://www.gnu.org/licenses/> 

# Process this file with autoconf to produce a configure script.

AC_INIT([bounds], [0.4.2], [matthew.love@colorado.edu])
AC_CONFIG_AUX_DIR([./config])
AM_INIT_AUTOMAKE([-Wall -Werror])
AC_CONFIG_MACRO_DIR([m4])
AM_PROG_AR

AC_CHECK_FUNCS(dlopen)
if test $ac_cv_func_dlopen = yes
then DYNAMIC_AG=-export-dynamic
else DYNAMIC_AG=""
fi
AC_SUBST(DYNAMIC_AG)

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for libraries.

AC_CHECK_LIB([m], fabs, [], [
		  echo "Math library is required..."
		  exit -1])

# Checks for header files.
AC_CHECK_HEADERS([stdio.h stdlib.h string.h math.h float.h limits.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile src/Makefile doc/Makefile])
AC_OUTPUT
