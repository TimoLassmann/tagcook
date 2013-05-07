#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(Tagdust, 2.0, timolassmann@gmail.com)
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE([1.10 -Wall no-define])

CXXFLAGS="$CXXFLAGS"
CFLAGS="$CFLAGS"

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([ assert.h ctype.h float.h getopt.h limits.h math.h pthread.h stdio.h stdlib.h string.h sys/stat.h time.h unistd.h xmmintrin.h ])



# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

AC_CONFIG_FILES([Makefile])
AC_OUTPUT