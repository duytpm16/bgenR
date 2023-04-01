define([pkgversion], esyscmd([sh -c "grep Version: DESCRIPTION | cut -d' ' -f2 | tr -d '\n'"]))dnl
AC_INIT([bgenR],[pkgversion],[duytpham18@gmail.com])
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   echo "could not determine R_HOME"
   exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

dnl adopted from: https://github.com/cran/dfmta/blob/master/configure.ac
AC_DEFUN([AX_GET_R_CONF_VAR],[$1=`"${R_HOME}/bin/R" CMD config $1`])
AX_GET_R_CONF_VAR([CC])
AX_GET_R_CONF_VAR([CFLAGS])
AX_GET_R_CONF_VAR([CXX])
AX_GET_R_CONF_VAR([CXXFLAGS])
CXX="${CXX}"
CXXFLAG="${CXXFLAGS}"

dnl print out flags to ease debugging
AC_MSG_NOTICE([CC = $CC])
AC_MSG_NOTICE([CFLAGS = $CFLAGS])
AC_MSG_NOTICE([CXX = $CXX])
AC_MSG_NOTICE([CXXFLAGS = $CXXFLAGS])

dnl skip checking C compilers as we also need to check C++ compiler
dnl AC_LANG(C)
dnl AC_PROG_CC
AC_LANG(C++)
AC_PROG_CXX

dnl check libraries
dnl refer:
dnl   https://www.gnu.org/savannah-checkouts/gnu/autoconf/manual/autoconf-2.69/html_node/Libraries.html
dnl   AC_CHECK_LIB (library, function, [action-if-found], [action-if-not-found], [other-libraries])
LIB_ZSTD=""
AC_CHECK_LIB([zstd], [ZSTD_decompress], [
  AC_CHECK_HEADERS([zstd.h], [
    HAVE_ZSTD=1
    AC_DEFINE([HAVE_ZSTD], [1], [Define to 1 if you have zstd library])
    LIB_ZSTD="-lzstd"
    AC_MSG_NOTICE([Will use system-wide zstd library])    
  ])                     
], [AC_MSG_NOTICE([Will build zstd library from the package])])
             

dnl dnl compile zstd-1.5.4
dnl AC_MSG_NOTICE([Finished building zstd library])
dnl base_dir=$(pwd)
dnl ZSTD_CFLAGS=$("${R_HOME}/bin/R" CMD config CFLAGS)
dnl ZSTD_LIBS="$base_dir/src/zstd/lib/libzstd.a"
dnl AC_CONFIG_COMMANDS([zstd compiling],
dnl                    [make -C src/zstd/lib lib CFLAGS="-fPIC ${ZSTD_CFLAGS}"])


LIB_LIBDEFLATE=""
AC_CHECK_LIB([deflate], [libdeflate_zlib_decompress], [
  AC_CHECK_HEADERS([libdeflate.h], [
    HAVE_LIBDEFLATE=1
    AC_DEFINE([HAVE_LIBDEFLATE], [1], [Define to 1 if you have libdeflate library])
    LIB_LIBDEFLATE="-ldeflate"
    AC_MSG_NOTICE([Will use system-wide libdeflate library])    
  ])                     
], [AC_MSG_NOTICE([Will build libdeflate library from the package])])
             

dnl dnl compile libdeflate-1.18
dnl AC_MSG_NOTICE([Finished building libdeflate library])
dnl base_dir=$(pwd)
dnl LIBDEFLATE_CFLAGS=$("${R_HOME}/bin/R" CMD config CFLAGS)
dnl LIBDEFLATE_LIBS="$base_dir/src/libdeflate/lib/libdeflate.a"
dnl AC_CONFIG_COMMANDS([libdeflate compiling],
dnl                    [make -C src/libdeflate/lib lib CFLAGS="-fPIC ${LIBDEFLATE_CFLAGS}"])


AC_SUBST([PKG_CFLAGS], ["-I. ${CFLAGS}"])
AC_SUBST([PKG_CPPFLAGS], ["-I. ${CXXFLAGS}"])
AC_SUBST([PKG_LIBS], ["-lz -lm ${LIB_ZSTD} ${LIB_LIBDEFLATE}"])
AC_SUBST([HAVE_ZSTD])
AC_SUBST([HAVE_LIBDEFLATE])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
