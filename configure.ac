AC_INIT([main.cxx],[0.1])
AC_CONFIG_SRCDIR([./])
AM_INIT_AUTOMAKE
AC_PROG_CC
AC_CHECK_HEADERS(zlib.h bam.h)
AC_PROG_CXX
AC_CHECK_SIZEOF(void*, 8)
AC_CHECK_SIZEOF(int, 4)
AC_CHECK_LIB(pthread, pthread_create)
AC_CHECK_LIB(bam, bam_read)
if test -f /usr/local/include/bam.h
then
	AC_DEFINE(HAVE_BAM_H)
fi
AC_OUTPUT(Makefile)
