#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([liblsb], [0.2.2], [salvodg@gmail.com])
AC_CONFIG_SRCDIR([liblsb_fortran.cpp])
AC_CONFIG_FILES([liblsb.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR(config)

LT_INIT[disable-shared]
AC_CONFIG_MACRO_DIR([m4])

AM_MAINTAINER_MODE

AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_MAKE_SET
AC_PROG_CXX
AC_PROG_CC

NG_WITH_MPICC
NG_WITH_MPICXX
LSB_WITH_PAPI
LSB_ENABLE_SYNC


# Checks for libraries.
AC_CHECK_LIB([unwind], [unw_init_local])

# Checks for header files.
AC_CHECK_HEADERS([inttypes.h stdint.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_UINT64_T

# Checks for library functions.
AC_CHECK_FUNCS([atexit uname])



if test "x$have_mpi" = "x1"; then
    AC_SUBST(HAVE_MPI_DEF, ["#define HAVE_MPI_H"])
else
    AC_SUBST(HAVE_MPI_DEF, [""])
fi


AC_MSG_NOTICE([**** checking for available high precision timers ****])
#AC_MSG_CHECKING([high resolution timer architecture])

if test "x$HRT_ARCH" = "x"; then
# not enforced on command line

  hr_subdir="hrtimer"
  HTOR_TIMERARCH(x86_64,2, $hr_subdir)
  if test "x$hrt_arch" = "x"; then # if not found yet 
    HTOR_TIMERARCH(x86_32,1,$hr_subdir)
  fi
  if test "x$hrt_arch" = "x"; then # if not found yet
    HTOR_TIMERARCH(ppc,3, $hr_subdir)
  fi
  if test "x$hrt_arch" = "x"; then # if not found yet
    HTOR_TIMERARCH(ia64,4, $hr_subdir)
  fi
  if test "x$hrt_arch" = "x"; then # if not found yet
    HTOR_TIMERARCH(SiCortex,5, $hr_subdir)
  fi
  if test "x$hrt_arch" = "x"; then # if not found yet
    HTOR_TIMERARCH(wtime,6, $hr_subdir)
  fi

  if test "x$hrt_arch" != "x"; then
    ######################################################
    # getting ticks guess ... should be portable
    printf "compiling tick-counter benchmark ... "
    $CC $CFLAGS -DHRT_ARCH=$hrt_arch_num $hr_subdir/getres.c -o getres
    if test "$?" != "0"; then 
      exit 1; 
    fi;
    printf "done\n" 
    if test "$hrt_arch" != "wtime"; then
      printf "getting clock resolution (may take some seconds)"
      this_ticks=$(./getres)
      if test "$?" != "0"; then 
        exit 1; 
      fi;
      printf " $this_ticks ticks/sec\n"
      rm -f getres
    fi
  else
    printf "*** warning - no high resolution timer found! ***\n"
  fi
else
  printf "*** warning - HRT_RESOLUTION not set because timer was forced to $hrt_arch!\n"
  AC_DEFINE_UNQUOTED([HRT_ARCH], $HRT_ARCH, [highrestimer architecture])
  AC_DEFINE_UNQUOTED([HAVE_HRTIMER], 1, [highrestimer enabled])

fi

# Check for 'diff' and get full path.
AC_ARG_VAR([DIFF],[the 'diff' program to use for test output comparison])
AC_PATH_PROG([DIFF],[diff],[])
if test "x$DIFF" = "x"; then
	AC_MSG_NOTICE([==> diff command not found!])
	AC_MSG_NOTICE([==> Set DIFF variable if present in non-standard path!])
	AC_MSG_ERROR([diff is mandatory to run the tests : will stop here!])
fi

# Check for 'cat' and get full path.
AC_ARG_VAR([CAT],[the 'cat' program used for printing test output files])
AC_PATH_PROG([CAT],[cat],[])
if test "x$CAT" = "x"; then
	AC_MSG_NOTICE([==> cat command not found!])
	AC_MSG_NOTICE([==> Set CAT variable if present in non-standard path!])
	AC_MSG_ERROR([cat is mandatory to run the tests : will stop here!])
fi

# Check for 'printf' and get full path.
AC_ARG_VAR([PRINTF],[the 'printf' program used to print test information])
AC_PATH_PROG([PRINTF],[printf],[])
if test "x$PRINTF" = "x"; then
	AC_MSG_NOTICE([==> printf command not found!])
	AC_MSG_NOTICE([==> Set PRINTF variable if present in non-standard path!])
	AC_MSG_ERROR([printf is mandatory to run the tests : will stop here!])
fi

AC_CONFIG_FILES([tests/Makefile])
AC_CONFIG_FILES([Makefile liblsb.pc])
AC_OUTPUT

