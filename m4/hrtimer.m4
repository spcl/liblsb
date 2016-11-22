# $1 - arch, $2 - arch_num, $3 - subdir
AC_DEFUN([HTOR_TIMERARCH], [
thisarch=$1
AC_MSG_CHECKING([if $thisarch (HRT_ARCH $2) assembler portion compiles])
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[
#define HRT_ARCH $2
#include "$3/hrtimer.h"
]],
[[
static volatile HRT_TIMESTAMP_T t1, t2;
static volatile UINT64_T elapsed_ticks;
unsigned long long freq;
HRT_INIT(0,freq);
HRT_GET_TIMESTAMP(t1);
HRT_GET_TIMESTAMP(t2);
HRT_GET_ELAPSED_TICKS(t1, t2, &elapsed_ticks);
]])], 
  [
   AC_MSG_RESULT([yes])
   # no sanity check for MPI_Wtime() :)
   if test "$thisarch" != "wtime"; then
     # do sanity check ...
     printf "compiling sanity check ... "
     $CC $CFLAGS -DHRT_ARCH=$2 $hr_subdir/sanity-check.c -o sanity
     if test "$?" != "0"; then exit 1; fi;
     printf "done\n"
     if ./sanity; then 
       hrt_arch=$thisarch
       hrt_arch_num=$2 
       AC_DEFINE([HRT_ARCH], $2, [highrestimer architecture])
       AC_DEFINE([HAVE_HRTIMER], [], [can we use a high resolution timer])
       AC_MSG_NOTICE([found $thisarch])
     else
       printf "sanity check failed\n"
     fi;
     rm -f sanity
   else
     # this is copied from above -- nasty ... 
     hrt_arch=$thisarch
     hrt_arch_num=$2 
     AC_DEFINE([HRT_ARCH], $2, [highrestimer architecture])
     AC_DEFINE([HAVE_HRTIMER, [], [can we use a high resolution timer]])
     AC_MSG_NOTICE([found $thisarch])
   fi;
  ],
  [
   AC_MSG_RESULT([no])
  ])])


