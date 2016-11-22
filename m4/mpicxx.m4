AC_DEFUN([NG_WITH_MPICXX],
    [AC_MSG_NOTICE([*** checking for MPI C++ support ***])
        if test x"${with_mpi}" != xno; then
            if test x"${with_mpi}" != x -a x"${with_mpi}" != xyes; then
                AC_CHECK_PROG(ng_mpicxx_found, mpicxx, yes, no, ${with_mpi}/bin/)
                CXX=${with_mpi}/bin/mpicxx
                if test x${ng_mpicxx_found} = xno; then
                    AC_MSG_ERROR(${with_mpi}/mpicxx selected but not found)
                fi
            else
                # if the environment variable mpicxx is set
                if test x${MPICXX} != x; then
                    AC_CHECK_PROG(ng_mpicxx_found, $MPICXX, yes, no)
                    MYCXX=$MPICXX
                else  
                    AC_CHECK_PROG(ng_mpicxx_found, mpicxx, yes, no)
                    MYCXX=mpicxx
                fi;

                if test x${ng_mpicxx_found} = xno; then
                    AC_MSG_NOTICE(compiling without MPI support)
                else 
                    CXX=$MYCXX
                    AC_MSG_NOTICE(using $CXX for MPI C++ support)
                fi
            fi
        else
            AC_MSG_NOTICE(compiling without MPI support)
        fi
    ]
)

