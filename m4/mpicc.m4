AC_DEFUN([NG_WITH_MPICC],
    [AC_MSG_NOTICE([*** checking for MPI support ***])
    AC_ARG_WITH(mpi,
        AC_HELP_STRING([--with-mpi], [compile with MPI support (ARG can be the path to the root MPI directory, if mpicc is not in PATH)]))
    
        printf "${with_mpi}"
        if test x"${with_mpi}" != xno; then
            if test x"${with_mpi}" != x -a x"${with_mpi}" != xyes; then
                AC_CHECK_PROG(ng_mpicc_found, mpicxx, yes, no, ${with_mpi}/bin/)
                CC=${with_mpi}/bin/mpicc
                if test x${ng_mpicc_found} = xno; then
                    AC_MSG_ERROR(${with_mpi}/mpicc selected but not found)
                fi
                AC_DEFINE(HAVE_MPI, 1, enables the MPI specific code)
                have_mpi=1
            else
                # if the environment variable MPICC is set
                if test x${MPICC} != x; then
                    AC_CHECK_PROG(ng_mpicc_found, $MPICC, yes, no)
                    MYCC=$MPICC
                else  
                    AC_CHECK_PROG(ng_mpicc_found, mpicc, yes, no)
                    MYCC=mpicc
                fi;

                if test x${ng_mpicc_found} = xno; then
                    AC_MSG_NOTICE(B compiling without MPI support)
                else 
                    CC=$MYCC
                    AC_MSG_NOTICE(using $CC for MPI support)
                    AC_DEFINE(HAVE_MPI, 1, enables the MPI specific code)
                    have_mpi=1
                fi
            fi
        else
            AC_MSG_NOTICE(A compiling without MPI support) 
        fi
    ]
)

