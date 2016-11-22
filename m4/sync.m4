AC_DEFUN([LSB_ENABLE_SYNC],
    [AC_ARG_ENABLE(sync, AC_HELP_STRING([--enable-sync], [enables time-window synchronization support.]))


    #echo "####### x${enable_sync}"
    if test "x${enable_sync}" == "xyes"; then

        AC_DEFINE(HAVE_SYNC, 1, enables time-window synch)
        AC_MSG_NOTICE([SYNC support enabled])

    fi
    ]
)

