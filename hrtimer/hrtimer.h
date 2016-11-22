/*
 * Copyright (c) 2009 The Trustees of Indiana University and Indiana
 *                    University Research and Technology
 *                    Corporation.  All rights reserved.
 *
 * Author(s): Torsten Hoefler <htor@cs.indiana.edu>
 *             Timo Schneider <timos@perlplexity.org>
 *
 */

#include <math.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdint.h>
#include <inttypes.h>
#include <unistd.h>

#ifdef HRT_ARCH

#if HRT_ARCH==1
#include "x86_32-gcc-rdtsc.h"
#endif
#if HRT_ARCH==2
#include "x86_64-gcc-rdtsc.h"
#endif
#if HRT_ARCH==3
#include "ppc-gcc-tb.h"
#endif
#if HRT_ARCH==4
#include "ia64-gcc-itc.h"
#endif
#if HRT_ARCH==5
#include "mips64-sicortex-gcc.h"
#endif
#if HRT_ARCH==6
#include "mpi-wtime.h"
#endif

/* global timer frequency in Hz */
extern uint64_t liblsb_g_timerfreq;

static double HRT_GET_USEC(uint64_t ticks) {
  return 1e6*(double)ticks/(double)liblsb_g_timerfreq;
}


static int sanity_check(int print) {

	HRT_TIMESTAMP_T t1, t2;
	uint64_t s, s2, s3;

	int sanity = 1;

	if (print < 0) return 0; // this is so that we can call sanity_check from HRT_INIT without it doing anything --- otherwise compilers might warn about unused functions

	if (print == 1) printf("# Sanity check of the timer\n");
	HRT_GET_TIMESTAMP(t1);
	sleep(1);
	HRT_GET_TIMESTAMP(t2);
	HRT_GET_ELAPSED_TICKS(t1, t2, &s);
	if (print == 1) printf("# sleep(1) needed %llu ticks.\n", (long long unsigned int) s);
	HRT_GET_TIMESTAMP(t1);
	sleep(2);
	HRT_GET_TIMESTAMP(t2);
	HRT_GET_ELAPSED_TICKS(t1, t2, &s2);
	if (print == 1) printf("# sleep(2)/2 needed %llu ticks.\n", (long long unsigned int) (s2/2));
	if (fabs((double)s - (double)s2/2) > (s*0.05) || s < 1) {
		sanity = 0;
		printf("# The high performance timer gives bogus results on this system!\n");
	}
	s = s2/2;
	HRT_GET_TIMESTAMP(t1);
	sleep(3);
	HRT_GET_TIMESTAMP(t2);
	HRT_GET_ELAPSED_TICKS(t1, t2, &s3);
	if (print == 1) printf("# sleep(3)/3 needed %llu ticks.\n", (long long unsigned int) s3/3);
	if (fabs((double)s - (double)s3/3) > (s*0.05) || s < 1) {
		sanity = 0;
		printf("# The high performance timer gives bogus results on this system!\n");
	}
	return sanity;
}

#endif

