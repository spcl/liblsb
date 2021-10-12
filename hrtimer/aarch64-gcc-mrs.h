/*
 * Copyright (c) 2020 ETH-Zurich. All rights reserved.
 * Use of this source code is governed by a BSD-style license that can be
 *  found in the LICENSE file.
 * Author(s): Salvatore Di Girolamo <digirols@inf.ethz.ch>
 *
 */

#include <inttypes.h>
#include <stdio.h>
#include "calibrate.h"
#define UINT32_T uint32_t
#define UINT64_T uint64_t

#define HRT_INIT(print, freq) do {\
  if(print) printf("# initializing aarch64 timer (takes some seconds)\n"); \
  HRT_CALIBRATE(freq); \
} while(0) 


#define HRT_TIMESTAMP_T UINT64_T

#define HRT_GET_TIMESTAMP(t1)  __asm__ __volatile__ ("mrs %0, cntvct_el0" : "=r"(t1));

#define HRT_GET_ELAPSED_TICKS(t1, t2, numptr) *numptr = (t2 - t1)*8;

#define HRT_GET_TIME(t1, time) time = t1;

