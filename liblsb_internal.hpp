/*
* Copyright (c) 2015 ETH-Zurich. All rights reserved.
* Use of this source code is governed by a BSD-style license that can be
* found in the LICENSE file.
*/

#include "config.h" 
#include "liblsb_lvector.hpp"
#include "uthash/uthash.h"
#include "liblsb.h"


#define MAX_DOUBLE 1e99

//#define HAVE_MPI 1 
#define USE_LVECTOR

#ifdef USE_LVECTOR
#define VECTOR_TYPE lvector
#else
#define VECTOR_TYPE std::vector
#endif

#ifdef HAVE_PAPI
#include <papi.h> /* This needs to be included every time you use PAPI */
#endif

#ifdef HAVE_UNWIND
#define UNW_LOCAL_ONLY
#include <libunwind.h>
#endif


#include <stdlib.h>
#include <vector>
#include <signal.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/utsname.h>
#include <stdarg.h>



#ifdef HAVE_LIKWID
extern "C" {
#include <bstrlib.h>
#include <cpuid.h>
#include <msr.h>
#include <numa.h>
#include <perfmon.h>
#include <timer.h>
#include <registers.h>
#include <hashTable.h>
#include <likwid.h>
}
#endif

#ifdef HAVE_PAPI
// todo: make this a configurable parameter
#define MAX_PAPI_COUNTERS 16
#endif


/**
 * a single record field to be inserted during each record call
 */
typedef struct {
#ifdef HAVE_PAPI
  long long papictrs[MAX_PAPI_COUNTERS]; //< PAPI counters values
#endif
#ifdef HAVE_LIKWID
  long long PMcounter0; ///< likwid MSR counter values
  long long PMcounter1; ///< likwid MSR counter values
  long long PMcounter2; ///< likwid MSR counter values
  long long PMcounter3; ///< likwid MSR counter values
  long long PMcounter4; ///< likwid MSR counter values
  long long PMcounter5; ///< likwid MSR counter values
  long long PMcounter6; ///< likwid MSR counter values
  long long PMcounter7; ///< likwid MSR counter values
  long long PMcounter8; ///< likwid MSR counter values
  long long PMcounter9; ///< likwid MSR counter values
#endif
  double tduration; ///< the duration of the current epoch
  unsigned long toverhead; ///< the duration (overhead) of profiling for this epoch
  unsigned int id; ///< the user-specified kernel id
} t_rec;


typedef enum{
    LSB_INT=0, LSB_DBL, LSB_STR, LSB_LNG
} lsb_type_t;

typedef struct{
  unsigned int start_idx;
  const char * pindex;
  lsb_type_t type;
  union{
    int value;
    const char * str;
    double dbl;
    int64_t lng;
  };
} t_recparam;


/*hash table*/
typedef struct{
  const char * pname;
  int value;
  UT_hash_handle hh;
} t_hparam;

/**
 * the global data for the profiler
 * this holds all state related with the profiling session
 */
typedef struct {
  VECTOR_TYPE<t_rec> recs; ///< records
  VECTOR_TYPE<t_recparam> recparams;

#ifdef HAVE_PAPI
  int eventset; ///< PAPI eventset config
  VECTOR_TYPE<int> papi_ctrs_ids;
  int num_papi_ctrs; //< len(papi_ctrs_ids)
#endif
#ifdef HAVE_UNWIND
  VECTOR_TYPE< std::vector<long> > iptrace; ///< the instruction pointer values for the stack for each epoch
#endif

  int r; ///< MPI rank
  int p; ///< MPI processes in world

  FILE *outfd; ///< the output file
  double tstart; ///< the program start time
  double tlast; ///< lasttime starttime of the current epoch
  char print; ///< internal flag if this process prints or not
  char write_file; ///< internal flag if this process writes a file or not
  char rec_enabled; ///< true if recording is enabled, false if it is not
  char lsb_disabled; ///< true if lsb is completely disabled
  int write_header;
  int next;
  int group_ptr; //< starting index of the current group
 
#if defined(SYNC_WINDOW) && defined(HAVE_MPI)
  double sync_window;
  MPI_Comm sync_comm;
#endif


} LSB_Data;

#define CHK_DISABLED if(lsb_data->lsb_disabled) return;
void lsb_sort (double *data, int count);
double lsb_median_from_sorted_data (const double sorted_data[], const size_t stride, const size_t n);
double lsb_quantile_from_sorted_data (const double sorted_data[], const size_t stride, const size_t n, const double f);
int lsb_fit_linear (const double *x, const size_t xstride, const double *y, const size_t ystride, const size_t n, double *c0, double *c1, double *cov_00, double *cov_01, double *cov_11, double *sumsq);
double lsb_mean (const double data[], const size_t stride, const size_t size);

double min(double a, double b);
double max(double a, double b);

// prototypes
#include "liblsb.h"
