/*
* Copyright (c) 2015 ETH-Zurich. All rights reserved.
* Use of this source code is governed by a BSD-style license that can be
* found in the LICENSE file.
*/

#include "config.h"
#include "liblsb_internal.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif
#ifdef HAVE_HRTIMER
#include "hrtimer/hrtimer.h"
#endif


#include <algorithm>

LSB_Data *lsb_data;

#ifdef HAVE_LIKWID
std::vector<long long> PMvalues(10);
#endif

#ifdef HAVE_HRTIMER
uint64_t liblsb_g_timerfreq;
double liblsb_g_hrtimer_startvalue;
#endif

/**
 * end the current epoch and start a new one
 * same as LSB_Rec() but doesn't record anything -> all measurements will be lost
 */
void LSB_Res() {
  CHK_DISABLED;
#ifdef _OPENMP
#pragma omp master
{
#endif
#if defined HAVE_PAPI 
  std::vector<long long> values(2);
  PAPI_stop(lsb_data->eventset, &values[0]);
  lsb_data->tlast = PAPI_get_real_usec();
  PAPI_start(lsb_data->eventset);
#elif defined HAVE_LIKWID
    int cpu_id;
    cpu_id = likwid_getProcessorId();
    PMvalues[0] = msr_read(cpu_id, MSR_AMD15_PMC0);
    PMvalues[1] = msr_read(cpu_id, MSR_AMD15_PMC1);
    PMvalues[2] = msr_read(cpu_id, MSR_AMD15_PMC2);
    PMvalues[3] = msr_read(cpu_id, MSR_AMD15_PMC3);
    PMvalues[4] = msr_read(cpu_id, MSR_AMD15_PMC4);
    PMvalues[5] = msr_read(cpu_id, MSR_AMD15_PMC5);
    PMvalues[6] = msr_read(cpu_id, MSR_AMD15_NB_PMC0);
    PMvalues[7] = msr_read(cpu_id, MSR_AMD15_NB_PMC1);
    PMvalues[8] = msr_read(cpu_id, MSR_AMD15_NB_PMC2);
    PMvalues[9] = msr_read(cpu_id, MSR_AMD15_NB_PMC3);
    unsigned long long ticks;
    HRT_TIMESTAMP_T t;
    HRT_GET_TIMESTAMP(t);
    HRT_GET_TIME(t, ticks);
    lsb_data->tlast =  HRT_GET_USEC(ticks);
#elif defined HAVE_HRTIMER
	unsigned long long ticks;
	HRT_TIMESTAMP_T t;
	HRT_GET_TIMESTAMP(t);
	HRT_GET_TIME(t, ticks);
	lsb_data->tlast = HRT_GET_USEC(ticks);
#else
#error "No possibility to do accurate timing!\n"
#endif
#ifdef _OPENMP
}
#endif
}

double LSB_Wait(double microseconds){

    double start, current;
	unsigned long long ticks;

#if defined HAVE_PAPI
    start = PAPI_get_real_usec();
#else
    HRT_TIMESTAMP_T t;
    HRT_GET_TIMESTAMP(t);
    HRT_GET_TIME(t, ticks);
    start = HRT_GET_USEC(ticks);
#endif
    
    current = start;

    while (current-start<microseconds){
#if defined HAVE_PAPI
        current = PAPI_get_real_usec();
#else
        HRT_TIMESTAMP_T t;
        HRT_GET_TIMESTAMP(t);
        HRT_GET_TIME(t, ticks);
        current = HRT_GET_USEC(ticks);
#endif       
    } 
   
    return current - start;
}


void LSB_Group_Begin(LSB_Group_t *group){
    group->start = lsb_data->recs.size();
    group->stop = -1;
}

void LSB_Group_End(LSB_Group_t *group){
    group->stop = lsb_data->recs.size() - 1;
}

void _LSB_Fold(unsigned int id, int start, int stop, lsb_op_t op, double * result){


    double res=-1;
    std::vector<double> measures;
    int j=0;
       

    for (int i=start; i<stop; i++){

        if (id!=LSB_ANY && lsb_data->recs[i].id != id) continue;

        switch (op){
            case LSB_SUM:
                res+=lsb_data->recs[i].tduration;            
                break;
            case LSB_COUNT:
                res++;
                break;        
            case LSB_MEDIAN:
                measures.push_back(lsb_data->recs[i].tduration);
                break;
            case LSB_MAX:
                if (i==0 || res<lsb_data->recs[i].tduration) 
                    res = lsb_data->recs[i].tduration;
                break;
            case LSB_MIN:
                if (i==0 || res>lsb_data->recs[i].tduration) 
                    res = lsb_data->recs[i].tduration;
                break;
        } 
    }


    if (op==LSB_MEDIAN && measures.size()>0){
        std::sort(measures.begin(), measures.end());

        //for (int j=0; j<measures.size(); j++) printf("%i: %f\n", j, measures[j]);
        res = measures[measures.size()/2];
        //printf("returning res: %lf; size: %i; measures[0]: %lf; size/2: %i; measures[size-1]: %lf\n", res, measures.size(), measures[0], measures.size()/2, measures[measures.size()-1]);
    }

    *result = res;

}

void LSB_Fold(unsigned int id, lsb_op_t op, double * result){
    _LSB_Fold(id, 0, lsb_data->recs.size(), op, result);
}

void LSB_Group_Fold(LSB_Group g, unsigned int id, lsb_op_t op, double * result){
    if (g.start < 0 || g.start > lsb_data->recs.size()) { *result=-1; return; }
    if (g.stop!=-1 && (g.stop<0 || g.stop > lsb_data->recs.size())) { *result=-1; return; }
    
    _LSB_Fold(id, g.start, (g.stop!=-1) ? g.stop : lsb_data->recs.size(), op, result);
}


double LSB_Rec(unsigned int id){
   double res = LSB_Stop(id, 1);
   lsb_data->next = lsb_data->recs.size();
   return res;
}

void LSB_Next(){
   lsb_data->next = lsb_data->recs.size();
}

double LSB_Check(unsigned int id){
   double res = LSB_Stop(id, 0);
   lsb_data->next = lsb_data->recs.size();
   return res;
}

/**
 * write a tracing record
 * this ends the epoch with the specified id and opens a new epoch
 * @param id user-defined kernel id
 */
double LSB_Stop(unsigned int id, unsigned int reset) {
  if (lsb_data->lsb_disabled) return -1;
  //CHK_DISABLED;
#ifdef _OPENMP
#pragma omp master
{
#endif
  double tstart, tlast; // epoch duration and function start time
#if defined HAVE_PAPI
  tstart  = PAPI_get_real_usec();

  std::vector<long long> values(2);
  PAPI_stop(lsb_data->eventset, &values[0]);
#elif defined HAVE_LIKWID
    //std::vector<long long> values(10);
    int cpu_id;
    cpu_id = likwid_getProcessorId();
    PMvalues[0] = msr_read(cpu_id, MSR_AMD15_PMC0) - PMvalues[0];
    PMvalues[1] = msr_read(cpu_id, MSR_AMD15_PMC1) - PMvalues[1];
    PMvalues[2] = msr_read(cpu_id, MSR_AMD15_PMC2) - PMvalues[2];
    PMvalues[3] = msr_read(cpu_id, MSR_AMD15_PMC3) - PMvalues[3];
    PMvalues[4] = msr_read(cpu_id, MSR_AMD15_PMC4) - PMvalues[4];
    PMvalues[5] = msr_read(cpu_id, MSR_AMD15_PMC5) - PMvalues[5];
    PMvalues[6] = msr_read(cpu_id, MSR_AMD15_NB_PMC0) - PMvalues[6];
    PMvalues[7] = msr_read(cpu_id, MSR_AMD15_NB_PMC1) - PMvalues[7];
    PMvalues[8] = msr_read(cpu_id, MSR_AMD15_NB_PMC2) - PMvalues[8];
    PMvalues[9] = msr_read(cpu_id, MSR_AMD15_NB_PMC3) - PMvalues[9];
    unsigned long long ticks;
    HRT_TIMESTAMP_T t;
    HRT_GET_TIMESTAMP(t);
    HRT_GET_TIME(t, ticks);
    tstart = HRT_GET_USEC(ticks);
#elif defined HAVE_HRTIMER
	unsigned long long ticks;
	HRT_TIMESTAMP_T t;
	HRT_GET_TIMESTAMP(t);
	HRT_GET_TIME(t, ticks);
	tstart = HRT_GET_USEC(ticks);
#else
#error "No possibility to do accurate timing!\n"
#endif
  
  double measure = tstart-lsb_data->tlast;
#if defined HAVE_PAPI
  t_rec rec = {values[0], values[1], measure, 0, id};
#elif defined HAVE_LIKWID
  t_rec rec = {PMvalues[0], PMvalues[1], PMvalues[2], PMvalues[3], PMvalues[4], \
               PMvalues[5], PMvalues[6], PMvalues[7], PMvalues[8], PMvalues[9], \
               measure, 0, id};
#else
  t_rec rec = {measure, 0, id};
#endif
  if(lsb_data->rec_enabled) lsb_data->recs.push_back(rec);

#ifdef HAVE_UNWIND
  unw_cursor_t cursor; unw_context_t uc;
  unw_word_t ip;
  unw_getcontext(&uc);

  if(lsb_data->rec_enabled) {
    lsb_data->iptrace.resize(lsb_data->iptrace.size()+1);
    unw_init_local(&cursor, &uc);
    while (unw_step(&cursor) > 0) {
      unw_get_reg(&cursor, UNW_REG_IP, &ip);
      //unw_get_reg(&cursor, UNW_REG_SP, &sp);
      lsb_data->iptrace[lsb_data->iptrace.size()-1].push_back(ip);
      //printf ("ip = %lx, sp = %lx\n", (long) ip, (long) sp);
    }
  }
#endif

#if defined HAVE_PAPI
  PAPI_start(lsb_data->eventset);

  tlast = PAPI_get_real_usec();
#elif  HAVE_LIKWID
    PMvalues[0] = msr_read(cpu_id, MSR_AMD15_PMC0);
    PMvalues[1] = msr_read(cpu_id, MSR_AMD15_PMC1);
    PMvalues[2] = msr_read(cpu_id, MSR_AMD15_PMC2);
    PMvalues[3] = msr_read(cpu_id, MSR_AMD15_PMC3);
    PMvalues[4] = msr_read(cpu_id, MSR_AMD15_PMC4);
    PMvalues[5] = msr_read(cpu_id, MSR_AMD15_PMC5);
    PMvalues[6] = msr_read(cpu_id, MSR_AMD15_NB_PMC0);
    PMvalues[7] = msr_read(cpu_id, MSR_AMD15_NB_PMC1);
    PMvalues[8] = msr_read(cpu_id, MSR_AMD15_NB_PMC2);
    PMvalues[9] = msr_read(cpu_id, MSR_AMD15_NB_PMC3);
    HRT_GET_TIMESTAMP(t);
    HRT_GET_TIME(t, ticks);
    tlast =  HRT_GET_USEC(ticks);
#elif defined HAVE_HRTIMER
	HRT_GET_TIMESTAMP(t);
	HRT_GET_TIME(t, ticks);
	tlast = HRT_GET_USEC(ticks);
#else
#error "No possibility to do accurate timing!\n"
#endif
  if(lsb_data->rec_enabled) lsb_data->recs.back().toverhead = tlast-tstart;
  if (reset) lsb_data->tlast = tlast;
  return measure;
#ifdef _OPENMP
}
#endif
}

/**
 * simple handler for SIGALRM callback
 * simply calls LSB_Rec to record an epoch
 * @param x epoch identifier
 */
void LSB_Handler(int x) {
  LSB_Rec(0);
}

/**
 *  flush data to disk and end current measurement
 *  this is only used if an application is started multiple times with different parameters
 */
void LSB_Flush() {
  CHK_DISABLED;
#ifdef _OPENMP
#pragma omp master
{
#endif
  //printf("starting flushing\n");
  double tend;
#if defined HAVE_PAPI
  tend = PAPI_get_real_usec();
#elif defined HAVE_HRTIMER
  	unsigned long long ticks;
	HRT_TIMESTAMP_T t;
	HRT_GET_TIMESTAMP(t);
	HRT_GET_TIME(t, ticks);
	tend = HRT_GET_USEC(ticks);
#else
#error "No possibility to do accurate timing!\n"
#endif

  
  FILE *fp=lsb_data->outfd;

  enum {PRETTY, EFFICIENT, ACCUMULATED} style = PRETTY; // printing style
  char *env = getenv("LSB_OUTPUT_FORMAT"); // print format
  if(env != NULL) {
    if(strstr(env, "efficient")) style = EFFICIENT;
    if(strstr(env, "accumulated")) style = ACCUMULATED;
  }

  //printf("flushing: style: %i\n", style);

  /* Create hashtable for the parameters */
  t_hparam * ptable = NULL;
  t_hparam * found;
  t_hparam * current;
  t_hparam * toadd;
  int pnum=0;
  //printf("creating hashtable\n");
  for (int i=0; i<lsb_data->recparams.size(); i++){    
    HASH_FIND_STR(ptable, lsb_data->recparams[i].pindex, found);
    if (found==NULL) {
        toadd =  (t_hparam *) malloc(sizeof(t_hparam));
        //printf("hashtable add: %s\n", lsb_data->recparams[i].pindex);
        toadd->pname = lsb_data->recparams[i].pindex;
        toadd->value = pnum++;
        HASH_ADD_KEYPTR(hh, ptable, toadd->pname, strlen(toadd->pname), toadd);
    }
  }

  int rparams_ptr = 0;
  t_recparam * rparams = (t_recparam *) malloc(sizeof(t_recparam)*(pnum)); 
  for (int i=0; i<pnum; i++) rparams[i].start_idx=-1; 
  const char ** pnames = (const char **) malloc(sizeof(char *)*(pnum));

  for (current=ptable; current!=NULL; current= (t_hparam *) current->hh.next){
     pnames[current->value] = current->pname;
  }
  if (lsb_data->write_header){  
#ifdef HAVE_MPI
    if(lsb_data->write_file==1) fprintf(fp, "# MPI execution on rank %i with %i processes in world\n", lsb_data->r, lsb_data->p);
#endif
  fprintf(fp, "# Reported time measurements are in microseconds\n");
  }
  unsigned long long dur = 0; for(unsigned int i=0; i<lsb_data->recs.size(); ++i) dur += lsb_data->recs[i].toverhead;
  //printf("tend: %lf; tstart: %lf\n", tend, lsb_data->tstart);
  //if(lsb_data->write_file==1) fprintf(fp, "# Runtime: %lf s (overhead: %lf %%) %i records\n", (double)(tend-lsb_data->tstart)/1e6, (double)dur/(tend-lsb_data->tstart)*100, (int)lsb_data->recs.size());
  

  if((style == PRETTY) && (lsb_data->write_file==1)) {
    //printf("pretty printing\n");
    if (lsb_data->write_header) {
      fprintf(fp, "# pretty output format\n");
      for (int j=0; j<pnum; j++) fprintf(fp, " %12s", pnames[j]);
      fprintf(fp, " %8s %12s %12s ", "id", "time", "overhead");
    
#ifdef HAVE_PAPI
      char p1[1024], p2[1024];
      PAPI_event_code_to_name(lsb_data->papi1, p1);
      PAPI_event_code_to_name(lsb_data->papi2, p2);
      fprintf(fp, "%16s %16s", p1, p2);
#endif
#ifdef HAVE_UNWIND
      fprintf(fp, " #ips ip-list");
#endif
      fprintf(fp, "\n");
    }
    // print all records
    for(unsigned int i=0; i < lsb_data->recs.size(); ++i) {
      //printf("aa\n");
      //printf("printing records; rparams_ptr: %i; recparams.size(): %i; i: %i; \n", rparams_ptr, lsb_data->recparams.size(), i);

      //printf("printing records; rparams_ptr: %i; recparams.size(): %i; i: %i; start_idx: %i; \n", rparams_ptr, lsb_data->recparams.size(), i, lsb_data->recparams[rparams_ptr].start_idx);
      /* search for rparam update */
      while (rparams_ptr < lsb_data->recparams.size() && lsb_data->recparams[rparams_ptr].start_idx <= i){
        //printf("advancing values\n");
        HASH_FIND_STR(ptable, lsb_data->recparams[rparams_ptr].pindex, current);
        rparams[current->value] = (lsb_data->recparams[rparams_ptr]);
        rparams_ptr++;
      }
     
      //printf("printing parameters\n");   
      for (int j=0; j<pnum; j++) {
        if (rparams[j].start_idx==-1) fprintf(fp, "NA");
        else if (rparams[j].type==LSB_STR) fprintf(fp, " %12s ", rparams[j].str);
        else if (rparams[j].type==LSB_INT) fprintf(fp, " %12i ", rparams[j].value);
        else if (rparams[j].type==LSB_DBL) fprintf(fp, " %12f ", rparams[j].dbl);
        else if (rparams[j].type==LSB_LNG) fprintf(fp, " %12ld ", rparams[j].lng);
      }
      //printf("end printing paramaters\nstart printing data\n");
    
      fprintf(fp, "  %8i ", lsb_data->recs[i].id);
      fprintf(fp, "%f ", lsb_data->recs[i].tduration);
      fprintf(fp, "%12lu ", lsb_data->recs[i].toverhead);
      //printf("end printing data111\n");
  #ifdef HAVE_PAPI
      //printf("have papi\n");
      fprintf(fp, "%16lli ", lsb_data->recs[i].papictr1);
      fprintf(fp, "%16lli ", lsb_data->recs[i].papictr2);
  #endif
  #ifdef HAVE_LIKWID
      //printf("have likwid\n");
      fprintf(fp, "%16lli ", lsb_data->recs[i].PMcounter0);
      fprintf(fp, "%16lli ", lsb_data->recs[i].PMcounter1);
      fprintf(fp, "%16lli ", lsb_data->recs[i].PMcounter2);
      fprintf(fp, "%16lli ", lsb_data->recs[i].PMcounter3);
      fprintf(fp, "%16lli ", lsb_data->recs[i].PMcounter4);
      fprintf(fp, "%16lli ", lsb_data->recs[i].PMcounter5);
      fprintf(fp, "%16lli ", lsb_data->recs[i].PMcounter6);
      fprintf(fp, "%16lli ", lsb_data->recs[i].PMcounter7);
      fprintf(fp, "%16lli ", lsb_data->recs[i].PMcounter8);
      fprintf(fp, "%16lli ", lsb_data->recs[i].PMcounter9);
  #endif
  #ifdef HAVE_UNWIND
      //printf("have unwind\n");
      fprintf(fp, "%i ", (int)lsb_data->iptrace[i].size());
      for(unsigned int j=0; j<lsb_data->iptrace[i].size(); ++j) fprintf(fp, "%lx ", lsb_data->iptrace[i][j]);
  #endif
      fprintf(fp, "\n");
    }

    //printf("printing ints\n");
    if(lsb_data->recints.size()) {
      // print all ints
      fprintf(fp, "# printing integer list, %i records\n  %8s %8s %8s\n", (int)lsb_data->recints.size(), "id", "int1", "int2");
      for(unsigned int i=0; i < lsb_data->recints.size(); ++i) {
        fprintf(fp, "  %8i ", lsb_data->recints[i].id);
        fprintf(fp, "%8i ", lsb_data->recints[i].int1);
        fprintf(fp, "%8i ", lsb_data->recints[i].int2);
        fprintf(fp, "\n");
      }
    }
    //printf("end printing ints\nprinting dbls");

    if(lsb_data->recintdbl.size()) {
      // print all ints
      fprintf(fp, "# printing integer/double list, %i records\n  %8s %8s %16s\n", (int)lsb_data->recintdbl.size(), "id", "int1", "double");
      for(unsigned int i=0; i < lsb_data->recintdbl.size(); ++i) {
        fprintf(fp, "  %8i ", lsb_data->recintdbl[i].id);
        fprintf(fp, "%8i ", lsb_data->recintdbl[i].int1);
        fprintf(fp, "%16f ", lsb_data->recintdbl[i].dbl);
        fprintf(fp, "\n");
      }
    }
    
    //printf("end pretty printing\n");
  } else if ((style == EFFICIENT) && (lsb_data->write_file==1)) {
    if (lsb_data->write_header){
      fprintf(fp, "# efficient output format\n");
      for (int j=0; j<pnum; j++) fprintf(fp, " P%i", j);
      fprintf(fp, "  %s %s %s ", "id", "time", "overhead");
    }

#ifdef HAVE_PAPI
    char p1[1024], p2[1024];
    PAPI_event_code_to_name(lsb_data->papi1, p1);
    PAPI_event_code_to_name(lsb_data->papi2, p2);
    fprintf(fp, "%s %s", p1, p2);
#endif
#ifdef HAVE_UNWIND
    fprintf(fp, " #ips ip-list");
#endif
    fprintf(fp, "\n");

    // print all records
    for(unsigned int i=0; i < lsb_data->recs.size(); ++i) {


      while (rparams_ptr < lsb_data->recparams.size() && lsb_data->recparams[rparams_ptr].start_idx <= i){
        HASH_FIND_STR(ptable, lsb_data->recparams[rparams_ptr].pindex, current);
        rparams[current->value] = (lsb_data->recparams[rparams_ptr]);
        rparams_ptr++;
      }

      for (int j=0; j<pnum; j++) {
        if (rparams[j].start_idx==-1) fprintf(fp, "NA");
        else if (rparams[j].type==LSB_STR) fprintf(fp, " %12s ", rparams[j].str);
        else if (rparams[j].type==LSB_INT) fprintf(fp, " %12i ", rparams[j].value);
        else if (rparams[j].type==LSB_DBL) fprintf(fp, " %12f ", rparams[j].dbl);
        else if (rparams[j].type==LSB_LNG) fprintf(fp, " %12ld ", rparams[j].lng);
      }


      fprintf(fp, "%i ", lsb_data->recs[i].id);
      fprintf(fp, "%f ", lsb_data->recs[i].tduration);
      fprintf(fp, "%lu ", lsb_data->recs[i].toverhead);
     

  #ifdef HAVE_PAPI
      fprintf(fp, "%lli ", lsb_data->recs[i].papictr1);
      fprintf(fp, "%lli ", lsb_data->recs[i].papictr2);
  #endif
  #ifdef HAVE_LIKWID
      fprintf(fp, "%lli ", lsb_data->recs[i].PMcounter0);
      fprintf(fp, "%lli ", lsb_data->recs[i].PMcounter1);
      fprintf(fp, "%lli ", lsb_data->recs[i].PMcounter2);
      fprintf(fp, "%lli ", lsb_data->recs[i].PMcounter3);
      fprintf(fp, "%lli ", lsb_data->recs[i].PMcounter4);
      fprintf(fp, "%lli ", lsb_data->recs[i].PMcounter5);
      fprintf(fp, "%lli ", lsb_data->recs[i].PMcounter6);
      fprintf(fp, "%lli ", lsb_data->recs[i].PMcounter7);
      fprintf(fp, "%lli ", lsb_data->recs[i].PMcounter8);
      fprintf(fp, "%lli ", lsb_data->recs[i].PMcounter9);
  #endif
  #ifdef HAVE_UNWIND
      fprintf(fp, "%i ", (int)lsb_data->iptrace[i].size());
      for(unsigned int j=0; j<lsb_data->iptrace[i].size(); ++j) fprintf(fp, "%lx ", lsb_data->iptrace[i][j]);
  #endif
      fprintf(fp, "\n");
    }

    if(lsb_data->recints.size()) {
      // print all ints
      fprintf(fp, "# printing integer list, %i records\n  %8s %8s %8s\n", (int)lsb_data->recints.size(), "id", "int1", "int2");
      for(unsigned int i=0; i < lsb_data->recints.size(); ++i) {
        fprintf(fp, "%i ", lsb_data->recints[i].id);
        fprintf(fp, "%i ", lsb_data->recints[i].int1);
        fprintf(fp, "%i ", lsb_data->recints[i].int2);
        fprintf(fp, "\n");
      }
    }

    if(lsb_data->recintdbl.size()) {
      // print all ints
      fprintf(fp, "# printing integer/double list, %i records\n  %8s %8s %16s\n", (int)lsb_data->recintdbl.size(), "id", "int1", "double");
      for(unsigned int i=0; i < lsb_data->recintdbl.size(); ++i) {
        fprintf(fp, "%i ", lsb_data->recintdbl[i].id);
        fprintf(fp, "%i ", lsb_data->recintdbl[i].int1);
        fprintf(fp, "%f ", lsb_data->recintdbl[i].dbl);
        fprintf(fp, "\n");
      }
    }
  } else if (style == ACCUMULATED) { /* TODO: extend with rparams */
#ifdef HAVE_MPI
    MPI_Comm comm = MPI_COMM_WORLD;
    int size = lsb_data->recs.size();
    int max, min;
    MPI_Allreduce(&size,&max,1,MPI_INT,MPI_MAX,comm);
    MPI_Allreduce(&size,&min,1,MPI_INT,MPI_MIN,comm);
    if(max != min) {
      if(lsb_data->r == 0) printf("not the same number of records on all processes (min: %i, max: %i)! No output is generated!", min, max);
    } else {
      std::vector<long long> my_duration(size); // convert to long long for MPI :-(
      std::vector<long> my_overhead(size); // convert to long for MPI :-(
#ifdef HAVE_PAPI
      std::vector<long long> my_papictr1(size), my_papictr2(size); // convert to long long for MPI :-(
#endif
      for(int i=0; i<size; ++i) {
        my_duration[i] = lsb_data->recs[i].tduration;
        my_overhead[i] = lsb_data->recs[i].toverhead;
#ifdef HAVE_PAPI
        my_papictr1[i] = lsb_data->recs[i].papictr1;
        my_papictr2[i] = lsb_data->recs[i].papictr2;
#endif
      }

      std::vector<long long> max_duration(size), min_duration(size), avg_duration(size);
      std::vector<long> max_overhead(size), min_overhead(size), avg_overhead(size);
#ifdef HAVE_PAPI
      std::vector<long long> max_papictr1(size), max_papictr2(size), min_papictr1(size),
                             min_papictr2(size), avg_papictr1(size), avg_papictr2(size);
#endif
      MPI_Allreduce(&my_overhead[0],&min_overhead[0],size,MPI_LONG,MPI_MIN,comm);
      MPI_Allreduce(&my_overhead[0],&max_overhead[0],size,MPI_LONG,MPI_MAX,comm);
      MPI_Allreduce(&my_overhead[0],&avg_overhead[0],size,MPI_LONG,MPI_SUM,comm);
      for(int i=0; i<size; ++i) avg_overhead[i]/= lsb_data->p;

      MPI_Allreduce(&my_duration[0],&min_duration[0],size,MPI_LONG_LONG,MPI_MIN,comm);
      MPI_Allreduce(&my_duration[0],&max_duration[0],size,MPI_LONG_LONG,MPI_MAX,comm);
      MPI_Allreduce(&my_duration[0],&avg_duration[0],size,MPI_LONG_LONG,MPI_SUM,comm);
      for(int i=0; i<size; ++i) avg_duration[i]/= lsb_data->p;

#ifdef HAVE_PAPI
      MPI_Allreduce(&my_papictr1[0],&min_papictr1[0],size,MPI_LONG_LONG,MPI_MIN,comm);
      MPI_Allreduce(&my_papictr1[0],&max_papictr1[0],size,MPI_LONG_LONG,MPI_MAX,comm);
      MPI_Allreduce(&my_papictr1[0],&avg_papictr1[0],size,MPI_LONG_LONG,MPI_SUM,comm);
      for(int i=0; i<size; ++i) avg_papictr1[i]/= lsb_data->p;

      MPI_Allreduce(&my_papictr2[0],&min_papictr2[0],size,MPI_LONG_LONG,MPI_MIN,comm);
      MPI_Allreduce(&my_papictr2[0],&max_papictr2[0],size,MPI_LONG_LONG,MPI_MAX,comm);
      MPI_Allreduce(&my_papictr2[0],&avg_papictr2[0],size,MPI_LONG_LONG,MPI_SUM,comm);
      for(int i=0; i<size; ++i) avg_papictr2[i]/= lsb_data->p;
#endif

      if(lsb_data->write_file==1) {
        fprintf(fp, "# accumulated output format\n");
        fprintf(fp, "%s %s - %s - ", "id", "time (min, avg, max)", "overhead (min, avg, max)");
#ifdef HAVE_PAPI
        char p1[1024], p2[1024];
        PAPI_event_code_to_name(lsb_data->papi1, p1);
        PAPI_event_code_to_name(lsb_data->papi2, p2);
        fprintf(fp, "%s - %s", p1, p2);
#endif
#ifdef HAVE_UNWIND
        fprintf(fp, " #ips ip-list");
#endif
        fprintf(fp, "\n");
      }

      // print all records
      if(lsb_data->write_file==1) for(unsigned int i=0; i < size; ++i) {
        fprintf(fp, "%i ", lsb_data->recs[i].id);
        fprintf(fp, "%lli %lli %lli - ", min_duration[i], avg_duration[i], max_duration[i] );
        fprintf(fp, "%li %li %li - ", min_overhead[i], avg_overhead[i], max_overhead[i]);
#ifdef HAVE_PAPI
        fprintf(fp, "%lli %lli %lli - ", min_papictr1[i], avg_papictr1[i], max_papictr1[i]);
        fprintf(fp, "%lli %lli %lli - ", min_papictr2[i], avg_papictr2[i], max_papictr2[i]);
#endif
#ifdef HAVE_UNWIND
        fprintf(fp, "%i ", (int)lsb_data->iptrace[i].size());
        for(unsigned int j=0; j<lsb_data->iptrace[i].size(); ++j) fprintf(fp, "%lx ", lsb_data->iptrace[i][j]);
#endif
        fprintf(fp, "\n");
      }
    }
#else
    printf("MPI support required for the ACCUMULATED output format\n");
#endif

  }

  if(lsb_data->write_file==1) fprintf(fp, "# Runtime: %lf s (overhead: %lf %%) %i records\n", (double)(tend-lsb_data->tstart)/1e6, (double)dur/(tend-lsb_data->tstart)*100, (int)lsb_data->recs.size());

  
  /* Update rparams to the last values. It covers the corner case in which there are 
     parameters but not measures. */
  while (rparams_ptr < lsb_data->recparams.size()){
    //printf("advancing values\n");
    HASH_FIND_STR(ptable, lsb_data->recparams[rparams_ptr].pindex, current);
    rparams[current->value] = (lsb_data->recparams[rparams_ptr]);
    rparams_ptr++;
  }


  lsb_data->next=0;
  lsb_data->group_ptr=0;
  //printf("flush: start clear\n");
  lsb_data->recs.clear();
  lsb_data->recints.clear();
  lsb_data->recintdbl.clear();
  lsb_data->recparams.clear();
    
  /* Add last known value of the parameters */ 
  for (int i=0; i<pnum; i++){
    rparams[i].start_idx=0;
    lsb_data->recparams.push_back(rparams[i]);
  }
  /******/

  //printf("flush: end clear\n");
  free(rparams);

  lsb_data->write_header = 0;
#ifdef HAVE_UNWIND
  lsb_data->ipctrace.clear();
#endif

#if defined HAVE_PAPI
  //lsb_data->tstart = PAPI_get_real_usec();
#elif defined HAVE_LIKWID
    int cpu_id;
    cpu_id = likwid_getProcessorId();
    PMvalues[0] = msr_read(cpu_id, MSR_AMD15_PMC0);
    PMvalues[1] = msr_read(cpu_id, MSR_AMD15_PMC1);
    PMvalues[2] = msr_read(cpu_id, MSR_AMD15_PMC2);
    PMvalues[3] = msr_read(cpu_id, MSR_AMD15_PMC3);
    PMvalues[4] = msr_read(cpu_id, MSR_AMD15_PMC4);
    PMvalues[5] = msr_read(cpu_id, MSR_AMD15_PMC5);
    PMvalues[6] = msr_read(cpu_id, MSR_AMD15_NB_PMC0);
    PMvalues[7] = msr_read(cpu_id, MSR_AMD15_NB_PMC1);
    PMvalues[8] = msr_read(cpu_id, MSR_AMD15_NB_PMC2);
    PMvalues[9] = msr_read(cpu_id, MSR_AMD15_NB_PMC3);
    HRT_GET_TIMESTAMP(t);
    HRT_GET_TIME(t, ticks);
    //lsb_data->tstart =  HRT_GET_USEC(ticks);
#elif defined HAVE_HRTIMER
	HRT_GET_TIMESTAMP(t);
	HRT_GET_TIME(t, ticks);
	//lsb_data->tstart = HRT_GET_USEC(ticks);
#else
#error "No possibility to do accurate timing!\n"
#endif
#ifdef _OPENMP
}
#endif
}

/**
 * ends profiling
 * last function to be called in a profiling session
 */
static int LSB_finalized;
void LSB_Finalize() {
  CHK_DISABLED;
  if(LSB_finalized) return; // so that it can be called before atexit() if atexit doesn't work (e.g., BG/P)
#ifdef _OPENMP
#pragma omp master
{
#endif

  LSB_finalized = 1;

  if(lsb_data->print) printf("******* LSB_Finalize *******\n");

  if(lsb_data->recs.size() > 0 || lsb_data->recints.size() || lsb_data->recintdbl.size()) LSB_Flush();

  if(lsb_data->write_file==1) fclose(lsb_data->outfd);
#ifdef _OPENMP
}
#endif
#ifdef HAVE_LIKWID
  perfmon_stopCounters();
  msr_finalize();
#endif
}

/**
 * internal helper function to write benchmark information in header
 * @param fd file descriptor to write to
 */
static void write_host_information(FILE *fd) {
#ifdef _OPENMP
#pragma omp master
{
#endif
  struct utsname uninfo;

  if (uname (&uninfo) >= 0) {
    fprintf(fd, "# Sysname : %s\n", uninfo.sysname);
    fprintf(fd, "# Nodename: %s\n", uninfo.nodename);
    fprintf(fd, "# Release : %s\n", uninfo.release);
    fprintf(fd, "# Version : %s\n", uninfo.version);
    fprintf(fd, "# Machine : %s\n", uninfo.machine);
  }

  time_t t = time(NULL);
  struct tm *local = gmtime(&t);
  fprintf(fd, "# Execution time and date (UTC): %s", asctime(local));
  local = localtime(&t);
  fprintf(fd, "# Execution time and date (local): %s", asctime(local));
#ifdef _OPENMP
}
#endif

}


/**
 * first function to be called to open a tracing session
 * @param projname name of the project
 * @param autoprof_interval automatic (SIGALRM) profiling interval, set 0 to disable automatic profiling
 */
void LSB_Init(const char* projname, int autoprof_interval /* in ms, off if 0 */) {
  const char *fname = projname;
  char name[1024];
  int rev;
  struct stat st;
  char *env;

  lsb_data = new(LSB_Data);
  lsb_data->rec_enabled = 1;
  lsb_data->lsb_disabled = 0;
  lsb_data->print=1;
  lsb_data->next = 0;
  lsb_data->group_ptr=0;

  env = getenv("LSB_OUTFILE");
  if(env != NULL) fname = env;
  snprintf(name, 1023, "lsb.%s", fname);
  rev = 0;
  while(stat(name, &st) != -1) { // path exists!
    snprintf(name, 1023, "lsb.%s-%i", fname, ++rev);
  }



  int r;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &r);
  MPI_Comm_size(MPI_COMM_WORLD, &lsb_data->p);
#else
  r=0;
  lsb_data->p = 1;
#endif

  lsb_data->r=r;

  lsb_data->write_header=1;

  if(r > 0) lsb_data->print = 0;
  snprintf(name, 1023, "lsb.%s.r%i", fname, r);
  rev = 0;
  while(stat(name, &st) != -1) { // path exists!
    snprintf(name, 1023, "lsb.%s.r%i-%i", fname, r, ++rev);
  }

  lsb_data->write_file=1;
  env = getenv("LSB_OUTPUT_FORMAT"); // if accumulated is specified, then all processes but rank 0 are disabled!
  if((env != NULL) && (strstr(env, "accumulated") != NULL)) {
#ifndef HAVE_MPI
    printf("Compiled without MPI: accumulated format disabled (check LSB_OUTPUT_FORMAT envirorment variable)!\n");
    exit(-1);
#else
    lsb_data->write_file=0;
    if(lsb_data->r == 0) lsb_data->write_file=1;
#endif
  }

  env = getenv("LSB_ENABLE_PROCMASK"); // list of enabled processes
  if(env != NULL) {
    lsb_data->lsb_disabled=1;
    // parse comma-separated list of ints
    for( char *s2 = env; s2; ){
      while( *s2 == ' ' || *s2 == '\t' ) s2++;
      char *s1 = strsep( &s2, "," );
      if( !*s1 ){
        //printf("val: (empty)\n" );
      } else {
        int val;
        char ch;
        int ret = sscanf( s1, " %i %c", &val, &ch );
        if( ret != 1 ){
          //printf("val: (syntax error)\n" );
        }
        else{
          //printf("val: %i\n", val );
          if(r==val) lsb_data->lsb_disabled=0;
        }
      }
    }
  }

  CHK_DISABLED;

  //printf("printing enabled: %i\n", lsb_data->write_file);
  if(lsb_data->write_file==1) {
    lsb_data->outfd = fopen(name, "w");
    write_host_information(lsb_data->outfd);
  }

#ifdef HAVE_PAPI
  if(PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT ) {
    printf("PAPI initialization error! \n");
    exit(1);
  }
  lsb_data->eventset = PAPI_NULL;
  assert(PAPI_create_eventset(&lsb_data->eventset) == PAPI_OK);

  // define what to count! -- check "papi_event_chooser PRESET event1 event2"
  env = getenv("LSB_PAPI1");
  if(env != NULL) PAPI_event_name_to_code(env, &lsb_data->papi1);
  else lsb_data->papi1 = PAPI_TOT_INS;
  env = getenv("LSB_PAPI2");
  if(env != NULL) PAPI_event_name_to_code(env, &lsb_data->papi2);
  else lsb_data->papi2 = PAPI_TOT_CYC;

  if(PAPI_add_event(lsb_data->eventset, lsb_data->papi1) != PAPI_OK) {
    char p1[1024];
    PAPI_event_code_to_name(lsb_data->papi1, p1);
    if(lsb_data->print) printf("Adding PAPI counter 1 (%s) failed!\n", p1);
    exit(11);
  }
  if(PAPI_add_event(lsb_data->eventset, lsb_data->papi2) != PAPI_OK) {
    char p1[1024];
    PAPI_event_code_to_name(lsb_data->papi1, p1);
    if(lsb_data->print) printf("Adding PAPI counter 2 (%s) failed!\n", p1);
    exit(11);
  }
#endif

#ifdef HAVE_LIKWID
// init likwid using the low level API 
  FILE *OUTSTREAM = stdout;
  int numThreads=1;
  int threads[1];
  char *myeventstring;


// start likwid things
  cpuid_init();
  // printf("cpuid_topology.numHWThreads %d\n",cpuid_topology.numHWThreads);
  // numa_init(); // do we need this one ?
  // affinity_init(); // skip this, let aprun pin
  msr_init();

  threads[0]=likwid_getProcessorId();
  perfmon_init(numThreads,threads,OUTSTREAM);
  timer_init();

// Set the event set string
  myeventstring=(char *)malloc(60*sizeof(char));
  sprintf(myeventstring,"UNC_L3_CACHE_MISS_CORE_%d:UPMC%d",likwid_getProcessorId()%8,lsb_data->r%4);
  printf("LIKWID: rank %d event %s\n",lsb_data->r,myeventstring);
  bstring eventString = bfromcstr(myeventstring);
  perfmon_setupEventSet(eventString);
  perfmon_startCounters();
  free(myeventstring);
#endif

  if(lsb_data->print) {
    printf("***** LSB_Init >%s< writing to >%s< ", projname, name);
    if(autoprof_interval) printf("prof int: %i ms ", autoprof_interval);
    printf("*****\n");
  }

  if(autoprof_interval) {
    // register signal handler
    struct sigaction newact;
    newact.sa_handler = LSB_Handler;
    newact.sa_flags = 0;
    sigemptyset (&newact.sa_mask);

    sigaction(SIGVTALRM, &newact, NULL);

    // set timer
    struct itimerval timer;
    timer.it_interval.tv_sec=timer.it_value.tv_sec=0;
    timer.it_interval.tv_usec=timer.it_value.tv_usec=autoprof_interval*1000; // 1ms

    setitimer(ITIMER_VIRTUAL, &timer, NULL);
    
    // this is a workaround for autoprof malloc non-reentrant issues --
    // see USAGE for details
    env = getenv("LSB_AUTOPROF_WORKAROUND");
    if(env != NULL) {
      lsb_data->recs.reserve(atoi(env));
      if(lsb_data->print) printf("***** autoprof workaround - preallocating %i elements\n", atoi(env));
    }
  }

  LSB_finalized = 0;
  atexit(LSB_Finalize);

#if defined HAVE_PAPI
  // start PAPI
  PAPI_start(lsb_data->eventset);
  lsb_data->tlast = lsb_data->tstart = PAPI_get_real_usec();
#elif defined HAVE_LIKWID
    int cpu_id;
    cpu_id = likwid_getProcessorId();
    PMvalues[0] = msr_read(cpu_id, MSR_AMD15_PMC0);
    PMvalues[1] = msr_read(cpu_id, MSR_AMD15_PMC1);
    PMvalues[2] = msr_read(cpu_id, MSR_AMD15_PMC2);
    PMvalues[3] = msr_read(cpu_id, MSR_AMD15_PMC3);
    PMvalues[4] = msr_read(cpu_id, MSR_AMD15_PMC4);
    PMvalues[5] = msr_read(cpu_id, MSR_AMD15_PMC5);
    PMvalues[6] = msr_read(cpu_id, MSR_AMD15_NB_PMC0);
    PMvalues[7] = msr_read(cpu_id, MSR_AMD15_NB_PMC1);
    PMvalues[8] = msr_read(cpu_id, MSR_AMD15_NB_PMC2);
    PMvalues[9] = msr_read(cpu_id, MSR_AMD15_NB_PMC3);
    HRT_INIT(0, liblsb_g_timerfreq);
    unsigned long long ticks;
    HRT_TIMESTAMP_T t;
    HRT_GET_TIMESTAMP(t);
    HRT_GET_TIME(t, ticks);
    lsb_data->tstart = lsb_data->tlast =  HRT_GET_USEC(ticks);
#elif defined HAVE_HRTIMER
	HRT_INIT(0, liblsb_g_timerfreq);
  	unsigned long long ticks;
	HRT_TIMESTAMP_T t;
	HRT_GET_TIMESTAMP(t);
	HRT_GET_TIME(t, ticks);
	lsb_data->tlast = lsb_data->tstart = HRT_GET_USEC(ticks);
#else
#error "No possibility to do accurate timing!\n" 
#endif 
  
}

/**
 * registers a program parameter
 * @param format free float text (in printf-style)
 */
void LSB_Reg_param(const char *format, ...) {
  CHK_DISABLED;
#ifdef _OPENMP
#pragma omp master
{
#endif
  va_list val;
  char msg[1024];
  va_start(val, format);
  vsnprintf(msg, 1024, format, val);
  va_end(val);
  if(lsb_data->write_file) fprintf(lsb_data->outfd, "# Param: %s\n", msg);
#ifdef _OPENMP
}
#endif
}


/**
 * Set a record parameter
 * @id parameter identifier 
 * @val parameter new value
 */
void LSB_Set_Rparam_int(const char * id, int32_t val){
#ifdef _OPENMP
#pragma omp master
{
#endif
    t_recparam t;
    t.start_idx = lsb_data->next; //lsb_data->recs.size();
    t.pindex = id;
    t.type = LSB_INT;
    t.value = val;
    lsb_data->recparams.push_back(t);
    //lsb_data->recparams.push_back((t_recparam) {lsb_data->recs.size(), id, 0, {val}});
#ifdef _OPENMP
}
#endif
}

void LSB_Set_Rparam_long(const char * id, int64_t val){
#ifdef _OPENMP
#pragma omp master
{
#endif
    t_recparam t;
    t.start_idx = lsb_data->next; //lsb_data->recs.size();
    t.pindex = id;
    t.type = LSB_LNG;
    t.lng = val;
    lsb_data->recparams.push_back(t);
    //lsb_data->recparams.push_back((t_recparam) {lsb_data->recs.size(), id, 0, {val}});
#ifdef _OPENMP
}
#endif
}



void LSB_Set_Rparam_string(const char * id, const char * val){
#ifdef _OPENMP
#pragma omp master
{
#endif
    t_recparam t;
    //printf("setting %s to %s; starting from: %i\n", id, val, lsb_data->next);
    t.start_idx = lsb_data->next; //lsb_data->recs.size();
    t.pindex = id;
    t.type = LSB_STR;
    t.str = val;
    lsb_data->recparams.push_back(t);
    //lsb_data->recparams.push_back((t_recparam) {lsb_data->recs.size(), id, 1, {.str=val}});
#ifdef _OPENMP
}
#endif
}

void LSB_Set_Rparam_double(const char * id, double val){
#ifdef _OPENMP
#pragma omp master
{
#endif
    t_recparam t;
    t.start_idx = lsb_data->next; //lsb_data->recs.size();
    t.pindex = id;
    t.type = LSB_DBL;
    t.dbl = val;
    lsb_data->recparams.push_back(t);
    //lsb_data->recparams.push_back((t_recparam) {lsb_data->recs.size(), id, 1, {.str=val}});
#ifdef _OPENMP
}
#endif
}

/**
 * registers an identifier for output in the trace
 * @param format free float text (in printf-style)
 */
void LSB_Reg_id(const char *format, ...) {
  CHK_DISABLED;
#ifdef _OPENMP
#pragma omp master
{
#endif

  va_list val;
  char msg[1024];
  va_start(val, format);
  vsnprintf(msg, 1024, format, val);
  va_end(val);
  if(lsb_data->write_file) fprintf(lsb_data->outfd, "# Id: %s\n", msg);

#ifdef _OPENMP
}
#endif
}




void LSB_Rec_disable() {
  CHK_DISABLED;
#ifdef _OPENMP
#pragma omp master
{
#endif
  lsb_data->rec_enabled = 0;
#ifdef _OPENMP
}
#endif
}

void LSB_Rec_enable() {
  CHK_DISABLED;
#ifdef _OPENMP
#pragma omp master
{
#endif
  lsb_data->rec_enabled = 1;
#ifdef _OPENMP
}
#endif
}

