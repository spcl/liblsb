/* Copyright (c) 2009 The Trustees of Indiana University and Indiana
 *                    University Research and Technology
 *                    Corporation.  All rights reserved.
 *
 * Author(s): Torsten Hoefler <htor@inf.ethz.ch>
 *            Salvatore Di Girolamo <digirols@inf.ethz.ch>
 */

#include <stdio.h>
#include <stdlib.h>
#include "liblsb_internal.hpp"
#include "config.h"
#include "sync/hca_sync.h"

#define MAX_DOUBLE 1e99


#if defined(HAVE_SYNC) && defined(HAVE_MPI)

extern LSB_Data *lsb_data;



void LSB_Sync_init(MPI_Comm comm, double window){
   lsb_data->sync_comm = comm;
   hca_init_synchronization_module_def(comm, window/1e6);
   LSB_Sync_reset(window);
   //hca_synchronize_clocks(); 
}

void LSB_Sync_reset(double window){
   double nw;
   MPI_Allreduce(&window, &nw, 1, MPI_DOUBLE, MPI_MAX, lsb_data->sync_comm);
   hca_set_local_window(nw/1e6); //hca wants the window in secs

   /*TODO: here the full synchronization can be avoided. 
     We just need to broadcast the new synctime. However, resyncrhonizing 
     should be safer since it relies less on the model. */
   hca_synchronize_clocks(); 
}




/* for every single measurement */
double LSB_Sync() {
  return hca_start_synchronization()*1e6;
}
#else

void LSB_Sync_init(MPI_Comm comm, double window) { printf("Error: This feature needs MPI support\n"); exit(-1); }
void LSB_Sync_reset(double window) { printf("Error: This feature needs MPI support\n"); exit(-1); } 
double LSB_Sync() { printf("Error: This feature needs MPI support\n"); exit(-1); } 
 

#endif
