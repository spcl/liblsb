/*
* Copyright (c) 2015 ETH-Zurich. All rights reserved.
* Use of this source code is governed by a BSD-style license that can be
* found in the LICENSE file.
*/

#include <stdint.h>

#define HAVE_MPI_H

#ifndef LIBLSB_H_
#define LIBLSB_H_
#define SYNC_WINDOW

#define LSB_ANY -1

#ifdef HAVE_MPI_H
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef enum { LSB_SUM=0, LSB_COUNT, LSB_MEDIAN, LSB_MAX, LSB_MIN} lsb_op_t;

struct LSB_Group{
    int start, stop;
};

typedef struct LSB_Group LSB_Group_t;

void LSB_Init(const char *projname, int autoprof_interval);
void LSB_Finalize();
void LSB_Reg_param(const char *format, ...);
void LSB_Set_Rparam_int(const char * id, int val);
void LSB_Set_Rparam_long(const char * id, int64_t val);
void LSB_Set_Rparam_string(const char * id, const char * val);
void LSB_Set_Rparam_double(const char * id, double val);
void LSB_Reg_id(const char *format, ...);
double LSB_Rec(unsigned int id);
double LSB_Check(unsigned int id);
double LSB_Stop(unsigned int id, unsigned int reset);
void LSB_Next();
//void LSB_Rec_ints(unsigned int id, int int1, int int2);
//void LSB_Rec_intdbl(unsigned int id, int int1, double dbl);
void LSB_Flush();
void LSB_Res();
void LSB_Rec_enable();
void LSB_Rec_disable();
void LSB_Fold(unsigned int id, lsb_op_t op, double * result);
double LSB_Wait(double microseconds);


#ifdef HAVE_MPI_H

void LSB_Sync_init(MPI_Comm comm, double window);
void LSB_Sync_reset(double window);
double LSB_Sync();

#endif



#ifdef __cplusplus
}
#endif
#endif /* LIBLSB_H_ */
