#include <stdio.h>
#include <mpi.h>
#include "liblsb.h"
#include <unistd.h>
#include <stdlib.h>

#define N 4000
#define MSGSIZE 8192
#define SAMPLES 4000


#define BINSIZE 100

#define WAIT_TIME 500.0
#define EXP_TIME 10.0



static volatile int dummyvar=0;

//helper function: waits for usec usec
double bwait(double usec){   
    double curr;
    double start = MPI_Wtime();
    double end = MPI_Wtime() + usec/1e6;
    //printf("%lf %lf\n", MPI_Wtime(), end);
    while ((curr=MPI_Wtime()) < end){ dummyvar++; }
    return curr - start;
}

// the function to benchmark
inline void tobench(void * sndbuff, void * rcvbuff, int size){
    MPI_Bcast(sndbuff, size, MPI_CHAR, 0, MPI_COMM_WORLD);
}

int main(int argc, char * argv[]){
    
    
    int rank, i;
    double err;
    double win;
 
    /* init */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    LSB_Init("t1", 0);

    /* set output params */
    LSB_Set_Rparam_int("rank", rank);
    LSB_Set_Rparam_double("err", 0); // meaningless here
    LSB_Set_Rparam_string("type", "WARMUP");

    /* allocate buff */
    void * sndbuff = malloc(sizeof(char)*MSGSIZE);


    /* Warmup */
    for (int i=0; i<SAMPLES; i++){
        tobench(sndbuff, NULL, MSGSIZE);
    }


    /* First step: determine local win size */
    for (int i=0; i<SAMPLES; i++){
        LSB_Res();
        tobench(sndbuff, NULL, MSGSIZE);
        LSB_Rec(0);
    }

    /* Compute the maximum measurement taken until now with id==0 */
    LSB_Fold(0, LSB_MAX, &win);

    
    
    /* Bench 1: sync with liblsb */
    LSB_Set_Rparam_string("type", "LSB");

    /* init sync */
    LSB_Sync_init(MPI_COMM_WORLD, win*4);
    printf("win: %lf\n", win);   
    
    /* Take N measurements */
    for (int i=0; i<N; i++){
        /* Sync */
        err=LSB_Sync();

        /* Print if the synchronization was successful */
        LSB_Set_Rparam_double("err", err);

        /* Take the time */
        LSB_Res();
        tobench(sndbuff, NULL, MSGSIZE);
        LSB_Rec(i/BINSIZE);

    }

    /* Bench 2: sync with MPI_BARRIER */
    LSB_Set_Rparam_string("type", "MPI");
    LSB_Set_Rparam_double("err", 0); //meaningless from now on
    
    for (int i=0; i<N; i++){
        /* sYNC */
        MPI_Barrier(MPI_COMM_WORLD);

        /* Take the time */
        LSB_Res();
        tobench(sndbuff, NULL, MSGSIZE);
        LSB_Rec(i/BINSIZE);
    }

    /* Finalize */
    LSB_Finalize();
    MPI_Finalize();
    free(sndbuff);
}

