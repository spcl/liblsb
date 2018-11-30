
#include "../liblsb.h"
#include "../config.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define N 100

int main(int argc, char * argv[]){

    srand(time(NULL));

    MPI_Init(&argc, &argv);
    LSB_Init("test_group_fold", 0);

    LSB_Res();
    LSB_Wait(100);
    LSB_Rec(0);


    LSB_Group_t mgroup; //group of measurments

    LSB_Group_Begin(&mgroup);
    for (int i=0; i<N; i++){
        LSB_Res();
        LSB_Wait((rand() % 9) + 1);
        LSB_Rec(i+1);
    }
    LSB_Group_End(&mgroup);

    double min, max, median;
    //one can have multiple measurements with the same ID, so you could
    //fold all the measurements with a given ID in the specified group.
    LSB_Group_Fold(mgroup, LSB_ANY, LSB_MIN, &min); 
    LSB_Group_Fold(mgroup, LSB_ANY, LSB_MAX, &max); 
    LSB_Group_Fold(mgroup, LSB_ANY, LSB_MEDIAN, &median); 
    
    //these statistics shouldn't be influenced by the initial LSB_Wait
    printf("max: %lf; min: %lf; median: %lf\n", max, min, median);

    LSB_Finalize();
    MPI_Finalize();

    return 0;
}
