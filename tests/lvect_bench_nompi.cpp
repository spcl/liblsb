
#include <liblsb.h>


#include <stdio.h>
#include <vector>

//#include <mpi.h>
#include "../liblsb_lvector.hpp"

#define RUNS 100


/* Benchmark for lvector */
int bench_lvector(int n){
    lvector<int> vec;
    int acc=0; 

    /* Set the size of the current test test */
    LSB_Set_Rparam_int("size", n);

    /* Start the lvector_write region */
    LSB_Set_Rparam_string("region", "lvector_write");
    
    for (int i=0; i<RUNS; i++){
        /* Reset */
        LSB_Res();

        /* Push n integers */
        for (int j=0; j<n; j++) {
            vec.push_back(j);
        }

        /* Record and complete the current measure */
        LSB_Rec(i);
    }

    /* Start the lvector_read region */    
    LSB_Set_Rparam_string("region", "lvector_read");

  
    for (int i=0; i<RUNS; i++){
        /* Reset */
        LSB_Res();

        /* Read n integers */
        for(int j=0; j<n; j++) {
            acc += vec[j];
            asm("");
        }

        /* Register and complete the current measure */
        LSB_Rec(i);
    }
}


/* Benchmark for std::vector */
int bench_vector(int n){
    std::vector<int> vec;
    int acc=0; 

    /* Update the "size" parameter */
    LSB_Set_Rparam_int("size", n);

    /* The read section for std::vector is starting */
    LSB_Set_Rparam_string("region", "vector_write");
    

    for (int i=0; i<RUNS; i++){
        /* Reset */
        LSB_Res();

        /* Insert n integers */
        for (int j=0; j<n; j++) {
            vec.push_back(j);
        }
        
        /* Record and terminate the current measure */
        LSB_Rec(i);
    }

    /* Start of the write section for std::vector */
    LSB_Set_Rparam_string("region", "vector_read");

    for (int i=0; i<RUNS; i++){
        /* Reset */
        LSB_Res();

        /* Read n integers */
        for(int j=0; j<n; j++) {
            acc += vec[j];
            asm("");
        }

        /* Record and terminate the current measure */
        LSB_Rec(i);
    }
}



int main(){

    //MPI_Init(NULL, NULL);

    /* Initialize lsblib */
    LSB_Init("test_lvect", 0);
      

    /* Perform the measurements */
    for (int i=1000; i<5000; i += 50) bench_lvector(i);
    for (int i=1000; i<5000; i += 50) bench_vector(i);

    /* Finalize lsblib */
    LSB_Finalize();

    //MPI_Finalize();

    return 0;
}

