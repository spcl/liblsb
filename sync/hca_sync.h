/*  ReproMPI Benchmark
 *
 *  Copyright 2015 Alexandra Carpen-Amarie, Sascha Hunold
    Research Group for Parallel Computing
    Faculty of Informatics
    Vienna University of Technology, Austria

<license>
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
</license>
*/

#ifndef HCA_SYNC_H_
#define HCA_SYNC_H_

typedef struct {
    long n_rep; /* --repetitions */
    double window_size_sec; /* --window-size */

    int n_fitpoints; /* --fitpoints */
    int n_exchanges; /* --exchanges */
    MPI_Comm comm;

} hca_options_t;


int hca_init_synchronization_module(int argc, char* argv[]);
int hca_init_synchronization_module_def(MPI_Comm comm, double window);
void hca_set_local_window(double window);
void hca_synchronize_clocks(void);

double hca_start_synchronization(void);
double hca_stop_synchronization(void);
void hca_cleanup_synchronization_module(void);

int* hca_get_local_sync_errorcodes(void);

double hca_get_normalized_time(double local_time);
double hca_get_adjusted_time(void);

void hca_print_sync_parameters(void);

int my_pow_2(int exp);

#endif /* HCA_SYNC_H_ */
