#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <mpi.h>

#define sqr(x) ((x) * (x))
#define DEFAULT_NUMBER_OF_POINTS "1000000000"
#define DEFAULT_A "2"
#define DEFAULT_B "1"
#define DEFAULT_RANDOM_SEED "1"
#define ROOT 0

uint c_const = (uint)RAND_MAX + (uint)1;
inline double get_random_coordinate(uint *random_seed) {
  return ((double)rand_r(random_seed)) / c_const;  // thread-safe random number generator
}

struct GetPointsResult {    
	unsigned long curve_point;
	double time_taken;
};

GetPointsResult get_points_in_curve(unsigned long n, uint random_seed, float a, float b) {
	timer timer;
	timer.start();
  unsigned long curve_count = 0;
  double x_coord, y_coord;
  for (unsigned long i = 0; i < n; i++) {
    x_coord = ((2.0 * get_random_coordinate(&random_seed)) - 1.0);
    y_coord = ((2.0 * get_random_coordinate(&random_seed)) - 1.0);
    if ((a*sqr(x_coord) + b*sqr(sqr(y_coord))) <= 1.0)
      curve_count++;
  }
	double time_taken = timer.stop();
	GetPointsResult res = {curve_count, time_taken};
  return res;
}

void curve_area_calculation_parallel(unsigned long n, float a, float b, uint r_seed) {
  timer overall_timer;
  double time_taken = 0.0;
  uint random_seed = r_seed;

  overall_timer.start();

	// Dividing up n vertices on P processes. 
	// Total number of processes is world_size. This process rank is world_rank
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	unsigned long min_points_per_process = n / world_size;
	unsigned long excess_points = n % world_size;
	unsigned long points_to_be_generated;
	if (world_rank < excess_points) {
		points_to_be_generated = min_points_per_process + 1;
	}	else {
		points_to_be_generated = min_points_per_process;
	}
	// Each process will work on points_to_be_generated and estimate curve_points.

	GetPointsResult res = get_points_in_curve(points_to_be_generated, random_seed + world_rank, a, b);
	unsigned long local_count = res.curve_point;
	double local_time_taken = res.time_taken;

	unsigned long* global_curve_points = nullptr;
	double* global_time_taken = nullptr;
	if (world_rank == ROOT) {
		global_curve_points = static_cast<unsigned long*>(malloc(sizeof(unsigned long) * world_size));
		global_time_taken = static_cast<double*>(malloc(sizeof(double) * world_size));
	}

	MPI_Gather(&local_count, 1, MPI_UNSIGNED_LONG, global_curve_points, 1, MPI_UNSIGNED_LONG, ROOT, MPI_COMM_WORLD);
	MPI_Gather(&local_time_taken, 1, MPI_DOUBLE, global_time_taken, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

  //*------------------------------------------------------------------------
  time_taken = overall_timer.stop();

	if (world_rank != ROOT) return;

	std::cout << "rank, points_generated, curve_points, time_taken\n";
	for (int i = 0; i < world_size; i++) {
		std::cout << i << ", "
							<< points_to_be_generated << ", "
							<< global_curve_points[i] << ", " << std::setprecision(TIME_PRECISION)
							<< global_time_taken[i] << "\n";
	}

	std::cout << "Total points generated : " << points_to_be_generated * world_size << "\n";
	unsigned long total_points_count = 0;
	for(int i = 0; i < world_size; i++) {
		total_points_count += global_curve_points[i];
	}
	std::cout << "Total points in curve : " << total_points_count << "\n";
	std::cout << "Area : " << std::setprecision(VAL_PRECISION) << 4.0 * (double) total_points_count / (double) n << "\n";
	std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
	
	free(global_time_taken);
	free(global_curve_points);
}

int main(int argc, char *argv[]) {
	MPI_Init(NULL, NULL);

  // Initialize command line arguments
  cxxopts::Options options("Curve_area_calculation",
                           "Calculate area inside curve a x^2 + b y ^4 = 1 using serial and parallel execution");
  options.add_options(
      "custom",
      {
          {"nPoints", "Number of points",         
           cxxopts::value<unsigned long>()->default_value(DEFAULT_NUMBER_OF_POINTS)},
	        {"coeffA", "Coefficient a",
	         cxxopts::value<float>()->default_value(DEFAULT_A)},
          {"coeffB", "Coefficient b",
           cxxopts::value<float>()->default_value(DEFAULT_B)},
          {"rSeed", "Random Seed",
           cxxopts::value<uint>()->default_value(DEFAULT_RANDOM_SEED)}
      });
  auto cl_options = options.parse(argc, argv);
  unsigned long n_points = cl_options["nPoints"].as<unsigned long>();
  float a = cl_options["coeffA"].as<float>();
  float b = cl_options["coeffB"].as<float>();
  uint r_seed = cl_options["rSeed"].as<uint>();

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	if (world_rank == ROOT) {
		std::cout << "Number of processes : " << world_size << "\n";
		std::cout << "Number of points : " << n_points << "\n";
		std::cout << "A : " << a << "\n" << "B : " << b << "\n";
		std::cout << "Random Seed : " << r_seed << "\n";
	}
  curve_area_calculation_parallel(n_points, a, b, r_seed);
	
	MPI_Finalize();
  return 0;
}
