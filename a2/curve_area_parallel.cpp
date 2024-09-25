// CMPT431 Assignment 2 Siu Yu Chau 301604828 2024-9-24

#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>

#define sqr(x) ((x) * (x))
#define DEFAULT_NUMBER_OF_THREADS "4"
#define DEFAULT_NUMBER_OF_POINTS "1000000000"
#define DEFAULT_A "2"
#define DEFAULT_B "1"
#define DEFAULT_RANDOM_SEED "1"

uint c_const = (uint)RAND_MAX + (uint)1;
inline double get_random_coordinate(uint *random_seed) {
  return ((double)rand_r(random_seed)) / c_const;  // thread-safe random number generator
}

void get_points_in_curve(unsigned long n, uint random_seed,
  std::vector<unsigned long>& curve_points, unsigned int thread_id,
  std::vector<double>& jobs_time_taken, float a, float b) {

  timer job_timer;

  // Start job
  job_timer.start();
  unsigned long curve_count = 0;
  double x_coord, y_coord;
  for (unsigned long i = 0; i < n; i++) {
    x_coord = ((2.0 * get_random_coordinate(&random_seed)) - 1.0);
    y_coord = ((2.0 * get_random_coordinate(&random_seed)) - 1.0);
    if ((a * sqr(x_coord) + b * sqr(sqr(y_coord))) <= 1.0)
      curve_count++;
  }

  // End job
  jobs_time_taken[thread_id] = job_timer.stop();
  curve_points[thread_id] = curve_count;
}

void curve_area_calculation_parallel(unsigned int n_threads, unsigned long n_points, float a, float b, uint r_seed) {
  timer overall_timer;
  double overall_time_taken = 0.0;
  uint random_seed = r_seed;
  std::vector<std::thread> threads;
  std::vector<unsigned long> curve_points(n_threads, 0);
  std::vector<double> jobs_time_taken(n_threads, 0.0);
  

  overall_timer.start();
  
  // Calculate points per thread, if not divisible first thread take extra
  unsigned long itPerThreads = n_points / n_threads;
  unsigned long itInFirstThread = itPerThreads + n_points % n_threads;

  // Assign job to threads
  threads.push_back(std::thread(get_points_in_curve, itInFirstThread, r_seed,
    std::ref(curve_points), 0, std::ref(jobs_time_taken), a, b));
  for (unsigned int i = 1; i < n_threads; i++) {
    threads.push_back(std::thread(get_points_in_curve, itInFirstThread, r_seed + i, // different seed
      std::ref(curve_points), i, std::ref(jobs_time_taken), a, b));
  }

  // Join threads
  for (auto& t : threads) {
    t.join();
  }  

  // Calculate area
  unsigned long total_curve_points = 0;
  for (const unsigned long& p : curve_points) {
    total_curve_points += p;
  }
  double area_value = 4.0 * (double)total_curve_points / (double)n_points;

  //*------------------------------------------------------------------------
  overall_time_taken = overall_timer.stop();

  std::cout << "thread_id, points_generated, curve_points, time_taken\n";
  for (int i = 0; i < n_threads; ++i) {
    unsigned long n_it = (i == 0) ? itInFirstThread : itPerThreads;
    std::cout << i << ", " << n_it << ", "
                << curve_points[i] << ", " << std::setprecision(TIME_PRECISION)
                << jobs_time_taken[i] << "\n";
  }

  std::cout << "Total points generated : " << n_points << "\n";
  std::cout << "Total points in curve : " << total_curve_points << "\n";
  std::cout << "Area : " << std::setprecision(VAL_PRECISION) << area_value
            << "\n";
  std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
            << overall_time_taken << "\n";
}

int main(int argc, char *argv[]) {
  // Initialize command line arguments
  cxxopts::Options options("Curve_area_calculation",
                           "Calculate area inside curve a x^2 + b y ^4 = 1 using serial and parallel execution");
  options.add_options(
      "custom",
      {
          {"nThreads", "Number of threads",         
           cxxopts::value<unsigned int>()->default_value(DEFAULT_NUMBER_OF_THREADS)},
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
  unsigned int n_threads = cl_options["nThreads"].as<unsigned int>();
  unsigned long n_points = cl_options["nPoints"].as<unsigned long>();
  float a = cl_options["coeffA"].as<float>();
  float b = cl_options["coeffB"].as<float>();
  uint r_seed = cl_options["rSeed"].as<uint>();
  std::cout << "Number of points : " << n_points << "\n";;
  std::cout << "A : " << a << "\n" << "B : " << b << "\n";
  std::cout << "Random Seed : " << r_seed << "\n";

  curve_area_calculation_parallel(n_threads, n_points, a, b, r_seed);
  return 0;
}
