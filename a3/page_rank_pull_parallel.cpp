#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>

#ifdef USE_INT
#define INIT_PAGE_RANK 100000
#define EPSILON 1000
#define PAGE_RANK(x) (15000 + (5 * x) / 6)
#define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
typedef int64_t PageRankType;
#else
#define INIT_PAGE_RANK 1.0
#define EPSILON 0.01
#define DAMPING 0.85
#define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
#define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
typedef double PageRankType;
#endif

// Include start_index but exclude end_index
void pageRankAlgo(
	Graph &g, std::vector<double> &threads_time_taken, CustomBarrier& barrier,
	PageRankType *pr_curr, PageRankType *pr_next,
	int max_iters, int thread_id, int start_index, int end_index
) {
	timer timer;
	timer.start();

	for (int iter = 0; iter < max_iters; iter++) {
		// for each vertex 'v', process all its inNeighbors 'u'
		for (uintV v = start_index; v < end_index; v++) {
			uintE in_degree = g.vertices_[v].getInDegree();
			for (uintE i = 0; i < in_degree; i++) {
				uintV u = g.vertices_[v].getInNeighbor(i);
				uintE u_out_degree = g.vertices_[u].getOutDegree();
				if (u_out_degree > 0) {
					pr_next[v] += (pr_curr[u] / (PageRankType) u_out_degree);
				}
			}
		}
		barrier.wait();
		for (uintV v = start_index; v < end_index; v++) {
			pr_next[v] = PAGE_RANK(pr_next[v]);
			// reset pr_curr for the next iteration
			pr_curr[v] = pr_next[v];
			pr_next[v] = 0.0;
		}
		barrier.wait();
	}
	threads_time_taken[thread_id] = timer.stop();
}

void pageRankParallel(Graph &g, int max_iters, int n_threads) {
  uintV n = g.n_;

  PageRankType *pr_curr = new PageRankType[n];
  PageRankType *pr_next = new PageRankType[n];

  for (uintV i = 0; i < n; i++) {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = 0.0;
  }

  // Pull based pagerank
  timer t1;
  double time_taken = 0.0;
  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  t1.start();

	std::vector<std::thread> threads;
	std::vector<double> threads_time_taken(n_threads, 0.0);

  // Calculate vertices per thread, if not divisible first thread take extra
  unsigned long vertices_per_thread = n / n_threads;
  unsigned long vertices_in_first_thread = vertices_per_thread + n % n_threads;

  // Assign job to threads
	CustomBarrier barrier{n_threads};

  threads.emplace_back(pageRankAlgo,
		std::ref(g), std::ref(threads_time_taken), std::ref(barrier), 
		pr_curr, pr_next,
		max_iters, 0, 0, vertices_in_first_thread
	);

  for (int i = 1; i < n_threads; i++) {
		int start_index = vertices_in_first_thread + (i - 1) * vertices_per_thread;
		int end_index = vertices_in_first_thread + i * vertices_per_thread;

    threads.emplace_back(pageRankAlgo,
			std::ref(g), std::ref(threads_time_taken), std::ref(barrier), 
			pr_curr, pr_next,
			max_iters, i, start_index, end_index
		);
  }

  // Join threads
  for (auto& t : threads) {
    t.join();
  }  

  time_taken = t1.stop();
  // -------------------------------------------------------------------
  std::cout << "thread_id, time_taken" << std::endl;
  // Print the above statistics for each thread
  // Example output for 2 threads:
  // thread_id, time_taken
  // 0, 0.12
  // 1, 0.12

	// Let thread 0 be main thread and created thread be 1 to n_threads
	for (int i = 0; i < n_threads; ++i) {
		std::cout << i << ", " << std::setprecision(6) << threads_time_taken[i] << "\n";
	}

  PageRankType sum_of_page_ranks = 0;
  for (uintV u = 0; u < n; u++) {
    sum_of_page_ranks += pr_curr[u];
  }
  std::cout << "Sum of page ranks : " << sum_of_page_ranks << "\n";
  std::cout << "Time taken (in seconds) : " << time_taken << "\n";
  delete[] pr_curr;
  delete[] pr_next;
}

int main(int argc, char *argv[]) {
  cxxopts::Options options(
      "page_rank_pull",
      "Calculate page_rank using serial and parallel execution");
  options.add_options(
      "",
      {
          {"nThreads", "Number of Threads",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_THREADS)},
          {"nIterations", "Maximum number of iterations",
           cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_threads = cl_options["nThreads"].as<uint>();
  uint max_iterations = cl_options["nIterations"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();

#ifdef USE_INT
  std::cout << "Using INT\n";
#else
  std::cout << "Using DOUBLE\n";
#endif
  std::cout << std::fixed;
  std::cout << "Number of Threads : " << n_threads << std::endl;
  std::cout << "Number of Iterations: " << max_iterations << std::endl;

  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";
  pageRankParallel(g, max_iterations, n_threads);

  return 0;
}
