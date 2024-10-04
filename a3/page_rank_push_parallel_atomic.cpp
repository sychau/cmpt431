#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <atomic>

#ifdef USE_INT
#define INIT_PAGE_RANK 100000
#define EPSILON 1000
#define PAGE_RANK(x) (15000 + (5 * x) / 6)
#define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
typedef int64_t PageRankType;
typedef std::atomic<int64_t> PageRankTypeAtomic;
#else
#define INIT_PAGE_RANK 1.0
#define EPSILON 0.01
#define DAMPING 0.85
#define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
#define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
typedef double PageRankType;
typedef std::atomic<double> PageRankTypeAtomic;
#endif

// Include start_index but exclude end_index
void pageRankAlgo(
	Graph &g, std::vector<double> &threads_time_taken, CustomBarrier& barrier,
	PageRankTypeAtomic *pr_curr, PageRankTypeAtomic *pr_next,
	int max_iters, int thread_id, int start_index, int end_index
) {
	timer timer;
	timer.start();

  for (int iter = 0; iter < max_iters; iter++) {
    // for each vertex 'u', process all its outNeighbors 'v'
    for (uintV u = start_index; u < end_index; u++) {
			// u here is only used by this thread
      uintE out_degree = g.vertices_[u].getOutDegree();
      for (uintE i = 0; i < out_degree; i++) {
				// v here could be an index used by other thread in the same time
        uintV v = g.vertices_[u].getOutNeighbor(i);
        // Only one thread should read and modify pr_next[v] at a time

  			PageRankType curr = pr_next[v].load(std::memory_order_relaxed);
				PageRankType newWeight = pr_curr[u] / (PageRankType) out_degree;
				while(!pr_next[v].compare_exchange_weak(curr, curr + newWeight,
                                        	std::memory_order_release,
                                          std::memory_order_relaxed));   
      }
    }
    barrier.wait();
    for (uintV v = start_index; v < end_index; v++) {
			// v here is only used by this thread
      pr_next[v] = PAGE_RANK(pr_next[v]);
      // reset pr_curr for the next iteration
      pr_curr[v] = pr_next[v].load();
      pr_next[v] = 0.0;
    }
    barrier.wait();
  }

	threads_time_taken[thread_id] = timer.stop();
}

void pageRankParallel(Graph &g, int max_iters, int n_threads) {
  uintV n = g.n_;

  PageRankTypeAtomic *pr_curr = new PageRankTypeAtomic[n];
  PageRankTypeAtomic *pr_next = new PageRankTypeAtomic[n];

  for (uintV i = 0; i < n; i++) {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = 0.0;
  }

  // Push based pagerank
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
    sum_of_page_ranks += pr_curr[u].load();
  }
  std::cout << "Sum of page ranks : " << sum_of_page_ranks << "\n";
  std::cout << "Time taken (in seconds) : " << time_taken << "\n";
  delete[] pr_curr;
  delete[] pr_next;
}

int main(int argc, char *argv[]) {
  cxxopts::Options options(
      "page_rank_push",
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
  std::cout << "Using INT" << std::endl;
#else
  std::cout << "Using DOUBLE" << std::endl;
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
