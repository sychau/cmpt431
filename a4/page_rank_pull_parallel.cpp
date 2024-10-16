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
void pageRankVertexBasedDecompStatic(
	Graph &g, CustomBarrier& barrier,
  std::vector<uint> &threads_vertices, std::vector<uint> &threads_edges,
  std::vector<double> &threads_barrier1_time, std::vector<double> &threads_barrier2_time,
  std::vector<double> &threads_get_next_vertex_time, std::vector<double> &threads_total_time,
	PageRankType *pr_curr, PageRankType *pr_next,
	int max_iters, int thread_id, int start_index, int end_index
) {
	timer total_timer;
  timer barrier1_timer;
  timer barrier2_timer;

  double barrier1_total_time = 0.0;
  double barrier2_total_time = 0.0;
  uintV vertices_processed = 0;
  uintE edges_processed = 0;

	total_timer.start();

	for (int iter = 0; iter < max_iters; iter++) {
		// for each vertex 'v', process all its inNeighbors 'u'
    vertices_processed += end_index - start_index;
		for (uintV v = start_index; v < end_index; v++) {
			uintE in_degree = g.vertices_[v].getInDegree();
      edges_processed += in_degree;
			for (uintE i = 0; i < in_degree; i++) {
				uintV u = g.vertices_[v].getInNeighbor(i);
				uintE u_out_degree = g.vertices_[u].getOutDegree();
				if (u_out_degree > 0) {
					pr_next[v] += (pr_curr[u] / (PageRankType) u_out_degree);
				}
			}
		}
    barrier1_timer.start();
		barrier.wait();
    barrier1_total_time += barrier1_timer.stop();

		for (uintV v = start_index; v < end_index; v++) {
			pr_next[v] = PAGE_RANK(pr_next[v]);
			// reset pr_curr for the next iteration
			pr_curr[v] = pr_next[v];
			pr_next[v] = 0.0;
		}
    barrier2_timer.start();
		barrier.wait();
    barrier2_total_time += barrier2_timer.stop();
	}
  threads_vertices[thread_id] = vertices_processed;
  threads_edges[thread_id] = edges_processed;
  threads_barrier1_time[thread_id] = barrier1_total_time;
  threads_barrier2_time[thread_id] = barrier2_total_time;
  threads_total_time[thread_id] = total_timer.stop();
}

void pageRankEdgeBasedDecomp(
	Graph &g, CustomBarrier& barrier,
  std::vector<uint> &threads_vertices, std::vector<uint> &threads_edges,
  std::vector<double> &threads_barrier1_time, std::vector<double> &threads_barrier2_time,
  std::vector<double> &threads_get_next_vertex_time, std::vector<double> &threads_total_time,
	PageRankType *pr_curr, PageRankType *pr_next,
	int max_iters, int thread_id, const std::vector<uint> &verticesGroup
) {
	timer total_timer;
  timer barrier1_timer;
  timer barrier2_timer;

  double barrier1_total_time = 0.0;
  double barrier2_total_time = 0.0;
  uintV vertices_processed = 0;
  uintE edges_processed = 0;

	total_timer.start();

  for (int iter = 0; iter < max_iters; iter++) {
		// for each vertex 'v', process all its inNeighbors 'u'
    vertices_processed += verticesGroup.size();
		for (const uint &v: verticesGroup) {
			uintE in_degree = g.vertices_[v].getInDegree();
      edges_processed += in_degree;
			for (uintE i = 0; i < in_degree; i++) {
				uintV u = g.vertices_[v].getInNeighbor(i);
				uintE u_out_degree = g.vertices_[u].getOutDegree();
				if (u_out_degree > 0) {
					pr_next[v] += (pr_curr[u] / (PageRankType) u_out_degree);
				}
			}
		}
    barrier1_timer.start();
		barrier.wait();
    barrier1_total_time += barrier1_timer.stop();

		for (const uint &v: verticesGroup) {
			pr_next[v] = PAGE_RANK(pr_next[v]);
			// reset pr_curr for the next iteration
			pr_curr[v] = pr_next[v];
			pr_next[v] = 0.0;
		}

    barrier2_timer.start();
		barrier.wait();
    barrier2_total_time += barrier2_timer.stop();
	}
  threads_vertices[thread_id] = vertices_processed;
  threads_edges[thread_id] = edges_processed;
  threads_barrier1_time[thread_id] = barrier1_total_time;
  threads_barrier2_time[thread_id] = barrier2_total_time;
  threads_total_time[thread_id] = total_timer.stop();
}

void pageRankVertexBasedDecompDynamic(
	Graph &g, CustomBarrier& barrier,
  std::vector<uint> &threads_vertices, std::vector<uint> &threads_edges,
  std::vector<double> &threads_barrier1_time, std::vector<double> &threads_barrier2_time,
  std::vector<double> &threads_get_next_vertex_time, std::vector<double> &threads_total_time,
	PageRankType *pr_curr, PageRankType *pr_next,
	int max_iters, int thread_id, int granularity, std::atomic<uint> &vertex_pool, uint n_vertices
) {
	timer total_timer;
  timer barrier1_timer;
  timer barrier2_timer;
  timer get_next_vertex_timer;

  double barrier1_total_time = 0.0;
  double barrier2_total_time = 0.0;
  double get_next_vertex_total_time = 0.0;
  uintV vertices_processed = 0;
  uintE edges_processed = 0;

	total_timer.start();

  std::vector<uint> v_list; // vertices taken by this thread
  for(int iter = 0; iter < max_iters; iter++) {
    while(true) {
      get_next_vertex_timer.start();
      uint curr_v_index = vertex_pool.load();
      while (!vertex_pool.compare_exchange_weak(curr_v_index, curr_v_index + granularity));
      uint v = curr_v_index - granularity;
      get_next_vertex_total_time += get_next_vertex_timer.stop();
      if (v >= n_vertices) break;
      v_list.push_back(v);

      vertices_processed += (v + granularity >= n_vertices) ? (n_vertices - v) : granularity;
      for (int gi = 0; gi < granularity; ++gi) {
        if (v + gi >= n_vertices) break;
        uintE in_degree = g.vertices_[v + gi].getInDegree();
        edges_processed += in_degree;
        for (uintE i = 0; i < in_degree; i++) {
          uintV u = g.vertices_[v + gi].getInNeighbor(i);
          uintE u_out_degree = g.vertices_[u].getOutDegree();
          if (u_out_degree > 0) {
            pr_next[v + gi] += (pr_curr[u] / (PageRankType) u_out_degree);
          }
        }
      }
    }
    barrier1_timer.start();
		barrier.wait();
    barrier1_total_time += barrier1_timer.stop();

    for (uint& v : v_list) {
      for (int gi = 0; gi < granularity; ++gi) {
        if (v + gi >= n_vertices) break;
        pr_next[v + gi] = PAGE_RANK(pr_next[v + gi]);
        // reset pr_curr for the next iteration
        pr_curr[v + gi] = pr_next[v + gi];
        pr_next[v + gi] = 0.0;
      }
    }
    if (thread_id == 0) vertex_pool.store(0);
    v_list.clear();
    barrier2_timer.start();
		barrier.wait();
    barrier2_total_time += barrier2_timer.stop();
  }
  threads_vertices[thread_id] = vertices_processed;
  threads_edges[thread_id] = edges_processed;
  threads_barrier1_time[thread_id] = barrier1_total_time;
  threads_barrier2_time[thread_id] = barrier2_total_time;
  threads_get_next_vertex_time[thread_id] = get_next_vertex_total_time;
  threads_total_time[thread_id] = total_timer.stop();
}

void pageRankParallel(Graph &g, int max_iters, int n_threads, int strategy, int granularity) {
  uintV n = g.n_;
  uintE m = g.m_;

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

  // Statistics
  std::vector<uint> threads_vertices(n_threads, 0); 
  std::vector<uint> threads_edges(n_threads, 0); 
  std::vector<double> threads_barrier1_time(n_threads, 0.0);
  std::vector<double> threads_barrier2_time(n_threads, 0.0);
  std::vector<double> threads_get_next_vertex_time(n_threads, 0.0);
  std::vector<double> threads_total_time(n_threads, 0.0);

  CustomBarrier barrier{n_threads};
  std::vector<std::vector<uint>> vertices_groups; // all groups will have similar edge count for strategy 2
  std::atomic<uint> vertex_pool(0); // distribute vertex dynamically for strategy 3 & 4

  switch (strategy) {
    // Vertex-based decomposition and static mapping, the same strategy as your Assignment 3 programs.
    // Vertices are statically distributed among threads such that each thread performs all computations on n/T vertices.
    case 1: {
      std::cout << "one\n";
      // Calculate vertices per thread, if not divisible first thread take extra
      unsigned long vertices_per_thread = n / n_threads;
      unsigned long vertices_in_first_thread = vertices_per_thread + n % n_threads;

      // Assign job to threads
      threads.emplace_back(pageRankVertexBasedDecompStatic,
          std::ref(g), std::ref(barrier),
          std::ref(threads_vertices), std::ref(threads_edges),
          std::ref(threads_barrier1_time), std::ref(threads_barrier2_time),
          std::ref(threads_get_next_vertex_time), std::ref(threads_total_time),
        pr_curr, pr_next,
        max_iters, 0, 0, vertices_in_first_thread
      );

      for (int thread_i = 1; thread_i < n_threads; thread_i++) {
        int start_index = vertices_in_first_thread + (thread_i - 1) * vertices_per_thread;
        int end_index = vertices_in_first_thread + thread_i * vertices_per_thread;

        threads.emplace_back(pageRankVertexBasedDecompStatic,
          std::ref(g), std::ref(barrier),
          std::ref(threads_vertices), std::ref(threads_edges),
          std::ref(threads_barrier1_time), std::ref(threads_barrier2_time),
          std::ref(threads_get_next_vertex_time), std::ref(threads_total_time),
          pr_curr, pr_next,
          max_iters, thread_i, start_index, end_index
        );
      }
      break;
    }
    case 2: {
      // Vertices are distributed among threads such that each thread performs computations on m/T edges.
      // Pull-based PageRank distributes vertices based on their in-degrees, while Push-based PageRank
      // distributes vertices based on their out-degrees
      vertices_groups.reserve(n_threads);
      int vertex_index = 0;
      int remaining_edges = m;
      for(int thread_i = 0; thread_i < n_threads; ++thread_i) {
        vertices_groups.push_back(std::vector<uint>());
        int target_edges = remaining_edges / (n_threads - thread_i);
        int current_edge_count = 0;
        while (current_edge_count < target_edges && vertex_index < n) {
          vertices_groups[thread_i].push_back(vertex_index);
          current_edge_count += g.vertices_[vertex_index].getInDegree();
          vertex_index++;
        }
        remaining_edges -= current_edge_count;
        
        threads.emplace_back(pageRankEdgeBasedDecomp,
          std::ref(g), std::ref(barrier),
          std::ref(threads_vertices), std::ref(threads_edges),
          std::ref(threads_barrier1_time), std::ref(threads_barrier2_time),
          std::ref(threads_get_next_vertex_time), std::ref(threads_total_time),
          pr_curr, pr_next,
          max_iters, thread_i, std::cref(vertices_groups[thread_i])
        );
      }
      break;
    }
    case 3: case 4: {
      for(int thread_i = 0; thread_i < n_threads; ++thread_i) {
        threads.emplace_back(pageRankVertexBasedDecompDynamic,
          std::ref(g), std::ref(barrier),
          std::ref(threads_vertices), std::ref(threads_edges),
          std::ref(threads_barrier1_time), std::ref(threads_barrier2_time),
          std::ref(threads_get_next_vertex_time), std::ref(threads_total_time),
          pr_curr, pr_next,
          max_iters, thread_i, (strategy == 3) ? 1 : granularity, std::ref(vertex_pool), n
        );
      }
      break;
    }
    default: {
      std::cout << "strategy " << strategy << " invalid";
      break;
    }
  }

  // Join threads
  for (auto& t : threads) {
    t.join();
  }  

  time_taken = t1.stop();
  // -------------------------------------------------------------------
  std::cout << "thread_id, num_vertices, num_edges, barrier1_time, "
    "barrier2_time, getNextVertex_time, total_time" << std::endl;

	// Let thread 0 be main thread and created thread be 1 to n_threads
	for (int i = 0; i < n_threads; ++i) {
		std::cout
    << i << ", "
    << threads_vertices[i] << ", "
    << threads_edges[i] << ", "
    << std::setprecision(6)
    << threads_barrier1_time[i] << ", "
    << threads_barrier2_time[i] << ", "
    << threads_get_next_vertex_time[i] << ", "
    << threads_total_time[i]
    << "\n";
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
          {"strategy", "Strategy",
           cxxopts::value<uint>()->default_value("1")},
          {"granularity", "Granularity",
           cxxopts::value<uint>()->default_value("1")}
      });

  auto cl_options = options.parse(argc, argv);
  uint n_threads = cl_options["nThreads"].as<uint>();
  uint max_iterations = cl_options["nIterations"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();
  uint strategy = cl_options["strategy"].as<uint>();
  if (strategy < 1 || strategy > 4) {
    std::cout << "Strategy must be between 1 to 4\n";
  }
  uint granularity = cl_options["granularity"].as<uint>();
  if (granularity <= 0) {
    std::cout << "Strategy must be positive integer\n";
  }
#ifdef USE_INT
  std::cout << "Using INT\n";
#else
  std::cout << "Using DOUBLE\n";
#endif
  std::cout << std::fixed;
  std::cout << "Number of Threads : " << n_threads << std::endl;
  std::cout << "Strategy : " << strategy << std::endl;
  std::cout << "Granularity : " << granularity << std::endl;
  std::cout << "Iterations : " << max_iterations << std::endl;

  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";
  pageRankParallel(g, max_iterations, n_threads, strategy, granularity);

  return 0;
}
