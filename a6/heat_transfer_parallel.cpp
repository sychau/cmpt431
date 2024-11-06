#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <mpi.h>

#define DEFAULT_GRID_SIZE "1000"
#define DEFAULT_CX "1"
#define DEFAULT_CY "1"
#define DEFAULT_TIME_STEPS "1000"
#define DEFAULT_MIDDLE_TEMP "600"
#define ROOT 0

class TemperatureArray {
private:
	uint size;
	uint step;
	double Cx;
	double Cy;
	double *CurrArray;
	double *PrevArray;
  void assign(double *A, uint x, uint y, double newvalue) {
    A[x*size+y] = newvalue;
  };
  double read(double *A, uint x, uint y) {
    return A[x*size+y];
  }; 
public:	
	TemperatureArray(uint input_size, double iCx, double iCy, double init_temp) { // create array of dimension sizexsize
		size = input_size;
 		Cx = iCx;
    		Cy = iCy;
		step = 0;
		CurrArray = (double *)malloc(size*size*sizeof(double));
		PrevArray = (double *)malloc(size*size*sizeof(double));
		for (uint i = 0; i < size; i++)
			for (uint j = 0; j < size; j++) {
				if ((i > size/3) && (i < 2*size/3) && (j > size/3) && (j < 2*size/3)) {
					assign(PrevArray, i, j, init_temp); assign(CurrArray, i, j, init_temp);
				}
				else {
					assign(PrevArray, i, j, 0); assign (CurrArray, i, j, 0);
				}	
			}
	};
 
	~TemperatureArray() {
		free (PrevArray);   free (CurrArray);
	};

	void IncrementStepCount() { step ++; };

	uint ReadStepCount() { return(step); };

	void ComputeNewTemp(uint x, uint y) {
		if ((x > 0) && (x < size-1) && (y > 0) && (y < size-1))
			assign(CurrArray, x, y , read(PrevArray,x,y)	
				+ Cx * (read(PrevArray, x-1, y) + read(PrevArray, x+1, y) - 2*read(PrevArray, x, y)) 
				+ Cy * (read(PrevArray, x, y-1) + read(PrevArray, x, y+1) - 2*read(PrevArray, x, y)));
	};

	void SwapArrays() {
		double *temp = PrevArray;
		PrevArray = CurrArray;
		CurrArray = temp;
	};	
 
	double temp(uint x, uint y) {
		return read(CurrArray, x, y);
	};

	void write(uint x, uint y, double newvalue) {
		assign(CurrArray, x, y, newvalue);
	};
	
	double* GetCurrArray() {
		return CurrArray;
	}

	void SetCurrArray(double* newCurrArray) {
		CurrArray = newCurrArray;
	}
};

inline double heat_transfer_calculation(uint size, uint start, uint end, TemperatureArray* T, uint steps) {
  timer t1;
  t1.start();

	// // Prepare arguments for Allgatherv
	// int world_size;
	// MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // int world_rank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	// int ele_per_process = size * (end - start + 1);
	// double* recv_buffer = static_cast<double*>(malloc(size * size * sizeof(double)));
	// double* send_buffer = static_cast<double*>(malloc(ele_per_process * sizeof(double)));

	// int* recv_counts = static_cast<int*>(malloc(ele_per_process * sizeof(int)));
	// for (int i = 0; i < world_size; i++) {
	// 	recv_counts[i] = ele_per_process;
	// }
	// int* displs = static_cast<int*>(malloc(ele_per_process * sizeof(int)));
	// for (int i = 0; i < world_size; i++) {
	// 	displs[i] = ele_per_process * i;
	// }

  // for (uint stepcount = 1; stepcount <= steps; stepcount ++) {
	// 	// x: row index, y: col index
	//   for (uint x = start; x <= end; x++) {
	// 	  for (uint y = 0; y < size; y++) {
	// 		  T->ComputeNewTemp(x, y);
	// 	  }
	//   }
	// 	MPI_Barrier(MPI_COMM_WORLD);

	// 	MPI_Allgatherv(T->GetCurrArray() + world_rank * ele_per_process, ele_per_process, MPI_DOUBLE,
	// 			recv_buffer, recv_counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

  //   MPI_Barrier(MPI_COMM_WORLD);

	// 	// Swap array buffer to avoid copying
	// 	double* temp_pointer = T->GetCurrArray();
	// 	T->SetCurrArray(recv_buffer);
	// 	recv_buffer = temp_pointer;

  //   T->SwapArrays();
  //   T->IncrementStepCount();
  // }  // end of current step

	// free(recv_buffer);
	// free(recv_counts);
	// free(displs);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	// --- synchronization: Send and Receive boundary columns from neighbors
	// Even processes communicate with right proces first
	// Odd  processes communicate with left process first
	for (uint stepcount = 1; stepcount <= steps; stepcount ++) {
		// x: row index, y: col index
	  for (uint x = start; x <= end; x++) {
		  for (uint y = 0; y < size; y++) {
			  T->ComputeNewTemp(x, y);
		  }
	  }

		// Wait till calculation in all processes finish
		MPI_Barrier(MPI_COMM_WORLD);

		if (world_rank % 2 == 0)  {   // even rank
			if (world_rank < world_size - 1)  {  // not last process
				// Send my column "end" to the right process world_rank+1
				MPI_Send(T->GetCurrArray()+end*size, size, MPI_DOUBLE,  world_rank+1, end, MPI_COMM_WORLD);
				// Receive column "end+1" from the right process world_rank+1, populate local Curr Array
				MPI_Recv(T->GetCurrArray()+(end+1)*size, size, MPI_DOUBLE, world_rank+1, end+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			if (world_rank > 0)   {  // not first process
				// Receive column "start-1" from the left process world_rank-1, populate local Curr Array
				MPI_Recv(T->GetCurrArray()+(start-1)*size, size, MPI_DOUBLE, world_rank-1, start-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// Send my column "start" to the left process world_rank-1
				MPI_Send(T->GetCurrArray()+start*size, size, MPI_DOUBLE, world_rank-1, start, MPI_COMM_WORLD);
			}
		} else {  // odd rank
			if (world_rank > 0)   {  // not first process
				// Receive column "start-1" from the left process world_rank-1, populate local Curr Array
				MPI_Recv(T->GetCurrArray()+(start-1)*size, size, MPI_DOUBLE, world_rank-1, start-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// Send my column "start" to the left process world_rank-1
				MPI_Send(T->GetCurrArray()+start*size, size, MPI_DOUBLE, world_rank-1, start, MPI_COMM_WORLD);
			}
			if (world_rank < world_size - 1)  {  // not last process
				// Send my column "end" to the right process world_rank+1
				MPI_Send(T->GetCurrArray()+end*size, size, MPI_DOUBLE, world_rank+1, end, MPI_COMM_WORLD);
				// Receive column "end+1" from the right process world_rank+1, populate local Curr Array
				MPI_Recv(T->GetCurrArray()+(end+1)*size, size, MPI_DOUBLE, world_rank+1, end+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}

		// Wait till all communication across processes finish
		MPI_Barrier(MPI_COMM_WORLD);

		T->SwapArrays(); // CurrArray is only partially true
		T->IncrementStepCount();
	}

	return t1.stop();
}

void heat_transfer_calculation_parallel(uint size, TemperatureArray* T, uint steps) {
  timer overall_timer;
  overall_timer.start();
  //*------------------------------------------------------------------------
  int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  int min_columns = size / world_size;
  int excess_columns = size % world_size;

	auto getStartEnd = [min_columns, excess_columns](int world_rank) {
			int startx, endx;
			if (world_rank < excess_columns) {
					startx = world_rank * (min_columns + 1);
					endx = startx + min_columns;
			} else {
					startx = (excess_columns * (min_columns + 1)) + ((world_rank - excess_columns) * min_columns);
					endx = startx + min_columns - 1;
			}
			return std::make_pair(startx, endx);
	};
	std::pair<int, int> res = getStartEnd(world_rank);
	int startx = res.first;
	int endx = res.second;

	double local_time_taken = heat_transfer_calculation(size, startx, endx, T, steps);

	double* global_time_taken = nullptr;
	if (world_rank == ROOT) {
		global_time_taken = static_cast<double*>(malloc(sizeof(double) * world_size));
	}
	MPI_Gather(&local_time_taken, 1, MPI_DOUBLE, global_time_taken, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

	// Gather all Temp Array to root process
	if (world_rank == ROOT) {
		for (int wr = ROOT+1; wr < world_size; wr++) {
			std::pair<int, int> wr_res = getStartEnd(wr);
			int wr_startx = wr_res.first;
			int wr_endx = wr_res.second;
			int wr_col = wr_endx - wr_startx + 1;
    	MPI_Recv(T->GetCurrArray()+wr_startx*size, wr_col*size, MPI_DOUBLE, wr, wr, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	} else {
		MPI_Send(T->GetCurrArray()+startx*size, (endx-startx+1)*size, MPI_DOUBLE, ROOT, world_rank, MPI_COMM_WORLD);
	}
  //*------------------------------------------------------------------------
	if (world_rank != ROOT) return;

  // Print these statistics for each thread 
  std::cout << "rank, start_column, end_column, time_taken\n";

	int c = endx + 1;
	for(int i = 0; i < world_size; i++) {
  	std::cout << i << ", "
							<< i * c + startx << ", "
							<< i * c + endx << ", " << std::setprecision(TIME_PRECISION)
							<< global_time_taken[i] << "\n";
	}
	std::vector<uint> points_of_interest;
  uint step = size/6;
  uint position = 0;
  for (uint x = 0; x < 6; x++) {
			points_of_interest.push_back(position);
      position += step;
  } 
  // Print temparature at select boundary points;  
	for (uint i = 0; i < world_size; i++) {
		points_of_interest.push_back(endx + i * c);
  	
	}
	std::sort(points_of_interest.begin(), points_of_interest.end());
	for (auto& p : points_of_interest) {
		std::cout<< std::setprecision(TIME_PRECISION) << "Temp[" << p << "," << p << "]=" << T->temp(p, p) << "\n";
	}
  //*------------------------------------------------------------------------
  double overall_time_taken = overall_timer.stop();
  std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
            << overall_time_taken << "\n";

	free(global_time_taken);
}

int main(int argc, char *argv[]) {
  MPI_Init(NULL, NULL);
  
  // Initialize command line arguments
  cxxopts::Options options("Heat_transfer_calculation",
                           "Model heat transfer in a grid using serial and parallel execution");
  options.add_options(
      "custom",
      {
          {"gSize", "Grid Size",         
           cxxopts::value<uint>()->default_value(DEFAULT_GRID_SIZE)},
          {"mTemp", "Temperature in middle of array",         
           cxxopts::value<double>()->default_value(DEFAULT_MIDDLE_TEMP)},
	        {"iCX", "Coefficient of horizontal heat transfer",
           cxxopts::value<double>()->default_value(DEFAULT_CX)},
          {"iCY", "Coefficient of vertical heat transfer",
           cxxopts::value<double>()->default_value(DEFAULT_CY)},
          {"tSteps", "Time Steps",
           cxxopts::value<uint>()->default_value(DEFAULT_TIME_STEPS)}
      });
  auto cl_options = options.parse(argc, argv);
  uint grid_size = cl_options["gSize"].as<uint>();
  double init_temp = cl_options["mTemp"].as<double>();
  double Cx = cl_options["iCX"].as<double>();
  double Cy = cl_options["iCY"].as<double>();
  uint steps = cl_options["tSteps"].as<uint>();

	int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (world_rank == ROOT) {
    std::cout << "Number of processes : " << world_size << "\n";
    std::cout << "Grid Size : " << grid_size << "x" << grid_size << "\n";
    std::cout << "Cx : " << Cx << "\n" << "Cy : " << Cy << "\n";
    std::cout << "Temperature in the middle of grid : " << init_temp << "\n";
    std::cout << "Time Steps : " << steps << "\n";

    std::cout << "Initializing Temperature Array..." << "\n";
  }

  TemperatureArray *T = new TemperatureArray(grid_size, Cx, Cy, init_temp);
  if (!T) {
      std::cout << "Cannot Initialize Temperature Array...Terminating" << "\n";
      return 2;
  }
  heat_transfer_calculation_parallel (grid_size, T, steps);

  delete T;

  MPI_Finalize();
  return 0;
}
