#include "Array.hh"
#include "CGSolver.hh"
#include "HeatSolver.hh"
#include <iostream>
#include "Timer.h"
#include <mpi.h>

#ifndef PI
#define PI (3.14159265358979323846264338327)
#endif

#define ROOT_THREAD (0)

int main(int argc, char** args)
{

	MPI_Init(&argc, &args);

	int size;
	int rank;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (argc != 10 && rank == ROOT_THREAD)
	{
		std::cout << "Usage: ./cg <grid_points_x> <grid_points_y> <iterations> <eps> <timesteps> <dt> <k> <alpha> <vtk_spacing>" << std::endl;
		MPI_Finalize();
		return 1;
	}
	else if (argc != 10 && rank != ROOT_THREAD)
	{
		MPI_Finalize();
		return 1;
	}

	int nx;
	int ny;
	int maxIter;
	double eps;
	double k = 2.0 * PI;

#if 1
	double alpha;
	double kk;
	double dt;
	int timesteps;
	int vtk_spacing;
#endif

#if 0
	double alpha = .5;
	double kk = 0.8;
	double dt = 0.001;
	int timesteps = 50;
#endif

	siwir::Timer t;
	double elapsedTime;

	std::istringstream iss(args[1]);
	if (!(iss >> nx))
	{
		std::cerr << "Could not parse first argument: " << args[1] << std::endl;
		return 1;
	}
	iss.str("");
	iss.clear();

	iss.str(args[2]);
	if (!(iss >> ny))
	{
		std::cerr << "Could not parse second argument: " << args[2] << std::endl;
		return 1;
	}
	iss.str("");
	iss.clear();

	iss.str(args[3]);
	if (!(iss >> maxIter))
	{
		std::cerr << "Could not parse third argument: " << args[3] << std::endl;
		return 1;
	}

	iss.str("");
	iss.clear();

	iss.str(args[4]);
	if (!(iss >> eps))
	{
		std::cerr << "Could not parse 4th argument: " << args[4] << std::endl;
		return 1;
	}

	iss.str("");
	iss.clear();

	iss.str(args[5]);
	if (!(iss >> timesteps))
	{
		std::cerr << "Could not parse 5th argument: " << args[5] << std::endl;
		return 1;
	}

	iss.str("");
	iss.clear();

	iss.str(args[6]);
	if (!(iss >> dt))
	{
		std::cerr << "Could not parse 6th argument: " << args[6] << std::endl;
		return 1;
	}

	iss.str("");
	iss.clear();

	iss.str(args[7]);
	if (!(iss >> kk))
	{
		std::cerr << "Could not parse 7th argument: " << args[7] << std::endl;
		return 1;
	}

	iss.str("");
	iss.clear();

	iss.str(args[8]);
	if (!(iss >> alpha))
	{
		std::cerr << "Could not parse 8th argument: " << args[8] << std::endl;
		return 1;
	}

	iss.str("");
	iss.clear();

	iss.str(args[9]);
	if (!(iss >> vtk_spacing))
	{
		std::cerr << "Could not parse 9th argument: " << args[9] << std::endl;
		return 1;
	}


	if (maxIter < 1)
	{
		if (rank == ROOT_THREAD)
		{
			std::cout << "Please set the number of iterations to a value higher than 0.\nExiting." << std::endl;
			MPI_Finalize();
			return 0;
		}
		else
		{
			MPI_Finalize();
			return 0;
		}
	}


	CGSolver c(nx, ny, k, maxIter, eps);
	HeatSolver h(c);
	t.reset();
	h.solve(alpha, kk, dt, timesteps, vtk_spacing);
	elapsedTime = t.elapsed();
	if (rank == ROOT_THREAD)
	{
		std::cout << "Elapsed time: " << elapsedTime << " seconds" << std::endl;

		//std::cout << "Saving solution to solution.gnuplot ..." << std::endl;
	}



	MPI_Finalize();

	return 0;
}
