#include "Array.hh"
#include "CGSolver.hh"
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

	if (argc != 5 && rank == ROOT_THREAD)
	{
		std::cout << "Usage: ./cg <grid_points_x> <grid_points_y> <iterations> <eps>" << std::endl;
		MPI_Finalize();
		return 1;
	}
	else if (argc != 5 && rank != ROOT_THREAD)
	{
		MPI_Finalize();
		return 1;
	}

	int nx;
	int ny;
	int maxIter;
	double eps;
	double k = 2.0 * PI;

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
		std::cerr << "Could not parse third argument: " << args[3] << std::endl;
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
	t.reset();
	c.solve();
	elapsedTime = t.elapsed();
	if (rank == ROOT_THREAD)
	{
		std::cout << "Elapsed time: " << elapsedTime << " seconds" << std::endl;

		std::cout << "Saving solution to solution.gnuplot ..." << std::endl;
		c.saveToFile("solution.gnuplot", c.getU());
	}



	MPI_Finalize();

	return 0;
}
