#include "Array.hh"
#include "CGSolver.hh"
#include <iostream>
#include "Timer.h"

#ifndef PI
#define PI (3.1415)
#endif

int main(int argc, char** args)
{

	if (argc != 5)
	{
		std::cout << "Usage: ./rgbs <grid_points_x> <grid_points_y> <iterations>" << std::endl;
		return 0;
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


	// suppress unused warning
	(void) elapsedTime;

	CGSolver c(nx, ny, k, maxIter, eps);
	c.solve();

	c.saveToFile("solution.gnuplot", c.getU());

	return 0;
}
