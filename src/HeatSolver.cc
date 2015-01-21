#include "HeatSolver.hh"
#include "Array.hh"
#include "CGSolver.hh"
#include <cmath>
#include <mpi.h>
#include <string>

#ifndef PI
#define PI (3.14159265358979323846264338327)
#endif

HeatSolver::HeatSolver(CGSolver & cgSolver)
	: cgSolver_(cgSolver)
{

}

void HeatSolver::initU()
{
	Array &  u_ = cgSolver_.getU();
	double dx_ = cgSolver_.getDX();
	double dy_ = cgSolver_.getDY();

	for (int j = 0; j < u_.getSize(DIM_2D); j++)
	{
		for (int i = 0; i < u_.getSize(DIM_1D); i++)
		{
			u_(i,j) = sin(PI * (i+1)*dx_) * sin(PI * (j+1)*dy_);
		}
	}
}

void HeatSolver::solve(double alpha, double k, double dt, int timesteps)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	initU();

	for (int step = 0; step < timesteps; step++)
	{
		updateF(alpha, k, dt);
		cgSolver_.solve(dt, alpha, k);

		if (rank == 0)
		{
			std::string name = "solution_" + std::to_string(step);
			cgSolver_.saveToFile(name, cgSolver_.getU());
		}
	}

}

void HeatSolver::updateF(double alpha, double k, double dt)
{
	Array & u = cgSolver_.getU();
	Array & f = cgSolver_.getF();
	double dx_ = cgSolver_.getDX();
	double dy_ = cgSolver_.getDY();

	double hx =  ((1.0-alpha) * k) / (dx_*dx_);
	double hy =  ((1.0-alpha) * k) / (dy_*dy_);
	double hxy = - (1.0/dt - 2.0*hx - 2.0*hy);



	for (int j = 0; j < f.getSize(DIM_2D); j++)
	{
		for (int i = 0; i < f.getSize(DIM_1D); i++)
		{
			// inner domain
			if (i > 0 && j > 0 && i < u.getSize(DIM_1D)-1 && j < u.getSize(DIM_2D)-1)
			{
				f(i, j) = hx * (u(i-1,j)+u(i+1,j)) + hy * (u(i,j+1)+u(i,j-1)) - hxy * u(i,j);
			}
#if 1
			// bottom left corner
			else if (i == 0 && j == 0)
			{
				f(i, j) = hx * u(i+1,j) + hy * u(i,j+1) - hxy * u(i,j);
			}

			// top left corner
			else if (i == 0 && j == u.getSize(DIM_2D)-1)
			{
				f(i, j) = hx * u(i+1,j) + hy * u(i,j-1) - hxy * u(i,j);
			}

			// top right corner
			else if (i == u.getSize(DIM_1D)-1 && j == u.getSize(DIM_2D)-1)
			{
				f(i, j) = hx * u(i-1,j) + hy * u(i,j-1) - hxy * u(i,j);
			}

			// bottom right corner
			else if (i == u.getSize(DIM_1D)-1 && j == 0)
			{
				f(i, j) = hx * u(i-1,j) + hy * u(i,j+1) - hxy * u(i,j);
			}

			// right border
			else if (i == 0 && j != 0 && j != u.getSize(DIM_2D)-1)
			{
				f(i, j) = hx * u(i+1,j) + hy * (u(i,j+1)+u(i,j-1)) - hxy * u(i,j);
			}

			// left border
			else if (i == u.getSize(DIM_1D)-1 && j != 0 && j != u.getSize(DIM_2D)-1)
			{
				f(i, j) = hx * u(i-1,j) + hy * (u(i,j+1)+u(i,j-1)) - hxy * u(i,j);
			}

			// top border
			else if (j == u.getSize(DIM_2D)-1 && i != 0 && i != u.getSize(DIM_1D)-1)
			{
				f(i, j) = hx * (u(i-1,j)+u(i+1,j)) + hy * u(i,j-1) - hxy * u(i,j);
			}

			// bottom border
			else if (j == 0 && i != 0 && i != u.getSize(DIM_1D)-1)
			{
				f(i, j) = hx * (u(i-1,j)+u(i+1,j)) + hy * u(i,j+1) - hxy * u(i,j);
			}
#endif
		}
	}
}





































