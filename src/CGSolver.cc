#include "CGSolver.hh"
#include <cmath>
#include <iostream>
#include <fstream>
#include <mpi.h>

#ifndef PI
#define PI (3.1415)
#endif

#define VERBOSE (1)
#define ROOT_THREAD (0)

CGSolver::CGSolver
	(const int nx
	,const int ny
	,const double k
	,const int maxIter
	,const double eps
	)
	: u_(nx-1, ny-1)
	, f_(nx-1, ny-1)
	, r_(nx-1, ny-1)
	, nx_(nx)
	, ny_(ny)
	, dx_(2.0 / (double)nx)
	, dy_(1.0 / (double)ny)
	, k_(k)
	, maxIter_(maxIter)
	, eps_(eps)
{

	// init Arrays
	u_.fill(0.0);
	f_.fill(0.0);
	r_.fill(0.0);

	for (int j = 0; j < f_.getSize(DIM_2D); j++)
	{
		for (int i = 0; i < f_.getSize(DIM_1D); i++)
		{
			f_(i,j) = 4.0 * PI * PI * sin(2.0 * PI * (i+1)*dx_) * sinh(2.0 * PI * (j+1)*dy_);
		}
	}

	// subtract upper border from f
	double hy = 1.0 / (dy_*dy_);

	for (int i = 0; i < f_.getSize(DIM_1D); i++)
	{
		f_(i, f_.getSize(DIM_2D)-1) -=  hy * sin(2.0 * PI * (i+1)*dx_) * sinh(2.0 * PI);
	}

}

void CGSolver::solve()
{

	int size;
	int rank;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double delta0;
	double delta1;
	double alpha;
	double beta;

	Array z(u_);
	Array d(u_);

	z.fill(0.0);
	d.fill(0.0);

	// 1: r = f - Au
	computeResidual();

	// 2: delta0 = (r,r)
	delta0 = r_.scalProdSelf();

	// 3: if ||r||_2 <= eps then stop
	if (sqrt(delta0/((nx_-1)*(ny_-1))) <= eps_) return;

	// 4: d = r
	d.copyFrom(r_);

	// 5: for number of iterations do
	for (int it = 0; it < maxIter_; it++)
	{
		// 6: z= Ad
		applyOperator(d, z);

		// 7: alpha = delta0 / (d,z)
		alpha = delta0 / d.scalProd(z);

		// 8: u := u + alpha * d
		u_.addScaleOperand(d, alpha);

		// 9: r := r - alpha * z
		r_.addScaleOperand(z, -alpha);

		// 10: delta1 = (r,r)
		delta1 = r_.scalProdSelf();

		// 11: if ||r||_2 <= eps then stop
		if (VERBOSE && it%1000 == 0 && rank == ROOT_THREAD)
		{
			std::cout << "Step: " << it << "\tResidual: " << sqrt(delta1/((nx_-1)*(ny_-1))) << std::endl;
		}
		if (sqrt(delta1/((nx_-1)*(ny_-1))) <= eps_) return;

		// 12: beta = delta1 / delta0
		beta = delta1 / delta0;

		// 13: d = r + beta * d
		d.addScaleTarget(r_, -beta);


		// 14: delta0 = delta1
		delta0 = delta1;
	}


}

int CGSolver::saveToFile(std::string filename, Array & u) const
{
	std::ofstream gnuFile(filename);
	if (gnuFile.is_open())
	{
		for (int j = 0; j < u.getSize(DIM_2D); j++)
		{
			for (int i = 0; i < u.getSize(DIM_1D); i++)
			{
				gnuFile << i << " " << j << " " << u(i,j) << "\n\n";
			}
		}
		gnuFile.close();
		return 0;
	}
	else
	{
		return 1;
	}
}

void CGSolver::applyOperator(Array & u, Array & target)
{
	double hx = 1.0 / (dx_*dx_);
	double hy = 1.0 / (dy_*dy_);
	double hxy = 2.0 * hx + 2.0 * hy + k_ * k_;

	// TODO: MPI here
	// dont know to handle boundaries the fastest way
	// also dont know if that matters at all
	// so here is the easy solution
	for (int j = 0; j < u.getSize(DIM_2D); j++)
	{
		for (int i = 0; i < u.getSize(DIM_1D); i++)
		{
			// inner domain
			if (i > 0 && j > 0 && i < u.getSize(DIM_1D)-1 && j < u.getSize(DIM_2D)-1)
			{
				target(i, j) = hx * (u(i-1,j)+u(i+1,j)) + hy * (u(i,j+1)+u(i,j-1)) - hxy * u(i,j);
			}

			// bottom left corner
			else if (i == 0 && j == 0)
			{
				target(i, j) = hx * u(i+1,j) + hy * u(i,j+1) - hxy * u(i,j);
			}

			// top left corner
			else if (i == 0 && j == u.getSize(DIM_2D)-1)
			{
				target(i, j) = hx * u(i+1,j) + hy * u(i,j-1) - hxy * u(i,j);
			}

			// top right corner
			else if (i == u.getSize(DIM_1D)-1 && j == u.getSize(DIM_2D)-1)
			{
				target(i, j) = hx * u(i-1,j) + hy * u(i,j-1) - hxy * u(i,j);
			}

			// bottom right corner
			else if (i == u.getSize(DIM_1D)-1 && j == 0)
			{
				target(i, j) = hx * u(i-1,j) + hy * u(i,j+1) - hxy * u(i,j);
			}

			// right border
			else if (i == 0 && j != 0 && j != u.getSize(DIM_2D)-1)
			{
				target(i, j) = hx * u(i+1,j) + hy * (u(i,j+1)+u(i,j-1)) - hxy * u(i,j);
			}

			// left border
			else if (i == u.getSize(DIM_1D)-1 && j != 0 && j != u.getSize(DIM_2D)-1)
			{
				target(i, j) = hx * u(i-1,j) + hy * (u(i,j+1)+u(i,j-1)) - hxy * u(i,j);
			}

			// top border
			else if (j == u.getSize(DIM_2D)-1 && i != 0 && i != u.getSize(DIM_1D)-1)
			{
				target(i, j) = hx * (u(i-1,j)+u(i+1,j)) + hy * u(i,j-1) - hxy * u(i,j);
			}

			// bottom border
			else if (j == u.getSize(DIM_2D)-1 && i != 0 && i != u.getSize(DIM_1D)-1)
			{
				target(i, j) = hx * (u(i-1,j)+u(i+1,j)) + hy * u(i,j+1) - hxy * u(i,j);
			}
		}
	}
}

void CGSolver::computeResidual()
{
	// compute A*u and store in r_
	applyOperator(u_, r_);

	// there are nicer solutions to this... cannot be called with __restrict__ yet
	f_.sub(r_, r_);

}


