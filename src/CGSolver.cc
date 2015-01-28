#include "CGSolver.hh"
#include <cmath>
#include <iostream>
#include <fstream>
#include <mpi.h>

#ifndef PI
#define PI (3.14159265358979323846264338327)
#endif

#define VERBOSE (0)
#define ROOT_THREAD (0)
#define SEND_TAG (1000)

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
	, dx_(1.0 / (double)nx)
	, dy_(1.0 / (double)ny)
	, k_(k)
	, maxIter_(maxIter)
	, eps_(eps)
{

	// init Arrays
	u_.fill(0.0);
	f_.fill(0.0);
	r_.fill(0.0);
#if 0
	for (int j = 0; j < f_.getSize(DIM_2D); j++)
	{
		for (int i = 0; i < f_.getSize(DIM_1D); i++)
		{
			f_(i,j) = 4.0 * PI * PI * sin(2.0 * PI * (i+1)*dx_) * sinh(2.0 * PI * (j+1)*dy_);
		}
	}

	// substract upper border from f
	double hy =  - 1.0 / (dy_*dy_);

	for (int i = 0; i < f_.getSize(DIM_1D); i++)
	{
		f_(i, f_.getSize(DIM_2D)-1) -=  hy * sin(2.0 * PI * (i+1)*dx_) * sinh(2.0 * PI);
	}
#endif

}

void CGSolver::solve(double dt, double alphaCrank, double k)
{



	int size;
	int rank;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double delta0 = 0.0;
	double delta1;
	double alpha;
	double beta;

	double residual = 0.0;

	Array z(u_);
	Array d(u_);

	z.fill(0.0);
	d.fill(0.0);

	// 1: r = f - Au
	computeResidual(dt, alphaCrank, k);

	// 2: delta0 = (r,r)
	delta0 = r_.scalProdSelf();

	// 3: if ||r||_2 <= eps then stop
	if (sqrt(delta0/((nx_-1)*(ny_-1))) <= eps_) 
	{
		if (rank == ROOT_THREAD) std::cout << "CG iterations: 0, Residual: "  << sqrt(delta0/((nx_-1)*(ny_-1))) << std::endl;
		return;
	}

	// 4: d = r
	d.copyFrom(r_);

	// 5: for number of iterations do
	for (int it = 0; it < maxIter_; it++)
	{

		// 6: z= Ad
		applyOperatorHeat(d, z, dt, alphaCrank, k);

		// 7: alpha = delta0 / (d,z)
		alpha = delta0 / d.scalProd(z);

		// 8: u := u + alpha * d
		u_.addScaleOperand(d, alpha);


		// 9: r := r - alpha * z
		r_.addScaleOperand(z, -alpha);

		// 10: delta1 = (r,r)
		delta1 = r_.scalProdSelf();

		// 11: if ||r||_2 <= eps then stop
		residual = sqrt(delta1/((nx_-1)*(ny_-1)));
		if (VERBOSE && it%1000 == 0 && rank == ROOT_THREAD)
		{
		//	std::cout << "Step: " << it << "\tResidual: " << residual << std::endl;
		}
		if (residual <= eps_)
		{
			// if (rank == ROOT_THREAD) std::cout << "Residual dropped below eps after " << it << " iterations. Residual = " << residual << std::endl;
			if (rank == ROOT_THREAD) std::cout << "CG iterations: " << it << " Residual: "  << residual << std::endl;
			return;
		}

		// 12: beta = delta1 / delta0
		beta = delta1 / delta0;

		// 13: d = r + beta * d
		d.addScaleTarget(r_, beta);


		// 14: delta0 = delta1
		delta0 = delta1;
	}

	// if (rank == ROOT_THREAD) std::cout << "Solver ran through all " << maxIter_ << " iterations. Residual = " << residual << std::endl;
	if (rank == ROOT_THREAD) std::cout << "CG iterations: " << maxIter_ << " Residual: "  << residual << std::endl;
}

int CGSolver::saveToPvd(std::string filename, int time_steps) const
{
	std::ofstream gnuFile(filename);
	if (gnuFile.is_open())
	{
		gnuFile << "<?xml version=\"1.0\"?>\n";		
		gnuFile << "<VTKFile type=\"Collection\" version=\"" << 0.1 << "\">\n";		
		gnuFile << "<Collection>\n";		
		for(int i = 0; i < time_steps; i++)
		{
			gnuFile << "<DataSet timestep=\"" << i << "\" file=\"time_step_" << i << ".pvtu\"" << "/>\n";		
		}
		gnuFile << "</Collection>\n";		
		gnuFile << "</VTKFile>\n";		
		
		gnuFile.close();
		return 0;
	}
	else
	{
		return 1;
	}
}

int CGSolver::saveToPvtu(std::string filename, int time_steps, int size) const
{
	std::ofstream gnuFile(filename);
	if (gnuFile.is_open())
	{
		gnuFile << "<?xml version=\"1.0\"?>\n";		
		gnuFile << "<VTKFile type=\"PUnstructuredGrid\" version=\"" << 0.1 << "\">\n";		
		gnuFile << "<PUnstructuredGrid>\n";		
		gnuFile << "<PPoints>\n";
		gnuFile << "<PDataArray type=\"Float64\" NumberOfComponents=\"" << 3 << "\" format=\"ascii\"/>\n";
		gnuFile << "</PPoints>\n";		
		gnuFile << "<PPointData>\n";		
		gnuFile << "<PDataArray type=\"Float64\" Name=\"Temperature\" NumberOfComponents=\"" << 1 << "\" format=\"ascii\"/>\n";
		gnuFile << "</PPointData>\n";
		// std::cout << "size: " << size << std::endl;
		for(int i = 0; i < size; i++)
		{
			gnuFile << "<Piece Source=\"time_step_" << time_steps << "_" << i << ".vtu\"" << "/>\n";		
		}
		gnuFile << "</PUnstructuredGrid>\n";		
		gnuFile << "</VTKFile>\n";		
		
		gnuFile.close();
		return 0;
	}
	else
	{
		return 1;
	}
}

int CGSolver::saveToVtu(std::string filename, Array & u, int firstRow, int lastRow) const
{
//std::cout << "data test" << " u: " << u(1,1) << " DIM_1D: " << u.getSize(DIM_1D) << std::endl;
	std::ofstream gnuFile(filename);
	if (gnuFile.is_open())
	{
		gnuFile << "<?xml version=\"1.0\"?>\n";		
		gnuFile << "<VTKFile type=\"UnstructuredGrid\" version=\"" << 0.1 << "\">\n";		
		gnuFile << "<UnstructuredGrid>\n";		
		gnuFile << "<Piece NumberOfPoints=\"" << (lastRow-firstRow)*u.getSize(DIM_1D) << "\" NumberOfCells=\"" << (lastRow-firstRow)*u.getSize(DIM_1D) << "\">\n";
		gnuFile << "<Points>\n";		
		gnuFile << "<DataArray type=\"Float64\" NumberOfComponents=\"" << 3 << "\" format=\"ascii\">\n";
//KOORD ARRAY
		for (int j = firstRow; j < lastRow; j++)
		{
			for (int i = 0; i < u.getSize(DIM_1D); i++)
			{
				if(i == 0 || i == u.getSize(DIM_1D)-1)
				{
					if(i == u.getSize(DIM_1D)-1)
					{
						gnuFile << "1" << "." << "00" << " ";	
					}else{
						gnuFile << "0" << "." << "00" << " ";
					}
				}else{
					double x = ((double)i)/(u.getSize(DIM_1D)-1);
					x*=100;
					x = (int)x;
					x/=100;
					gnuFile << x << " "; 
				}
				if(j == 0 || j == u.getSize(DIM_2D)-1)
				{
					if(j == u.getSize(DIM_2D)-1)
					{
						gnuFile << "1" << "." << "00" << " ";	
					}else{
						gnuFile << "0" << "." << "00" << " ";
					}
				}else{
					double x = ((double)j)/(u.getSize(DIM_2D)-1); 
					x*=100;
					x = (int)x;
					x/=100;
					gnuFile << x << " "; 
				}
				gnuFile << "0.00" << "\n";
			}
		}
		
		gnuFile << "\n</DataArray>";
		gnuFile << "\n</Points>";
		gnuFile << "\n<Cells>";
		gnuFile << "\n<DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">";
		for(int i = 0; i < (lastRow - firstRow)*u.getSize(DIM_1D); i++)
		{
			gnuFile << " " << i << " ";
		}
		gnuFile << "</DataArray>\n";
		gnuFile << "<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">";
		for(int i = 0; i < (lastRow - firstRow)*u.getSize(DIM_1D); i++)
		{
			gnuFile << " " << i+1 << " ";
		}
		gnuFile << "</DataArray>\n";
		gnuFile << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
		for(int i = 0; i < (lastRow - firstRow)*u.getSize(DIM_1D); i++)
		{
			gnuFile << " " << 1 << " ";
		}
		gnuFile << "</DataArray>\n";
		gnuFile << "</Cells>\n";
		gnuFile << "<PointData>\n";
		gnuFile << "<DataArray type=\"Float64\" Name=\"Temperature\" NumberOfComponents=\"1\" format=\"ascii\">\n";
//DATA ARRAY
		for (int j = firstRow; j < lastRow; j++)
		{
			for(int i = 0; i < u.getSize(DIM_1D); i++)
                	{
                		gnuFile << u(i,j) << " ";
                	}
		}

		gnuFile << "\n</DataArray>";
		gnuFile << "\n</PointData>";
		gnuFile << "\n</Piece>";
		gnuFile << "\n</UnstructuredGrid>";
		gnuFile << "\n</VTKFile>";
		
		return 0;
	}
	else
	{
		return 1;
	}
}

int CGSolver::saveToFile(std::string filename, Array & u) const
{
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// ATTENTION: PLOT INCLUDES HACK TO ADD BOUNDARY VALUES (which normally do not belong to the vector u_) //
	//////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::ofstream gnuFile(filename);
	if (gnuFile.is_open())
	{
#if 0
		// BOUNDARY HACK ///////////////////////////////////////////////////////////////
		// bottom row full of zeros
		gnuFile << 0 << " " << 0 << " " << "0" << "\n\n";
		for (int i = 0; i < u.getSize(DIM_1D); i++)
		{
			gnuFile << i+1 << " " << 0 << " " << u(i,0) << "\n\n";
		}
		gnuFile << u.getSize(DIM_1D) + 1 << " " << 0 << " " << "0" << "\n\n";
		////////////////////////////////////////////////////////////////////////////////
#endif
		// NORMAL PLOT
		// middle rows
		for (int j = 0; j < u.getSize(DIM_2D); j++)
		{
			gnuFile << 0 << " " << j+1 << " " << "0" << "\n\n";
			for (int i = 0; i < u.getSize(DIM_1D); i++)
			{
				gnuFile << i+1 << " " << j+1 << " " << u(i,j) << "\n\n";
			}
			gnuFile << u.getSize(DIM_1D) + 1 << " " << j+1 << " " << "0" << "\n\n";
		}
#if 0
		// BOUNDARY HACK ///////////////////////////////////////////////////////////////
		// top row
		gnuFile << 0 << " " << u.getSize(DIM_2D) + 1 << " " << "0" << "\n\n";
		for (int i = 0; i < u_.getSize(DIM_1D); i++)
		{
			gnuFile << i+1 << " " << u.getSize(DIM_2D) + 1 << " " <<   sin(2.0 * PI * (i+1)*dx_) * sinh(2.0 * PI) << "\n\n";
		}
		gnuFile << u.getSize(DIM_1D) + 1 << " " << u.getSize(DIM_2D) + 1 << " " << "0" << "\n\n";
		////////////////////////////////////////////////////////////////////////////////
#endif


		gnuFile.close();
		return 0;
	}
	else
	{
		return 1;
	}
}

void CGSolver::applyOperatorHeat(Array & u, Array & target, double dt, double alpha, double k)
{
	int size;
	int rank;
	MPI_Status status;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double hx = - (dt* alpha * k) / (dx_*dx_);
	double hy = - (dt* alpha * k) / (dy_*dy_);
	double hxy = - (1.0 - 2.0*hx - 2.0*hy);


	int num_cols = u_.getSize(DIM_1D);
	int num_rows = u_.getSize(DIM_2D);
	int num_rows_to_compute;

	// first row this thread has to compute
	int first_row;
	// last row this thread has to compute (exclusive)
	int last_row;

	// each thread but the master will now get num_rows_to_send
	// rows of the array to do the calculation
	//
	// the master thread will additionally also compute the rest rows
	// (it is not guaranteed that the number of rows % size == 0)
	/*
	 * Example: 3 processes, 4 rows
	 *
	 * =========== -> rank == 2
	 * =========== -> rank == 1
	 * =========== -> master
	 * =========== -> master
	 */
	int num_rows_to_send = num_rows / size;
	int rest_rows = num_rows % size;

	int num_elements_to_send = num_rows_to_send * num_cols;
	int rest_elements = rest_rows * num_cols;

	if (rank == ROOT_THREAD)
	{
		num_rows_to_compute = num_rows_to_send + rest_rows;
		first_row = 0;
		last_row = first_row + num_rows_to_compute;
	}
	else
	{
		num_rows_to_compute = num_rows_to_send;
		first_row = rank * num_rows_to_compute + rest_rows;
		last_row = first_row + num_rows_to_compute;
	}

	// dont know to handle boundaries the fastest way
	// also dont know if that matters at all
	// so here is the easy solution
	for (int j = first_row; j < last_row; j++)
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
			else if (j == 0 && i != 0 && i != u.getSize(DIM_1D)-1)
			{
				target(i, j) = hx * (u(i-1,j)+u(i+1,j)) + hy * u(i,j+1) - hxy * u(i,j);
			}
		}
	}
/*
	//if steps -1 calculating residual
	if(steps >= 0)
	{
		std::cout << target(0,0) << " " << target(1,1) << " " << target(1,0) << std::endl;
		std::string name = "time_step_" + std::to_string(steps) + "_" + std::to_string(rank) + ".vtu";
		saveToVtu(name, target, first_row, last_row);
	}

MPI_Barrier(MPI_COMM_WORLD);
*/
	///////////////////////////////////////////
	// MERGE target and broadcast afterwards //
	///////////////////////////////////////////

	// get the double array from our u Array-object
	double * ar = target.getArray();

	// send it all back to root thread
	if (rank == ROOT_THREAD)
	{


		// receive from all slaves
		for (int i = 1; i < size; i++)
		{
			MPI_Recv(&ar[num_elements_to_send*i + rest_elements], num_elements_to_send, MPI_DOUBLE, i, SEND_TAG, MPI_COMM_WORLD, &status);
		}
	}
	// only for slaves
	else
	{
		MPI_Send(&ar[num_elements_to_send*rank + rest_elements], num_elements_to_send, MPI_DOUBLE, ROOT_THREAD, SEND_TAG, MPI_COMM_WORLD);
	}

	// broadcast to have the right array in every thread
	MPI_Bcast(&ar[0], target.getSize(), MPI_DOUBLE, ROOT_THREAD, MPI_COMM_WORLD);

}

void CGSolver::applyOperator(Array & u, Array & target)
{
	int size;
	int rank;
	MPI_Status status;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double hx = -1.0 / (dx_*dx_);
	double hy = -1.0 / (dy_*dy_);
	double hxy = 2.0 * hx + 2.0 * hy - k_ * k_;

	int num_cols = u_.getSize(DIM_1D);
	int num_rows = u_.getSize(DIM_2D);
	int num_rows_to_compute;

	// first row this thread has to compute
	int first_row;
	// last row this thread has to compute (exclusive)
	int last_row;

	// each thread but the master will now get num_rows_to_send
	// rows of the array to do the calculation
	//
	// the master thread will additionally also compute the rest rows
	// (it is not guaranteed that the number of rows % size == 0)
	/*
	 * Example: 3 processes, 4 rows
	 *
	 * =========== -> rank == 2
	 * =========== -> rank == 1
	 * =========== -> master
	 * =========== -> master
	 */
	int num_rows_to_send = num_rows / size;
	int rest_rows = num_rows % size;

	int num_elements_to_send = num_rows_to_send * num_cols;
	int rest_elements = rest_rows * num_cols;

	if (rank == ROOT_THREAD)
	{
		num_rows_to_compute = num_rows_to_send + rest_rows;
		first_row = 0;
		last_row = first_row + num_rows_to_compute;
	}
	else
	{
		num_rows_to_compute = num_rows_to_send;
		first_row = rank * num_rows_to_compute + rest_rows;
		last_row = first_row + num_rows_to_compute;
	}

	// dont know to handle boundaries the fastest way
	// also dont know if that matters at all
	// so here is the easy solution
	for (int j = first_row; j < last_row; j++)
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
			else if (j == 0 && i != 0 && i != u.getSize(DIM_1D)-1)
			{
				target(i, j) = hx * (u(i-1,j)+u(i+1,j)) + hy * u(i,j+1) - hxy * u(i,j);
			}
		}
	}

	///////////////////////////////////////////
	// MERGE target and broadcast afterwards //
	///////////////////////////////////////////

	// get the double array from our u Array-object
	double * ar = target.getArray();

	// send it all back to root thread
	if (rank == ROOT_THREAD)
	{


		// receive from all slaves
		for (int i = 1; i < size; i++)
		{
			MPI_Recv(&ar[num_elements_to_send*i + rest_elements], num_elements_to_send, MPI_DOUBLE, i, SEND_TAG, MPI_COMM_WORLD, &status);
		}
	}
	// only for slaves
	else
	{
		MPI_Send(&ar[num_elements_to_send*rank + rest_elements], num_elements_to_send, MPI_DOUBLE, ROOT_THREAD, SEND_TAG, MPI_COMM_WORLD);
	}

	// broadcast to have the right array in every thread
	MPI_Bcast(&ar[0], target.getSize(), MPI_DOUBLE, ROOT_THREAD, MPI_COMM_WORLD);

}

void CGSolver::computeResidual(double dt, double alphaCrank, double k)
{
	// compute A*u and store in r_
	applyOperatorHeat(u_, r_, dt,  alphaCrank,  k);

	// there are nicer solutions to this... cannot be called with __restrict__ yet
	r_.addScaleTarget(f_, -1.0);


}


