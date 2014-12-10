#include "Array.hh"
#include <cstdlib>
#include <iostream>
#include "Debug.hh"
#include <cmath>
#include <mpi.h>

#ifndef ROOT_THREAD
#define ROOT_THREAD (0)
#endif

#define SEND_TAG (1000)

//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================


Array::Array( int xSize )
	: xSize_(xSize)
	, ySize_(1)
	, zSize_(1)
	, dimension_(DIM_1D)
{
   // TODO construct 1D array here
	CHECK_MSG(xSize > 0, "xSize must be positive");
	size_ = xSize;
	ar = new double[size_];
}

Array::Array( int xSize, int ySize )
	: xSize_(xSize)
	, ySize_(ySize)
	, zSize_(1)
	, dimension_(DIM_2D)
{
   // TODO construct 2D array here
	CHECK_MSG(xSize > 0 && ySize > 0, "xSize and ySize must be positive");
	size_ = xSize * ySize;
	ar = new double[size_];
}

Array::Array( int xSize, int ySize, int zSize )
	: xSize_(xSize)
	, ySize_(ySize)
	, zSize_(zSize)
	, dimension_(DIM_3D)
{
   // TODO construct 3D array here
	CHECK_MSG(xSize > 0 && ySize > 0 && zSize > 0, "xSize, ySize and zSize must be positive");
	size_ = xSize * ySize * zSize;
	ar = new double[size_];
}

Array::Array(const Array & s)
{
	dimension_ = s.getDimension();
	switch (dimension_)
	{
	case DIM_1D:
		xSize_ = s.getSize(DIM_1D);
		ySize_ = 1;
		zSize_ = 1;
		break;
	case DIM_2D:
		xSize_ = s.getSize(DIM_1D);
		ySize_ = s.getSize(DIM_2D);
		zSize_ = 1;
		break;
	case DIM_3D:
		xSize_ = s.getSize(DIM_1D);
		ySize_ = s.getSize(DIM_2D);
		zSize_ = s.getSize(DIM_3D);
		break;
	}
	size_ = s.getSize();
	ar = new double[size_];

	double *p_s = s.getArray();
	for (int i = 0; i < size_; i++)
	{
		ar[i] = p_s[i];
	}

}

Array& Array::operator=(const Array& s)
{
	double *local_ar = new double[s.getSize()];
	double *s_ar = s.getArray();

	for (int i = 0; i < s.getSize(); i++)
	{
		local_ar[i] = s_ar[i];
	}

	this->ar = local_ar;
	this->xSize_ = s.getSize(DIM_1D);
	this->ySize_ = s.getSize(DIM_2D);
	this->zSize_ = s.getSize(DIM_3D);
	this->dimension_ = s.getDimension();

	return *this;
}

//	Destructor
Array::~Array()
{
	if (NULL != ar)
	{
		delete[] ar;
	}
}




//===================================================================================================================
//
//  Convenience Functions
//
//===================================================================================================================


//initialize the whole array with a constant value
void Array::fill( double value )
{
   // TODO
   // you might want to use std::fill() here
	for (int i = 0; i < getSize(); i++)
	{
		ar[i] = value;
	}

}

void Array::fillRandom()
{
	for (int i = 0; i < getSize(); i++)
		{
			ar[i] = (double)rand() / RAND_MAX;;
		}
}

// returns dimension
int Array::getDimension() const
{
	return dimension_;
}

double * Array::getArray() const
{
	return ar;
}

// Print the whole array (for debugging purposes)
void Array::print()
{
   // TODO
   // For 2D Arrays the positive x-coordinate goes to the right
   //                   positive y-coordinate goes upwards
   //      -> the line with highest y-value should be printed first
	switch (dimension_)
	{
	case DIM_1D:
		for (int i = 0; i < xSize_; i++)
		{
			std::cout << Array::operator ()(i) << " ";
		}
		std::cout << std::endl;
		break;
	case DIM_2D:
		for (int i = ySize_ - 1; i >= 0; i--)
		{
			for (int j = 0; j < xSize_; j++)
			{
				std::cout << Array::operator ()(j,i) << " ";
			}
			std::cout << std::endl;
		}
		break;
	case DIM_3D:
		for (int k = 0; k < zSize_; k++)
		{
			std::cout << "z = " << k << std::endl;
			for (int i = ySize_ - 1; i >= 0; i--)
				{
					for (int j = 0; j < xSize_; j++)
					{
						std::cout << Array::operator ()(j,i,k) << " ";
					}
					std::cout << std::endl;
				}
			std::cout << std::endl;
		}
		break;

	}

}

int Array::getSize( int dimension ) const
{
   //TODO
	int ret = 0;
	switch (dimension)
	{
	case DIM_1D:
		ret = xSize_;
		break;
	case DIM_2D:
		ret = ySize_;
		break;
	case DIM_3D:
		ret = zSize_;
		break;
	}
	return ret;
}

//return total size of the array
int Array::getSize() const
{
   //TODO
   return size_;
}

double Array::getAbsMax()
{
	double max = fabs(ar[0]);

	for (int i = 1; i < size_; i++)
	{
		if (fabs(ar[i]) > max)
		{
			max = fabs(ar[i]);
		}
	}

	return max;
}

void Array::copyFrom(Array & d)
{
	double * d_ar = d.getArray();

	for (int i = 0; i < size_; i++)
	{
		ar[i] = d_ar[i];
	}
}


double Array::scalProdSelf()
{
	double ret = 0.0;

	int size;
	int rank;

	int num_elements_to_compute;

	MPI_Status status;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// only for master thread
	if (rank == ROOT_THREAD)
	{

		// size_ is the size of the array
		// whereas size is the number of threads
		//
		// each thread but the master will now get num_elements_to_send
		// elements of the array to do the calculation
		//
		// the master thread will additionally also compute the rest elements
		// (it is not guaranteed that size_ % size == 0)
		int num_elements_to_send = size_ / size;
		int rest = size_ % size;

		// excluding the master thread here
		for (int i = 1; i < size; i++)
		{
			// sending number of array elements to process
			MPI_Send(&num_elements_to_send, 1, MPI_INT, i, SEND_TAG, MPI_COMM_WORLD);
			// sending array address of the first element to be processed by the receiver
			MPI_Send(&ar[num_elements_to_send*i + rest], num_elements_to_send, MPI_DOUBLE, i, SEND_TAG, MPI_COMM_WORLD);
		}

		// master computes one part plus the rest
		num_elements_to_compute = num_elements_to_send + rest;

	}
	// only for "slaves"
	else
	{
		// receive number of elements to compute / receiving array length
		MPI_Recv(&num_elements_to_compute, 1, MPI_INT, ROOT_THREAD, SEND_TAG, MPI_COMM_WORLD, &status);
		// receive array and store at pos 0 of local copy
		// THIS DESTROYS THE ARRAY AND MAKES IT NOT USABLE FOR THIS THREAD AGAIN
		MPI_Recv(&ar[0], num_elements_to_compute, MPI_DOUBLE, ROOT_THREAD, SEND_TAG, MPI_COMM_WORLD, &status);

	}

	// compute here
	for (int i = 0; i < num_elements_to_compute; i++)
	{
		ret += ar[i] * ar[i];
	}

	// use reduce function to add all ret values to one big suḿ in ROOT_THREAD
	double global_sum = 0.0;
	MPI_Reduce(&ret, &global_sum, 1, MPI_DOUBLE, MPI_SUM, ROOT_THREAD, MPI_COMM_WORLD);

	// broadcast global sum, so every thread returns the right value
	MPI_Bcast(&global_sum, 1, MPI_DOUBLE, ROOT_THREAD, MPI_COMM_WORLD);
	return global_sum;


}

double Array::scalProd(Array & d)
{
	// check if sizes fit
	CHECK
		(xSize_ == d.xSize_ &&
		 ySize_ == d.ySize_ &&
		 zSize_ == d.zSize_);

	// get array of d
	double * d_ar = d.getArray();

	double ret = 0.0;

	int size;
	int rank;

	int num_elements_to_compute;

	MPI_Status status;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// only for master thread
	if (rank == ROOT_THREAD)
	{

		// size_ is the size of the array
		// whereas size is the number of threads
		//
		// each thread but the master will now get num_elements_to_send
		// elements of the array to do the calculation
		//
		// the master thread will additionally also compute the rest elements
		// (it is not guaranteed that size_ % size == 0)
		int num_elements_to_send = size_ / size;
		int rest = size_ % size;

		// excluding the master thread here
		for (int i = 1; i < size; i++)
		{
			// sending number of array elements to process
			MPI_Send(&num_elements_to_send, 1, MPI_INT, i, SEND_TAG, MPI_COMM_WORLD);
			// sending array address of the first element to be processed by the receiver for both arrays
			MPI_Send(&ar[num_elements_to_send*i + rest], num_elements_to_send, MPI_DOUBLE, i, SEND_TAG, MPI_COMM_WORLD);
			MPI_Send(&d_ar[num_elements_to_send*i + rest], num_elements_to_send, MPI_DOUBLE, i, SEND_TAG, MPI_COMM_WORLD);
		}

		// master computes one part plus the rest
		num_elements_to_compute = num_elements_to_send + rest;

	}
	// only for "slaves"
	else
	{
		// receive number of elements to compute / receiving array length
		MPI_Recv(&num_elements_to_compute, 1, MPI_INT, ROOT_THREAD, SEND_TAG, MPI_COMM_WORLD, &status);
		// receive arrays and store at pos 0 of local copies
		// THIS DESTROYS THE ARRAYs AND MAKES IT NOT USABLE FOR THIS THREAD AGAIN
		MPI_Recv(&ar[0], num_elements_to_compute, MPI_DOUBLE, ROOT_THREAD, SEND_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(&d_ar[0], num_elements_to_compute, MPI_DOUBLE, ROOT_THREAD, SEND_TAG, MPI_COMM_WORLD, &status);

	}

	// compute here
	for (int i = 0; i < num_elements_to_compute; i++)
	{
		ret += ar[i] * d_ar[i];
	}

	// use reduce function to add all ret values to one big suḿ in ROOT_THREAD
	double global_sum = 0.0;
	MPI_Reduce(&ret, &global_sum, 1, MPI_DOUBLE, MPI_SUM, ROOT_THREAD, MPI_COMM_WORLD);

	// broadcast global sum, so every thread returns the right value
	MPI_Bcast(&global_sum, 1, MPI_DOUBLE, ROOT_THREAD, MPI_COMM_WORLD);

	return global_sum;

}

void Array::addScaleOperand(Array & d, double scaleValue)
{

	double * d_ar = d.getArray();

	int size;
	int rank;

	int num_elements_to_compute;

	MPI_Status status;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// size_ is the size of the array
	// whereas size is the number of threads
	//
	// each thread but the master will now get num_elements_to_send
	// elements of the array to do the calculation
	//
	// the master thread will additionally also compute the rest elements
	// (it is not guaranteed that size_ % size == 0)
	int num_elements_to_send = size_ / size;
	int rest = size_ % size;

	// only for master thread
	if (rank == ROOT_THREAD)
	{



		// excluding the master thread here
		for (int i = 1; i < size; i++)
		{
			// sending number of array elements to process
			MPI_Send(&num_elements_to_send, 1, MPI_INT, i, SEND_TAG, MPI_COMM_WORLD);
			// sending array address of the first element to be processed by the receiver for both arrays
			MPI_Send(&ar[num_elements_to_send*i + rest], num_elements_to_send, MPI_DOUBLE, i, SEND_TAG, MPI_COMM_WORLD);
			MPI_Send(&d_ar[num_elements_to_send*i + rest], num_elements_to_send, MPI_DOUBLE, i, SEND_TAG, MPI_COMM_WORLD);
		}

		// master computes one part plus the rest
		num_elements_to_compute = num_elements_to_send + rest;

	}
	// only for "slaves"
	else
	{
		// receive number of elements to compute / receiving array length
		MPI_Recv(&num_elements_to_compute, 1, MPI_INT, ROOT_THREAD, SEND_TAG, MPI_COMM_WORLD, &status);
		// receive arrays and store at pos 0 of local copies
		// THIS DESTROYS THE ARRAYs AND MAKES IT NOT USABLE FOR THIS THREAD AGAIN
		MPI_Recv(&ar[0], num_elements_to_compute, MPI_DOUBLE, ROOT_THREAD, SEND_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(&d_ar[0], num_elements_to_compute, MPI_DOUBLE, ROOT_THREAD, SEND_TAG, MPI_COMM_WORLD, &status);

	}

	// computation here
	for (int i = 0; i < num_elements_to_compute; i++)
	{
		ar[i] = ar[i] + scaleValue * d_ar[i];
	}

		// send it all back to root thread
	if (rank == ROOT_THREAD)
	{
		for (int i = 1; i < size; i++)
		{
			MPI_Recv(&ar[num_elements_to_send*i + rest], num_elements_to_send, MPI_DOUBLE, i, SEND_TAG, MPI_COMM_WORLD, &status);
		}
	}
	// only for slaves
	else
	{
		MPI_Send(&ar[0], num_elements_to_compute, MPI_DOUBLE, ROOT_THREAD, SEND_TAG, MPI_COMM_WORLD);
	}

	// broadcast to have the right array in every thread
	MPI_Bcast(&ar[0], size_, MPI_DOUBLE, ROOT_THREAD, MPI_COMM_WORLD);
}

void Array::addScaleTarget(Array & d, double scaleValue)
{

	double * d_ar = d.getArray();

	int size;
	int rank;

	int num_elements_to_compute;

	MPI_Status status;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// size_ is the size of the array
	// whereas size is the number of threads
	//
	// each thread but the master will now get num_elements_to_send
	// elements of the array to do the calculation
	//
	// the master thread will additionally also compute the rest elements
	// (it is not guaranteed that size_ % size == 0)
	int num_elements_to_send = size_ / size;
	int rest = size_ % size;

	// only for master thread
	if (rank == ROOT_THREAD)
	{



		// excluding the master thread here
		for (int i = 1; i < size; i++)
		{
			// sending number of array elements to process
			MPI_Send(&num_elements_to_send, 1, MPI_INT, i, SEND_TAG, MPI_COMM_WORLD);
			// sending array address of the first element to be processed by the receiver for both arrays
			MPI_Send(&ar[num_elements_to_send*i + rest], num_elements_to_send, MPI_DOUBLE, i, SEND_TAG, MPI_COMM_WORLD);
			MPI_Send(&d_ar[num_elements_to_send*i + rest], num_elements_to_send, MPI_DOUBLE, i, SEND_TAG, MPI_COMM_WORLD);
		}

		// master computes one part plus the rest
		num_elements_to_compute = num_elements_to_send + rest;

	}
	// only for "slaves"
	else
	{
		// receive number of elements to compute / receiving array length
		MPI_Recv(&num_elements_to_compute, 1, MPI_INT, ROOT_THREAD, SEND_TAG, MPI_COMM_WORLD, &status);
		// receive arrays and store at pos 0 of local copies
		// THIS DESTROYS THE ARRAYs AND MAKES IT NOT USABLE FOR THIS THREAD AGAIN
		MPI_Recv(&ar[0], num_elements_to_compute, MPI_DOUBLE, ROOT_THREAD, SEND_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(&d_ar[0], num_elements_to_compute, MPI_DOUBLE, ROOT_THREAD, SEND_TAG, MPI_COMM_WORLD, &status);

	}

	// computation here
	for (int i = 0; i < num_elements_to_compute; i++)
	{
		ar[i] = d_ar[i] + scaleValue * ar[i] ;
	}

		// send it all back to root thread
	if (rank == ROOT_THREAD)
	{
		for (int i = 1; i < size; i++)
		{
			MPI_Recv(&ar[num_elements_to_send*i + rest], num_elements_to_send, MPI_DOUBLE, i, SEND_TAG, MPI_COMM_WORLD, &status);
		}
	}
	// only for slaves
	else
	{
		MPI_Send(&ar[0], num_elements_to_compute, MPI_DOUBLE, ROOT_THREAD, SEND_TAG, MPI_COMM_WORLD);
	}

	// broadcast to have the right array in every thread
	MPI_Bcast(&ar[0], size_, MPI_DOUBLE, ROOT_THREAD, MPI_COMM_WORLD);
}

void Array::add(Array & d, Array & target)
{

	double * d_ar = d.getArray();
	double * t_ar = target.getArray();

	// TODO: MPI here
	for (int i = 0; i < size_; i++)
	{
		t_ar[i] = ar[i] + d_ar[i];
	}
}

void Array::sub(Array & d)
{

	double * d_ar = d.getArray();

	// TODO: MPI here
	for (int i = 0; i < size_; i++)
	{
		ar[i] = ar[i] - d_ar[i];
	}
}

void Array::sub(Array & d, Array & target)
{

	double * d_ar = d.getArray();
	double * t_ar = target.getArray();

	// TODO: MPI here
	for (int i = 0; i < size_; i++)
	{
		t_ar[i] = ar[i] - d_ar[i];
	}
}


#if 0
void Array::normalize()
{
	// calc avg
	double avg;
	double sum = 0.0;

	for (int i = 0; i < size_; i++)
	{
		sum += ar[i];
	}

	avg = sum / size_;

	for (int i = 0; i < size_; i++)
	{
		ar[i] -= avg;
	}

}
#endif
