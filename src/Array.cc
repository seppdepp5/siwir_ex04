#include "Array.hh"
#include <cstdlib>
#include <iostream>
#include "Debug.hh"
#include <cmath>

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
	size = xSize;
	ar = new double[size];
}

Array::Array( int xSize, int ySize )
	: xSize_(xSize)
	, ySize_(ySize)
	, zSize_(1)
	, dimension_(DIM_2D)
{
   // TODO construct 2D array here
	CHECK_MSG(xSize > 0 && ySize > 0, "xSize and ySize must be positive");
	size = xSize * ySize;
	ar = new double[size];
}

Array::Array( int xSize, int ySize, int zSize )
	: xSize_(xSize)
	, ySize_(ySize)
	, zSize_(zSize)
	, dimension_(DIM_3D)
{
   // TODO construct 3D array here
	CHECK_MSG(xSize > 0 && ySize > 0 && zSize > 0, "xSize, ySize and zSize must be positive");
	size = xSize * ySize * zSize;
	ar = new double[size];
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
	size = s.getSize();
	ar = new double[size];

	double *p_s = s.getArray();
	for (int i = 0; i < size; i++)
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
   return size;
}

double Array::getAbsMax()
{
	double max = fabs(ar[0]);

	for (int i = 1; i < size; i++)
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

	for (int i = 0; i < size; i++)
	{
		ar[i] = d_ar[i];
	}
}


double Array::scalProdSelf()
{
	double ret = 0.0;

	// TODO MPI here
	for (int i = 0; i < this->size; i++)
	{
		ret += ar[i] * ar[i];
	}

	return ret;
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

	// TODO MPI here
	for (int i = 0; i < this->size; i++)
	{
		ret += ar[i] * d_ar[i];
	}

	return ret;
}

void Array::add(Array & d)
{

	double * d_ar = d.getArray();

	// TODO: MPI here
	for (int i = 0; i < size; i++)
	{
		ar[i] = ar[i] + d_ar[i];
	}
}

void Array::add(Array & d, Array & target)
{

	double * d_ar = d.getArray();
	double * t_ar = target.getArray();

	// TODO: MPI here
	for (int i = 0; i < size; i++)
	{
		t_ar[i] = ar[i] + d_ar[i];
	}
}

void Array::sub(Array & d)
{

	double * d_ar = d.getArray();

	// TODO: MPI here
	for (int i = 0; i < size; i++)
	{
		ar[i] = ar[i] - d_ar[i];
	}
}

void Array::sub(Array & d, Array & target)
{

	double * d_ar = d.getArray();
	double * t_ar = target.getArray();

	// TODO: MPI here
	for (int i = 0; i < size; i++)
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

	for (int i = 0; i < size; i++)
	{
		sum += ar[i];
	}

	avg = sum / size;

	for (int i = 0; i < size; i++)
	{
		ar[i] -= avg;
	}

}
#endif
