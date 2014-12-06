#ifndef ARRAY_HH
#define ARRAY_HH

#include <iostream>
#include "Debug.hh"

#define DIM_1D (0x0)
#define DIM_2D (0x1)
#define DIM_3D (0x2)

//*******************************************************************************************************************
/*!  Array class for 1,2 and 3 dimensions
*
*    - all elements should be stored in a contiguous chunk of memory ( no vector<vector> ! )
*/
//*******************************************************************************************************************
class Array
{
public:
   // Constructors for 1D,2D and 3D
   Array( int xSize );
   Array( int xSize, int ySize );
   Array( int xSize, int ySize, int zSize );


   // Depending on your implementation you might need the following:
   ~Array();
   Array(const Array& s);
   Array& operator= (const Array& s);


   // Access Operators for 1D, 2D and 3D
   inline double & operator () ( int i );
   inline double & operator () ( int i ,int j );
   inline double & operator () ( int i, int j, int k );

   // for const Arrays the following access operators are required
   inline const double & operator () ( int i ) const;
   inline const double & operator () ( int i ,int j ) const;
   inline const double & operator () ( int i, int j, int k ) const;



   // initialize the whole array with a constant value
   void fill( double value );

   // initialize with random doubles between 0.0 and 1.0
   void fillRandom();

   // return total size of the array
   int getSize() const;

   // return xSize for dimension==0, ySize for dimension==1 and zSize for dimension==2
   // other dimension values are not allowed
   int getSize(int dimension ) const;

   // added to skeleton
   // return dimension of the array
   int getDimension() const;

   // added to skeleton
   // return address of first element of the array
   double * getArray() const;

   // Print the whole array ( for debugging purposes )
   void print();

   // returns the maximum (absolute) value of the array
   double getAbsMax();

   // normalizes the array around zero (subtracting by avg)
   void normalize();

   // copies every value in d to this array
   void copyFrom(Array & d);

   // scalar prod with itself, stored in d
   double scalProdSelf();

   // scalar prod with other Array
   double scalProd(Array & d);

   // add to this
   void add(Array & d);

   // add to this, store in target
   void add(Array & d, Array & target);

   // subtract from this
   void sub(Array & d);

   // subtract from this, store in target
   void sub(Array & d, Array & target);

private:

   // Sizes
   int xSize_;
   int ySize_;
   int zSize_;
   int size;

   // Dimension of the array
   int dimension_;

   // Array pointer holding the actual elements
   double *ar;


};


//===================================================================================================================
//
//  Inline Access Operators and Sizes
//
//===================================================================================================================



// Operator() 1D
inline double& Array::operator ()(int i)
{
   //TODO
   //static double dummy;
	CHECK_MSG(i >= 0 && i < size, "Index i out of bounds");
	CHECK_MSG(dimension_ == DIM_1D, "Wrong dimension.");
   return ar[i];
}

// Operator() 2D
inline double& Array::operator ()(int i,int j)
{
   //TODO
   //static double dummy;
	CHECK_MSG(i >= 0 && i < xSize_, "Index i out of bounds");
	CHECK_MSG(j >= 0 && j < ySize_, "Index j out of bounds");
	CHECK_MSG(dimension_ == DIM_2D, "Wrong dimension.");
   return ar[xSize_*j + i];
}

// Operator() 3D
inline double& Array::operator ()(int i, int j, int k)
{
   //TODO
   //static double dummy;
	CHECK_MSG(i >= 0 && i < xSize_, "Index i out of bounds");
	CHECK_MSG(j >= 0 && j < ySize_, "Index j out of bounds");
	CHECK_MSG(k >= 0 && k < zSize_, "Index k out of bounds");
	CHECK_MSG(dimension_ == DIM_3D, "Wrong dimension.");
   return ar[(xSize_ * ySize_)*k + xSize_*j + i];
}

///////////
// CONST //
///////////

// Operator() 1D
inline const double& Array::operator ()(int i) const
{
   //TODO
   //static double dummy;
	CHECK_MSG(i >= 0 && i < size, "Index i out of bounds");
	CHECK_MSG(dimension_ == DIM_1D, "Wrong dimension.");
   return ar[i];
}

// Operator() 2D
inline const double& Array::operator ()(int i,int j) const
{
   //TODO
   //static double dummy;
	CHECK_MSG(i >= 0 && i < xSize_, "Index i out of bounds");
	CHECK_MSG(j >= 0 && j < ySize_, "Index j out of bounds");
	CHECK_MSG(dimension_ == DIM_2D, "Wrong dimension.");
   return ar[xSize_*j + i];
}

// Operator() 3D
inline const double& Array::operator ()(int i, int j, int k) const
{
   //TODO
   //static double dummy;
	CHECK_MSG(i >= 0 && i < xSize_, "Index i out of bounds");
	CHECK_MSG(j >= 0 && j < ySize_, "Index j out of bounds");
	CHECK_MSG(k >= 0 && k < zSize_, "Index k out of bounds");
	CHECK_MSG(dimension_ == DIM_3D, "Wrong dimension.");
   return ar[(xSize_ * ySize_)*k + xSize_*j + i];
}



#endif //ARRAY_HH

