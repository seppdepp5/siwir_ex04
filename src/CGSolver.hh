#ifndef CGSOLVER_HH
#define CGSOLVER_HH

#include "Array.hh"

class CGSolver
{
public:

	// constructors
	CGSolver
		(const int nx
		,const int ny
		,const double k
		,const int maxIter
		,const double eps	// for residual comparison
		);

	// yay
	void solve();

	int saveToFile(std::string filename, Array & u) const;

	Array & getU() {return u_;}
	Array & getF() {return f_;}
	Array & getR() {return r_;}



private:

	// does target = A*u where A is the system matrix of our problem
	// assuming u does NOT contain a boundary layer (thus this is plain matrix-vector multiplication)
	void applyOperator(Array & u, Array & target);

	// computes f - A*u and stores the result in r_
	void computeResidual();



	Array u_;
	Array f_;
	Array r_;

	int nx_;
	int ny_;

	double dx_;
	double dy_;
	double k_;

	int maxIter_;
	double eps_;



};

#endif
