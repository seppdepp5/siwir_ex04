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
	void solve(double dt, double alphaCrank, double k);
 	int saveToVtu(std::string filename, Array & u, int firstRow, int lastRow) const;
	int saveToPvtu(std::string filename, int time_steps, int size) const;
	int saveToPvd(std::string filename, int time_steps) const;
	int saveToFile(std::string filename, Array & u) const;

	Array & getU() {return u_;}
	Array & getF() {return f_;}
	Array & getR() {return r_;}

	double getDX() {return dx_;}
	double getDY() {return dy_;}


private:

	// does target = A*u where A is the system matrix of our problem
	// assuming u does NOT contain a boundary layer (thus this is plain matrix-vector multiplication)
	void applyOperator(Array & u, Array & target);

	void applyOperatorHeat(Array & u, Array & target, double dt, double alpha, double k);

	// computes f - A*u and stores the result in r_
	void computeResidual(double dt, double alphaCrank, double k);



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
