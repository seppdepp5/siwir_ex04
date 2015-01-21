#include "CGSolver.hh"

class HeatSolver
{
public:

	HeatSolver(CGSolver & cgSolver);

	void initU();
	void updateF(double alpha, double k, double dt);
	void solve(double alpha, double k, double dt, int timesteps);

private:

	CGSolver cgSolver_;


};
