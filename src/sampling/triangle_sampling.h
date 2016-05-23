#include "storage/mgraph.h"
#include "solver/solver.h"
#include <cstring>
class Triangle_Sampling//:public sae::Solver<int> {
{
public:
	Triangle_Sampling(std::string path);
	~Triangle_Sampling();
    double solve(double p, double q);
    double solve(double q);
	std::string filePath;
private:
};

