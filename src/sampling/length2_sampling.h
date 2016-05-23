
#include "storage/mgraph.h"
#include "solver/solver.h"
#include <cstring>

class length2_sampling//:public sae::Solver<int>
{
public:
    length2_sampling(std::string path);
    ~length2_sampling();
    double solve(double p , double q);
    double solve(double q);

	std::string filePath;
};
