#include "storage/mgraph.h"
#include "solver/solver.h"

class Triangle_Count:public sae::Solver<int> {
public:
	Triangle_Count(sae::io::MappedGraph *graph);
	~Triangle_Count();
	int solve();
    int solve(double k);
private:
};
