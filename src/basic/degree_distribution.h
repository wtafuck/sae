#include "storage/mgraph.h"
#include "solver/solver.h"

class Degree_Distribution:public sae::Solver<std::vector<std::pair<int, double>>>
{
public:
	Degree_Distribution(sae::io::MappedGraph *graph);
	~Degree_Distribution();
	std::vector<std::pair<int, double>> solve();
};
