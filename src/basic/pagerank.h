#include "storage/mgraph.h"
#include "solver/solver.h"

class PageRank:public sae::Solver<std::vector<std::pair<sae::io::vid_t , double>>>{
public:
	PageRank(sae::io::MappedGraph *graph);
	~PageRank();
	std::vector<std::pair<sae::io::vid_t , double>> solve();
};
