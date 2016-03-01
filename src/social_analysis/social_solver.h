#include "storage/mgraph.h"
#include "solver/solver.h"
#include <utility>
#include <cassert>
class Social_Solver: public sae::Solver<void>{
    void solver(){};
public:
    Social_Solver(sae::io::MappedGraph *graph): Solver(graph){}
    ~Social_Solver(){};
    double Density();
    std::vector< std::vector<int> > Components();
    std::vector<int> Degree_Centrality();
    std::vector<std::pair<double, int> > Unweighted_Betweenness();
    std::vector<double> Effective_Size();
    std::vector<int> K_Core();
};