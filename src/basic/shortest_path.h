#include "storage/mgraph.h"
#include "solver/solver.h"

class Shortest_Path:public sae::Solver<double> {
public:
    Shortest_Path(sae::io::MappedGraph *graph);
    ~Shortest_Path();
    double solve(sae::io::vid_t, sae::io::vid_t);
    double solve(sae::io::vid_t, sae::io::vid_t, std::vector<sae::io::EdgeIteratorPtr>&);
private:
};
