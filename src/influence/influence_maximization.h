#include "storage/mgraph.h"
#include "solver/solver.h"

class Influence_Maximization:public sae::Solver<std::pair<std::vector<sae::io::vid_t>, double>> {
public:
	Influence_Maximization(sae::io::MappedGraph *graph);
	~Influence_Maximization();
	std::pair<std::vector<sae::io::vid_t>, double> solve(int K, int R, double W);
    bool* active_flag;
private:
    double* rand_weight;
    double simulation(sae::io::vid_t seed, int R, double W); 
};

