#include "storage/mgraph.h"
#include "solver/solver.h"
#include "../statistic_analysis/statistic_analysis.h"

class Propensity_Score_Matching:public sae::Solver<double> {
public:
    Propensity_Score_Matching(sae::io::MappedGraph *graph);
    ~Propensity_Score_Matching();
    std::vector<int> solve(double);
    double compare(int);

private:
};
