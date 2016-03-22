#include "storage/mgraph.h"
#include "solver/solver.h"

class Propensity_Score_Matching:public sae::Solver<double> {
public:
    Propensity_Score_Matching(sae::io::MappedGraph *graph);
    ~Propensity_Score_Matching();
    double solve(int influence_neighbours);
    std::vector<double> logistic_regression(const std::vector< std::vector<double> >& data
        , const std::vector<bool> is_treated);
private:
};
