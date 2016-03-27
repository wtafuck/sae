#include "storage/mgraph.h"
#include "solver/solver.h"

class Community_detection_sampling:public sae::Solver<int > {
public:
    Community_detection_sampling(sae::io::MappedGraph *graph);
    ~Community_detection_sampling();
    int solve();
private:
};
