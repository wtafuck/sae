#include "storage/mgraph.h"
#include "solver/solver.h"

class K_Core_Decomposition:public sae::Solver<std::vector<std::pair<sae::io::vid_t, sae::io::vid_t> >> {
public:
    K_Core_Decomposition(sae::io::MappedGraph *graph);
    ~K_Core_Decomposition();
    std::vector<std::pair<sae::io::vid_t, sae::io::vid_t> > solve();  //return ID and its k-shell
private:
};
