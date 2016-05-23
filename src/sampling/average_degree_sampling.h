#include "storage/mgraph.h"
#include "solver/solver.h"
#include "cstring"

class average_degree_sampling//:public sae::Solver<int>
{
public:
     // average_degree_sampling(sae::io::MappedGraph *graph);
       average_degree_sampling(std::string path);
       ~average_degree_sampling();
      double solve(double p , double q);
      double solve(double q);
     
     std::string filePath;	
};
