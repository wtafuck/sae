#ifndef SOLVERFORSTREAMING_HPP
#define SOLVERFORSTREAMING_HPP

#include <string>

namespace sae {
namespace io {

typedef uint64_t vid_t;
typedef uint64_t eid_t;

}
}

namespace sae {
/*
 * specific algorithm applied to graph data
 */
template <class return_t>
class SolverForStreaming {
    public:
		SolverForStreaming(std::string file_path) {
			this->file_path = file_path;
		}
		~SolverForStreaming() {
			//delete graph;
		}
		virtual return_t solve() {} // main process
		virtual return_t solve(sae::io::vid_t) {} // online method (vertex)
		virtual return_t solve(sae::io::vid_t, sae::io::vid_t) {} // online method (edge)
	private:
		SolverForStreaming(const SolverForStreaming&) {}
	protected:
		std::string file_path;
};
}

#endif
