#include "storage/mgraph.h"
#include "solver/solver.h"
#include <vector>

typedef long long LL;

class lctNode
{
	public:
    	lctNode(int,int,lctNode*,lctNode*);
    	void Reverse();
    	void PushUp();
    	void PushDown();
    	lctNode *fa,*ls,*rs,*node1,*node2;
    	int num, maxEdge,weight;
    	bool revMark;
};

class dynamicMinimumSpanningTree:public sae::Solver<std::vector<std::pair<std::pair<sae::io::vid_t, sae::io::vid_t>,int> >>
{
	public:
		void Zig(lctNode*);
		void Zag(lctNode*);
		void Splay(lctNode*);
		void Access(lctNode*);
		lctNode *findRoot(lctNode*);
		void moveToRoot(lctNode*);
		void Link(lctNode*,lctNode*);
		void Cut(lctNode*,lctNode*);
		int Query(lctNode*,lctNode*);
		int Insert(lctNode*,lctNode*,int);
		dynamicMinimumSpanningTree(sae::io::MappedGraph *graph);
    	~dynamicMinimumSpanningTree();
    	std::vector<std::pair<std::pair<sae::io::vid_t, sae::io::vid_t>,int> > solve();  //return edge's two vertexs and its weight in MST
};
