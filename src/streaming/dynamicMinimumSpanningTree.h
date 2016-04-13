#include "solver/solverForStreaming.h"
#include <vector>
#include <string>

typedef long long LL;
typedef std::uint64_t vid_t;
typedef std::uint64_t eid_t;

class resultMST
{
	public:
		LL mstValue;
		vid_t n;
		eid_t m;
		std::vector<std::pair<std::pair<vid_t,vid_t>,int> > edge;
		resultMST()
		{
			n=m=mstValue=0;edge.clear();
		}
};

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

class dynamicMinimumSpanningTree:public sae::SolverForStreaming<resultMST*>
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
		dynamicMinimumSpanningTree(std::string file_path);
    	~dynamicMinimumSpanningTree();
    	resultMST* solve();  //return edge's two vertexs and its weight in MST
};
