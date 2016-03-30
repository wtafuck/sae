#include "k_core_decomposition.h"
#include <vector>
#include <queue>
using namespace std;
using namespace sae::io;

K_Core_Decomposition::K_Core_Decomposition(MappedGraph *graph)
    :Solver(graph)
{
}

K_Core_Decomposition::~K_Core_Decomposition()
{
}

vector<pair<vid_t, vid_t> > K_Core_Decomposition::solve()
{
	vid_t n=graph->VertexCount(),total=n,k_core=-1,x,i;
	vector<queue<int> > nodes;
	vector<pair<vid_t, vid_t> > ans;
	vector<int> degree,res;
	ans.clear();
	nodes.clear();
	res.clear();
	for (i=1;i<=n;++i)
	{
		nodes.push_back(queue<int>());
		res.push_back(-1);
		degree.push_back(-1);
	}
	for (auto v=graph->Vertices();v->Alive();v->Next())
	{
		nodes[v->OutEdgeCount()].push(v->GlobalId());
		degree[v->GlobalId()]=v->OutEdgeCount();
	}
	while (total)
	{
		++k_core;
		while (!nodes[k_core].empty())
		{
			x=nodes[k_core].front();
			nodes[k_core].pop();
			if (res[x]>=0)
				continue;
			res[x]=k_core;
			auto v=graph->Vertices();
			v->MoveTo(x);
			--total;
			for(auto iter=v->OutEdges();iter->Alive();iter->Next())
			{
				int target=iter->TargetId();
				if (res[target]==-1&&degree[target]>k_core)
				{
					i=--degree[target];
					nodes[i].push(iter->TargetId());
				}
			}
		}
	}
	for (i=0;i<n-1;++i)
		if (res[i]>=0)
			printf("%d %d\n",i,res[i]),
			ans.push_back(make_pair(i,res[i]));
	printf("\n");
	return ans;
}
