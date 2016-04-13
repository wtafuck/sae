#include <cstdio>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <set>
#include <utility>
#include <fstream>
#include "dynamicMinimumSpanningTree.h"

#define INFI 2147483647

using namespace std;

lctNode *null;
int tot;
LL mstValue;
vector<lctNode*> tree,edges;

dynamicMinimumSpanningTree::dynamicMinimumSpanningTree(string file_path)
:SolverForStreaming(file_path)
{}

dynamicMinimumSpanningTree::~dynamicMinimumSpanningTree()
{}

lctNode::lctNode(int edgeNum,int _weight,lctNode* _node1,lctNode* _node2)
{
	weight=_weight;
    fa = ls = rs = null;
    num = maxEdge = edgeNum;
    node1=_node1;node2=_node2;
    revMark = 0;
}

void lctNode::Reverse()
{
    revMark ^= 1;
    swap(ls , rs);
}

void lctNode::PushUp()
{
    int maxEdgeWeight = weight,o;
    maxEdge=num;
    o=edges[ls->maxEdge]->weight;
    if (o > maxEdgeWeight)
    {
        maxEdgeWeight=o;
        maxEdge = ls->maxEdge;
    }
    o=edges[rs->maxEdge]->weight;
    if (o > maxEdgeWeight)
    {
        maxEdgeWeight=o;
        maxEdge = rs->maxEdge;
    }
}

void lctNode::PushDown()
{
    if (fa->ls==this || fa->rs==this)
        fa->PushDown();
    if (!revMark) return;
    ls->Reverse();
    rs->Reverse();
    revMark=0;
}

void dynamicMinimumSpanningTree::Zig(lctNode *x)
{
    lctNode *y = x->fa;
	y->ls = x->rs;
	x->rs->fa = y;
	x->rs = y;
	x->fa = y->fa;
	if (y == y->fa->ls)
		y->fa->ls = x;
	else if (y == y->fa->rs)
		y->fa->rs = x;
	y->fa = x;
	y->PushUp();
}

void dynamicMinimumSpanningTree::Zag(lctNode *x)
{
    lctNode *y = x->fa;
	y->rs = x->ls;
	x->ls->fa = y;
	x->ls = y;
	x->fa = y->fa;
	if (y == y->fa->ls)
		y->fa->ls = x;
	else if (y == y->fa->rs)
		y->fa->rs = x;
	y->fa = x;
	y->PushUp();
}

void dynamicMinimumSpanningTree::Splay(lctNode *x)
{
	x->PushDown();
	while (x->fa->ls == x || x->fa->rs == x)
	{
		lctNode *y = x->fa;
		lctNode *z = y->fa;
		if (x == y->ls)
		{
			if (y == z->ls) Zig(y);
			Zig(x);
		}
		else
		{
			if (y == z->rs)	Zag(y);
			Zag(x);
		}
	}
	x->PushUp();
}


void dynamicMinimumSpanningTree::Access(lctNode *x)
{
	lctNode *y = null;
	while (x != null)
	{
		Splay(x);
		x->rs = y;
		x->PushUp();
		y = x;
		x = x->fa;
	}
}

lctNode* dynamicMinimumSpanningTree::findRoot(lctNode *x)
{
    while (x->fa != null)
        x = x->fa;
    return x;
}

void dynamicMinimumSpanningTree::moveToRoot(lctNode *x)
{
    Access(x);
    Splay(x);
    x->Reverse();
}

void dynamicMinimumSpanningTree::Link(lctNode *x , lctNode *y)
{
    moveToRoot(x);
    x->fa = y;
}

void dynamicMinimumSpanningTree::Cut(lctNode *x , lctNode *y)
{
    moveToRoot(x);
    Access(y);
    Splay(y);
    x->fa=y->ls=null;
    y->PushUp();
}

int dynamicMinimumSpanningTree::Query(lctNode *x , lctNode *y)
{
    moveToRoot(x);
    Access(y);
    Splay(y);
    return y->maxEdge;
}

int dynamicMinimumSpanningTree::Insert(lctNode *x , lctNode *y , int weight)
{
	int temp=-1;
	if (x == y)
		return -1;
	if (findRoot(x) == findRoot(y))
	{
		temp = Query(x, y);
		if (edges[temp]->weight <= weight) return -1;
		Cut(edges[temp], edges[temp]->node1);
		Cut(edges[temp], edges[temp]->node2);
        mstValue -= edges[temp]->weight;
	}
	if (temp==-1) temp=++tot;
	edges[temp]=new lctNode(temp,weight,x,y);
	Link(x, edges[temp]);
	Link(y, edges[temp]);
    mstValue += weight;
    return temp;
}

/*vector<pair<pair<vid_t, vid_t>,int> > dynamicMinimumSpanningTree::solve()
{
	vector<pair<pair<vid_t, vid_t>,int> > ans;
	ans.clear();tot=0;
	tree.clear();edges.clear();
	null=new lctNode(0,-1,0,0);
    null->fa=null->ls=null->rs=null->node1=null->node2=null;
    edges.push_back(null);
    vid_t i,n=graph->VertexCount(),x,y;
    int weight;
    mstValue = 0;
    tree.push_back(null);
    for (i = 1; i <= n; ++ i)
    {
        tree.push_back(new lctNode(0,-1,null,null));
        edges.push_back(new lctNode(0,-1,null,null));
        ans.push_back(make_pair(make_pair(-1,-1),-1));
    }
    for (auto v=graph->Vertices();v->Alive();v->Next())
    {
    	x=v->GlobalId();++x;
        for(auto iter=v->OutEdges();iter->Alive();iter->Next())
        {
        	y=iter->TargetId();++y;
        	weight=int(sae::serialization::convert_from_string<int>(iter->Data()));
        	int temp=Insert(tree[x],tree[y],weight);
        	if (temp>0)
        		ans[temp]=make_pair(make_pair(x-1,y-1),weight);
        }
    }
    for (i=0;i<ans.size();++i)
    	if (ans[i].first.first==-1)
    		ans.erase(ans.begin()+i),--i;
    cout << "mstValue:" << mstValue << endl;
    return ans;
}*/

resultMST* dynamicMinimumSpanningTree::solve()
{
	freopen((file_path).c_str(),"r",stdin);
	vid_t n,m;
	scanf("%llu%llu",&n,&m);
	cout<<n<<' '<<m<<endl;
	resultMST*ans=new resultMST;
	tot=0;ans->n=n;ans->m=m;
	tree.clear();edges.clear();
	null=new lctNode(0,-1,0,0);
    null->fa=null->ls=null->rs=null->node1=null->node2=null;
    edges.push_back(null);
    vid_t i,x,y;
    int weight;
    mstValue = 0;
    tree.push_back(null);
    for (i = 1; i <= n; ++ i)
    {
        tree.push_back(new lctNode(0,-1,null,null));
        edges.push_back(new lctNode(0,-1,null,null));
        ans->edge.push_back(make_pair(make_pair(-1,-1),-1));
    }
    for (i=1;i<=m;++i)
    {
    	scanf("%llu%llu%llu",&x,&y,&weight);
    	++x;++y;
    	if (i%1000000==0) cout<<i<<endl;
        int temp=Insert(tree[x],tree[y],weight);
        if (temp>0)
        	ans->edge[temp]=make_pair(make_pair(x-1,y-1),weight);
    }
    ans->mstValue=mstValue;
    return ans;
}
