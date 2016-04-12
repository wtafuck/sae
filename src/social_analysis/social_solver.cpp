#include "social_solver.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <queue>
using namespace std;
using namespace sae::io;

double Social_Solver::Density(){
    return (double)graph->EdgeCount() / (graph->VertexCount() * (graph->VertexCount() - 1));
}

std::vector< std::vector<int> > Social_Solver::Components()
{
    vector< vector<int> > ret;
    int n = graph->VertexCount();
    vector<bool> visited(n, false);
    vector<int> component;
    auto viter = graph->Vertices();
    
    for(int i = 0; i < n; i++)
    if(!visited[i])
    {
        component.clear();
        component.push_back(i);
        visited[i] = true;
        int head = -1;
        while(head + 1< component.size()){
            viter->MoveTo(component[++head]);
            for(auto eiter = viter->OutEdges(); eiter->Alive(); eiter->Next())
            {
                int target = eiter-> TargetId();
                if(!visited[target]){
                    visited[target] = true;
                    component.push_back(target);
                }
            }
        }
        ret.push_back(component);
    }
    return ret;
}

vector<int> Social_Solver::Degree_Centrality()
{
    vector<int> ids,outdegrees;
    auto viter = graph->Vertices();
    for(int i = 0; i < graph->VertexCount();i++){
        ids.push_back(i);
        viter -> MoveTo(i);
        outdegrees.push_back(viter->OutEdgeCount());
    }
    sort(ids.begin(), ids.end(), [&outdegrees](int x,int y){ return outdegrees[x] > outdegrees[y]; });
    // return ids in social_main
    return outdegrees;
}

std::vector<pair<double, int> > Social_Solver::Unweighted_Betweenness()
{
    vector<pair<double, int> > ret;
    int n = graph->VertexCount();
    for(int i = 0;i < n; i++)
        ret.push_back(make_pair(0,i));
    for(int s = 0; s < n; s++)
    {
        vector< int > stack;
        vector< int > path_count(n,0), dis(n,-1);
        vector< vector<int> > p(n);// preceding vertices
        dis[s] = 0;
        path_count[s] = 1;
        stack.push_back(s);
        int head = -1;
        auto viter = graph->Vertices();
        while(head + 1 < stack.size())
        {
            viter->MoveTo(stack[++head]);
            for(auto eiter = viter->OutEdges(); eiter->Alive(); eiter->Next())
            {
                int v = eiter->SourceId(), w = eiter->TargetId();
                if(dis[w] < 0){
                    stack.push_back(w);
                    dis[w] = dis[v] + 1;
                }
                if(dis[w] == dis[v] + 1){
                    path_count[w] += path_count[v];
                    p[w].push_back(v);
                }
            }
        }
        vector<double> part_value(n,0);
        for(int i = stack.size()-1; i >= 0; i--)
        {
            int w = stack[i];
            for(int j = 0;j < p[w].size(); j++)
                part_value[p[w][j]] += (double)path_count[p[w][j]] / path_count[w] * (1 + part_value[w]);
            if(w != s) ret[w].first += part_value[w];
        }
    }
    sort(ret.begin(), ret.end());
    return ret;
}

vector<double> Social_Solver::Effective_Size()
{
    std::vector<double> ret;
    auto viter = graph->Vertices();
    for(int i = 0; i < graph->VertexCount(); i++)
    {
        viter->MoveTo(i);
        ret.push_back(graph->VertexCount() 
            - ((double)(graph->EdgeCount() - viter->InEdgeCount() - viter->OutEdgeCount())) / graph->VertexCount());
    }
    return ret;
}

/* K - Core */

struct Node{
    int deg, id;
    Node* pre, *next;
    Node* Delete(){
        pre -> next = this -> next;
        next -> pre = this -> pre;
        return this;
    }
    Node* Insert(Node * behind){
        Node* before = behind -> pre;
        before -> next = this;
        this -> pre = before;
        behind -> pre = this;
        this -> next = behind;
        return this;
    }
    Node(int _deg = -1):deg(_deg), id(-1), pre(NULL), next(NULL){}
};

std::vector<int> Social_Solver::K_Core()
{
    int n = graph->VertexCount();
    vector< Node > sentinel(n * 2), index(n);
    for(int i = 0;i < n;i++)
    {
        sentinel[i].next = &sentinel[i + n];
        sentinel[i + n].pre = &sentinel[i];
    }
    auto viter = graph->Vertices();
    for(int i = 0;i < n;i++)
    {
        viter->MoveTo(i);
        assert(viter->InEdgeCount() == viter->OutEdgeCount());// must be non-direction graph
        assert(viter->InEdgeCount() < n); // must be simple graph
        index[i].deg = viter->OutEdgeCount();
        index[i].id = i;
        index[i].Insert(&sentinel[index[i].deg + n]);
    }
    std::vector<int> ret(n,-1);
    for(int core = 0; core < n; core++){
        for( Node* iter = sentinel[core].next; iter->deg != -1; iter = iter->next){
            int x = iter->id;
            ret[x] = core;
            viter-> MoveTo(x);
            for(auto eiter = viter->OutEdges(); eiter->Alive(); eiter->Next())
            {
                int y = eiter->TargetId();
                if(index[y].deg > core){
                    index[y].deg--;
                    index[y].Delete()->Insert(&sentinel[index[y].deg + n]);
                }
            }
        }
    }
    return ret;
}
