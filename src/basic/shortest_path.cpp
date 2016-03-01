#include "shortest_path.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <queue>
using namespace std;
using namespace sae::io;
struct PV
{
    vid_t point;
    double value;
    PV(){}
    PV(vid_t x,double y):point(x),value(y){}
    friend bool operator < (const PV & a, const PV & b)
    {return a.value > b.value;}
};


Shortest_Path::Shortest_Path(MappedGraph *graph)
    :Solver(graph) {
}

Shortest_Path::~Shortest_Path() {
}

double Shortest_Path::solve(vid_t start , vid_t end, vector< EdgeIteratorPtr >& path)
{
    int n = graph->VertexCount();
    vector<double> dist(n,1e9);
    vector<EdgeIteratorPtr> last(n);
    priority_queue <PV> h;  //heap
    dist[start] = 0;
    h.push(PV(start,0));
    while(!h.empty()){
        PV k = h.top();
        h.pop();
        if(k.value == dist[k.point]){
            if(k.point== end) break;
            auto v = graph -> Vertices() ;
            v-> MoveTo(k.point);
            for(auto iter = v->OutEdges(); iter -> Alive(); iter-> Next()){
                int j = iter -> TargetId();
                int edgeValue = double(sae::serialization::convert_from_string<int>(iter->Data()));
                if(k.value + edgeValue < dist[j])
                {
                    dist[j] = k.value + edgeValue;
                    last[j] = iter->Clone();
                    h.push( PV(j,dist[j]));
                }
            }
        }
    }
    path.clear();
    if(dist[end] > 1e9 -1)
        return -1;
    vid_t now = end;
    while( now != start ){
        int tmp = last[now] -> SourceId();
        path.push_back(std::move(last[now]));
        now = tmp;
    }
    reverse(path.begin(), path.end());
    return dist[end];
}

double Shortest_Path::solve(vid_t start, vid_t end)
{
    vector< EdgeIteratorPtr > path;
    return solve(start,end,path);
}