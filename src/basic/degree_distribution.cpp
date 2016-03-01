#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "degree_distribution.h"
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <cmath>

using namespace std;
using namespace sae::io;

Degree_Distribution::Degree_Distribution(MappedGraph *graph):
Solver(graph)
{}

Degree_Distribution::~Degree_Distribution(){}

bool equal(double p , double q ){
	return fabs(p - q) < 1e-5;
}

// return a vector with <degree, percentage>
vector<pair<int, double>> Degree_Distribution::solve(){
	set<pair<vid_t , vid_t>> edges;
	map<vid_t , size_t> node_flag;
	vector<double> node_degree;
	for (auto iter = graph->Edges(); iter->Alive(); iter->Next()){
	     vid_t a = iter -> Source() -> GlobalId();
	     vid_t b = iter -> Target() -> GlobalId();
	     if (a == b) continue;
	     if (a > b) swap(a , b);
	     if (edges.find(make_pair(a , b)) != edges.end()) continue;
	     edges.insert(make_pair(a , b));
	     if (node_flag.find(a) == node_flag.end()) {
			node_flag[a] = node_degree.size();
			node_degree.push_back(1.0);
		}	else node_degree[node_flag[a]] += 1.0;
	     if (node_flag.find(b) == node_flag.end()) {
			node_flag[b] = node_degree.size();
			node_degree.push_back(1.0);
		}	else node_degree[node_flag[b]] += 1.0;
	}
    // first: degree, second: percentage
	vector<pair<int, double> > Ans; 
    sort(node_degree.begin() , node_degree.end() );
	for (size_t start = 0; start < node_degree.size(); ){
		size_t end = start + 1;
		for ( ; end < node_degree.size() && equal(node_degree[end],node_degree[start]) ; ++ end);
		Ans.push_back(make_pair((int) node_degree[start] ,1.0 * (end - start) / node_degree.size() ));
		start = end;
	}
	return Ans;
}
