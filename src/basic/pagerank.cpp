#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>
#include "pagerank.h"
#include <algorithm>
#include <set>
#include <map>
#include <utility>
using namespace std;
using namespace sae::io;

typedef pair<sae::io::vid_t , double> Pair;

PageRank::PageRank(MappedGraph *graph):Solver(graph){
}

PageRank::~PageRank(){
}

vector<double> bf_vec;

void bf_matrix_mult(vector<double> & a , const vector<vector<double>> & b){
	bf_vec.clear();
	double q = 0.85;
	double sum;
	// for (int i = 0; i < b.size(); ++ i){
	// 	sum = (1-q)*1.0 / b.size();
	// 	for (int j = 0; j < b[i].size(); ++ j )
	// 		sum += q * b[j][i] * a[i];
	// 	bf_vec.push_back(sum);
	// }
	for(int i = 0; i < b.size(); ++ i){
	 	sum = (1-q);
		for(int j = 0; j < b.size(); ++ j){
			sum += q * a[j] * b[j][i];
		}
		bf_vec.push_back(sum);
	}
	a = bf_vec;
}

vector<double> bf_matrix_run(vector<vector<double>> & matr , int iter_num) {
	vector<double> pr_value;
	pr_value.clear();
	for (int i = 0; i < matr.size(); ++ i) pr_value.push_back(1.0);// / matr.size());
	for (; iter_num --; ) {
		bf_matrix_mult(pr_value , matr);
	}
	return pr_value;
}

bool bf_cmp(const Pair &a , const Pair &b){
	return a.second > b.second;
}

// output a vector with <vertex_id, pagerank_score> and sorted by pagerank_score
vector<Pair> PageRank::solve() {
	map<pair<vid_t , vid_t> , bool > edges;
	for (auto iter = graph->Edges(); iter -> Alive(); iter -> Next()){
		vid_t x = iter -> Source() -> GlobalId();
		vid_t y = iter -> Target() -> GlobalId();
		if (x == y) continue;
		//if (x > y) swap(x , y);
		if (edges.find(make_pair(x , y)) == edges.end() ){
			edges[make_pair(x , y)] = true;
		}
	}
	map<vid_t , size_t> node_map;
	vector<vector<vid_t>> node_edge;
	vector<vid_t> row_vec;

	for (auto iter = edges.begin(); iter != edges.end(); ++ iter){
		vid_t x = iter->first.first;
		vid_t y = iter->first.second;
		if (node_map.find(x) == node_map.end()){
			node_map[x] = node_edge.size();
			row_vec.clear();
			row_vec.push_back(y);
			node_edge.push_back(row_vec);
		}  else node_edge[node_map[x]].push_back(y);
		if (node_map.find(y) == node_map.end()){
			node_map[y] = node_edge.size();
			row_vec.clear();
			//row_vec.push_back(x);
			node_edge.push_back(row_vec);
		}  // node_edge[node_map[y]].push_back(x);
	}

	vector<vector<double> > Matr;
	vector<double> row;
	Matr.clear();
	row.clear();
	for (int i = 0; i < node_edge.size(); ++ i) row.push_back(0.0);
	for (int i = 0; i < node_edge.size(); ++ i){
		Matr.push_back(row);
		for (int j = 0; j < node_edge[i].size(); ++ j){
			vid_t x = node_edge[i][j];
			Matr[i][node_map[x]] = 1.0 / node_edge[i].size();
		}
	}
	auto value = bf_matrix_run(Matr, 200);
	vector<Pair> Ans;
	for (auto iter = node_map.begin(); iter != node_map.end(); ++ iter){
		// first is the vertex index while second is the PageRank value
        Ans.push_back(make_pair(iter->first , value[iter->second]));
	}
	sort(Ans.begin(), Ans.end(), bf_cmp);
	return Ans;
}
