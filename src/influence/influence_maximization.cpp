#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <queue>
#include "influence_maximization.h"
#include <algorithm>
#include <set>
#include <map>
#include <utility>
using namespace std;
using namespace sae::io;

typedef pair<sae::io::vid_t , double> Pair;

Influence_Maximization::Influence_Maximization(MappedGraph *graph):Solver(graph){
}

Influence_Maximization::~Influence_Maximization() {
}

//  parameters: 
//  K:  size of seed size
//  R:  number of simulations
pair<vector<vid_t>, double> Influence_Maximization::solve(int K, int R, double W) {
    priority_queue<pair<double, pair<int, int> > > Q;
    vector<vid_t> seed_list;
    active_flag = new bool[graph -> VertexCount()];
    rand_weight = new double[graph -> VertexCount()];
    for (unsigned int i = 0; i < graph -> VertexCount(); i ++) {
        rand_weight[i] = (double) rand() / RAND_MAX;
    }
    for (auto iter = graph -> Vertices(); iter -> Alive(); iter -> Next()) {
        vid_t v = iter -> GlobalId(); 
        pair<int, int> tp = make_pair(v, 1);
		double act_num = simulation(v, R, W);
		Q.push(make_pair(act_num, tp));
    }
    double total_num = 0;
    for (int i = 1; i <= K; i ++) {
        while (1) {
			pair<double, pair<int, int> > tp = Q.top();
			Q.pop();
			if (tp.second.second == i) {
				total_num += tp.first;
                seed_list.push_back(tp.second.first);
				break;
			} else {
				double act_num = simulation(tp.second.first, R, W);
				tp.second.second = i;
				tp.first = act_num - total_num;
				Q.push(tp);
			}
		}
    }
    delete rand_weight;
    delete active_flag;
    return make_pair(seed_list, total_num);
}

double  Influence_Maximization::simulation(vid_t seed, int R, double W) {
    double res = 0.0;
    vector<vid_t> active_list;
    for (int i = 0; i < R; i ++) 
    {
        for (unsigned int j = 0; j < graph -> VertexCount(); j ++)
            active_flag[j] = false;
        active_flag[seed] = true;
        active_list.clear();
        active_list.push_back(seed);
        
        for (unsigned int j = 0; j < active_list.size(); j ++) {
            vid_t v = active_list[j];
            auto vertex = graph -> Vertices();
            vertex -> MoveTo(v);
            for (auto iter = vertex -> OutEdges(); iter -> Alive(); iter -> Next()) {
                vid_t target = iter -> TargetId();
                double prob = W;
                if (W == 0)
                    prob = rand_weight[v] * rand_weight[target];
                if (W == -1)
                    prob = 1.0 / (vertex -> OutEdgeCount() + 0.0);
                if (! active_flag[target] && (double) rand() / RAND_MAX <= prob) {
                    active_flag[target] = true;
                    active_list.push_back(target);
                }
            }
        }
        res += active_list.size() - 1;
    }
    res /= R + 0.0;
    return res;
    
}

