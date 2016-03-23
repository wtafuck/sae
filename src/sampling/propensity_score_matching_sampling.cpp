#include "propensity_score_matching.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <queue>
#include <cmath>
#include <cassert>
using namespace std;
using namespace sae::io;

Propensity_Score_Matching_Sampling::Propensity_Score_Matching_Sampling(MappedGraph *graph)
    :Solver(graph) {
}

Propensity_Score_Matching_Sampling::~Propensity_Score_Matching_Sampling() {
}

double Propensity_Score_Matching_Sampling::solve(int influence_neighbours)
{
    int n = graph->VertexCount(); assert(n >= 3);
    vector< vector<double> > properties;
    int cnt = 0;
    for(auto v = graph->Vertices(); v-> Alive(); v->Next())
    {
        assert(cnt ++ == v -> GlobalId());
        //Here I assume the id of vertices equal to its order number.
        properties.push_back(sae::serialization::convert_from_string< vector<double> >(v -> Data()) );
    }
    //properties[i][0], whether it adopted the product (-1/1).
    vector<double> adopted_friends(n);
    vector<bool> is_treated(n);
    // check whether the node has enough adopted friends
    for(auto v = graph->Vertices(); v-> Alive(); v->Next())
    {
        for(auto iter = v->OutEdges(); iter -> Alive(); iter-> Next())
        {
            if(properties[iter -> TargetId()][0] > 0.1 )
                adopted_friends[v -> GlobalId()] ++ ;
        }
        is_treated[v -> GlobalId()] = (adopted_friends[v -> GlobalId()] >= influence_neighbours );
    }
    // logistic regression
    vector<double> covarieties = logistic_regression(properties, is_treated), possibilities;
    // pair the treated and untreated
    int m = properties[0].size();
    cnt = 0;
    for(auto v = graph->Vertices(); v-> Alive(); v->Next())
    {
        assert(cnt ++ == v -> GlobalId());
        double value = 0;
        for(int i = 1;i < m;i++){
            value += covarieties[i] * properties[v -> GlobalId()][i];
        }
        value = 1.0 / (1 + exp(-value));// sigmoid
        possibilities.push_back(value);
    }
    vector<int> pair_node(n);
    double mean = 0, treated_num = 0, sqr_deviation = 0;
    //calculate average
    for(int i = 0;i < n; i++) if(is_treated[i])
    {
        double min_dis = 1e7;
        for(int j = 0;j< n;j++) if(!is_treated[j])
            if(min_dis > abs(possibilities[i] - possibilities[j]))
                min_dis = abs(possibilities[i] - possibilities[j]), pair_node[i] = j;
         mean += min_dis;
         treated_num ++;
    }
    mean /= treated_num;
    //calculate sqr_deviation
    for(int i = 0; i < n; i++) if(is_treated[i])
    {
        double min_dis = abs(possibilities[i] - possibilities[ pair_node[i] ]);
        if(properties[i][0] > 0)
        //cout<<"property1["<<i<<"]="<<properties[i][1]<<" property1["<<pair_node[i]<<"]="<<properties[pair_node[i]][1]<<endl;
        sqr_deviation += (min_dis - mean)*(min_dis - mean);
    }
    sqr_deviation /= treated_num;
    //delete pair_nodes which are not far enough
    vector<bool> used(n, false);
    vector<int> treated_nodes, untreated_nodes;
    for(int i = 0; i < n; i++) if(is_treated[i])
    {
        double min_dis = abs(possibilities[i] - possibilities[ pair_node[i] ]);
        if(min_dis * min_dis <= 4 * sqr_deviation &&  !used[i] && !used[pair_node[i]])
        {
            treated_nodes.push_back(i);
            untreated_nodes.push_back(pair_node[i]);
            //unique the vector, I don't know how to understand the paper's algorithm when it is overlapped.
            used[i] = used[pair_node[i]] = true;
        }
    }

    //count adopt-situation
    int n1 = 0, n2 = 0;
    for(int i = 0; i < treated_nodes.size(); i++)
        if(properties[treated_nodes[i]][0] > 0.1) n1++;
    for(int i = 0; i < untreated_nodes.size(); i++)
        if(properties[untreated_nodes[i]][0] > 0.1) n2++;
    cout<< n1 << " "<< n2 <<" sum = "<<treated_nodes.size()<<endl;
    return 1.0 * n1 / n2;

}

vector<double> Propensity_Score_Matching_Sampling::logistic_regression
    (const std::vector< std::vector<double> >& data, const std::vector<bool> is_treated)
{
    double alpha = 0.01;
    int maxCycles = 500;
    int m = data[0].size(), n = data.size();
    vector<double> weights(m), delta(m);
    for(int k = 0; k < maxCycles; k++){
        vector<double> h(n, 0);
        for(int i = 0;i<n;i++){
            for(int j = 1; j<m; j++)
                h[i] += data[i][j]* weights[j];
            h[i] = 1.0 / (1 + exp(-h[i]));
            h[i] = is_treated[i] - h[i];
        }
        for(int i = 1;i<m;i++){
            delta[i] = 0;
            for(int j = 0;j < n;j++)
                delta[i] += h[j]* data[j][i];
            weights[i] += delta[i] * alpha;
        }
    }
    //for(int i =1; i< m ;i++) cout<<weights[i] << " "<<endl;
    return weights;
}

double Propensity_Score_Matching_Sampling::compare(int influence_neighbours)
{
    int n = graph->VertexCount(); assert(n >= 3);
    vector< vector<double> > properties;
    int cnt = 0;
    for(auto v = graph->Vertices(); v-> Alive(); v->Next())
    {
        assert(cnt ++ == v -> GlobalId());
        //Here I assume the id of vertices equal to its order number.
        properties.push_back(sae::serialization::convert_from_string< vector<double> >(v -> Data()) );
    }
    //properties[i][0], whether it adopted the product (-1/1).
    vector<double> adopted_friends(n);
    vector<bool> is_treated(n);
    // check whether the node has enough adopted friends
    vector<int> treated_nodes, untreated_nodes;
    for(auto v = graph->Vertices(); v-> Alive(); v->Next())
    {
        for(auto iter = v->OutEdges(); iter -> Alive(); iter-> Next())
        {
            if(properties[iter -> TargetId()][0] > 0.1 )
                adopted_friends[v -> GlobalId()] ++ ;
        }
        if(adopted_friends[v -> GlobalId()] >= influence_neighbours )
            treated_nodes.push_back(v -> GlobalId());
        else untreated_nodes.push_back(v -> GlobalId());
    }

    //count adopt-situation
    int n1 = 0, n2 = 0;
    for(int i = 0; i < treated_nodes.size(); i++)
        if(properties[treated_nodes[i]][0] > 0.1) n1++;
    for(int i = 0; i < untreated_nodes.size(); i++)
        if(properties[untreated_nodes[i]][0] > 0.1) n2++;
    cout<< "compare"<<endl<<treated_nodes.size() << " " << n1<<endl;
    cout<< untreated_nodes.size()<<" "<<n2<<endl;
    return (1.0 * n1 / treated_nodes.size())  / (1.0 * n2 / untreated_nodes.size());

}