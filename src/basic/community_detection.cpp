#include "community_detection.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <queue>
#include <map>
#include <set>
#include <cstdlib>
#include <ctime>

using namespace std;
using namespace sae::io;


Community_detection::Community_detection(MappedGraph *graph)
    :Solver(graph){
}

Community_detection::~Community_detection() {
}


vector<pair<double, int> > compute_betweenness(MappedGraph *graph,vector< bool> is_deleted)
{
    vector<pair<double, int> > edge_betweenness;
    unsigned int n = graph->VertexCount();
     for(auto eiter = graph->Edges(); eiter->Alive(); eiter->Next())
            edge_betweenness.push_back(make_pair(0,eiter->GlobalId()));
    for(unsigned int s = 0; s < n; s++)
    {
        vector< vector<eid_t> > pre_vertex(n);
        vector< int >  search_queue;
        vector< int > point_stack;
        vector< int > path_count(n,0);
        vector< int > dis(n,-1);
        vector<double> part_value(n,0);
        dis[s] = 0;
        path_count[s] = 1;
        search_queue.push_back(s);
        auto viter = graph->Vertices();
        unsigned int head=-1;
        while(head+1<search_queue.size())
        {
            int v=search_queue[++head];
            point_stack.push_back(v);
            viter->MoveTo(v);
            for(auto eiter = viter->OutEdges(); eiter->Alive(); eiter->Next())
            {
                if(is_deleted[eiter->GlobalId()]) continue;
                int w = eiter->TargetId();
                if(dis[w] < 0)
                {
                    search_queue.push_back(w);
                    dis[w] = dis[v] + 1;
                }
                if(dis[w] == dis[v] + 1)
                {
                    path_count[w] += path_count[v];
                    pre_vertex[w].push_back(v);
                }
            }
        }

        for(int i = point_stack.size()-1; i >= 0; i--)
        {
            unsigned int w = point_stack[i];
            if(w==s) continue;
            for(unsigned int j = 0;j < pre_vertex[w].size(); j++)
            {
                unsigned int wp=pre_vertex[w][j];
                part_value[wp] += (double)path_count[wp] / path_count[w] * (1 + part_value[w]);
                if(wp==s||w==s) continue;
                auto vex = graph->Vertices();
                vex->MoveTo(wp);
                for(auto eiter = vex->OutEdges(); eiter->Alive(); eiter->Next())
                {
                    if(is_deleted[eiter->GlobalId()]) continue;
                    if(eiter->TargetId()==w)
                    {
                        edge_betweenness[eiter->GlobalId()].first +=(double)path_count[wp] / path_count[w] * part_value[w];
                        break;
                    }
                }
            }
        }
    }
    sort(edge_betweenness.begin(), edge_betweenness.end());
    reverse(edge_betweenness.begin(), edge_betweenness.end());
    return edge_betweenness;
}

pair<vector<vid_t>,int >  split_community(MappedGraph *graph,vector< bool> is_deleted)
{
    unsigned int n=graph->VertexCount();
    vector <vid_t> community(n,0);
    int k=0,s=0;
    bool is_success=false;
    while(!is_success)
    {
        is_success=true;
        for (unsigned int i=0;i<n;i++)
            if(community[i]==0) {is_success=false;s=i; k++;break;}
         if(is_success) break;
         queue < int >  search_queue;
         search_queue.push(s);
        auto viter = graph->Vertices();
        while(!search_queue.empty())
        {
                int v=search_queue.front();
                search_queue.pop();
                community[v]=k;
                viter->MoveTo(v);
                for(auto eiter = viter->OutEdges(); eiter->Alive(); eiter->Next())
                {
                    if(is_deleted[eiter->GlobalId()]) continue;
                    int w = eiter->TargetId();
                    if(community[w]== 0)
                    {
                        search_queue.push(w);
                        community[w]=k;
                    }
                }
        }
    }
    return make_pair(community,k);
}


double  Community_detection::compute_modularity(MappedGraph *graph, vector<eid_t> community,int k)
{
    int m = graph->EdgeCount();
    vector< vector<eid_t> > community_set(k);
    for (unsigned int i=0;i<community.size();i++)
            community_set[community[i]-1].push_back(i);
    double ans=0,eii=0,ai=0;
    auto viter = graph->Vertices();
    for (int i=0;i<k;i++)
    {
        double eii_temp=0,ai_temp=0;
        for(unsigned int j=0;j<community_set[i].size();j++)
        {
                viter->MoveTo(community_set[i][j]);
                for(auto eiter = viter->OutEdges(); eiter->Alive(); eiter->Next())
                {
                        if (find(community_set[i].begin(), community_set[i].end(), eiter->TargetId()) != community_set[i].end())
                            eii_temp+=1;
                        ai_temp+=1;
                }
        }
        eii=eii_temp/m;
        ai=ai_temp/m;
        ans+=eii-ai*ai;
    }
    return ans;
}


pair<vector<vid_t>,double>   Community_detection::run_Girvan_NewMan(MappedGraph *graph,int K)
{
    int m=graph->EdgeCount();
    int n= graph->VertexCount();
    vector< vid_t> community(n,1);
    vector< bool> is_deleted(m,false);
    vector<double> modularity(m,0);
    int k=1,j=0;
    for (;j<m;j++)
    {
        vector<pair<double, int> >edge_betweenness =compute_betweenness(graph,is_deleted);
        for (unsigned int c=0;c<edge_betweenness.size();c++)
            if(edge_betweenness[c].first==edge_betweenness[0].first)
            {
                int first =edge_betweenness[c].second;
                int second=edge_betweenness[c].second%2==0?edge_betweenness[c].second+1:edge_betweenness[c].second-1;
                is_deleted[first]=true;
                is_deleted[second]=true;
            }

        if(edge_betweenness[0].first==0)  break;
        modularity[j]=modularity[j-1];
        pair<vector<vid_t>,int > split_com=split_community(graph,is_deleted);
        community=split_com.first;
        k=split_com.second;
        modularity[j]=compute_modularity(graph,community,k);
        if(k>=K) break;
    }
    return make_pair(community,modularity[j]);
}

pair<vector<vid_t>,double>   Community_detection::run_label_propagation(MappedGraph *graph)
{
        cout<<"\tRun community detection algorithm with lable propagation"<<endl;
        srand(time(NULL));
        unsigned int  n=graph->VertexCount();
        vector<vector<vid_t>> neighbor(n);
        vector<vid_t> C(n);
        vector<vid_t> temp(n);
        vector<vid_t> max_community(n);
        double max_modularity=0;
        for (unsigned int i=0;i<n;i++)
        {
                C[i]=i;temp[i]=i;
                auto viter = graph->Vertices();
                viter->MoveTo(i);
                for(auto eiter = viter->OutEdges(); eiter->Alive(); eiter->Next())
                    neighbor[i].push_back(eiter->TargetId());
        }
        bool is_success=false;
        while(!is_success)
        {
                is_success=true;
                vector<vid_t> rand_seq=temp;
                for (int j =n-1;j>=0;j--)
                {
                        int x = rand() % (n - 1);
                        unsigned tem =rand_seq[j];
                        rand_seq[j]=rand_seq[x];
                        rand_seq[x]=tem;
                }
                for(unsigned int j=0;j<n;j++)
                {
                        unsigned v=rand_seq[j];
                        vector <vid_t> vertex_neigh=neighbor[v];
                        for (unsigned int k=0;k<vertex_neigh.size();k++)
                            vertex_neigh[k]=C[vertex_neigh[k]];
                        sort(vertex_neigh.begin(),vertex_neigh.end());
                        int max_count=1,new_count=1,max_vertex=vertex_neigh[0],new_vertex=vertex_neigh[0];
                        for (unsigned int k=1;k<vertex_neigh.size();k++)
                        {
                                if (vertex_neigh[k]==new_vertex) new_count++;
                                if(vertex_neigh[k]!=max_vertex||k==vertex_neigh.size()-1)
                                {
                                    if(new_count>max_count)
                                    {
                                            max_count=new_count;
                                            max_vertex=new_vertex;
                                    }
                                    new_vertex=vertex_neigh[k];
                                    new_count=1;
                                }
                        }
                        if(max_count==1) max_vertex=vertex_neigh[rand()%vertex_neigh.size()];
                        if (C[v]!=max_vertex)
                        {
                                C[v]=max_vertex;
                                is_success=false;
                        }
                }
            set<vid_t> community;
            set<vid_t>::reverse_iterator rit;
            for(unsigned int i=0;i<n;i++)
                community.insert(C[i]);
            map<vid_t,vid_t> map_seq;
            int seq=1,community_count=community.size();
            for (rit = community.rbegin(); rit != community.rend(); rit++)
                    map_seq[*rit]=seq++;
            vector <vid_t> final_community=C;
            for(unsigned int i=0;i<C.size();i++)
                    final_community[i]=map_seq[C[i]];
            double modularity=compute_modularity(graph,final_community,community_count);
            if (modularity>max_modularity) {
                max_modularity=modularity;max_community=final_community;
            }
        }
            return make_pair(max_community,max_modularity);
}


pair<vector<vid_t>,double> Community_detection::solve(int K)
{

    //pair<vector<vid_t>,double> ans=run_label_propagation(graph);
    pair<vector<vid_t>,double> ans=run_Girvan_NewMan(graph,K);
    return ans;
}
