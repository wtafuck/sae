#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>
#include "simrank.h"
#include <algorithm>

using namespace std;
using namespace sae::io;

SimRank::SimRank(MappedGraph *graph):Solver(graph){
}

SimRank::~SimRank(){
}


vector<vector<double> > simrankwithsum(vector<vector<double> > s,vector<vector<int> > indegree,double C,int K)
{
    vector<vector<double> > temp=s;
    vector<vector<double> > part=s;
    while (K--)
    {
        for (int i=0;i<s.size();i++)
            for(int j=0;j<s.size();j++)
            {
                part[i][j]=0;
                for (int k=0;k<indegree[i].size();k++)
                    part[i][j]+=s[indegree[i][k]][j];
            }

        for (int i=0;i<s.size();i++)
            for(int j=0;j<s.size();j++)
            {
                if (i==j) {
                    temp[i][j]=1;continue;
                }
                double simsum=0;
                for(int l=0;l<indegree[j].size();l++)
                    simsum+=part[i][indegree[j][l]];
                temp[i][j]=C*simsum/(indegree[i].size()*indegree[j].size());
            }

       s=temp;
    }
    return temp;
}

vector<pair<double, vid_t> > run_simrank(MappedGraph *graph,vid_t v) {
    int n=graph->VertexCount();
    vector <double> temp(n,0);
    vector<vector<double> > s(n,temp);
    for (int i=0;i<n;i++) s[i][i]=1.0;
    vector <int> temp2;
    vector<vector<int> > indegree(n,temp2);
    auto viter = graph->Vertices();
    for (int i=0;i<n;i++)
    {
        viter->MoveTo(i);
        for(auto eiter = viter->InEdges(); eiter->Alive(); eiter->Next())
            indegree[i].push_back(eiter->SourceId());
    }
    s=simrankwithsum(s,indegree,0.8,10);
    vector<pair< double,vid_t> > similarity;
    for(int i=0;i<s.size();i++)
        similarity.push_back(make_pair(s[v][i],i));
    return similarity;
}

pair<vector<vector<vid_t>>,vector<vector<vid_t>>> generate_random_path(MappedGraph *graph,int  R,int T)
{
    int n=graph->VertexCount();
    vector<vid_t> temp(T+1,0);
    vector<vector<vid_t>> path(R,temp);
    vector<vector<vid_t>> outdegree(n);
    vector<vector<vid_t> >index(n);
    map<vid_t,vid_t> map_seq;
    map<vid_t,vid_t>::iterator it;
    auto viter = graph->Vertices();
    for(int i=0;i<n;i++)
    {
        viter->MoveTo(i);
        for(auto eiter = viter->OutEdges(); eiter->Alive(); eiter->Next())
                    outdegree[i].push_back(eiter->TargetId());
    }
    for (int i=0;i<R;i++)
    {
        int s=rand() % n;
        path[i][0]=s;
        index[s].push_back(i);
        viter->MoveTo(s);
        for(int j=1;j<=T;j++)
        {
            int m=outdegree[s].size();
            if(m==0) break;
            int  rand_num=rand() % m ;
            int v=outdegree[s][rand_num];
            s=v;
            path[i][j]=s;
            index[s].push_back(i);
            viter->MoveTo(s);
        }
        map_seq.clear();
        vector<vid_t> temp;
        for(int j=0;j<path[i].size();j++)
            map_seq[path[i][j]]=path[i][j];
        for(it=map_seq.begin();it!=map_seq.end();++it)
            temp.push_back(it->second);
        path[i]=temp;
    }
	for(int i=0;i<n;i++)
	{
        map_seq.clear();
        vector<vid_t> temp;
         for(int j=0;j<index[i].size();j++)
            map_seq[index[i][j]]=index[i][j];
        for(it=map_seq.begin();it!=map_seq.end();++it)
            temp.push_back(it->second);
        index[i]=temp;
	}
	return make_pair(path,index);
}

vector<pair<double, vid_t> > query_simrank(pair<vector<vector<vid_t>>,vector<vector<vid_t>>>data,vid_t v)
{
    vector<vector<vid_t>> path=data.first;
    vector<vector<vid_t> >index=data.second;
    vector<pair<double, vid_t>>num;
    vector<vid_t> path_list=index[v];
    int n=index.size(),R=path.size();
    for(int i=0;i<n;i++)
        num.push_back(make_pair(0.0,i));
    for (int i = 0; i < path_list.size(); ++i)
    {
            int path_number=path_list[i];
            for(int j=0;j<path[path_number].size();j++)
                if(path[path_number][j]!=v)
                    num[path[path_number][j]].first+=1.0/R;
    }
    num[v].first=1.0;
    return num;
}

vector<pair<double, vid_t> > run_simrank_rand(MappedGraph *graph,vid_t v)
{
    int T=5,m=graph->EdgeCount()/2;
    double c=0.5,eita=0.1,ep=sqrt(1.0/m);
    int R=c/(ep*ep)*(log(T*(T-1)/2)/log(2)+1+log(1/eita));
    pair<vector<vector<vid_t>>,vector<vector<vid_t>>>data=generate_random_path(graph,R,T);
    vector<pair<double, vid_t> > ans=query_simrank(data,v);
    return ans;
}

vector<pair<double, vid_t> > SimRank::solve(vid_t v,bool is_accurate) {
    vector<pair< double,vid_t> > similarity;
    if(is_accurate)
        similarity=run_simrank(graph,v);
    else
        similarity=run_simrank_rand(graph,v);
    sort(similarity.begin(), similarity.end());
    reverse(similarity.begin(), similarity.end());
    return similarity;
}
