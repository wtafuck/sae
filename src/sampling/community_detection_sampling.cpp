#include "community_detection_sampling.h"
#include "../basic/community_detection.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <queue>
#include <set>
#include <cstdlib>
#include <ctime>

using namespace std;
using namespace sae::io;

Community_detection_sampling::Community_detection_sampling(MappedGraph *graph)
    :Solver(graph){
}

Community_detection_sampling::~Community_detection_sampling() {
}


MappedGraph * select_sub_graph(MappedGraph *graph,int num)
{
    GraphBuilder<int> sub_graph_builder;
    map<vid_t,bool> is_part;
    auto viter = graph->Vertices();
    for(int i=0;i<num;i++)
    {
        is_part[i]=true;
        sub_graph_builder.AddVertex(i,0);
    }
    for(int i=0;i<num;i++)
    {
        viter->MoveTo(i);
        for(auto eiter = viter->OutEdges(); eiter->Alive(); eiter->Next())
            if(is_part.find(eiter->TargetId())!=is_part.end())
                sub_graph_builder.AddEdge(i,eiter->TargetId(),0);
    }
    sub_graph_builder.Save("./data/temp");
    MappedGraph *sub_graph = MappedGraph::Open("./data/temp");
    return sub_graph;
}

vector<vid_t> allocate_vertex_with_shortest_path(MappedGraph *graph,vector<vid_t> community,int sample_num)
{
        int n=graph->VertexCount();
        vector<vid_t> communityChange(n,1);
        for(int i=0;i<sample_num;i++)
            communityChange[i]=community[i];
        for(int i=sample_num;i<n;i++)
        {
            queue < int >  search_queue;
            vector< int > dis(n,-1);
            dis[i] = 0;
            search_queue.push(i);
            auto viter = graph->Vertices();
            while(!search_queue.empty())
            {
                int v=search_queue.front();
                search_queue.pop();
                viter->MoveTo(v);
                for(auto eiter = viter->OutEdges(); eiter->Alive(); eiter->Next())
                {
                    int w = eiter->TargetId();
                    if(dis[w] < 0)
                    {
                        search_queue.push(w);
                        dis[w] = dis[v] + 1;
                    }
                }
            }
            int min_value=10000,min_index=0;
            for(int j=0;j<sample_num;j++)
                if(dis[j]>0&&dis[j]<min_value)
                {
                    min_value=dis[j];
                    min_index=j;
                }
            communityChange[i]=community[min_index];
        }
        return communityChange;
}

vector<vid_t> allocate_vertex_with_shortest_path2(MappedGraph *graph,vector<vid_t> community,int sample_num,int K)
{
        int n=graph->VertexCount();
        vector<vid_t> communityChange(n,1);
        for(int i=0;i<sample_num;i++)
            communityChange[i]=community[i];
        vector <int> t(n,-1);
        vector <vector<int>> dis(K+1,t);
        for(int k=1;k<=K;k++)
        for(int i=0;i<sample_num;i++)
        {
            queue < int >  search_queue;
            if (community[i]==k)
            {
                    dis[k][i] = 0;
                    search_queue.push(i);
            }
            auto viter = graph->Vertices();
            while(!search_queue.empty())
            {
                int v=search_queue.front();
                search_queue.pop();
                viter->MoveTo(v);
                for(auto eiter = viter->OutEdges(); eiter->Alive(); eiter->Next())
                {
                    int w = eiter->TargetId();
                    if(dis[k][w] < 0)
                    {
                        search_queue.push(w);
                        dis[k][w] = dis[k][v] + 1;
                    }
                }
            }
        }
        for(int j=sample_num;j<n;j++)
        {
                int min_value=10000,min_index=0;
                 for(int k=1;k<=K;k++)
                if(dis[k][j]>0&&dis[k][j]<min_value)
                {
                    min_value=dis[k][j];
                    min_index=k;
                }
            communityChange[j]=community[min_index];
        }
        return communityChange;

}

pair<vector<vid_t>,double> allocate_vertex_with_label_propagation(MappedGraph *graph,vector<vid_t> community,int sample_num,int K)
{
        srand(time(NULL));
        unsigned int  n=graph->VertexCount();
        vector<vector<vid_t>> neighbor(n);
        vector<vid_t> C(n,1);
        vector<vid_t> temp(n);
        vector<vid_t> max_community=C;
        double max_modularity=0;
        vid_t max_community_count=K;
        auto viter = graph->Vertices();
        for (unsigned int i=0;i<n;i++)
        {
                temp[i]=i;
                if(i<sample_num)
                    C[i]=community[i];
                else{
                    C[i]=rand()%max_community_count+1;
                    viter->MoveTo(i);
                    for(auto eiter = viter->OutEdges(); eiter->Alive(); eiter->Next())
                        neighbor[i].push_back(eiter->TargetId());
                }
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
                        if (v<sample_num) continue;
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
            Community_detection cc(graph);
            double modularity=cc.compute_modularity(graph,C,max_community_count);
            if (modularity>max_modularity) {
                max_modularity=modularity;max_community=C;
            }
        }
            return make_pair(max_community,max_modularity);
}

void recurPermutation(vector<vid_t> c,vector<vid_t> c2,double * max_overlap,vector<int> arr, int n, int i)
{
    if(i==n-1) {
        int overlapping =0;
        for (int i=0;i<c.size();i++)
            if(arr[c[i]-1]==c2[i])
                overlapping++;
        double overlap=overlapping*1.0/c.size();
        if (overlap>*max_overlap)
            *max_overlap=overlap;
    }
    for(int j=i; j<n; j++) {
        swap(arr[i], arr[j]);
        recurPermutation(c,c2,max_overlap,arr, n, i+1);
        swap(arr[i], arr[j]);
    }
}

void  test_community_sampling(MappedGraph *graph,int min_rate,int max_rate,int min_k,int max_k)
{
    Community_detection cd(graph);
    cout<<"\tRun community detection sampling algorithm"<<endl;
    double overlap[(max_rate-min_rate)/10+1][max_k-min_k+1];
    double mod1[(max_rate-min_rate)/10+1][max_k-min_k+1];
    double mod2[(max_rate-min_rate)/10+1][max_k-min_k+1];
    double t1[(max_rate-min_rate)/10+1][max_k-min_k+1];
    double t2[(max_rate-min_rate)/10+1][max_k-min_k+1];
    double rate[(max_rate-min_rate)/10+1][max_k-min_k+1];
    for (int p=min_rate;p<=max_rate;p+=10)
        for (int K=min_k;K<=max_k;K++)
        {
            cout<<"\n\n\sampling rate:"<<p<<"%"<<"\tcommunity number:"<<K<<endl<<endl;
            time_t start_time = clock();
            int n=graph->VertexCount();
            int sample_num=n*p/100;
            pair<vector<vid_t>,double> ans=cd.run_Girvan_NewMan(graph,K);
            time_t end_time = clock();
            t1[int(p-min_rate)/10][K-min_k]=(end_time- start_time+ 0.0) / CLOCKS_PER_SEC ;;
            cout<<endl<<"\tmodularity is "<<ans.second<<endl;
            time_t start_time2 = clock();
            MappedGraph * sub_graph= select_sub_graph(graph,sample_num);
            pair<vector<vid_t>,double> ans2=cd.run_Girvan_NewMan(sub_graph,K);
            cout<<endl<<"\tmodularity is "<<ans2.second<<endl;
            time_t end_time2 = clock();
            t2[int(p-min_rate)/10][K-min_k]=(end_time2- start_time2 + 0.0) / CLOCKS_PER_SEC ;
            rate[int(p-min_rate)/10][K-min_k]=t1[int(p-min_rate)/10][K-min_k]/t2[int(p-min_rate)/10][K-min_k];

            double  overlapping=0;
            double *point_overlap=&overlapping;
            vector<int> a(K);
            for (int i=0;i<K;i++)
                a[i]=i+1;
            recurPermutation(ans.first,ans2.first,point_overlap,a,K,0);
            overlap[int(p-min_rate)/10][K-min_k]=*point_overlap;
            mod1[int(p-min_rate)/10][K-min_k]=ans.second;
            mod2[int(p-min_rate)/10][K-min_k]=ans2.second;
        }

    cout<<"\n\n\n\tmod1";
        for (int K=min_k;K<=max_k;K++)
            cout<<"\t"<<K;
        cout<<endl;
        for (int p=min_rate;p<=max_rate;p+=10)
        {
            cout<<"\t"<<p<<"%";
            for (int K=min_k;K<=max_k;K++)
                cout<<"\t"<<mod1[int(p-min_rate)/10][K-min_k];
            cout<<endl;
        }

  cout<<"\n\n\n\tmod2";
        for (int K=min_k;K<=max_k;K++)
            cout<<"\t"<<K;
        cout<<endl;
        for (int p=min_rate;p<=max_rate;p+=10)
        {
            cout<<"\t"<<p<<"%";
            for (int K=min_k;K<=max_k;K++)
                cout<<"\t"<<mod2[int(p-min_rate)/10][K-min_k];
            cout<<endl;
        }
}

pair<vector<vid_t>,double>  Community_detection_sampling::solve(double p,int K)
{
    //test_community_sampling(graph,10,90,2,4);
    int n=graph->VertexCount();
    int sample_num=n*p;
    MappedGraph * sub_graph= select_sub_graph(graph,sample_num);

    cout << "===== sub graph information =====" << endl;
    cout << "#vertices: " << sub_graph -> VertexCount() << endl;
    cout << "#edges: " << sub_graph -> EdgeCount() << endl;
    cout << "=============================" << endl;
    Community_detection cd(sub_graph);

    pair<vector<vid_t>,double> ans=cd.run_Girvan_NewMan(sub_graph,K);
    cout<<endl<<endl<<"\tmodularity is "<<ans.second<<endl<<endl;
    for(unsigned int i=0;i<ans.first.size();i++)
        cout<<"\t"<<ans.first[i];
    int k=*max_element(ans.first.begin(),ans.first.end());
    cout << endl<<"after assign other vertex to community" << endl;
    //vector<vid_t> final_community=allocate_vertex_with_shortest_path(graph,ans.first,sample_num);
    vector<vid_t> final_community=allocate_vertex_with_shortest_path2(graph,ans.first,sample_num,k);
    double modularity=cd.compute_modularity(graph,final_community,k);
    return make_pair(final_community,modularity);

	//pair<vector<vid_t>,double> final_community=allocate_vertex_with_label_propagation(graph,ans.first,sample_num,k);
	//return final_community;
}
