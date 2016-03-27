#include "community_detection_sampling.h"
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

vector<pair<double, int> > compute_betweenness_sampling(MappedGraph *graph,vector< bool> is_deleted,unsigned int n)
{
    vector<pair<double, int> > edge_betweenness;
     for(auto eiter = graph->Edges(); eiter->Alive(); eiter->Next())
            edge_betweenness.push_back(make_pair(0,eiter->GlobalId()));

    for(unsigned int s = 0; s < n; s++)
    {
        vector< int > stack;
        queue < int >  search_queue;
        vector< int > path_count(n,0);
        vector< int > dis(n,-1);
        vector< vector<eid_t> > pre_vertex(n);
        vector<double> part_value(n,0);
        dis[s] = 0;
        path_count[s] = 1;
        search_queue.push(s);
        auto viter = graph->Vertices();
        while(!search_queue.empty())
        {
            int v=search_queue.front();
            search_queue.pop();
            stack.push_back(v);
            viter->MoveTo(v);
            for(auto eiter = viter->OutEdges(); eiter->Alive(); eiter->Next())
            {
                if(is_deleted[eiter->GlobalId()]) continue;
                int w = eiter->TargetId();
                if(dis[w] < 0)
                {
                    search_queue.push(w);
                    dis[w] = dis[v] + 1;
                }
                if(dis[w] == dis[v] + 1)
                {
                    path_count[w] += path_count[v];
                    pre_vertex[w].push_back(v);
                }
            }
        }

        for(int i = stack.size()-1; i >= 0; i--)
        {
            unsigned int w = stack[i];
            if(w==s) continue;
            for(unsigned int j = 0;j < pre_vertex[w].size(); j++)
            {
                part_value[pre_vertex[w][j]] += (double)path_count[pre_vertex[w][j]] / path_count[w] * (1 + part_value[w]);
                if(pre_vertex[w][j]==s||w==s) continue;
                auto vex = graph->Vertices();
                vex->MoveTo(pre_vertex[w][j]);
                for(auto eiter = vex->OutEdges(); eiter->Alive(); eiter->Next())
                {
                    if(is_deleted[eiter->GlobalId()]) continue;
                    if(eiter->TargetId()==w)
                    {
                        edge_betweenness[eiter->GlobalId()].first +=(double)path_count[pre_vertex[w][j]] / path_count[w] * part_value[w];
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

pair<vector<vid_t>,int >  split_community_sampling(MappedGraph *graph,vector< bool> is_deleted,bool is_sampling,int sample_num)
{
    unsigned int n=is_sampling==false?graph->VertexCount():sample_num;
    vector <vid_t> community(graph->VertexCount(),0);
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

double  compute_modularity_sampling(MappedGraph *graph, vector<eid_t> community,int k,vector<bool> is_deleted,bool is_sampling,int sample_num)
{
    int m = graph->EdgeCount();
    vector< vector<eid_t> > comset(k);
    double ans=0,eii=0,ai=0;
    auto viter = graph->Vertices();
    int  number=is_sampling==true?sample_num:community.size();
    if (is_sampling)
    {
        m=0;
        auto eiter = graph->Edges();
        for(; eiter->Alive(); eiter->Next())
        {
            if(is_deleted[eiter->GlobalId()]) continue;
                m++;
        }
    }

    for (unsigned int i=0;i<number;i++)
            comset[community[i]-1].push_back(i);
    for (int i=0;i<k;i++)
    {
        double eii_temp=0,ai_temp=0;
        for(unsigned int j=0;j<comset[i].size();j++)
        {
                viter->MoveTo(comset[i][j]);
                for(auto eiter = viter->OutEdges(); eiter->Alive(); eiter->Next())
                {
                        if(is_deleted[eiter->GlobalId()]&&is_sampling) continue;
                        if (find(comset[i].begin(), comset[i].end(), eiter->TargetId()) != comset[i].end())
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


pair<vector< bool>, vector< bool>>select_sub_graph(MappedGraph *graph,int num)
{
    vector< bool> is_part(graph->VertexCount(),false);
    vector< bool> is_deleted(graph->EdgeCount(),true);
    auto viter = graph->Vertices();
    for (int i=0;i<num;i++)
        is_part[i]=true;

    for(int i=0;i<num;i++)
    {
        viter->MoveTo(i);
        for(auto eiter = viter->OutEdges(); eiter->Alive(); eiter->Next())
            if(is_part[eiter->TargetId()]&&is_part[i])
                is_deleted[eiter->GlobalId()]=false;
    }

    return make_pair(is_part,is_deleted);
}

vector<vid_t> allocate_vertex_with_shortest_path(MappedGraph *graph,vector<vid_t> community,int sample_num)
{
        int n=graph->VertexCount();
        vector<vid_t> communityChange=community;
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
            int min_value=10000,min_index=i;
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


pair<vector<vid_t>,double>   allocate_vertex_with_label_propagation(MappedGraph *graph,vector<vid_t> community,int sample_num,vector<bool>is_deleted,bool is_sampling)
{
        srand(time(NULL));
        unsigned int  n=graph->VertexCount();
        vector<vector<vid_t>> neighbor(n);
        vector<vid_t> C(n,1);
        vector<vid_t> temp(n);
        vector<vid_t> max_community=C;
        double max_modularity=0;
        vid_t max_community_count=*max_element(community.begin(),community.end());
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
            double modularity=compute_modularity_sampling(graph,C,max_community_count,is_deleted,is_sampling,sample_num);
            if (modularity>max_modularity) {
                max_modularity=modularity;max_community=C;
            }
        }
            return make_pair(max_community,max_modularity);
}

pair<vector<vid_t>,double>   run_Girvan_NewMan_sampling(MappedGraph *graph,vector< bool>is_deleted,int m,int n,vector<bool> initialDeleted,bool is_sampling,int sample_num,int K)
{
    vector< vector<eid_t> > community(m);
    vector< eid_t> temp(n,1);
    vector<pair<double,int>> modularity(m);
    for(int i=0;i<n;i++)
        community[i]=temp;
    int k=1,j=0;
    for (;j<m;j++)
    {
        vector<pair<double, int> >edge_betweenness =compute_betweenness_sampling(graph,is_deleted,n);

        for (unsigned int c=0;c<edge_betweenness.size();c++)
            if(edge_betweenness[c].first==edge_betweenness[0].first)
            {
                int first =edge_betweenness[c].second;
                int second=edge_betweenness[c].second%2==0?edge_betweenness[c].second+1:edge_betweenness[c].second-1;
                is_deleted[first]=true;
                is_deleted[second]=true;
            }

        if(edge_betweenness[0].first==0)  break;

        if(j>0) community[j]=community[j-1];
        modularity[j]=modularity[j-1];

        pair<vector<vid_t>,int > splitAns=split_community_sampling(graph,is_deleted,is_sampling,sample_num);
        community[j]=splitAns.first;
        k=splitAns.second;
        modularity[j].second=k;
        modularity[j].first=compute_modularity_sampling(graph,community[j],k,initialDeleted,is_sampling,sample_num);
        if (k==K) break;
    }
    int maxIndex=0;
    if (k==K) maxIndex=j;
    else{
        for(int i=0;i<=j;i++)
        if (modularity[i].first>modularity[maxIndex].first) maxIndex=i;
    }

    if(is_sampling)
    {
        cout<<"\n\tthe sub graph's best community structure"<<endl;
        for(int i=0;i<community[maxIndex].size();i++)
            cout<<"\t"<<community[maxIndex][i];
        cout<<endl<<"\tmodularity is "<<modularity[maxIndex].first<<endl<<endl;

        cout<<"\tafter allocate other vertexs to these community"<<endl;
//        vector <vid_t> assignedCommunity = allocate_vertex_with_shortest_path(graph,community[maxIndex],sample_num);
//        double finalModularity=compute_modularity_sampling(graph,assignedCommunity,k,is_deleted,false,sample_num);
//         return make_pair(assignedCommunity,finalModularity);
        pair<vector <vid_t>,double> final_community = allocate_vertex_with_label_propagation(graph,community[maxIndex],sample_num,is_deleted,false);
        return make_pair(final_community .first,final_community.second);
    }
    return make_pair(community[maxIndex],modularity[maxIndex].first);
}

void  test_community_sampling(MappedGraph *graph,int min_rate,int max_rate,int min_k,int max_k)
{
    cout<<"\tRun community detection sampling algorithm"<<endl;
    double overlap[(max_rate-min_rate)/10+1][max_k-min_k+1];
    double mod1[(max_rate-min_rate)/10+1][max_k-min_k+1];
    double mod2[(max_rate-min_rate)/10+1][max_k-min_k+1];
    for (int p=min_rate;p<=max_rate;p+=10)
        for (int K=min_k;K<=max_k;K++)
        {
            cout<<"\n\n\sampling rate:"<<p<<"%"<<"\tcommunity number:"<<K<<endl<<endl;
            int n=graph->VertexCount();
            int sample_num=n*p/100;
            vector< bool> is_deleted(graph->EdgeCount(),false);
            pair<vector<vid_t>,double> ans=run_Girvan_NewMan_sampling(graph,is_deleted,graph->EdgeCount(),n,is_deleted,false,n,K);
            cout<<"\n\tthe whole graph's best community structure"<<endl;
            for(int i=0;i<ans.first.size();i++)
                cout<<"\t"<<ans.first[i];
            cout<<endl<<"\tmodularity is "<<ans.second<<endl;
            pair<vector< bool>, vector< bool>> isDeletForAll= select_sub_graph(graph,sample_num);
            vector< bool> is_part =isDeletForAll.first;
            is_deleted=isDeletForAll.second;
            pair<vector<vid_t>,double> ans2=run_Girvan_NewMan_sampling(graph,is_deleted,graph->EdgeCount(),is_part.size(),is_deleted,true,sample_num, K);
            cout<<"\n\tthe sampling  graph's best community structure"<<endl;
            for(int i=0;i<ans2.first.size();i++)
                cout<<"\t"<<ans2.first[i];
            cout<<endl<<"\tmodularity is "<<ans2.second<<endl;

            int overlapping=0;
            for(int i=0;i<n;i++)
                if (ans.first[i]==ans2.first[i])
                    overlapping++;
            overlap[int(p-min_rate)/10][K-min_k]=overlapping*1.0/n;
            mod1[int(p-min_rate)/10][K-min_k]=ans.second;
            mod2[int(p-min_rate)/10][K-min_k]=ans2.second;
        }

        cout<<"\n\n\n\tmodularity1";
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

        cout<<"\n\n\n\tmodularity2";
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

        cout<<"\n\n\n\toverlap";
        for (int K=min_k;K<=max_k;K++)
            cout<<"\t"<<K;
        cout<<endl;
        for (int p=min_rate;p<=max_rate;p+=10)
        {
            cout<<"\t"<<p<<"%";
            for (int K=min_k;K<=max_k;K++)
                cout<<"\t"<<overlap[int(p-min_rate)/10][K-min_k];
            cout<<endl;
        }
}


int Community_detection_sampling::solve()
{

    test_community_sampling(graph,10,90,2,4);
    return 0;
}
