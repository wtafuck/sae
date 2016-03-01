#pragma once
#include <iomanip>
#include "social_errors.h"
#include "social_solver.h"
int social_main(sae::io::MappedGraph *graph)
{
    using namespace std;
    vector<string> list;
    Social_Solver solver(graph);
    list.push_back("Density");
    list.push_back("Components");
    list.push_back("Degree Centrality");
    list.push_back("Betweenness(Unweighted)");
    list.push_back("StructureHole Effective Size");
    list.push_back("K-Core");
    cout << "List of algorithms in social network analysis:"<< endl;
    for(int i = 0 ;i < list.size(); i++)
    {
        cout << i<<". "<<list[i]<<endl;
    }
    while(1){
        int task = -1;
        cout << "Please input the number of task: ";
        try{
            cin >> task;
            if(!(task >= 0 && task < list.size()))
                throw Exception("WrongRange");
        }
        catch(...)
        {
            cout<<"Invalid! Please input again!"<<endl;
            continue;
        }
        // switch tasks
        if(task == 0){
            cout << solver.Density() << endl;
        }
        if(task == 1){
            vector< vector<int> > components = solver.Components();
            cout << "There are "<< components.size() <<" components"<<endl;
            cout <<"==========================================="<<endl;
            cout <<setw(8)<<"Component"<<setw(8)<<"Nodes"<<setw(12)<<"Proportion"<<endl;
            cout <<"==========================================="<<endl;           
            for(int i = 0; i < components.size();i++)
                cout <<setw(8)<<i+1<<setw(8)<<components[i].size()<<setw(12)
                     <<(double)components[i].size() / graph->VertexCount()<<endl;
        }
        if(task == 2){
            vector<int> ids = solver.Degree_Centrality();
            cout << setw(5) << "No."
                 << setw(10) << "OutDegree" << setw(10) << "InDegree"
                 << setw(12) << "NrmOutDeg" << setw(12) << "NrmInDeg" <<endl;
            cout <<"===================================================="<<endl;          
            auto viter = graph->Vertices(); 
            for(int i = 0;i < ids.size();i++)
            {
                viter->MoveTo(ids[i]);
                cout<< setw(5) << ids[i]
                    << setw(10) << viter->OutEdgeCount()
                    << setw(10) << viter->InEdgeCount()
                    << setw(12) << (double)viter->OutEdgeCount() / (ids.size()-1)
                    << setw(12) << (double)viter->InEdgeCount() / (ids.size()-1)
                    <<endl;
            }
        }
        if(task == 3){
            vector<pair<double, int> > betweenness = solver.Unweighted_Betweenness();
            reverse(betweenness.begin(), betweenness.end());
            cout << setw(8) << "No." << setw(15) << "Betweenness" << setw(15) << "nBetweenness" <<endl;
            cout <<"===================================================="<<endl;
            for(int i = 0; i< betweenness.size(); i++)
                cout << setw(8) << betweenness[i].second << setw(15) <<  betweenness[i].first
                     << setw(15) << betweenness[i].first/ (graph->VertexCount() - 1) / (graph->VertexCount() - 2)
                     <<endl;
        }
        if(task == 4){
            vector< double > effective_size = solver.Effective_Size();
            cout<< setw(5) << "No." 
                << setw(20) << "Effective Size" <<endl;
            cout <<"===================================================="<<endl;          
            for(int i = 0; i < graph->VertexCount(); i++)
                cout<< setw(5) << i << setw(20) << effective_size[i] <<endl;            
        }
        if(task == 5){
            vector<int> k_core = solver.K_Core();
            cout<< setw(5) << "No." << setw(10) << "K-Core" <<endl;
            cout <<"===================================================="<<endl;          
            for(int i = 0;i < k_core.size();i++)
                cout<< setw(5) << i << setw(10) << k_core[i] <<endl;           
        }
        break;
    }
}