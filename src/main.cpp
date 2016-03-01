#include "argparse/macro-argparse-jquery.hh"
#include "solver/solver.h"
#include "sampling/triangle_sampling.h"
#include "basic/triangle_count.h"
#include "basic/pagerank.h"
#include "basic/degree_distribution.h"
#include "storage/graph_builder.h"
#include "storage/mgraph.h"
#include "report/table_generator.h"
#include "influence/influence_maximization.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <map>
#include <cstdlib>

using namespace std;
using namespace sae::io;

DEF_ARGUMENT_CLASS(
	Argument,
	std::string,	input,		"",				REQUIRED,   OPT_SLH(-i, --input, "input data"),
    std::string,    output,     "./output",     OPTIONAL,   OPT_SLH(-o, --output, "output direction"),
    std::string,    task,       "",             OPTIONAL,   OPT_SLH(-t, --task, "declear task:\nim(inf-max),dd(degree),pr(PageRank),tr(triangle)\n"),
    int,            para_im_k,  "0",            OPTIONAL,   OPT_SLH(-k, --seed, "seed size"),
    std::string,    para_edge_w,    "rand",     OPTIONAL,   OPT_SLH(-w, --weight, "edge weights:rand,const,deg\n"),
    double,         para_const, 0.0,            OPTIONAL,   OPT_SLH(-c, -constant, "constant value")
);

string output_dir;

void makeFakeData(int numVertex=10, double p = 0.1) {
	//int numEdge = rand() % (numVertex * numVertex / 3) + numVertex;
    int numEdge = (numVertex * (numVertex - 1)) / 2 * p;
	GraphBuilder<int> graph;
	for (int i = 0; i < numVertex; ++i) graph.AddVertex(i, 0);
	map<pair<int, int>, bool> edges;
	for (int i = 0; i < numEdge; ++i) {
		pair<int, int> edge;
		do {
			int x = rand() % (numVertex - 1);
			int y = rand() % (numVertex - x - 1) + x + 1;
			edge = make_pair(x, y);
		} while (edges.find(edge) != edges.end());
		edges[edge] = true;
		graph.AddEdge(edge.first, edge.second, 0);
		//cout << edge.first << " " << edge.second << endl;
	}
	system("mkdir -p fake");
	graph.Save("./fake/graph");
}

void runPageRank(MappedGraph *graph) {
    PageRank pr(graph);
    time_t start_time = clock();
    vector<pair<sae::io::vid_t, double>> res = pr.solve();
    FILE *fout = fopen((output_dir + "/pagerank").c_str(), "w");
    fprintf(fout, "vid\tpage_rank_score\n");
    for (unsigned int i = 0; i < res.size(); i ++)
        fprintf(fout, "%llu\t%.5lf\n", res[i].first, res[i].second);
    fclose(fout);
    time_t end_time = clock();
    cout << "Running time of PageRank: " << (end_time - start_time + 0.0) / CLOCKS_PER_SEC << endl;
}

void runDegreeDistribution(MappedGraph *graph) {
    Degree_Distribution dd(graph);
    time_t start_time = clock();
    vector<pair<int, double>> res = dd.solve();
    FILE *fout = fopen((output_dir + "/degree_distribution").c_str(), "w");
    fprintf(fout, "degree\tproportion\n");
    for (unsigned int i = 0; i < res.size(); i ++) {
        fprintf(fout, "%d\t%.5lf\n", res[i].first, res[i].second);
    }
    fclose(fout);
    fout = fopen((output_dir + "/degree_of_each_vertex").c_str(), "w");
    for (auto iter = graph -> Vertices(); iter -> Alive(); iter -> Next()) {
        fprintf(fout, "%llu\t%d\n", iter -> GlobalId(), iter -> OutEdgeCount());
    }
    fclose(fout);
    time_t end_time = clock();
    cout << "Running time of Degree_Distribution: " << (end_time - start_time + 0.0) / CLOCKS_PER_SEC << endl;
}

void runInfluenceMaximization(MappedGraph *graph, int K, double W) {
    time_t start_time = clock();
    Influence_Maximization im(graph);
    pair<vector<sae::io::vid_t>, double> res = im.solve(K, 20, W);
    FILE* fout = fopen((output_dir + "/influence_maximization").c_str(), "w");
    fprintf(fout, "Expected #actived users: %.5lf\n", res.second);
    cout << "Expected #actived users: " << res.second << endl;
    for (unsigned int i = 0; i < res.first.size(); i ++) {
        fprintf(fout, "%llu\n", res.first[i]);
    }
    fclose(fout);
    time_t end_time = clock();
    cout << "Running time of Influence_Maximization: " << (end_time - start_time + 0.0) / CLOCKS_PER_SEC << endl;
}

int main(int argc, char **argv) {
    int vertexNum = 100;
    double edgeProb = 0.1;
    srand(time(NULL));
    //makeFakeData(vertexNum, edgeProb);
    
	// parse arguments
	Argument args;
	if (! args.parse_args(argc, argv)) 
        return 1;

    // declare input file
    if (args.input().length() == 0)
        return 1;
    MappedGraph *graph = MappedGraph::Open(args.input().c_str());
    cout << "===== graph information =====" << endl;
    cout << "#vertices: " << graph -> VertexCount() << endl;
    cout << "#edges: " << graph -> EdgeCount() << endl;
    cout << "=============================" << endl;

    // declare output file direction
    if (args.output().length() > 0) {
        output_dir = args.output().c_str();
    }
    system(("mkdir -p " + output_dir).c_str());
        
    //declare task
    string task = args.task();

    // call influence maximization
    if (task == "im") {
        int k = args.para_im_k();
        if (k <= 0) {
            cout << "The number of seeds must greater than zero!" << endl;
            return 1;
        }
        string weight_set = args.para_edge_w();
        if (weight_set == "const") {
            double w = args.para_const();
            if (w == 0) {
                cout << "Please define a constant value greater than zero!" << endl;
                return 1;
            }
            cout << "Running influence maximization with #seeds as " << k << endl;
            runInfluenceMaximization(graph, k, w);
        }
        if (weight_set == "rand") {
            cout << "Running influence maximization with #seeds as " << k << endl;
            runInfluenceMaximization(graph, k, 0.0);
        }
        if (weight_set == "deg") {
            cout << "Running influence maximization with #seeds as " << k << endl;
            runInfluenceMaximization(graph, k, -1.0);
        }
    }

    // call PageRank
    if (task == "pr") {
        runPageRank(graph);
    }

    // call degree distribution
    if (task == "dd") {
        runDegreeDistribution(graph);
    }
	
    // generate a graph
    if (task == "gg") {
        makeFakeData(vertexNum, edgeProb);
    }
	//testTable();
	//makeFakeData(vertexNum, edgeProb);
	//testTriangle();
   
    return 0;
}

