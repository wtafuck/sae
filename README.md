## SAE: Social Analytic Engine - MarkI
### Usage
/bin/sae –i STRING:INPUT [-o STRING:OUTPUT] [-t STRING:TASK] [TASK REPEVANT PARAMETERS]

#### Parameter Description
-i/--input: declare the direction of input graph

-o/--output: declare the direction of output files

-t/--task: declare the task, with the options as

md: make and save sae graph from file

mt: make and save sae graph from file for tencent weibo data

im: run influence maximization;

[-k INT:K]: the size of seed set;

[-w STRING:WEIGHT]: set the weights of edges, with options as

rand: assign weights randomly within (0, 1);

deg: for a directed edge (u, v), let its weight to be 1/[out-degree-of-u]);

const: assign a constant weight to all edges;

A parameter is followed to declare the constant value:

-c FLOAT:CONSTANT-WEIGHT

pr: run PageRank;

sp: shortest path;

dd: demonstrate the degree distribution;

tr: count the number of triangles;

cd: community detection;

[-k INT:K]: the number of community;

[-r INT:run]: choose a algorithm for community detection;

cs: community detection with sampling:

[-k INT:K]: the number of community;

[-p FLOAT:probability]: sample probability for community detection;

sr: run SimRank ;

[-r INT:run]: choose a algorithm for SimRank;

[-s INT:start]: the querying node,which should be exits in the dataset;

[-k INT:K]: the number of Top K;

psm: propensity score matching;

ec: expert classfication.

#### Example for make data
The input file is placed at ./resource/facebook,and output file will be placed at ./data/facebook

To transform data to SAE graph

./bin/sae –i ./resource/facebook -o ./data/facebook -t md

#### Example for make tencent data
The input file is placed at /tmp/tencent8.graph,and output file will be placed at ./data/tencent8

To transform data to SAE graph

./bin/sae –i /tmp/tencent8.graph -o ./data/tencent8 -t mt

#### Example for influence maximization
The input file is placed at ./data/facebook

To run influence maximization with constant edge weight as 0.5 and number of seed users as 10:

./bin/sae -i ./data/facebook -t im -k 10 -w const -c 0.5

#### Example for dynamicMST
The input file is placed at ./data.txt

case 1:
data.txt contains vertex number n, edge number m and all the edges of the graph:
n m
x1 y1 w1
x2 y2 w2
...

./bin/sae -i ./data.txt -t dm

case 2:
data.txt contains only all the edges of the graph:
x1 y1 w1
x2 y2 w2
...

./bin/sae -i ./data.txt -t dmraw

case 3:
data.txt contains the edges of the graph, but without weight. Program will give a random weight each edge:
x1 y1
x2 y2
...

./bin/sae -i ./data.txt -t dmrawnw




##### Example for Community Detection
The input file is placed at ./data/facebook


To run community detection with community number as 5 and aglorithm 2:

./bin/sae -i ./data/facebook -t cd -k 5 -r 2

-r 1 means using Girvan-Newman aglorithm

-r 2 means using label propagation aglorithm

-r 3 means using louvain method

-r 4 means using k community core 

To run community detection sampling method with community number as 5 and sample probability as 0.01:

./bin/sae -i ./data/facebook -t cs -k 5 -p 0.01

#### Example for SimRank
The input file is placed at ./data/facebook

To run simrank precisely with querying node 1's Top 50 similar nodes,which use Partial Sums Memoization algorithm

./bin/sae -i ./data/facebook -t sr -r 0 -s 1 -k 50

To run simrank approximately with querying node 2's Top 20 similar nodes,which use Random Walk algorithm

./bin/sae -i ./data/facebook -t sr -r 1 -s 2 -k 20

##### Example for Propensity Score Matching
./bin/sae -i ./data/expert -t psm

3
##### Example for Expert classfication
./bin/sae -i ./data/expert -t ec

./resource/dm-experts.txt
#### Basic Social Analysis Tools
./bin/sae -i ./data/fake -t social

#### Example for Sampling Algorithm
the input file of sampling algorithm is placed at ./data.txt

data.txt contains vertex number n, edge number m  and all edges of graph:

x1 y1 w1

x2 y2 w2

the sampling probability of sampling algorithm is 0.001 , this probability can be adjusted according to actual situation.  

To get the average degree of graph

./bin/sae -i ./data.txt -t sad -p 0.001

To get the length2 numbers of graph

./bin/sae -i ./data.txt -t sle -p 0.001

To get the triangle numbers of graph

./bin/sae -i ./data.txt -t str -p 0.001

then choose tasks listed on the terminate.
