## SAE: Social Analytic Engine - MarkI
### Usage
/bin/sae –i STRING:INPUT [-o STRING:OUTPUT] [-t STRING:TASK] [TASK REPEVANT PARAMETERS]

#### Parameter Description
-i/--input: declare the direction of input graph

-o/--output: declare the direction of output files

-t/--task: declare the task, with the options as

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

[-s INT:start]: the querying node;

[-k INT:K]: the number of Top K;

psm: propensity score matching;

ec: expert classfication.

#### Example
The input file is placed at ./data/facebook

To run influence maximization with constant edge weight as 0.5 and number of seed users as 10:

./bin/sae –i ./data/facebook –t im –k 10 –w const –c 0.5

#### Example for dynamicMST
The input file is placed at ./data.txt

data.txt contains vertex number n, edge number m and all the edges of the graph:
n m
x1 y1 w1
x2 y2 w2
...

./bin/sae –i ./data.txt –t dm

##### Example for Community Detection
The input file is placed at ./data/facebook

To run community detection with community number as 5 and aglorithm 2:

./bin/sae -i ./data/facebook -t cd -k 5 -r 2

To run community detection sampling method with community number as 5 and sample probability as 0.01:

./bin/sae -i ./data/facebook -t cs -k 5 -p 0.01

##### Example for SimRank
The input file is placed at ./data/facebook

To run simrank precisely with querying node 1's Top 50 similar nodes

./bin/sae -i ./data/facebook -t sr -r 0 -v 1 -k 50

To run simrank approximately with querying node 2's Top 20 similar nodes

./bin/sae -i ./data/facebook -t sr -r 1 -v 2 -k 20

##### Example for Propensity Score Matching
./bin/sae -i ./data/expert -t psm

3
##### Example for Expert classfication
./bin/sae -i ./data/expert -t ec

./resource/dm-experts.txt
#### Basic Social Analysis Tools
./bin/sae -i ./data/fake -t social

then choose tasks listed on the terminate.
