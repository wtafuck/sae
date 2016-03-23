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

dd: demonstrate the degree distribution;

tr: count the number of triangles.

#### Example
The input file is placed at ./data/facebook

To run influence maximization with constant edge weight as 0.5 and number of seed users as 10:

./bin/sae –i ./data/facebook –t im –k 10 –w const –c 0.5
