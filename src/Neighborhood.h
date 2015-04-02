// Brian Goldman

// Collection of utility functions for extracting graph based information.
// Used by the Improvement Harness to understand problem epistasis.
#ifndef NEIGHBORHOOD_H_
#define NEIGHBORHOOD_H_

#include "Evaluation.h"
#include <unordered_map>
using std::unordered_map;
#include <unordered_set>
using std::unordered_set;
#include "Util.h"

// Constructs a sparse graph from the epistasis tables of the evaluator
void build_graph(shared_ptr<GrayBox> evaluator,
                 vector<unordered_set<size_t>> & graph);

// Finds all connected induced subgraphs with "radius" or less vertices.
vector<vector<size_t>> k_order_subgraphs(
    const vector<unordered_set<size_t>>& graph, size_t radius);

// Recursive function used internally by k_order_subgraphs.
void recurse(const vector<unordered_set<size_t>>& graph,
             size_t v, unordered_set<size_t> & closed, vector<size_t> & prev,
             unordered_set<size_t> & prevopen, size_t radius,
             vector<vector<size_t>> & found);

// Find a random induced subgraph which contains "start" and has "k" vertices.
// Used by Gray Box P3's crossover operator
unordered_set<size_t> random_induced_subgraph(
    const vector<unordered_set<size_t>> & graph, size_t start,
    size_t k, Random& rand);
#endif /* NEIGHBORHOOD_H_ */
