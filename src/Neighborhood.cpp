// Brian Goldman

// Utility functions related to graphs,
// used by Improvement Harness to find connected induced
// sug graphs.
#include "Neighborhood.h"

// Convert the epistasis of a Gray Box problem into an adjacency list
// representation of a graph
void build_graph(shared_ptr<GrayBox> evaluator,
                 vector<unordered_set<size_t>> & graph) {
  graph.clear();
  graph.resize(evaluator->length());
  const auto& subfunctions = evaluator->epistasis();
  for (auto subfunction : subfunctions) {
    for (auto x : subfunction) {
      for (auto y : subfunction) {
        if (x != y) {
          graph[x].insert(y);
        }
      }
    }
  }
}

// Finds all possible connected induced subgraphs with k or less vertices
vector<vector<size_t>> k_order_subgraphs(
    const vector<unordered_set<size_t>>& graph, size_t radius) {
  // The list of all subgraphs found so far.
  vector<vector<size_t>> found;
  unordered_set<size_t> closed, prevopen;
  // Vertices already in the set, initially empty
  vector<size_t> prev;
  // Call the recursive function starting at each vertex
  for (size_t v=0; v < graph.size(); v++) {
    closed.insert(v);
    recurse(graph, v, closed, prev, prevopen, radius, found);
  }
  return found;
}

// Recursive function used by k_order_subgraphs to expand each subset by "v"
// and then find all possible subgraphs which contain "prev" and "v"
void recurse(const vector<unordered_set<size_t>>& graph,
             size_t v, unordered_set<size_t> & closed, vector<size_t> & prev,
             unordered_set<size_t> & prevopen, size_t radius,
             vector<vector<size_t>> & found) {
  // Create a new subset containing "prev" and "v"
  vector<size_t> inset(prev);
  inset.push_back(v);
  found.push_back(inset);
  // Add that subset to "found"
  if (inset.size() >= radius) {
    return;
  }
  vector<size_t> closed_here;
  // open everything adjacent to v
  unordered_set<size_t> openset(prevopen);
  const auto & adjacent = graph.at(v);
  openset.insert(adjacent.begin(), adjacent.end());
  // recurse using each not-closed option
  for (const auto& working : openset) {
    if (closed.count(working) == 0) {
      closed_here.push_back(working);
      closed.insert(working);
      recurse(graph, working, closed, inset, openset, radius, found);
    }
  }
  // reopen anything closed at this level
  for (const auto& working : closed_here) {
    closed.erase(working);
  }
}

// Finds a random induced subgraph given an adjacency graph and a starting point,
// such that "start" is in the graph. Vertices are added at random biased by
// how many different edges connect them to the current subgraph.
unordered_set<size_t> random_induced_subgraph(
    const vector<unordered_set<size_t>> & graph, size_t start,
    size_t k, Random& rand) {
  // The subgraph
  unordered_set<size_t> subset;
  // options initially contains 1 copy of "start"
  vector<size_t> options(1, start);
  // stop if you run out of options or the subset size is reached.
  while (options.size() > 0 and subset.size() < k) {
    // select a random option and remove it from the list
    size_t index = std::uniform_int_distribution<size_t>(0, options.size() - 1)(rand);
    size_t working = options[index];
    std::swap(options[index], options.back());
    options.pop_back();

    // attempt to add it to the result
    const auto& result = subset.insert(working);
    // If it was inserted
    if (result.second) {
      // add its adjacency to "open"
      const auto& adjacent = graph.at(working);
      options.insert(options.end(), adjacent.begin(), adjacent.end());
    }
  }
  return subset;
}
