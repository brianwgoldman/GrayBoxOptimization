// Brian Goldman

// Implementation of the Gray-Box P3 optimization method,
// which uses a pyramid of stored populations and linkage informed
// crossover to improve locally optimal solutions.
#include "Pyramid.h"

using namespace std;

Pyramid::Pyramid(Random& _rand, Configuration& _config,
                 ImprovementHarness& _harness)
    : Optimizer(_rand, _config, _harness) {
}

// Construct a collection of subsets (blocks) by generating
// random induced subgraphs. There are 2*N-2 of these subgraphs
// and the sum of their sizes is O(N log N)
void Pyramid::cis_tree(vector<vector<size_t>> & blocks) const {
  uniform_int_distribution<size_t> get_start(0, length - 1);
  // Initial sizes are the full size split into two random parts
  vector<size_t> sizes(2);
  sizes[0] = uniform_int_distribution<size_t>(1, length - 1)(rand);
  sizes[1] = length - sizes[0];

  // Works through "sizes", adding more as it goes.
  for (size_t sizes_index = 0; sizes_index < sizes.size(); sizes_index++) {
    size_t working_size = sizes[sizes_index];
    if (working_size > 1) {
      // split in half, add 2 new sizes
      sizes.push_back(uniform_int_distribution<size_t>(1, working_size - 1)(rand));
      sizes.push_back(working_size - sizes.back());
    }
    // Create a random connected induced subgraph with a random starting point
    // give the working size.
    const auto& subset = random_induced_subgraph(harness.adjacency(),
                                                 get_start(rand), working_size,
                                                 rand);
    blocks.emplace_back(subset.begin(), subset.end());
  }
}

// Adds the solution to the level if it doesn't already exist somewhere
// in the pyramid
void Pyramid::add_if_unique(const vector<bool>& candidate, size_t level) {
  if (seen.count(candidate) == 0) {
    // Create the new level if necessary
    if (solutions.size() == level) {
      solutions.push_back(vector<vector<bool>>(0));
      selector_tool.push_back(vector<size_t>(0));
    }
    // Add the solution to the population
    selector_tool[level].push_back(selector_tool[level].size());
    solutions[level].push_back(candidate);
    seen.insert(candidate);
  }
}

// Performs a full iteration of the pyramid, including
// random generation, r-bit hamming ball hill climbing,
// and crossover with each pyramid level.
int Pyramid::iterate() {
  rand_vector(rand, solution);
  harness.attach(&solution);
  auto fitness = harness.optimize(rand);

  bool improved = true;
  for (size_t level = 0; level < solutions.size(); level++) {
    // only add it if it has improved since the last time it was added
    if (improved) {
      add_if_unique(solution, level);
      improved = false;
    }
    // construct the crossover subsets and shuffle them
    vector<vector<size_t>> blocks;
    cis_tree(blocks);
    shuffle(blocks.begin(), blocks.end(), rand);

    auto& options = selector_tool[level];
    size_t limiter = blocks.size();
    for (size_t index = 0; index < limiter; index++) {
      size_t limit = options.size();
      harness.set_check_point();
      // choose donors at random until a crossover causes a modification
      while (limit > 0 and harness.modified() == 0) {
        // Select a random donor without replacement
        size_t choice = uniform_int_distribution<size_t>(0, limit - 1)(rand);
        size_t donor = options[choice];
        swap(options[choice], options[limit - 1]);
        limit--;

        // copy all bits for this block into solution
        for (size_t bit : blocks[index]) {
          if (solution[bit] != solutions[level][donor][bit]) {
            harness.modify_bit(bit);
          }
        }
      }

      // Perform hill climbing after the crossover
      auto new_fitness = harness.optimize(rand);

      // Keep it if the fitness is no worse.
      if (fitness <= new_fitness) {
        if (fitness < new_fitness) {
          improved = true;
        }
        fitness = new_fitness;
        harness.set_check_point();
        // drop out now if you reached the global optimum
        if (fitness >= harness.max_fitness()) {
          return fitness;
        }
      } else {
        // undo all changes made by this donation.
        harness.revert();
      }
    }
  }
  if (improved) {
    add_if_unique(solution, solutions.size());
  }
  return fitness;
}
