// Brian Goldman

// For Black Box P3, the population is in charge of storing solutions,
// performing linkage learning, and performing crossover.
// Included only for timing purposes.

#ifndef BLACKBOXPOPULATION_H_
#define BLACKBOXPOPULATION_H_

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <cmath>
#include <algorithm>
#include <numeric>

#include "Util.h"
#include "ImprovementHarness.h"

using std::vector;
using std::unordered_map;
using std::array;
using std::log;

// Stores all information related to a population
class Population {
 public:
  Population(Configuration& config);
  virtual ~Population() = default;

  // Adds a solution to the population.
  void add(const vector<bool> & solution);
  // Improves the passed in solution using crossover, as guided by the
  // entropy table and solutions in the population.  "fitness" should
  // start as the fitness of "solution", and both will be set to the
  // improved values after application.
  void improve(Random& rand, vector<bool> & solution, int & fitness,
               ImprovementHarness& evaluator);
  // Access the entropy table to construct a linkage tree
  void rebuild_tree(Random& rand);

 private:
  // Storage for all of the solutions in the population
  vector<vector<bool> > solutions;

  // Tools used for cluster creation
  vector<vector<int> > clusters;
  vector<int> cluster_ordering;

  // Keeps track of how often each possible pairing of bit values occurs in the population
  // for all pairs of genes
  vector<vector<array<int, 4>>> occurrences;
  // Scratch space for the entropy distance between clusters
  static vector<vector<float> > distances;

  // Given a list of occurrences, return the pairwise_distance
  float get_distance(const array<int, 4>& entry) const;

  // tool to access pairwise distance for two genes
  float get_distance(int x, int y) const;

  // Donates the genes specified by "cluster" from "source" into "solution".
  // If the solution was modified by this action, it is evaluated.  If the new
  // fitness is no worse the modification is kept and the fitness value is updated.
  // Returns true if an evaluation was performed.
  bool donate(vector<bool> & solution, int & fitness, vector<bool> & source,
      const vector<int> & cluster, ImprovementHarness& evaluator) const;

  // Sort clusters by size, smallest first.  Randomize order of equal sized clusters
  static void smallest_first(Random& rand, const vector<vector<int>>& clusters, vector<int>& cluster_ordering);

  // General configuration options
  // Number of genes in the genome
  size_t length;
};

#endif /* POPULATION_H_ */
