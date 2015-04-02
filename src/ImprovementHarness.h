// Brian Goldman

// An Improvement Harness is used to efficiently perform
// r-bit hamming ball hill climbing and partial solution evaluation.
// The harness is first attached to a solution, and then used to modify
// that solution.
// Based on the paper "Efficient Identification of Improving Moves in a Ball
// for Pseudo-Boolean Problems" by Francisco Chicano, Darrell Whitley,
// and Andrew M. Sutton.

#ifndef IMPROVEMENTHARNESS_H_
#define IMPROVEMENTHARNESS_H_
#include "Evaluation.h"
#include "Neighborhood.h"
#include "RandSet.hpp"
#include "Record.h"

#include <unordered_map>
using std::unordered_map;
#include <unordered_set>
using std::unordered_set;

class ImprovementHarness {
  // List of all moves, which are just collections of indices
  vector<vector<size_t>> moves_;
  // Lookup tables to find which moves effect what subfunctions, and vice versa
  vector<unordered_set<size_t>> move_to_sub, sub_to_move;
  vector<size_t> single_bit_moves;
  // Table storing the fitness effect of making a particular move
  vector<int> delta;
  // The actual problem being solved, used to perform partial evaluations.
  shared_ptr<GrayBox> evaluator;
  // The attached solution
  vector<bool> * solution;

  // Keeps track of which moves of each size are potential fitness improvements
  vector<RandSet<size_t, Random>> improvements;

  // Internal function to flip the bits related to a move
  void flip_move(size_t move_index);
  // Internal function which performs all updates caused by making a move
  int make_move(size_t move_index);

  // Track fitness improvement information
  Record& recording;

  // The current fitness of the attached solution
  int fitness, check_point_fitness;
  // Used for returning efficiently to previous states
  unordered_map<size_t, int> saved_delta;
  unordered_set<size_t> flipped;

  // The adjacency information in the epistasis graph
  vector<unordered_set<size_t>> graph;

 public:
  ImprovementHarness(shared_ptr<GrayBox> evaluator_, size_t radius,
                     Record& _recording);
  virtual ~ImprovementHarness() = default;
  // Connects a solution to the harness, Calculates initial delta information
  // in preparation for partial solution evaluation and hill climbing. Returns
  // current fitness
  int attach(vector<bool>* solution_);
  // Performs r-bit hamming ball hill climbing on the attached solution.
  // Returns the fitness of the result.
  int optimize(Random & rand);
  // Wrapper allowing for black box solution evaluation
  int evaluate(const vector<bool>& solution_) {
    auto fit = evaluator->evaluate(solution_);
    recording.record(fit);
    return fit;
  }

  // When "revert" is called, the solution and the harness's state are
  // reverted to how they looked when the last "set_check_point" was called.
  void set_check_point();

  // Flip a bit by index in the solution and update the state.
  // Returns the fitness result.
  int modify_bit(size_t bit);
  // Returns how many bits have been flipped since the last checkpoint call
  size_t modified() {
    return flipped.size();
  }
  int revert();

  // Provides access to the problem's epistasis
  const vector<vector<size_t>>& epistasis() const {
    return evaluator->epistasis();
  }
  // Forward max_fitness so algorithms can know if they have reached the global optimum.
  int max_fitness() const {
    return evaluator->max_fitness();
  }
  // Provides access to the problem's epistasis graph,
  // used by Gray Box P3's crossover operator.
  const vector<unordered_set<size_t>>& adjacency() const {
    return graph;
  }
};

#endif /* IMPROVEMENTHARNESS_H_ */
