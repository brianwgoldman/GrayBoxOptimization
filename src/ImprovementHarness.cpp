// Brian Goldman

// When attached to a solution, the Improvement Harness
// keeps track of the fitness effects of flipping bits
// and allows for partial solution evaluation.
#include "ImprovementHarness.h"

using namespace std;

ImprovementHarness::ImprovementHarness(shared_ptr<GrayBox> evaluator_,
                                       size_t radius, Record& _recording)
    : improvements(radius),
      recording(_recording) {
  recording.start_clock();
  evaluator = evaluator_;
  // Construct neighborhood
  build_graph(evaluator, graph);
  moves_ = k_order_subgraphs(graph, radius);

  const auto& subfunctions = evaluator->epistasis();

  // Efficiency tool for converting a bit to the subfunctions it participates in
  vector<vector<size_t>> bit_to_sub(evaluator->length());
  for (size_t sub = 0; sub < subfunctions.size(); sub++) {
    for (const auto& bit : subfunctions[sub]) {
      bit_to_sub[bit].push_back(sub);
    }
  }

  sub_to_move.resize(subfunctions.size());
  move_to_sub.resize(moves_.size());
  single_bit_moves.resize(evaluator->length(), -1);

  for (size_t m = 0; m < moves_.size(); m++) {
    improvements[moves_[m].size() - 1].insert(m);
    if (moves_[m].size() == 1) {
      single_bit_moves[moves_[m][0]] = m;
    }
    for (const auto& bit : moves_[m]) {
      for (const auto& sub : bit_to_sub[bit]) {
        move_to_sub[m].insert(sub);
        sub_to_move[sub].insert(m);
      }
    }
  }

  // Table for tracking the fitness effects of making each move
  delta.resize(moves_.size(), 0);

  fitness = 0;
  check_point_fitness = 0;
  solution = nullptr;
}

// Connect the harness to the solution and calculate initial move quality information
int ImprovementHarness::attach(vector<bool>* solution_) {
  solution = solution_;
  fitness = 0;
  // fill delta with 0s
  delta.assign(delta.size(), 0);
  const size_t number_of_subs = evaluator->epistasis().size();
  for (size_t sub = 0; sub < number_of_subs; sub++) {
    auto score = evaluator->evaluate(sub, *solution);
    fitness += score;
    // Update the effect each move has on this subfunction
    for (const auto& move : sub_to_move[sub]) {
      delta[move] -= score;
      flip_move(move);
      delta[move] += evaluator->evaluate(sub, *solution);
      flip_move(move);
    }
  }
  // Initially all moves are on
  for (auto& options : improvements) {
    options.all_on();
  }
  recording.record(fitness);
  return fitness;
}

// Performs r-bit hamming ball hill climbing on the attached solution
int ImprovementHarness::optimize(Random & rand) {
  bool improvement_found = true;
  while (improvement_found) {
    improvement_found = false;
    // scan improvements until a step size has an improving move
    for (auto & potential : improvements) {
      while (potential.size() and not improvement_found) {
        auto move = potential.random(rand);
        if (0 < delta[move]) {
          make_move(move);  // Put in move
          improvement_found = true;
        }
        potential.turn_off(move);
      }
      if (improvement_found) {
        break;
      }
    }
  }
  // no recording is necessary as that is done by make_move
  return fitness;
}

// Allows for fast return to this solution state
void ImprovementHarness::set_check_point() {
  saved_delta.clear();
  flipped.clear();
  check_point_fitness = fitness;
  for (auto& options : improvements) {
    options.all_off();
  }
}

// Quickly returns to a previously saved state
int ImprovementHarness::revert() {
  for (const auto& changed : saved_delta) {
    delta[changed.first] = changed.second;
  }
  for (const auto& bit : flipped) {
    (*solution)[bit] = not (*solution)[bit];
  }
  fitness = check_point_fitness;
  set_check_point();
  return fitness;
}

// Given the index of a move, apply it to the solution
// and update auxiliary information.
int ImprovementHarness::make_move(size_t move) {
  // update fitness and record it
  fitness += delta[move];
  recording.record(fitness);
  // For each subfunction effected by this move
  for (const auto& sub : move_to_sub[move]) {
    auto pre_move = evaluator->evaluate(sub, *solution);
    flip_move(move);  // Put in move
    auto just_move = evaluator->evaluate(sub, *solution);
    flip_move(move);  // Take out move
    // for each move that overlaps the effected subfunction
    for (const auto& next : sub_to_move[sub]) {
      flip_move(next);  // Put in next
      auto just_next = evaluator->evaluate(sub, *solution);
      flip_move(move);  // Put in move
      auto move_next = evaluator->evaluate(sub, *solution);
      flip_move(move);  // Take out move
      flip_move(next);  // Take out next
      // Check point the delta if it hasn't been added already
      saved_delta.insert({next, delta[next]});
      // Take out old information and add in new information
      delta[next] += (pre_move - just_next + move_next - just_move);
      // turn the move on
      improvements[moves_[next].size() - 1].turn_on(next);
    }
  }
  // keep track of which bits have been flipped since the last checkpoint
  for (const auto& bit : moves_[move]) {
    if (flipped.count(bit)) {
      // two flips cancel out
      flipped.erase(bit);
    } else {
      flipped.insert(bit);
    }
  }
  flip_move(move);  // Put in move
  return fitness;
}

// Converts a bit index to a move index
int ImprovementHarness::modify_bit(size_t bit) {
  return make_move(single_bit_moves[bit]);
}

// Internal function to flip all bits associated with a move
void ImprovementHarness::flip_move(size_t move_index) {
  for (const auto& bit : moves_[move_index]) {
    (*solution)[bit] = not (*solution)[bit];
  }
}
