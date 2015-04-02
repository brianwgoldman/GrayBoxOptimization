// Brian Goldman

// This is a collection of utilities used throughout programming.
// Basically just a bunch of misfits.

#ifndef UTIL_H_
#define UTIL_H_

#include <vector>
#include <iostream>
#include <random>
#include <algorithm>
#include <sstream>

using std::vector;

// Defines that the random number generator in use is the Mersenne Twister
using Random=std::mt19937;

// Generate a random bit vector of the requested length
void rand_vector(Random& rand, vector<bool>& solution);

// Print out a vector to a stream.
void print(const vector<bool> & vect, std::ostream & out = std::cout);

// Read in a vector from a stream.
void read(vector<bool> & vect, std::istream & in = std::cin);

// Returns the entropy given the list of counts and a total number,
// where total = sum(counts)
template<size_t T>
float entropy(const std::array<int, T>& counts, const float& total) {
  float sum = 0;
  float p;
  for (const auto& value : counts) {
    if (value) {
      p = value / total;
      sum -= (p * log(p));
    }
  }
  return sum;
}


#endif /* UTIL_H_ */
