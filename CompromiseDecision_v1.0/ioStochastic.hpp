//
//  ioStochastic.hpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/18/25.
//

#ifndef ioStochastic_hpp
#define ioStochastic_hpp

#include <stdio.h>

#include <random>
#include "DataStructure.hpp"

// read sto file
stoMap readStochasticMap(const std::string& stochasticPath);
// read ground truth
std::vector<stoPoint> readGroundTruth_discreteRV(const std::string& distPath);

// discrete random number generator, random_number must be in [0,1)
double discrete_random_number_generator(double random_number, const discreteDistribution dist);
// uniform random number generator
double uniform_random_number_generator(std::mt19937& generator, double lb, double ub);
// normal random number generator
double normal_random_number_generator(std::mt19937& generator, double mean, double stddev);
// truncated normal random number generator
double truncated_normal_random_number_generator(std::mt19937& generator, double mean, double stddev, double lb, double ub);
#endif /* ioStochastic_hpp */
