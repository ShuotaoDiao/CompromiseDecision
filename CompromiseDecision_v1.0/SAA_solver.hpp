//
//  SAA_solver.hpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/18/25.
//

#ifndef SAA_solver_hpp
#define SAA_solver_hpp

#include <stdio.h>
#include <stdlib.h> // rand
#include <ctime>
#include <cmath>
#include <random>
#include <list>

#include "ioModel.hpp"
#include "ioStochastic.hpp"
#include "ioDB.hpp"

// SAA problem per replication
solverOutput saa_2slp(const std::string& folder_path,
                      std::mt19937& generator,
                      int sample_size,
                      int idx_rep = 0);

solverOutput saa_sqqp(const std::string& folder_path, std::mt19937& generator, int sample_size);

// SAA problem after aggragation, compromise decision
// master problem for optimization
compromiseOutput saa_2slp_compromise_master(standardTwoStageParameters& model_parameters,
                                            const stoMap& stochastic_map,
                                            std::vector<solverOutput>& aggregate_outputs,
                                            int sample_size_per_replication,
                                            int num_replications,
                                            double regularizer_coefficient);

compromiseOutput saa_sqqp_compromise_master(standardTwoStageParameters& model_parameters,
                                            const stoMap& stochastic_map,
                                            std::vector<solverOutput>& aggregate_outputs,
                                            int sample_size_per_replication,
                                            int num_replications,
                                            double regularizer_coefficient);

// compromise decision process to coordinate replication and aggregation
compromiseOutput saa_2slp_compromise(const std::string& folder_path,
                                     std::mt19937& generator,
                                     int sample_size_per_replication,
                                     int num_replications,
                                     double regularizer_coefficient);

compromiseOutput saa_sqqp_compromise(const std::string& folder_path,
                                     std::mt19937& generator,
                                     int sample_size_per_replication,
                                     int num_replications,
                                     double regularizer_coefficient);

// ground truth (finite support)
solverOutput ground_truth_discrete(const std::string& folder_path);

solverOutput sqqp_ground_truth_discrete(const std::string& folder_path);

// validation
// two-stage SLP
// true dist is fine to enumerate
double validation_ground_truth_discrete(const std::string& folder_path,
                                              std::vector<double>& x);
// validation set is given
double validation_saa(const std::string& folder_path,
                          std::vector<double>& x);

// SQQP
double sqqp_validation_ground_truth_discrete(const std::string& folder_path,
                                        std::vector<double>& x);

double sqqp_validation_saa(const std::string& folder_path,
                           std::vector<double>& x);

// random number generator
stoPoint generator_random(standardTwoStageParameters& model_parameters,
                          const stoMap& stochastic_map,
                          std::mt19937& generator,
                          int dataPoint_idx,
                          const std::vector<std::vector<dataPoint>>& e_DB,
                          const std::vector<std::vector<dataPoint>>& C_DB,
                          int idx_rep = 0);


#endif /* SAA_solver_hpp */
