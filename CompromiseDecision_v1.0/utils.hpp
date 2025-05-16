//
//  utils.hpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 5/12/25.
//

#ifndef utils_hpp
#define utils_hpp

#include <stdio.h>

#include "CuttingPlane_solver.hpp"

// two-stage stochastic linear program
// sample average approximation
void interface_compromise_saa(const std::string& folder_path,
                              const std::string& validation_folder_path,
                              int random_seed,
                              int sample_size,
                              int num_replications,
                              double regularizer_coefficient,
                              bool saa_validation = false);

// Benders decomposition
void interface_compromise_bd(const std::string& folder_path,
                             const std::string& validation_folder_path,
                             int random_seed,
                             int sample_size,
                             int num_replications,
                             double regularizer_coefficient,
                             double error,
                             double error2,
                             bool saa_validation = false);

// SQQP
// sample average approximation
void interface_compromise_saa_sqqp(const std::string& folder_path,
                                   const std::string& validation_folder_path,
                                   int random_seed,
                                   int sample_size,
                                   int num_replications,
                                   double regularizer_coefficient,
                                   bool saa_validation = false);

// cutting plane algorithm
void interface_compromise_cp(const std::string& folder_path,
                             const std::string& validation_folder_path,
                             int random_seed,
                             int sample_size,
                             int num_replications,
                             double regularizer_coefficient,
                             double error,
                             double error2,
                             bool saa_validation = false);

// stochastic decomposition
void interface_compromise_sd(const std::string& folder_path,
                             const std::string& validation_folder_path,
                             int random_seed,
                             stoPoint& point_est,
                             int max_iteration,
                             int num_replications,
                             double f_lowerbound,
                             double sigma_init,
                             double regularizer_coefficient,
                             double lb_error,
                             bool saa_validation = false);

// generate saa set
void saa_set_generator(const std::string& folder_path,
                       int random_seed,
                       int sample_size);

void test_gurobi01();
void test_basics01();
void test_solver01();
#endif /* utils_hpp */
