//
//  BD_solver.hpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/21/25.
//

#ifndef BD_solver_hpp
#define BD_solver_hpp

#include <stdio.h>
#include "SAA_solver.hpp"

// master problem
bd_solution bd_master(standardTwoStageParameters& model_parameters,
                      bd_subproblem_cuts& subproblem_constraints);

// compromise decision master
compromise_bd_output bd_compromise_master(standardTwoStageParameters& model_parameters,
                                 std::vector<bd_output>& aggregate_outputs,
                                   double regularizer_coefficient);

// single-cut Benders decomposition (assume relative complete recourse)
bd_output bd_solver(const std::string& folder_path,
                    std::mt19937& generator,
                    int sample_size,
                    double error);

// compromise decision
compromise_bd_output bd_compromise_solver(const std::string& folder_path,
                                          std::mt19937& generator,
                                          int sample_size_per_replication,
                                          int num_replications,
                                          double regularizer_coefficient,
                                          double error,
                                          double error2);

#endif /* BD_solver_hpp */
