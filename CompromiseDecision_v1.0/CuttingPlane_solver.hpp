//
//  CuttingPlane_solver.hpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/27/25.
//

#ifndef CuttingPlane_solver_hpp
#define CuttingPlane_solver_hpp

#include <stdio.h>
#include "SD_solver.hpp"

// master problem
cp_solution cp_master(standardTwoStageParameters& model_parameters,
                      cp_subproblem_cuts& subproblem_constraints);

compromise_cp_output cp_compromise_master(standardTwoStageParameters& model_parameters,
                                         std::vector<cp_output>& aggregate_outputs,
                                         double regularizer_coefficient);

// cutting plane algorithm (assume relative complete recourse)
cp_output cp_solver(const std::string& folder_path,
                    std::mt19937& generator,
                    int sample_size,
                    double error);

// compromise decisioin
compromise_cp_output cp_compromise_solver(const std::string& folder_path,
                                          std::mt19937& generator,
                                          int sample_size_per_replication,
                                          int num_replications,
                                          double regularizer_coefficient,
                                          double error,
                                          double error2);

#endif /* CuttingPlane_solver_hpp */
