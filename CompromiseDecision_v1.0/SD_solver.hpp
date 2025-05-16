//
//  SD_solver.hpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/23/25.
//

#ifndef SD_solver_hpp
#define SD_solver_hpp

#include <stdio.h>
#include "BD_solver.hpp"

// supplemental fucntions
bool if_vec_equal(const std::vector<double>& vec1, const std::vector<double>& vec2);

bool if_vec_equal(const std::vector<int>& vec1, const std::vector<int>& vec2);

bool if_face_equal(const face& face1, const face& face2);

bool if_face_new(const std::vector<face>& face_collection, const face& face_candidate);

face compute_face(const std::vector<double> s);

double compute_obj(std::vector<double> x, standardTwoStageParameters& model_parameters, std::vector<minorant>& minorant_collection);

std::vector<double> sd_compromise_master(standardTwoStageParameters& model_parameters,
                                         const stoMap& stochastic_map,
                                         std::vector<sd_output> aggregate_outputs,
                                         double regularizer_coefficient,
                                         int num_replications,
                                         double lb_error);
// --- SQLP ---
// presolve
std::vector<double> sd_sqlp_presolve(standardTwoStageParameters& model_parameters,
                                 const stoMap& stochastic_map,
                                 stoPoint& xi);

masterOutput sd_master(standardTwoStageParameters& model_parameters,
                                   std::vector<minorant>& minorant_collection,
                                   std::vector<double>& x_incumbent,
                                   double sigma);

sd_output sd_sqlp_solver(const std::string& folder_path,
                            std::mt19937& generator,
                            stoPoint& point_est,
                            int max_iterations,
                            double f_lowerbound,
                            double sigma_init);

compromise_sd_output compromise_sd_sqlp_solver(const std::string& folder_path,
                                               std::mt19937& generator,
                                               stoPoint& point_est,
                                               int max_iterations,
                                               double f_lowerbound,
                                               double sigma_init,
                                               double regularizer_coefficient,
                                               int num_replications,
                                               double lb_error);

// --- SQQP ---
std::vector<double> sd_sqqp_presolve(standardTwoStageParameters& model_parameters,
                                     const stoMap& stochastic_map,
                                     stoPoint& xi);

dualMultipliersQP sd_sqqp_2ndStagePrimal(standardTwoStageParameters& model_parameters,
                                         const stoMap& stochastic_map,
                                         stoPoint& xi,
                                         const std::vector<double>& x);

dualMultipliersQP sd_sqqp_2ndStageDual(standardTwoStageParameters& model_parameters,
                                       const stoMap& stochastic_map,
                                       stoPoint& xi,
                                       const std::vector<double>& x);

dualMultipliersQP sd_sqqp_2ndStageDual(standardTwoStageParameters& model_parameters,
                                       const stoMap& stochastic_map,
                                       stoPoint& xi,
                                       const std::vector<double>& x,
                                       face curFace);

sd_output sd_sqqp_solver(const std::string& folder_path,
                            std::mt19937& generator,
                            stoPoint& point_est,
                            int max_iterations,
                            double f_lowerbound,
                            double sigma_init);

compromise_sd_output compromise_sd_sqqp_solver(const std::string& folder_path,
                                               std::mt19937& generator,
                                               stoPoint& point_est,
                                               int max_iterations,
                                               double f_lowerbound,
                                               double sigma_init,
                                               double regularizer_coefficient,
                                               int num_replications,
                                               double lb_error);

#endif /* SD_solver_hpp */
