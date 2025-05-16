//
//  DataStructure.hpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/17/25.
//

#ifndef DataStructure_hpp
#define DataStructure_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <ctime>
#include <stdlib.h>
#include <cassert>
#include <unordered_map>
#include <map>
#include <utility>

#include "gurobi_c++.h"

#include "SparseMatrix.hpp"

// ****************************************************
// Target: SD and Cutting Plane Method
// dual multipliers
struct dualMultipliers {
    std::vector<double> dual;
    std::vector<double> sto_pi_e; // store temporary stochastic pi e values
    double l1_norm;
    //double optimal_value = 0;
    bool feasible_flag;
};

// feasibility cut $a^\top x \leq b$
struct feasibilityCut {
    std::vector<double> A_newRow;
    double b_newRow;
};

struct masterOutput {
    std::vector<double> x;
    double eta; // value of the epigraphical representation at x_candidate
    std::vector<double> dual_minorant;
    unsigned int sol_flag;
};

// parameters of standard two stage stochastic linear program
// A x <= b
// Dy = e - Cx
struct standardTwoStageParameters {
    // first stage problem
    SparseVector c;
    SparseMatrix Q; // quadratic positive definite
    SparseMatrix AE; // A x = b
    SparseMatrix AL; // A x <= b
    SparseMatrix AG; // A x >= b
    SparseVector bE;
    SparseVector bL;
    SparseVector bG;
    // second stage subproblem
    // min d^\top y s.t. D y + C x + ISlack s = e , y >= 0, s >= 0
    // y and s share the indices
    SparseVector d;
    SparseMatrix P; // quadratic positive definite, include slack variables
    SparseMatrix D;
    SparseMatrix C;
    SparseMatrix ISlack; // s
    SparseVector e; // right hand side deterministic part
    // intermediate matrices
    SparseMatrix Pinv;
    SparseMatrix Pinv_Dt;
    SparseMatrix D_Pinv_Dt; // here, D' = [D Islack]
    SparseMatrix D_Pinv;
    // intermediate scalar
    double dt_Pinv_d;
    // intermediate vectors
    SparseVector dt_Pinv_Dt;
    SparseVector dt_Pinv;
    // extra paramters
    long num_E_1stStage = 0;
    long num_L_1stStage = 0;
    long num_G_1stStage = 0;
    long num_E_2ndStage = 0;
    long num_L_2ndStage = 0;
    long num_G_2ndStage = 0;
    long num_var_1stStage = 0;
    long num_var_2ndStage = 0;
    long num_slack_2ndStage = 0;
    double var_lb_1stStage = 0;
    double var_lb_2ndStage = 0;
    double var_ub_1stStage = 1e9;
};


// ****************************************************
// Target:ioStochastic
// data structure
struct randomVector {
    std::vector<double> component; // all the entries of a random vector
    std::vector<int> randomIndices; // indices of random entries
};

struct randomScalar {
    double component = 0;
    bool flag_random = false; // flag which tells whether this scalar is random
};

// random vector xi
struct random_vector_xi {
    // first stage
    std::vector<double> sto_log_c;
    // second stage
    std::vector<double> be;
    std::vector<double> bi;
    std::vector<double> Ce;
    std::vector<double> Ci;
};

// vectors on the right hand side of second stage problem
struct secondStageRHS {
    randomVector be;
    randomVector bi;
};
// database of vectors on the right hand side of second stage
struct secondStageRHSDB {
    std::vector<std::vector<std::vector<double>>> be_database;
    std::vector<std::vector<std::vector<double>>> bi_database;
    std::vector<std::vector<std::vector<double>>> Ce_database;
};

struct secondStageRHSpoint {
    std::vector<double> be;
    std::vector<double> bi;
    std::vector<double> Ce;
    std::vector<double> Ci;
};

// store the location of randomness
struct secondStageRHSmap {
    std::vector<int> be_map;
    std::vector<int> bi_map;
    std::vector<std::pair<int,int>> Ce_map;
    std::vector<std::pair<int,int>> Ci_map;
};

// Target: ioStochastic
// discrete distribution
struct discreteDistribution {
    std::vector<double> value;
    std::vector<double> probability;
    std::vector<double> cdf;
};

struct stoMap {
    std::vector<int> e_map;
    std::vector<std::pair<int,int>> C_map;
    // --- comment out for now ---
    // set to store the indices which are random, need them to perform priority check
    // when solving 2nd stage subproblem
    // std::set<int> e_set;
    // std::set<std::pair<int,int>> C_set;
    // --- comment out for now ---
    std::string e_flag_dist = "DATA"; // default setting is reading the realizations of random variables from the data file
    std::string C_flag_dist = "DATA"; // default setting is reading the realizations of random variables from the data file
    // second stage
    // e
    std::vector<double> e_normal_mean;
    std::vector<double> e_normal_stddev;
    std::vector<double> e_normal_lb;
    std::vector<double> e_normal_ub;
    std::vector<double> e_uniform_lower;
    std::vector<double> e_uniform_upper;
    std::vector<discreteDistribution> e_discrete;
    // C
    std::vector<double> C_normal_mean;
    std::vector<double> C_normal_stddev;
    std::vector<double> C_normal_lb;
    std::vector<double> C_normal_ub;
    std::vector<double> C_uniform_lower;
    std::vector<double> C_uniform_upper;
    std::vector<discreteDistribution> C_discrete;
};
// ****************************************************
// Target: SAA_solver
struct stoPoint {
    std::vector<double> e;
    std::vector<double> C;
    double prob = 0;
};

// ****************************************************
// Target: ioDB
// data structures for the dataPoint
struct dataPoint { // definition of dataPoint
    std::vector<double> response;
    double weight;
};

// ****************************************************
// Target: SAA, BD_solver, single-cut Benders decomposition, SD_solver
struct solverOutput {
    std::vector<double> sol;
    std::vector<stoPoint> Xi;
    double obj;
    int num_it = 0;
    double sol_time = 0;
};

struct compromiseOutput {
    std::vector<double> compromise_decision;
    std::vector<solverOutput> replications;
};

struct bd_solution {
    std::vector<double> x;
    double z;
    unsigned int sol_flag = 0; // 0 optimal solution exists; 1 infeasible or unbounded
};

struct bd_subproblem_cuts {
    std::vector<double> alpha_array;
    std::vector<std::vector<double>> beta_array;
};

struct bd_output {
    std::vector<double> x;
    std::vector<stoPoint> Xi;
    bd_subproblem_cuts subproblem_constraints;
    double max_gap;
    double sol_time;
    unsigned int sol_flag;
    int it_num;
    int cuts_count;
};

struct compromise_bd_output {
    std::vector<double> compromise_decision;
    std::vector<bd_output> replications;
};

// Target: CuttingPlane_solver
struct cp_solution {
    std::vector<double> x;
    double z;
    double obj;
    unsigned int sol_flag = 0; // 0 optimal solution exists; 1 infeasible or unbounded
};

struct cp_subproblem_cuts {
    std::vector<double> alpha_array;
    std::vector<std::vector<double>> beta_array;
};

struct cp_output {
    std::vector<double> x;
    std::vector<stoPoint> Xi;
    cp_subproblem_cuts subproblem_constraints;
    double max_gap;
    double sol_time;
    unsigned int sol_flag;
    int it_num;
    int cuts_count;
};

struct compromise_cp_output {
    std::vector<double> compromise_decision;
    std::vector<cp_output> replications;
};


// ****************************************************
// Target: SD_solver
// minorant
struct minorant {
    double alpha = 0;
    std::vector<double> beta;
    bool if_active = true;
};

struct sd_output {
    std::vector<double> sol;
    std::vector<minorant> minorant_collection;
    int it_num;
    double obj;
    double sol_time;
};

struct compromise_sd_output {
    std::vector<double> compromise_decision;
    std::vector<std::vector<double>> replication_decision;
};

// ****************************************************
// Target: SQQP Solver
// dual multipliers
struct dualMultipliersQP {
    std::vector<double> s;
    std::vector<double> t;
    double obj;
};

// faces
struct face {
    std::vector<int> axis; // store indices where the components of s are 0, e.g.,
                            // axis = [0,3], then s[0] = 0, s[3] = 0
    long dim = 0;
};

// functions
double operator*(const std::vector<double>& vec1, const std::vector<double>& vec2);

std::vector<double> operator*(double a, const std::vector<double>& vec1);

std::vector<double> operator+(const std::vector<double>& vec1, const std::vector<double>& vec2);

double max(double a, double b);

double min(double a, double b);

// Functions for testing
void SparseMatrix_demo_test1();

#endif /* DataStructure_hpp */
