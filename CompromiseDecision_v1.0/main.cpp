//
//  main.cpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/17/25.
//

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <random>
#include "utils.hpp"
#include "gurobi_c++.h"

// simpleQP
void test_simpleQP() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleQP";
    std::mt19937 engine(41); // Seed with a fixed value (e.g., 41)
    //sqqp_ground_truth_discrete(folderPath);
    // SAA
    solverOutput res = saa_sqqp(folderPath, engine, 100);
    //sqqp_validation_ground_truth_discrete(folderPath, res.sol);
    
    // test on cutting plane algorithm
    //cp_output res = cp_solver(folderPath, engine, 10, 1e-3);
    // validate
    //sqqp_validation_ground_truth_discrete(folderPath, res.x);
    
    // test on SD for solving SQQP
    //stoPoint point_est;
    //point_est.e.push_back(11.6);
    //sd_output sd_res = sd_sqqp_solver(folderPath, engine, point_est, 10, 0.0, 1.0);
    //sqqp_validation_ground_truth_discrete(folderPath, sd_res.sol);
}

// --- lands ---
void cd_saa_lands() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/lands";
    std::string validationFolderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/lands/validation";
    int random_seed = 41;
    int sample_size = 10;
    int num_replication = 5;
    double regularizer_coefficient = 1.0;
    interface_compromise_saa(folderPath, validationFolderPath, random_seed, sample_size, num_replication, regularizer_coefficient);
}

void cd_bd_lands() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/lands";
    std::string validationFolderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/lands/validation";
    int random_seed = 41;
    int sample_size = 10;
    int num_replication = 5;
    double regularizer_coefficient = 1.0;
    double error = 1e-3;
    double error2 = 1e-4;
    interface_compromise_bd(folderPath, validationFolderPath, random_seed, sample_size, num_replication, regularizer_coefficient, error, error2);
}

// --- simpleLands ---
void ground_truth_simpleLands() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleLands/validation";
    ground_truth_discrete(folderPath);
    
}

void cd_saa_simpleLands() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleLands";
    std::string validationFolderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleLands/validation";
    int random_seed = 41;
    int sample_size = 20;
    int num_replication = 5;
    double regularizer_coefficient = 1.0;
    interface_compromise_saa(folderPath, validationFolderPath, random_seed, sample_size, num_replication, regularizer_coefficient);
}

// simpleLands2
void ground_truth_simpleLands2() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleLands2/validation";
    ground_truth_discrete(folderPath);
    
}

void cd_saa_simpleLands2() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleLands2";
    std::string validationFolderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleLands2/validation";
    int random_seed = 41;
    int sample_size = 10;
    int num_replication = 5;
    double regularizer_coefficient = 1.0;
    interface_compromise_saa(folderPath, validationFolderPath, random_seed, sample_size, num_replication, regularizer_coefficient);
}

// --- simpleQP ---
void cd_saa_simpleQP() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleQP";
    std::string validationFolderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleQP/validation";
    int random_seed = 41;
    int sample_size = 10;
    int num_replication = 5;
    double regularizer_coefficient = 1.0;
    interface_compromise_saa_sqqp(folderPath, validationFolderPath, random_seed, sample_size, num_replication, regularizer_coefficient);
}

void cd_cp_simpleQP() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleQP";
    std::string validationFolderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleQP/validation";
    int random_seed = 41;
    int sample_size = 10;
    int num_replication = 5;
    double regularizer_coefficient = 1.0;
    double error = 1e-3;
    double error2 = 1e-4;
    interface_compromise_cp(folderPath, validationFolderPath, random_seed, sample_size, num_replication, regularizer_coefficient, error, error2);
}

void cd_sd_simpleQP() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleQP";
    std::string validationFolderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleQP/validation";
    int random_seed = 41;
    int max_iteration = 10;
    int num_replication = 5;
    double regularizer_coefficient = 1.0;
    double lb_error = 1e-3;
    stoPoint point_est;
    point_est.e.push_back(11.6);
    double f_lowerbound = 0.0;
    double sigma_init = 1.0;
    interface_compromise_sd(folderPath, validationFolderPath, random_seed, point_est, max_iteration, num_replication, f_lowerbound, sigma_init, regularizer_coefficient, lb_error);
}

// saa set generator
void saa_set() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/retail/saa_validation";
    int random_seed = 11279013;
    int sample_size = 3e4;
    saa_set_generator(folderPath, random_seed, sample_size);
}

// --- retail ---
// saa_retail
void saa_validation_retail() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/retail/saa_validation";
    ground_truth_discrete(folderPath);
}

void cd_saa_retail() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/retail";
    std::string validationFolderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/retail/saa_validation";
    int random_seed = 41;
    int sample_size = 10;
    int num_replication = 5;
    double regularizer_coefficient = 1.0;
    interface_compromise_saa(folderPath, validationFolderPath, random_seed, sample_size, num_replication, regularizer_coefficient);
}

void cd_bd_retail() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/retail";
    std::string validationFolderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/retail/saa_validation";
    int random_seed = 41;
    int sample_size = 10;
    int num_replication = 5;
    double regularizer_coefficient = 1.0;
    double error = 1e-3;
    double error2 = 1e-4;
    interface_compromise_bd(folderPath, validationFolderPath, random_seed, sample_size, num_replication, regularizer_coefficient, error, error2);
}

// --- pgp2 ---
// validation using true dist
void validation_pgp2() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/pgp2/validation";
    ground_truth_discrete(folderPath);
}

void cd_saa_pgp2() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/pgp2";
    std::string validationFolderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/pgp2/validation";
    int random_seed = 41;
    int sample_size = 10;
    int num_replication = 5;
    double regularizer_coefficient = 1.0;
    interface_compromise_saa(folderPath, validationFolderPath, random_seed, sample_size, num_replication, regularizer_coefficient);
}

void cd_saa_pgp2_batch() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/pgp2";
    std::string validationFolderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/pgp2/validation";
    
    int N = 100;
    int num_replication = 20;
    int sample_size = N / num_replication;
    
    double K = 0.1;
    double regularizer_coefficient = K * ((double) sample_size);
    
    for (int rng_seed = 1; rng_seed < 91; ++rng_seed) {
        interface_compromise_saa(folderPath, validationFolderPath, rng_seed, sample_size, num_replication, regularizer_coefficient);
    }
}


void cd_bd_pgp2() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/pgp2";
    std::string validationFolderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/pgp2/validation";
    int random_seed = 41;
    int sample_size = 10;
    int num_replication = 5;
    double regularizer_coefficient = 1.0;
    double error = 1e-3;
    double error2 = 1e-4;
    interface_compromise_bd(folderPath, validationFolderPath, random_seed, sample_size, num_replication, regularizer_coefficient, error, error2);
}

// --- baa99 ---
void validation_baa99() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/baa99/validation";
    ground_truth_discrete(folderPath);
}

void cd_saa_baa99() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/baa99";
    std::string validationFolderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/baa99/validation";
    int random_seed = 53; // 41 53 79 89 93 101 113 179 213 227 231 273 289 313 357
    int sample_size = 50;
    int num_replication = 1;
    double regularizer_coefficient = 1.0;
    interface_compromise_saa(folderPath, validationFolderPath, random_seed, sample_size, num_replication, regularizer_coefficient);
}

void cd_saa_baa99_batch() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/baa99";
    std::string validationFolderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/baa99/validation";
    
    int N = 100;
    int num_replication = 20;
    double K = 0.1;
    
    int sample_size = N / num_replication;
    double regularizer_coefficient = K * ((double) sample_size);
    
    /*
    int random_seed_list[] = {2, 3, 5, 7, 11, 13, 17, 19, 41, 53, 79, 89, 93, 101, 113, 179, 213, 227, 231, 273, 289, 313, 357, 92828, 2293, 255729, 611165, 728615, 122889, 360084, 375109, 616803, 458554, 492987, 123462, 50001, 50004, 50009,
        500010, 500014, 500018, 500022, 500029, 500052, 500083, 500088, 500091, 5000111,
        5000115, 5000167, 5000177, 5000215, 5000235, 5000272, 5000285, 5000361, 5000441, 5000447, 5000449, 5000479, 5000487, 5000491, 5000503, 5000507, 5000515};
    */
    for (int rng_seed = 1; rng_seed < 91; ++rng_seed) {
        interface_compromise_saa(folderPath, validationFolderPath, rng_seed, sample_size, num_replication, regularizer_coefficient);
    }
    
}

void cd_bd_baa99() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/baa99";
    std::string validationFolderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/baa99/validation";
    int random_seed = 41;
    int sample_size = 10;
    int num_replication = 5;
    double regularizer_coefficient = 1.0;
    double error = 1e-3;
    double error2 = 1e-4;
    interface_compromise_bd(folderPath, validationFolderPath, random_seed, sample_size, num_replication, regularizer_coefficient, error, error2);
}

void cd_bd_baa99_batch() {
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/baa99";
    std::string validationFolderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/baa99/validation";
    int N = 100;
    int num_replication = 20;
    double K = 0.1;
    
    int sample_size = N / num_replication;
    double regularizer_coefficient = K * ((double) sample_size);
    
    double error = 1e-3;
    double error2 = 1e-4;
    
    /*
    int random_seed_list[] = {2, 3, 5, 7, 11, 13, 17, 19, 41, 53, 79, 89, 93, 101, 113, 179, 213, 227, 231, 273, 289, 313, 357, 92828, 2293, 255729, 611165, 728615, 122889, 360084, 375109, 616803, 458554, 492987, 123462, 50001, 50004, 50009,
        500010, 500014, 500018, 500022, 500029, 500052, 500083, 500088, 500091, 5000111,
        5000115, 5000167, 5000177, 5000215, 5000235, 5000272, 5000285, 5000361, 5000441, 5000447, 5000449, 5000479, 5000487, 5000491, 5000503, 5000507, 5000515};
    */
    
    for (int rng_seed = 1; rng_seed < 91; ++rng_seed) {
        interface_compromise_bd(folderPath, validationFolderPath, rng_seed, sample_size, num_replication, regularizer_coefficient, error, error2);
    }
    
}

int main(int argc, const char * argv[]) {
    
    // test on simleQP, previous model parameters in lands are not suitable for QP
    //test_simpleQP();
    //cd_saa_lands();
    //cd_bd_lands();
    //cd_saa_simpleQP();
    //cd_cp_simpleQP();
    //cd_sd_simpleQP();
    //ground_truth_simpleLands();
    //cd_saa_simpleLands();
    
    //ground_truth_simpleLands2();
    //cd_saa_simpleLands2();
    
    //saa_set();
    //saa_validation_retail();
    //cd_saa_retail();
    //cd_bd_retail();
    
    // --- pgp2 ---
    //validation_pgp2();
    //cd_saa_pgp2();
    //cd_bd_pgp2();
    cd_saa_pgp2_batch();
    
    // --- baa99 ---
    //validation_baa99();
    //cd_saa_baa99();
    //cd_bd_baa99();
    //cd_saa_baa99_batch();
    //cd_bd_baa99_batch();
    return 0;
}
