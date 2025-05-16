//
//  utils.cpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 5/12/25.
//

#include "utils.hpp"

void interface_compromise_saa(const std::string& folder_path,
                              const std::string& validation_folder_path,
                              int random_seed,
                              int sample_size,
                              int num_replications,
                              double regularizer_coefficient,
                              bool saa_validation) {
    std::mt19937 engine(random_seed); // random seed with a fixed value (e.g., 41)
    
    // compute compromise decisions
    compromiseOutput compromise_saa = saa_2slp_compromise(folder_path,
                                                          engine,
                                                          sample_size,
                                                          num_replications,
                                                          regularizer_coefficient);
    
    std::string resultsOutput_path = folder_path + "/summary(cd_saa_2slpv1.0).txt";
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    // write results to CSV files
    std::string writeCD_path = folder_path + "/compromise_saa_" + std::to_string(num_replications) + "_" + std::to_string(sample_size) + ".csv";
    const char* writeCD_path_const = writeCD_path.c_str();
    std::fstream writeCD;
    writeCD.open(writeCD_path_const,std::fstream::app);
    //
    std::string writeRep_path = folder_path + "/replication_saa_" + std::to_string(num_replications) + "_" + std::to_string(sample_size) + ".csv";
    const char* writeRep_path_const = writeRep_path.c_str();
    std::fstream writeRep;
    writeRep.open(writeRep_path_const,std::fstream::app);
    
    std::cout << "=================================\n";
    writeFile << "=================================\n";
    std::cout << "SAA is used in the replication step.\n";
    writeFile << "SAA is used in the replication step.\n";
    std::cout << "Problem Type: Two-Stage Stochastic Linear Program.\n";
    writeFile << "Problem Type: Two-Stage Stochastic Linear Program.\n";
    std::cout << "[seed: " << random_seed << "]\n";
    writeFile << "[seed: " << random_seed << "]\n";
    std::cout << "Sample Size per Replication: " << sample_size << std::endl;
    writeFile << "Sample Size per Replication: " << sample_size << std::endl;
    std::cout << "Number of Replications: " << num_replications << std::endl;
    writeFile << "Number of Replications: " << num_replications << std::endl;
    std::cout << "Regularizer: " << regularizer_coefficient << std::endl;
    writeFile << "Regularizer: " << regularizer_coefficient << std::endl;
    std::cout << "*** Compromise Decision Results ***\n";
    writeFile << "*** Compromise Decision Results ***\n";
    std::cout << "Compromise Decision: ";
    writeFile << "Compromise Decision: ";
    for (int idx_x = 0; idx_x < compromise_saa.compromise_decision.size() - 1; ++idx_x) {
        std::cout << compromise_saa.compromise_decision[idx_x] << ", ";
        writeFile << compromise_saa.compromise_decision[idx_x] << ", ";
        writeCD << compromise_saa.compromise_decision[idx_x]  << ", ";
    }
    std::cout << compromise_saa.compromise_decision[compromise_saa.compromise_decision.size() - 1] << std::endl;
    writeFile << compromise_saa.compromise_decision[compromise_saa.compromise_decision.size() - 1] << std::endl;
    writeCD << compromise_saa.compromise_decision[compromise_saa.compromise_decision.size() - 1] << ", *, ";
    // validation
    double val = 0;
    if (saa_validation == true) {
        val = validation_saa(validation_folder_path, compromise_saa.compromise_decision);
    }
    else {
        val = validation_ground_truth_discrete(validation_folder_path, compromise_saa.compromise_decision);
    }
    std::cout << "Validation Cost: " << val << std::endl;
    writeFile << "Validation Cost: " << val << std::endl;
    writeCD << val << ", *, " << random_seed << std::endl;
    writeCD.close(); // close the csv file writer for compromise decisions
    std::cout << "*** Replication Results ***\n";
    writeFile << "*** Replication Results ***\n";
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        std::cout << "---------------------------------\n";
        writeFile << "---------------------------------\n";
        std::cout << "Replication " << idx_rep << " Decision : ";
        writeFile << "Replication " << idx_rep << " Decision : ";
        for (int idx_x = 0; idx_x < compromise_saa.replications[idx_rep].sol.size() - 1; ++idx_x) {
            std::cout << compromise_saa.replications[idx_rep].sol[idx_x] << ", ";
            writeFile << compromise_saa.replications[idx_rep].sol[idx_x] << ", ";
            writeRep << compromise_saa.replications[idx_rep].sol[idx_x] << ", ";
        }
        std::cout << compromise_saa.replications[idx_rep].sol[compromise_saa.replications[idx_rep].sol.size() - 1] << std::endl;
        writeFile << compromise_saa.replications[idx_rep].sol[compromise_saa.replications[idx_rep].sol.size() - 1] << std::endl;
        writeRep << compromise_saa.replications[idx_rep].sol[compromise_saa.replications[idx_rep].sol.size() - 1] << ", *, ";
        // validation
        double val_rep = 0;
        if (saa_validation == true) {
            val_rep = validation_saa(validation_folder_path, compromise_saa.replications[idx_rep].sol);
        }
        else {
            val_rep = validation_ground_truth_discrete(validation_folder_path, compromise_saa.replications[idx_rep].sol);
        }
        std::cout << "Replication " << idx_rep << " Validation Cost : " << val_rep << std::endl;
        writeFile << "Replication " << idx_rep << " Validation Cost : " << val_rep << std::endl;
        writeRep << val_rep << ", *, " << random_seed << ", " << idx_rep << std::endl;
        std::cout << "---------------------------------\n";
        writeFile << "---------------------------------\n";
    }
    std::cout << "=================================\n";
    writeFile << "=================================\n";
    writeFile.close();
    writeRep.close();
}


// Benders decomposition
void interface_compromise_bd(const std::string& folder_path,
                             const std::string& validation_folder_path,
                             int random_seed,
                             int sample_size,
                             int num_replications,
                             double regularizer_coefficient,
                             double error,
                             double error2,
                             bool saa_validation) {
    std::mt19937 engine(random_seed); // random seed with a fixed value (e.g., 41)
    
    // compute compromise decisions
    compromise_bd_output compromise_bd = bd_compromise_solver(folder_path, engine, sample_size, num_replications, regularizer_coefficient, error, error2);
    
    std::string resultsOutput_path = folder_path + "/summary(cd_bd_2slpv1.0).txt";
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    
    // write results to CSV files
    std::string writeCD_path = folder_path + "/compromise_bd_" + std::to_string(num_replications) + "_" + std::to_string(sample_size) + ".csv";
    const char* writeCD_path_const = writeCD_path.c_str();
    std::fstream writeCD;
    writeCD.open(writeCD_path_const,std::fstream::app);
    //
    std::string writeRep_path = folder_path + "/replication_bd_" + std::to_string(num_replications) + "_" + std::to_string(sample_size) + ".csv";
    const char* writeRep_path_const = writeRep_path.c_str();
    std::fstream writeRep;
    writeRep.open(writeRep_path_const,std::fstream::app);
    std::cout << "=================================\n";
    writeFile << "=================================\n";
    std::cout << "Benders Decomposition is used in the replication step.\n";
    writeFile << "Benders Decomposition is used in the replication step.\n";
    std::cout << "Problem Type: Two-Stage Stochastic Linear Program.\n";
    writeFile << "Problem Type: Two-Stage Stochastic Linear Program.\n";
    std::cout << "[seed: " << random_seed << "]\n";
    writeFile << "[seed: " << random_seed << "]\n";
    std::cout << "Sample Size per Replication: " << sample_size << std::endl;
    writeFile << "Sample Size per Replication: " << sample_size << std::endl;
    std::cout << "Number of Replications: " << num_replications << std::endl;
    writeFile << "Number of Replications: " << num_replications << std::endl;
    std::cout << "Regularizer: " << regularizer_coefficient << std::endl;
    writeFile << "Regularizer: " << regularizer_coefficient << std::endl;
    std::cout << "Error in the Termination Criterion of the Replication Step: " << error << std::endl;
    writeFile << "Error in the Termination Criterion of the Replication Step: " << error << std::endl;
    std::cout << "Error in the Epsilon-Subgradient Generation of the Function Augmentation: " << error2 << std::endl;
    writeFile << "Error in the Epsilon-Subgradient Generation of the Function Augmentation: " << error2 << std::endl;
    std::cout << "*** Compromise Decision Results ***\n";
    writeFile << "*** Compromise Decision Results ***\n";
    std::cout << "Compromise Decision: ";
    writeFile << "Compromise Decision: ";
    for (int idx_x = 0; idx_x < compromise_bd.compromise_decision.size() - 1; ++idx_x) {
        std::cout << compromise_bd.compromise_decision[idx_x] << ", ";
        writeFile << compromise_bd.compromise_decision[idx_x] << ", ";
        writeCD << compromise_bd.compromise_decision[idx_x]  << ", ";
    }
    std::cout << compromise_bd.compromise_decision[compromise_bd.compromise_decision.size() - 1] << std::endl;
    writeFile << compromise_bd.compromise_decision[compromise_bd.compromise_decision.size() - 1] << std::endl;
    writeCD << compromise_bd.compromise_decision[compromise_bd.compromise_decision.size() - 1] << ", *, ";
    // validation
    double val = 0;
    if (saa_validation == true) {
        val = validation_saa(validation_folder_path, compromise_bd.compromise_decision);
    }
    else {
        val = validation_ground_truth_discrete(validation_folder_path, compromise_bd.compromise_decision);
    }
    std::cout << "Validation Cost: " << val << std::endl;
    writeFile << "Validation Cost: " << val << std::endl;
    writeCD << val << ", *, " << random_seed << std::endl;
    writeCD.close(); // close the csv file writer for compromise decisions
    std::cout << "*** Replication Results ***\n";
    writeFile << "*** Replication Results ***\n";
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        std::cout << "---------------------------------\n";
        writeFile << "---------------------------------\n";
        std::cout << "Replication " << idx_rep << " Decision : ";
        writeFile << "Replication " << idx_rep << " Decision : ";
        for (int idx_x = 0; idx_x < compromise_bd.replications[idx_rep].x.size() - 1; ++idx_x) {
            std::cout << compromise_bd.replications[idx_rep].x[idx_x] << ", ";
            writeFile << compromise_bd.replications[idx_rep].x[idx_x] << ", ";
            writeRep << compromise_bd.replications[idx_rep].x[idx_x] << ", ";
        }
        std::cout << compromise_bd.replications[idx_rep].x[compromise_bd.replications[idx_rep].x.size() - 1] << std::endl;
        writeFile << compromise_bd.replications[idx_rep].x[compromise_bd.replications[idx_rep].x.size() - 1] << std::endl;
        writeRep << compromise_bd.replications[idx_rep].x[compromise_bd.replications[idx_rep].x.size() - 1] << ", *, ";
        // validation
        double val_rep = 0;
        if (saa_validation == true) {
            val_rep = validation_saa(validation_folder_path, compromise_bd.replications[idx_rep].x);
        }
        else {
            val_rep = validation_ground_truth_discrete(validation_folder_path, compromise_bd.replications[idx_rep].x);
        }
        std::cout << "Replication " << idx_rep << " Validation Cost : " << val_rep << std::endl;
        writeFile << "Replication " << idx_rep << " Validation Cost : " << val_rep << std::endl;
        writeRep << val_rep << ", *, " << random_seed << ", " << idx_rep << std::endl;
        std::cout << "---------------------------------\n";
        writeFile << "---------------------------------\n";
    }
    std::cout << "=================================\n";
    writeFile << "=================================\n";
    writeFile.close();
    writeRep.close();
}


void interface_compromise_saa_sqqp(const std::string& folder_path,
                                   const std::string& validation_folder_path,
                                   int random_seed,
                                   int sample_size,
                                   int num_replications,
                                   double regularizer_coefficient,
                                   bool saa_validation) {
    std::mt19937 engine(random_seed); // random seed with a fixed value (e.g., 41)
    
    // compute compromise decisions
    compromiseOutput compromise_saa = saa_sqqp_compromise(folder_path,
                                                          engine,
                                                          sample_size,
                                                          num_replications,
                                                          regularizer_coefficient);
    
    std::string resultsOutput_path = folder_path + "/summary(cd_saa_2sqqpv1.0).txt";
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    std::cout << "=================================\n";
    writeFile << "=================================\n";
    std::cout << "SAA is used in the replication step.\n";
    writeFile << "SAA is used in the replication step.\n";
    std::cout << "Problem Type: Two-Stage Stochastic Quadratic-Quadratic Program.\n";
    writeFile << "Problem Type: Two-Stage Stochastic Quadratic-Quadratic Program.\n";
    std::cout << "[seed: " << random_seed << "]\n";
    writeFile << "[seed: " << random_seed << "]\n";
    std::cout << "Sample Size per Replication: " << sample_size << std::endl;
    writeFile << "Sample Size per Replication: " << sample_size << std::endl;
    std::cout << "Number of Replications: " << num_replications << std::endl;
    writeFile << "Number of Replications: " << num_replications << std::endl;
    std::cout << "Regularizer: " << regularizer_coefficient << std::endl;
    writeFile << "Regularizer: " << regularizer_coefficient << std::endl;
    std::cout << "*** Compromise Decision Results ***\n";
    writeFile << "*** Compromise Decision Results ***\n";
    std::cout << "Compromise Decision: ";
    writeFile << "Compromise Decision: ";
    for (int idx_x = 0; idx_x < compromise_saa.compromise_decision.size() - 1; ++idx_x) {
        std::cout << compromise_saa.compromise_decision[idx_x] << ", ";
        writeFile << compromise_saa.compromise_decision[idx_x] << ", ";
    }
    std::cout << compromise_saa.compromise_decision[compromise_saa.compromise_decision.size() - 1] << std::endl;
    writeFile << compromise_saa.compromise_decision[compromise_saa.compromise_decision.size() - 1] << std::endl;
    // validation
    double val = 0;
    if (saa_validation == true) {
        val = sqqp_validation_saa(validation_folder_path, compromise_saa.compromise_decision);
    }
    else {
        val = sqqp_validation_ground_truth_discrete(validation_folder_path, compromise_saa.compromise_decision);
    }
    std::cout << "Validation Cost: " << val << std::endl;
    writeFile << "Validation Cost: " << val << std::endl;
    std::cout << "*** Replication Results ***\n";
    writeFile << "*** Replication Results ***\n";
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        std::cout << "---------------------------------\n";
        writeFile << "---------------------------------\n";
        std::cout << "Replication " << idx_rep << " Decision : ";
        writeFile << "Replication " << idx_rep << " Decision : ";
        for (int idx_x = 0; idx_x < compromise_saa.replications[idx_rep].sol.size() - 1; ++idx_x) {
            std::cout << compromise_saa.replications[idx_rep].sol[idx_x] << ", ";
            writeFile << compromise_saa.replications[idx_rep].sol[idx_x] << ", ";
        }
        std::cout << compromise_saa.replications[idx_rep].sol[compromise_saa.replications[idx_rep].sol.size() - 1] << std::endl;
        writeFile << compromise_saa.replications[idx_rep].sol[compromise_saa.replications[idx_rep].sol.size() - 1] << std::endl;
        // validation
        double val_rep = 0;
        if (saa_validation == true) {
            val_rep = sqqp_validation_saa(validation_folder_path, compromise_saa.replications[idx_rep].sol);
        }
        else {
            val_rep = sqqp_validation_ground_truth_discrete(validation_folder_path, compromise_saa.replications[idx_rep].sol);
        }
        std::cout << "Replication " << idx_rep << " Validation Cost : " << val_rep << std::endl;
        writeFile << "Replication " << idx_rep << " Validation Cost : " << val_rep << std::endl;
        std::cout << "---------------------------------\n";
        writeFile << "---------------------------------\n";
    }
    std::cout << "=================================\n";
    writeFile << "=================================\n";
    writeFile.close();
}


// cutting plane algorithm
void interface_compromise_cp(const std::string& folder_path,
                             const std::string& validation_folder_path,
                             int random_seed,
                             int sample_size,
                             int num_replications,
                             double regularizer_coefficient,
                             double error,
                             double error2,
                             bool saa_validation) {
    std::mt19937 engine(random_seed); // random seed with a fixed value (e.g., 41)
    
    // compute compromise decisions
    compromise_cp_output compromise_cp = cp_compromise_solver(folder_path, engine, sample_size, num_replications, regularizer_coefficient, error, error2);
    
    std::string resultsOutput_path = folder_path + "/summary(cd_cp_sqqpv1.0).txt";
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    std::cout << "=================================\n";
    writeFile << "=================================\n";
    std::cout << "Cutting Plane Algorithm is used in the replication step.\n";
    writeFile << "Cutting Plane Algorithm is used in the replication step.\n";
    std::cout << "Problem Type: Two-Stage Stochastic Quadratic-Quadratic Program.\n";
    writeFile << "Problem Type: Two-Stage Stochastic Quadratic-Quadratic Program.\n";
    std::cout << "[seed: " << random_seed << "]\n";
    writeFile << "[seed: " << random_seed << "]\n";
    std::cout << "Sample Size per Replication: " << sample_size << std::endl;
    writeFile << "Sample Size per Replication: " << sample_size << std::endl;
    std::cout << "Number of Replications: " << num_replications << std::endl;
    writeFile << "Number of Replications: " << num_replications << std::endl;
    std::cout << "Regularizer: " << regularizer_coefficient << std::endl;
    writeFile << "Regularizer: " << regularizer_coefficient << std::endl;
    std::cout << "Error in the Termination Criterion of the Replication Step: " << error << std::endl;
    writeFile << "Error in the Termination Criterion of the Replication Step: " << error << std::endl;
    std::cout << "Error in the Epsilon-Subgradient Generation of the Function Augmentation: " << error2 << std::endl;
    writeFile << "Error in the Epsilon-Subgradient Generation of the Function Augmentation: " << error2 << std::endl;
    std::cout << "*** Compromise Decision Results ***\n";
    writeFile << "*** Compromise Decision Results ***\n";
    std::cout << "Compromise Decision: ";
    writeFile << "Compromise Decision: ";
    for (int idx_x = 0; idx_x < compromise_cp.compromise_decision.size() - 1; ++idx_x) {
        std::cout << compromise_cp.compromise_decision[idx_x] << ", ";
        writeFile << compromise_cp.compromise_decision[idx_x] << ", ";
    }
    std::cout << compromise_cp.compromise_decision[compromise_cp.compromise_decision.size() - 1] << std::endl;
    writeFile << compromise_cp.compromise_decision[compromise_cp.compromise_decision.size() - 1] << std::endl;
    // validation
    double val = 0;
    if (saa_validation == true) {
        val = sqqp_validation_saa(validation_folder_path, compromise_cp.compromise_decision);
    }
    else {
        val = sqqp_validation_ground_truth_discrete(validation_folder_path, compromise_cp.compromise_decision);
    }
    std::cout << "Validation Cost: " << val << std::endl;
    writeFile << "Validation Cost: " << val << std::endl;
    std::cout << "*** Replication Results ***\n";
    writeFile << "*** Replication Results ***\n";
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        std::cout << "---------------------------------\n";
        writeFile << "---------------------------------\n";
        std::cout << "Replication " << idx_rep << " Decision : ";
        writeFile << "Replication " << idx_rep << " Decision : ";
        for (int idx_x = 0; idx_x < compromise_cp.replications[idx_rep].x.size() - 1; ++idx_x) {
            std::cout << compromise_cp.replications[idx_rep].x[idx_x] << ", ";
            writeFile << compromise_cp.replications[idx_rep].x[idx_x] << ", ";
        }
        std::cout << compromise_cp.replications[idx_rep].x[compromise_cp.replications[idx_rep].x.size() - 1] << std::endl;
        writeFile << compromise_cp.replications[idx_rep].x[compromise_cp.replications[idx_rep].x.size() - 1] << std::endl;
        // validation
        double val_rep = 0;
        if (saa_validation == true) {
            val_rep = sqqp_validation_saa(validation_folder_path, compromise_cp.replications[idx_rep].x);
        }
        else {
            val_rep = sqqp_validation_ground_truth_discrete(validation_folder_path, compromise_cp.replications[idx_rep].x);
        }
        std::cout << "Replication " << idx_rep << " Validation Cost : " << val_rep << std::endl;
        writeFile << "Replication " << idx_rep << " Validation Cost : " << val_rep << std::endl;
        std::cout << "---------------------------------\n";
        writeFile << "---------------------------------\n";
    }
    std::cout << "=================================\n";
    writeFile << "=================================\n";
    writeFile.close();
}


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
                             bool saa_validation) {
    std::mt19937 engine(random_seed); // random seed with a fixed value (e.g., 41)
    
    // compute compromise decisions
    compromise_sd_output compromise_sd = compromise_sd_sqqp_solver(folder_path, engine, point_est, max_iteration, f_lowerbound, sigma_init, regularizer_coefficient, num_replications, lb_error);
    std::string resultsOutput_path = folder_path + "/summary(cd_sd_sqqpv1.0).txt";
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    std::cout << "=================================\n";
    writeFile << "=================================\n";
    std::cout << "Cutting Plane Algorithm is used in the replication step.\n";
    writeFile << "Cutting Plane Algorithm is used in the replication step.\n";
    std::cout << "Problem Type: Two-Stage Stochastic Quadratic-Quadratic Program.\n";
    writeFile << "Problem Type: Two-Stage Stochastic Quadratic-Quadratic Program.\n";
    std::cout << "[seed: " << random_seed << "]\n";
    writeFile << "[seed: " << random_seed << "]\n";
    std::cout << "Number of Iterations per Replication: " << max_iteration << std::endl;
    writeFile << "Number of Iterations per Replication: " << max_iteration << std::endl;
    std::cout << "Number of Replications: " << num_replications << std::endl;
    writeFile << "Number of Replications: " << num_replications << std::endl;
    std::cout << "Regularizer: " << regularizer_coefficient << std::endl;
    writeFile << "Regularizer: " << regularizer_coefficient << std::endl;
    std::cout << "Error in the Generation of the Function Augmentation (lower bound): " << lb_error << std::endl;
    writeFile << "Error in the Generation of the Function Augmentation (lower bound): " << lb_error << std::endl;
    std::cout << "*** Compromise Decision Results ***\n";
    writeFile << "*** Compromise Decision Results ***\n";
    std::cout << "Compromise Decision: ";
    writeFile << "Compromise Decision: ";
    for (int idx_x = 0; idx_x < compromise_sd.compromise_decision.size() - 1; ++idx_x) {
        std::cout << compromise_sd.compromise_decision[idx_x] << ", ";
        writeFile << compromise_sd.compromise_decision[idx_x] << ", ";
    }
    std::cout << compromise_sd.compromise_decision[compromise_sd.compromise_decision.size() - 1] << std::endl;
    writeFile << compromise_sd.compromise_decision[compromise_sd.compromise_decision.size() - 1] << std::endl;
    // validation
    double val = 0;
    if (saa_validation == true) {
        val = sqqp_validation_saa(validation_folder_path, compromise_sd.compromise_decision);
    }
    else {
        val = sqqp_validation_ground_truth_discrete(validation_folder_path, compromise_sd.compromise_decision);
    }
    std::cout << "Validation Cost: " << val << std::endl;
    writeFile << "Validation Cost: " << val << std::endl;
    std::cout << "*** Replication Results ***\n";
    writeFile << "*** Replication Results ***\n";
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        std::cout << "---------------------------------\n";
        writeFile << "---------------------------------\n";
        std::cout << "Replication " << idx_rep << " Decision : ";
        writeFile << "Replication " << idx_rep << " Decision : ";
        for (int idx_x = 0; idx_x < compromise_sd.replication_decision[idx_rep].size() - 1; ++idx_x) {
            std::cout << compromise_sd.replication_decision[idx_rep][idx_x] << ", ";
            writeFile << compromise_sd.replication_decision[idx_rep][idx_x] << ", ";
        }
        std::cout << compromise_sd.replication_decision[idx_rep][compromise_sd.replication_decision[idx_rep].size() - 1] << std::endl;
        writeFile << compromise_sd.replication_decision[idx_rep][compromise_sd.replication_decision[idx_rep].size() - 1] << std::endl;
        // validation
        double val_rep = 0;
        if (saa_validation == true) {
            val_rep = sqqp_validation_saa(validation_folder_path, compromise_sd.replication_decision[idx_rep]);
        }
        else {
            val_rep = sqqp_validation_ground_truth_discrete(validation_folder_path, compromise_sd.replication_decision[idx_rep]);
        }
        std::cout << "Replication " << idx_rep << " Validation Cost : " << val_rep << std::endl;
        writeFile << "Replication " << idx_rep << " Validation Cost : " << val_rep << std::endl;
        std::cout << "---------------------------------\n";
        writeFile << "---------------------------------\n";
    }
    std::cout << "=================================\n";
    writeFile << "=================================\n";
    writeFile.close();
}
// -----------------------------------------------------
// generate saa set
void saa_set_generator(const std::string& folder_path,
                       int random_seed,
                       int sample_size) {
    std::mt19937 engine(random_seed); // random seed with a fixed value (e.g., 41)
    std::uniform_real_distribution<double> dist_uniform(0.0,1.0);
    bool flag_data_e; // tell if e stochastic is generated
    bool flag_data_C; // tell if C stochastic is generated
    // create directory paths for database and model
    std::string e_DB_path = folder_path + "/e_DB.txt";
    std::string C_DB_path = folder_path + "/C_DB.txt";
    // sto file, model file
    std::string sto_path = folder_path + "/sto.txt";
    std::string model_path = folder_path + "/model.txt";
    // output file path
    std::string output_path = folder_path + "/true_dist.txt";
    const char* writeFilePath = output_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    // convert all the paths into constant chars
    const char* e_DB_path_const = e_DB_path.c_str();
    const char* C_DB_path_const = C_DB_path.c_str();
    // create stream object
    std::ifstream readFile_e(e_DB_path_const);
    std::ifstream readFile_C(C_DB_path_const);
    // create database
    std::vector<std::vector<dataPoint>> e_DB;
    std::vector<std::vector<dataPoint>> C_DB;
    // create sto object
    stoMap stochastic_map = readStochasticMap(sto_path);
    // create model object
    standardTwoStageParameters model_parameters = readStandardTwoStageParameters(model_path);
    // read  be
    if (readFile_e.is_open()) {
        std::cout << "e_DB (RHS in the 2nd stage) data file is found." << std::endl;
        readFile_e.close(); // close the file
        // read be database
        e_DB = readDB(e_DB_path);
        flag_data_e = true;
    }
    else {
        readFile_e.close(); // close the file
        flag_data_e = false;
        std::cout << "e_DB (RHS in the 2nd stage) data file is not found!" << std::endl;
    }
    // read C
    if (readFile_C.is_open()) {
        std::cout << "C_DB (2nd stage) data file is found." << std::endl;
        readFile_C.close(); // close the file
        // Ce database
        C_DB = readDB(C_DB_path);
        flag_data_C = true;
    }
    else {
        readFile_C.close(); // close the file
        flag_data_C = false;
        std::cout << "C_DB (2nd stage) data file is not found!" << std::endl;
    }
    // generate saa file
    double one_over_sample_size = 1.0 / ((double) sample_size);
    writeFile << "<dist>\n";
    // generate data points
    for (int idx_sample = 0; idx_sample < sample_size ; ++idx_sample) {
        stoPoint xi_cur = generator_random(model_parameters, stochastic_map, engine, idx_sample, e_DB, C_DB);
        // e
        writeFile << "e:";
        if (xi_cur.e.size() > 0) {
            for (int idx_e = 0; idx_e < xi_cur.e.size() - 1; ++idx_e) {
                writeFile << xi_cur.e[idx_e] << ",";
            }
            writeFile << xi_cur.e[xi_cur.e.size() - 1] << ";";
        }
        else {
            writeFile << "*;";
        }
        // C
        writeFile << "C:";
        if (xi_cur.C.size() > 0) {
            for (int idx_C = 0; idx_C < xi_cur.C.size() - 1; ++idx_C) {
                writeFile << xi_cur.C[idx_C] << ",";
            }
            writeFile << xi_cur.C[xi_cur.C.size() - 1] << ";";
        }
        else {
            writeFile << "*;";
        }
        // prob
        writeFile << "prob:" << one_over_sample_size << ";\n";
    }
    writeFile << "</dist>\n";
    writeFile.close();
}


// -----------------------------------------------------
void test_gurobi01() {
    try {
            // Create the Gurobi environment
            GRBEnv env = GRBEnv();

            // Create an empty model
            GRBModel model = GRBModel(env);

            // Add decision variables (x and y)
            GRBVar x = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x");
            GRBVar y = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y");

            // Set the objective function: maximize 3x + 4y
            model.setObjective(3 * x + 4 * y, GRB_MAXIMIZE);

            // Add constraints
            GRBConstr c1 = model.addConstr(x + 2 * y <= 4, "c1");
            GRBConstr c2 = model.addConstr(4 * x + 3 * y <= 6, "c2");

            // Optimize the model
            model.optimize();

            // Display results
            std::cout << "Optimal solution:\n";
            std::cout << "x = " << x.get(GRB_DoubleAttr_X) << std::endl;
            std::cout << "y = " << y.get(GRB_DoubleAttr_X) << std::endl;

            // Get and print dual values (shadow prices)
            std::cout << "Dual multipliers (shadow prices):\n";
            std::cout << "Constraint c1: " << c1.get(GRB_DoubleAttr_Pi) << std::endl;
            std::cout << "Constraint c2: " << c2.get(GRB_DoubleAttr_Pi) << std::endl;
        }
        catch (GRBException e) {
            std::cerr << "Gurobi Error code: " << e.getErrorCode() << ": " << e.getMessage() << std::endl;
        }
        catch (...) {
            std::cerr << "Unexpected error occurred" << std::endl;
        }
}

void test_basics01() {
    std::mt19937 engine(41); // Seed with a fixed value (e.g., 41)
    std::uniform_int_distribution<int> dist(1, 100);
        
    int random_num = dist(engine); // Will always produce the same sequence
    std::cout << "Random number 1 (seed1): " << random_num << std::endl;
    random_num = dist(engine); // Will always produce the same sequence
    std::cout << "Random number 2 (seed1): " << random_num << std::endl;
    random_num = dist(engine); // Will always produce the same sequence
    std::cout << "Random number 3 (seed1): " << random_num << std::endl;
     
    std::mt19937 engine2(99); // Seed with a fixed value (e.g., 99)
    random_num = dist(engine2); // Will always produce the same sequence
    std::cout << "Random number 1 (seed2): " << random_num << std::endl;
    random_num = dist(engine2); // Will always produce the same sequence
    std::cout << "Random number 2 (seed2): " << random_num << std::endl;
    random_num = dist(engine2); // Will always produce the same sequence
    std::cout << "Random number 3 (seed2): " << random_num << std::endl;
    
    // case 1
    // M1 = | 1 3 0 |
    //      | 2 0 4 |
    SparseMatrix M1;
    M1.insert(0, 0, 1);
    M1.insert(1, 0, 2);
    M1.insert(0, 1, 3);
    M1.insert(1, 2, 4);
    M1.insert_end(2, 3);
    std::vector<std::vector<double>> M1_dense = M1.scatter();
    for (int row_idx = 0; row_idx < 2; ++row_idx) {
        for (int col_idx = 0; col_idx < 3; ++col_idx) {
            std::cout << M1_dense[row_idx][col_idx] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "***************" << std::endl;
    // case 2
    // M2 = |1 0 0|
    //      |2 0 4|
    SparseMatrix M2;
    M2.insert(0, 0, 1);
    M2.insert(1, 0, 2);
    M2.insert(1, 2, 4);
    M2.insert_end(2, 3);
    std::vector<std::vector<double>> M2_dense = M2.scatter();
    for (int row_idx = 0; row_idx < 2; ++row_idx) {
        for (int col_idx = 0; col_idx < 3; ++col_idx) {
            std::cout << M2_dense[row_idx][col_idx] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "***************" << std::endl;
    // case 3
    // M3 = |0 0 0|
    //      |0 0 4|
    SparseMatrix M3;
    M3.insert(1, 2, 4);
    M3.insert_end(2, 3);
    std::vector<std::vector<double>> M3_dense = M3.scatter();
    for (int row_idx = 0; row_idx < 2; ++row_idx) {
        for (int col_idx = 0; col_idx < 3; ++col_idx) {
            std::cout << M3_dense[row_idx][col_idx] << " ";
        }
        std::cout << std::endl;
    }
}

void test_solver01(){
    // test on SAA_solver
    //std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/lands";
    std::string folderPath = "/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/lands_sqqp";
    std::mt19937 engine(41); // Seed with a fixed value (e.g., 41)
    //saa_2slp(folderPath, engine, 500);
    //solverOutput res = saa_sqqp(folderPath, engine, 50);
    // validate
    //sqqp_validation_ground_truth_discrete(folderPath, res.sol);
    // test on the ground truth
    //ground_truth_discrete(folderPath);
    //sqqp_ground_truth_discrete(folderPath);
    // test on Benders Decomposition
    //bd_solver(folderPath, engine, 2, 1e-3);
    //compromise_bd_output cd_res = bd_compromise_solver(folderPath, engine, 10, 5, 1, 1e-3);
    // validaton
    //validation_ground_truth_discrete(folderPath, cd_res.compromise_decision);
    // test on cutting plane algorithm
    //cp_output res = cp_solver(folderPath, engine, 10, 1e-3);
    // validate
    //sqqp_validation_ground_truth_discrete(folderPath, res.x);
    // test on SD
    //stoPoint point_est;
    //point_est.e.push_back(5);
    //solverOutput sd_res = sd_sqlp_solver(folderPath, engine, point_est, 200, 0.0, 1.0);
    //validation_ground_truth_discrete(folderPath, sd_res.sol);
    // test on SD for solving SQQP
    //stoPoint point_est;
    //point_est.e.push_back(5);
    //solverOutput sd_res = sd_sqqp_solver(folderPath, engine, point_est, 10, 0.0, 1.0);
    //sqqp_validation_ground_truth_discrete(folderPath, sd_res.sol);
}
