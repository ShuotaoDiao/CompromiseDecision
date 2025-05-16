//
//  ioStochastic.cpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/18/25.
//

#include "ioStochastic.hpp"

// read sto file
stoMap readStochasticMap(const std::string& stochasticPath) {
    stoMap stochasticMapping; // randomness of sto_log_c, and be, bi, and Ce on the right hand side of constraints in the second stage problem
    const std::string nameBeginSto("<sto>");
    const std::string nameEndSto("</sto>");
    
    // second stage
    const std::string nameBeginParameter_e("<e>");
    const std::string nameEndParameter_e("</e>");
    const std::string nameBeginParameter_C("<C>");
    const std::string nameEndParameter_C("</C>");
    // DIST
    const std::string nameBeginParameter_DIST_e("<DIST_e>");
    const std::string nameEndParameter_DIST_e("</DIST_e>");
    const std::string nameBeginParameter_DIST_C("<DIST_C>");
    const std::string nameEndParameter_DIST_C("</DIST_C>");
    
    // e (RHS in the seond stage)
    std::vector<int> discrete_position_e;
    std::vector<double> discrete_value_e;
    std::vector<double> discrete_prob_e;
    std::vector<int> normal_position_e;
    std::vector<double> normal_mean_e;
    std::vector<double> normal_std_e;
    std::vector<double> normal_lb_e;
    std::vector<double> normal_ub_e;
    std::vector<int> uniform_position_e;
    std::vector<double> uniform_lb_e;
    std::vector<double> uniform_ub_e;
    // C (coefficient before x in the second stage)
    std::vector<std::pair<int,int>> discrete_position_C;
    std::vector<double> discrete_value_C;
    std::vector<double> discrete_prob_C;
    std::vector<std::pair<int,int>> normal_position_C;
    std::vector<double> normal_mean_C;
    std::vector<double> normal_std_C;
    std::vector<double> normal_lb_C;
    std::vector<double> normal_ub_C;
    std::vector<std::pair<int,int>> uniform_position_C;
    std::vector<double> uniform_lb_C;
    std::vector<double> uniform_ub_C;
    std::string readCondition("null"); // current condition of reading
    int DIST_line = 0;
    const char* stochasticPathConst = stochasticPath.c_str();
    std::ifstream readFile(stochasticPathConst);
    
    if (readFile.is_open()) {
        std::string line1;
        while (getline(readFile,line1)) {
            std::stringstream ss1(line1);
            if (nameBeginSto.compare(line1) != 0 && nameEndSto.compare(line1) != 0) {
                if (nameBeginParameter_e.compare(line1) == 0) { // start reading e
                    readCondition = "e";
                }
                else if (nameEndParameter_e.compare(line1) == 0) { // end reading e
                    readCondition = "null";
                }
                else if (nameBeginParameter_C.compare(line1) == 0) { // start reading e
                    readCondition = "C";
                }
                else if (nameEndParameter_C.compare(line1) == 0) { // end reading e
                    readCondition = "null";
                }
                else if (nameBeginParameter_DIST_e.compare(line1) == 0) {
                    DIST_line = 0;
                    readCondition = "DIST_e";
                }
                else if (nameEndParameter_DIST_e.compare(line1) == 0) {
                    readCondition = "null";
                }
                else if (nameBeginParameter_DIST_C.compare(line1) == 0) {
                    DIST_line = 0;
                    readCondition = "DIST_C";
                }
                else if (nameEndParameter_DIST_C.compare(line1) == 0) {
                    readCondition = "null";
                }
                else {
                    if (readCondition.compare("e") == 0) { // input syntax: pos;
                        while (getline(ss1, line1, ';')) {
                            int pos = 0;
                            std::stringstream ss2(line1);
                            ss2 >> pos;
                            stochasticMapping.e_map.push_back(pos);
                        }
                    } // end if (readCondition.compare("e"))
                    else if (readCondition.compare("C") == 0) { // input syntax: row,col;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            int row = 0;
                            int col = 0;
                            while (getline(ss2, line1, ',')) {
                                std::stringstream ss3(line1);
                                if (count_item == 0) {
                                    ss3 >> row;
                                }
                                else if (count_item == 1) {
                                    ss3 >> col;
                                }
                                count_item += 1;
                            }
                            std::pair<int,int> loc(row,col);
                            stochasticMapping.C_map.push_back(loc);
                        }
                    } // end else if (readCondition.compare("C"))
                    else if (readCondition.compare("DIST_e") == 0) { // input syntax:
                        // line 0: <<Distribution Name>>
                        // line 1: pos,val,prob; (DISCRETE)
                        // line 1: pos,mean,std; (NORMAL)
                        // line 1: pos,lb,ub; (UNIFORM)
                        // line 1: pos,mean,std,lb,ub; (TRUN_NORMAL)
                        if (DIST_line == 0) {
                            while (getline(ss1, line1, ';')) {
                                stochasticMapping.e_flag_dist = line1;
                            }
                        }
                        else if (DIST_line > 0) {
                            if (stochasticMapping.e_flag_dist.compare("DISCRETE") == 0) {
                                while (getline(ss1, line1, ';')) {
                                    if (line1.compare("\r") != 0) { // revisit later
                                        //std::cout << line1 << std::endl;
                                        //std::cout << "***\n";
                                        std::stringstream ss2(line1);
                                        int count_item = 0;
                                        while (getline(ss2, line1, ',')) {
                                            std::stringstream ss3(line1);
                                            if (count_item == 0) {
                                                int pos = 0;
                                                ss3 >> pos;
                                                discrete_position_e.push_back(pos);
                                            }
                                            else if (count_item == 1) {
                                                double val = 0;
                                                ss3 >> val;
                                                discrete_value_e.push_back(val);
                                            }
                                            else if (count_item == 2) {
                                                double prob = 0;
                                                ss3 >> prob;
                                                discrete_prob_e.push_back(prob);
                                            }
                                            count_item += 1;
                                        }
                                    }
                                }
                            } // end if (stochasticMapping.e_flag_dist.compare("DISCRETE"))
                            else if (stochasticMapping.e_flag_dist.compare("NORMAL") == 0) {
                                while (getline(ss1, line1, ';')) {
                                    std::stringstream ss2(line1);
                                    int count_item = 0;
                                    while (getline(ss2, line1, ',')) {
                                        std::stringstream ss3(line1);
                                        if (count_item == 0) {
                                            int pos = 0;
                                            ss3 >> pos;
                                            normal_position_e.push_back(pos);
                                        }
                                        else if (count_item == 1) {
                                            double mean = 0;
                                            ss3 >> mean;
                                            normal_mean_e.push_back(mean);
                                        }
                                        else if (count_item == 2) {
                                            double std_deviation = 0;
                                            ss3 >> std_deviation;
                                            normal_std_e.push_back(std_deviation);
                                        }
                                        count_item += 1;
                                    }
                                }
                            } // end else if (stochasticMapping.e_flag_dist.compare("NORMAL"))
                            else if (stochasticMapping.e_flag_dist.compare("UNIFORM") == 0) {
                                while (getline(ss1, line1, ';')) {
                                    std::stringstream ss2(line1);
                                    int count_item = 0;
                                    while (getline(ss2, line1, ',')) {
                                        std::stringstream ss3(line1);
                                        if (count_item == 0) {
                                            int pos = 0;
                                            ss3 >> pos;
                                            uniform_position_e.push_back(pos);
                                        }
                                        else if (count_item == 1) {
                                            double lb = 0;
                                            ss3 >> lb;
                                            uniform_lb_e.push_back(lb);
                                        }
                                        else if (count_item == 2) {
                                            double ub = 0;
                                            ss3 >> ub;
                                            uniform_ub_e.push_back(ub);
                                        }
                                        count_item += 1;
                                    }
                                }
                            } // end else if (stochasticMapping.e_flag_dist.compare("UNIFORM"))
                            else if (stochasticMapping.e_flag_dist.compare("TRUN_NORMAL") == 0) {
                                while (getline(ss1, line1, ';')) {
                                    std::stringstream ss2(line1);
                                    int count_item = 0;
                                    while (getline(ss2, line1, ',')) {
                                        std::stringstream ss3(line1);
                                        if (count_item == 0) {
                                            int pos = 0;
                                            ss3 >> pos;
                                            normal_position_e.push_back(pos);
                                        }
                                        else if (count_item == 1) {
                                            double mean = 0;
                                            ss3 >> mean;
                                            normal_mean_e.push_back(mean);
                                        }
                                        else if (count_item == 2) {
                                            double std_deviation = 0;
                                            ss3 >> std_deviation;
                                            normal_std_e.push_back(std_deviation);
                                        }
                                        else if (count_item == 3) {
                                            double lb = 0;
                                            ss3 >> lb;
                                            normal_lb_e.push_back(lb);
                                        }
                                        else if (count_item == 4) {
                                            double ub = 0;
                                            ss3 >> ub;
                                            normal_ub_e.push_back(ub);
                                        }
                                        count_item += 1;
                                    }
                                }
                            } // end else if (stochasticMapping.e_flag_dist.compare("TRUN_NORMAL"))
                        } // end else if (DIST_line > 0)
                        //std::cout << "line: " << DIST_line << std::endl;
                        DIST_line += 1;
                    }// end else if (readCondition.compare("DIST_e"))
                    else if (readCondition.compare("DIST_C") == 0) {
                        // input syntax:
                            // line 0: <<Distribution Name>>
                            // line 1: row,col,val,prob; (DISCRETE)
                            // line 1: row,col,mean,std; (NORMAL)
                            // line 1: row,col,lb,ub; (UNIFORM)
                            // line 1: row,col,mean,std,lb,ub; (TRUN_NORMAL)
                        if (DIST_line == 0) {
                            while (getline(ss1, line1, ';')) {
                                stochasticMapping.C_flag_dist = line1;
                            }
                        }
                        else if (stochasticMapping.C_flag_dist.compare("DISCRETE") == 0) {
                            while (getline(ss1, line1, ';')) {
                                std::stringstream ss2(line1);
                                int count_item = 0;
                                int row = 0;
                                int col = 0;
                                double val = 0;
                                double prob = 0;
                                while (getline(ss2, line1, ',')) {
                                    std::stringstream ss3(line1);
                                    if (count_item == 0) {
                                        ss3 >> row;
                                    }
                                    else if (count_item == 1) {
                                        ss3 >> col;
                                    }
                                    else if (count_item == 2) {
                                        ss3 >> val;
                                    }
                                    else if (count_item == 3) {
                                        ss3 >> prob;
                                    }
                                    count_item += 1;
                                }
                                std::pair<int, int> loc(row,col);
                                discrete_position_C.push_back(loc);
                                discrete_value_C.push_back(val);
                                discrete_prob_C.push_back(prob);
                            }
                        } // end if (stochasticMapping.C_flag_dist.compare("DISCRETE") == 0)
                        else if (stochasticMapping.C_flag_dist.compare("NORMAL") == 0) {
                            while (getline(ss1, line1, ';')) {
                                std::stringstream ss2(line1);
                                int count_item = 0;
                                int row = 0;
                                int col = 0;
                                double mean = 0;
                                double stddev = 0;
                                while (getline(ss2, line1, ',')) {
                                    std::stringstream ss3(line1);
                                    if (count_item == 0) {
                                        ss3 >> row;
                                    }
                                    else if (count_item == 1) {
                                        ss3 >> col;
                                    }
                                    else if (count_item == 2) {
                                        ss3 >> mean;
                                    }
                                    else if (count_item == 3) {
                                        ss3 >> stddev;
                                    }
                                    count_item += 1;
                                }
                                std::pair<int, int> loc(row,col);
                                normal_position_C.push_back(loc);
                                normal_mean_C.push_back(mean);
                                normal_std_C.push_back(stddev);
                            }
                        } // end else if (stochasticMapping.C_flag_dist.compare("DISCRETE") == 0)
                        else if (stochasticMapping.C_flag_dist.compare("UNIFORM") == 0) {
                            while (getline(ss1, line1, ';')) {
                                std::stringstream ss2(line1);
                                int count_item = 0;
                                int row = 0;
                                int col = 0;
                                double lb = 0;
                                double ub = 0;
                                while (getline(ss2, line1, ',')) {
                                    std::stringstream ss3(line1);
                                    if (count_item == 0) {
                                        ss3 >> row;
                                    }
                                    else if (count_item == 1) {
                                        ss3 >> col;
                                    }
                                    else if (count_item == 2) {
                                        ss3 >> lb;
                                    }
                                    else if (count_item == 3) {
                                        ss3 >> ub;
                                    }
                                    count_item += 1;
                                }
                                std::pair<int, int> loc(row,col);
                                uniform_position_C.push_back(loc);
                                uniform_lb_C.push_back(lb);
                                uniform_ub_C.push_back(ub);
                            }
                        } // end else if (stochasticMapping.C_flag_dist.compare("UNIFORM") == 0)
                        else if (stochasticMapping.C_flag_dist.compare("TRUN_NORMAL") == 0) {
                            while (getline(ss1, line1, ';')) {
                                std::stringstream ss2(line1);
                                int count_item = 0;
                                int row = 0;
                                int col = 0;
                                double mean = 0;
                                double stddev = 0;
                                double lb = 0;
                                double ub = 0;
                                while (getline(ss2, line1, ',')) {
                                    std::stringstream ss3(line1);
                                    if (count_item == 0) {
                                        ss3 >> row;
                                    }
                                    else if (count_item == 1) {
                                        ss3 >> col;
                                    }
                                    else if (count_item == 2) {
                                        ss3 >> mean;
                                    }
                                    else if (count_item == 3) {
                                        ss3 >> stddev;
                                    }
                                    else if (count_item == 4) {
                                        ss3 >> lb;
                                    }
                                    else if (count_item == 5) {
                                        ss3 >> ub;
                                    }
                                    count_item += 1;
                                }
                                std::pair<int, int> loc(row,col);
                                normal_position_C.push_back(loc);
                                normal_mean_C.push_back(mean);
                                normal_std_C.push_back(stddev);
                                normal_lb_C.push_back(lb);
                                normal_ub_C.push_back(ub);
                            }
                        } // end else if (stochasticMapping.C_flag_dist.compare("TRUN_NORMAL") == 0)
                        DIST_line += 1;
                    } // end else if (readCondition.compare("DIST_C") == 0)
                } // end reading distribution parameters
            } // end if (nameBeginSto.compare(line1) != 0 && nameEndSto.compare(line1) != 0)
        }
    } // end if (readFile.is_open())
    // e
    if (stochasticMapping.e_flag_dist.compare("DISCRETE") == 0) {
        int cur_pos = -1;
        int cur_idx = -1;
        double cur_cdf = 0;
        for (int idx = 0; idx < discrete_position_e.size(); ++idx) {
            if (cur_pos != discrete_position_e[idx]) { // new position
                cur_cdf = 0; // reset cdf
                cur_idx += 1; // go to the next component
                discreteDistribution temp_dist;
                stochasticMapping.e_discrete.push_back(temp_dist);
            }
            cur_pos = discrete_position_e[idx];
            cur_cdf += discrete_prob_e[idx];
            stochasticMapping.e_discrete[cur_idx].value.push_back(discrete_value_e[idx]);
            stochasticMapping.e_discrete[cur_idx].probability.push_back(discrete_prob_e[idx]);
            stochasticMapping.e_discrete[cur_idx].cdf.push_back(cur_cdf);
        }
    }
    stochasticMapping.e_normal_mean = normal_mean_e;
    stochasticMapping.e_normal_stddev = normal_std_e;
    stochasticMapping.e_normal_lb = normal_lb_e;
    stochasticMapping.e_normal_ub = normal_ub_e;
    stochasticMapping.e_uniform_lower = uniform_lb_e;
    stochasticMapping.e_uniform_upper = uniform_ub_e;
    // C
    if (stochasticMapping.e_flag_dist.compare("DISCRETE") == 0) {
        int cur_row = -1;
        int cur_col = -1;
        int cur_idx = -1;
        double cur_cdf = 0;
        for (int idx = 0; idx < discrete_position_C.size(); ++idx) {
            std::pair<int, int> pos = discrete_position_C[idx];
            if (cur_row != pos.first || cur_col != pos.second) { // new position
                cur_cdf = 0; // reset cdf
                cur_idx += 1; // go to the next component
                discreteDistribution temp_dist;
                stochasticMapping.C_discrete.push_back(temp_dist);
            }
            cur_row = pos.first; // get current position
            cur_col = pos.second;
            cur_cdf += discrete_prob_C[idx];
            stochasticMapping.C_discrete[cur_idx].value.push_back(discrete_value_C[idx]);
            stochasticMapping.C_discrete[cur_idx].probability.push_back(discrete_prob_C[idx]);
            stochasticMapping.C_discrete[cur_idx].cdf.push_back(cur_cdf);
        }
    }
    stochasticMapping.C_normal_mean = normal_mean_C;
    stochasticMapping.C_normal_stddev = normal_std_C;
    stochasticMapping.C_normal_lb = normal_lb_C;
    stochasticMapping.C_normal_ub = normal_ub_C;
    stochasticMapping.C_uniform_lower = uniform_lb_C;
    stochasticMapping.C_uniform_upper = uniform_ub_C;
    return stochasticMapping;
}

// read ground truth
std::vector<stoPoint> readGroundTruth_discreteRV(const std::string& distPath) {
    std::vector<stoPoint> true_dist;
    std::string nameBeginDist("<dist>");
    std::string nameEndDist("</dist>");
    std::string skip("*");
    const char* readPathConst = distPath.c_str(); // convert the string type path to constant
    std::ifstream readFile(readPathConst); // create a readFile object
    if (readFile.is_open()) {
        std::string line1;
        while (getline(readFile,line1)) {
            while (getline(readFile, line1)) { // get the whole line
                std::stringstream ss1(line1); // convert a string into stream
                stoPoint rv_realization;
                unsigned int index_position = 0; // 1 for response, 3 for weight
                if (nameBeginDist.compare(line1) != 0 && nameEndDist.compare(line1) != 0) { // main
                    // syntax:  e:(numbers);C:(numbers);prob:(number);
                    while (getline(ss1, line1, ';')) {
                        std::stringstream ss2(line1);
                        while (getline(ss2, line1, ':')) {
                            if (index_position == 1){ // e read vector
                                std::stringstream ss3(line1);
                                while (getline(ss3, line1, ',')) {
                                    double value;
                                    std::stringstream ss4(line1);
                                    ss4 >> value;
                                    rv_realization.e.push_back(value);
                                }
                            }
                            else if (index_position == 3) { // C
                                std::stringstream ss3(line1);
                                if (skip.compare(line1) != 0) { // no *, so do not skip
                                    while (getline(ss3, line1, ',')) {
                                        double value;
                                        std::stringstream ss4(line1);
                                        ss4 >> value;
                                        rv_realization.C.push_back(value);
                                    }
                                }
                            }
                            else if (index_position == 5) { // prob
                                double value;
                                std::stringstream ss3(line1);
                                ss3 >> value;
                                rv_realization.prob = value;
                            }
                            index_position++;
                        }
                    } // end while (getline(ss1, line1, ';'))
                    true_dist.push_back(rv_realization);
                }
            }
        } // end while (getline(readFile,line1))
    } // end if (readFile.is_open())
    readFile.close();
    return true_dist;
}

// discrete random number generator, random_number must be in [0,1)
double discrete_random_number_generator(double random_number, const discreteDistribution dist) {
    double res = 0;
    for (int idx = 0; idx < dist.cdf.size(); ++idx) {
        if (random_number < dist.cdf[idx]) {
            res = dist.value[idx];
            break;
        }
    }
    //std::cout << "DEBUG: random_number: " << random_number << std::endl;
    return res;
}

// uniform random number generator
double uniform_random_number_generator(std::mt19937& generator, double lb, double ub) {
    std::uniform_real_distribution<double> uniform(lb,ub);
    double val = uniform(generator);
    return val;
}

// normal random number generator
double normal_random_number_generator(std::mt19937& generator, double mean, double stddev) {
    std::normal_distribution<double> normal(mean, stddev);
    double val = normal(generator);
    return val;
}

// truncated normal random number generator
double truncated_normal_random_number_generator(std::mt19937& generator, double mean, double stddev, double lb, double ub) {
    std::normal_distribution<double> normal(mean, stddev);
    double val = normal(generator);
    while (val < lb || val > ub) {
        val = normal(generator);
    }
    return val;
}


