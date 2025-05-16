//
//  SAA_solver.cpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/18/25.
//

#include "SAA_solver.hpp"

solverOutput saa_2slp(const std::string& folder_path,
                      std::mt19937& generator,
                      int sample_size,
                      int idx_rep) {
    // initialize uniform (0,1)
    std::uniform_real_distribution<double> dist_uniform(0.0,1.0);
    std::vector<stoPoint> Xi; // set of historic data
    // data file
    bool flag_data_e; // tell if e stochastic is generated
    bool flag_data_C; // tell if C stochastic is generated
    // create directory paths for database and model
    std::string e_DB_path = folder_path + "/e_DB.txt";
    std::string C_DB_path = folder_path + "/C_DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(rep_SAAv1.0).txt";
    // convert all the paths into constant chars
    const char* e_DB_path_const = e_DB_path.c_str();
    const char* C_DB_path_const = C_DB_path.c_str();
    // create stream object
    std::ifstream readFile_e(e_DB_path_const);
    std::ifstream readFile_C(C_DB_path_const);
    // create database
    std::vector<std::vector<dataPoint>> e_DB;
    std::vector<std::vector<dataPoint>> C_DB;
    // create model structure
    standardTwoStageParameters model_parameters;
    // create sto object
    stoMap stochastic_map;
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
    // read model file
    model_parameters = readStandardTwoStageParameters(model_path);
    // read sto file
    stochastic_map = readStochasticMap(sto_path);
    // --- Finish Reading ---
    // --- gurobi solver environment ---
    // Create the Gurobi environment
    GRBEnv env = GRBEnv();
    
    // Create an empty model
    GRBModel model = GRBModel(env);
    
    // create decision variables
    // first-stage decision variable
    std::vector<GRBVar> x;
    x.reserve(model_parameters.num_var_1stStage); // Pre-allocate memory
    
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        x.push_back(model.addVar(model_parameters.var_lb_1stStage, model_parameters.var_ub_1stStage, 0.0, GRB_CONTINUOUS));
    }
    // second-stage decision variable
    std::vector<std::vector<GRBVar>> y;
    // slack variable
    std::vector<std::vector<GRBVar>> y_slack;
    long total_num_var_2ndStage = model_parameters.num_var_2ndStage + model_parameters.num_slack_2ndStage;
    for (int idx_scenario = 0; idx_scenario < sample_size; ++idx_scenario) {
        std::vector<GRBVar> y_scenario;
        y.reserve(model_parameters.num_var_2ndStage);
        for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
            y_scenario.push_back(model.addVar(model_parameters.var_lb_2ndStage,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
        }
        y.push_back(y_scenario);
        if (model_parameters.num_slack_2ndStage > 0) {
            std::vector<GRBVar> y_slack_scenario;
            for (int idx_y_slack = 0; idx_y_slack < model_parameters.num_slack_2ndStage; ++idx_y_slack) {
                y_slack_scenario.push_back(model.addVar(0,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
            }
            y_slack.push_back(y_slack_scenario);
        }
    }
    
    // objective
    GRBLinExpr expr_obj;
    // first stage objective
    for (int idx_x = 0; idx_x < model_parameters.c.getNzeroLen(); ++idx_x) {
        expr_obj += model_parameters.c.getVal(idx_x) * x[idx_x];
    }
    // second stage objective
    double one_over_N = 1.0 / ((double) sample_size);
    for (int idx_scenario = 0; idx_scenario < sample_size; ++idx_scenario) {
        for (int idx_y = 0; idx_y < model_parameters.d.getNzeroLen(); ++idx_y) {
            expr_obj += one_over_N * model_parameters.d.getVal(idx_y) * y[idx_scenario][model_parameters.d.getLoc(idx_y)];
        }
    }
    
    // add objective
    model.setObjective(expr_obj, GRB_MINIMIZE);
    
    
    // constrants
    // first stage constraints
    // <= constrants
    std::vector<GRBLinExpr> exprs_L;
    for (int index_cons = 0; index_cons < model_parameters.AL.getRowLength(); ++index_cons) {
        GRBLinExpr expr_L;
        exprs_L.push_back(expr_L);
    }
    for (int col_idx = 0; col_idx < model_parameters.AL.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AL.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AL.getClen(col_idx) + beg_idx; ++idx) {
            exprs_L[model_parameters.AL.getRow(idx)] += model_parameters.AL.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bL.getNzeroLen(); ++idx) {
        exprs_L[model_parameters.bL.getLoc(idx)] -= model_parameters.bL.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AL.getRowLength(); ++idx) {
        model.addConstr(exprs_L[idx] <= 0);
    }
    
    // == constraints
    std::vector<GRBLinExpr> exprs_E;
    for (int index_cons = 0; index_cons < model_parameters.AE.getRowLength(); ++index_cons) {
        GRBLinExpr expr_E;
        exprs_E.push_back(expr_E);
    }
    for (int col_idx = 0; col_idx < model_parameters.AE.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AE.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AE.getClen(col_idx) + beg_idx; ++idx) {
            exprs_E[model_parameters.AE.getRow(idx)] += model_parameters.AE.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bE.getNzeroLen(); ++idx) {
        exprs_E[model_parameters.bE.getLoc(idx)] -= model_parameters.bE.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AE.getRowLength(); ++idx) {
        model.addConstr(exprs_E[idx] == 0);
    }
    
    // >= constraints
    std::vector<GRBLinExpr> exprs_G;
    for (int index_cons = 0; index_cons < model_parameters.AG.getRowLength(); ++index_cons) {
        GRBLinExpr expr_G;
        exprs_G.push_back(expr_G);
    }
    for (int col_idx = 0; col_idx < model_parameters.AG.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AG.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AG.getClen(col_idx) + beg_idx; ++idx) {
            exprs_G[model_parameters.AG.getRow(idx)] += model_parameters.AG.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bG.getNzeroLen(); ++idx) {
        exprs_G[model_parameters.bG.getLoc(idx)] -= model_parameters.bG.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AG.getRowLength(); ++idx) {
        model.addConstr(exprs_G[idx] >= 0);
    }
    
    // second stage
    // equality constraints D y + C x + I_slack y_slack = e
    for (int idx_scenario = 0; idx_scenario < sample_size; ++idx_scenario) {
        std::vector<GRBLinExpr> exprs_eq;
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            GRBLinExpr expr;
            exprs_eq.push_back(expr);
        }
        // coefficients before y; Dy
        for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.D.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[idx_scenario][col_idx];
            }
        }
        // slack variables; I_slack y_slack
        if (model_parameters.num_slack_2ndStage > 0){
            for (int col_idx = 0; col_idx < model_parameters.ISlack.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.ISlack.getCbeg(col_idx);
                for (int idx = beg_idx; idx < model_parameters.ISlack.getClen(col_idx) + beg_idx; ++idx) {
                    exprs_eq[model_parameters.ISlack.getRow(idx)] += model_parameters.ISlack.getVal(idx) * y_slack[idx_scenario][col_idx];
                }
            }
        }
        // coefficients before x (deterministic part)
        for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.C.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.C.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.C.getRow(idx)] += model_parameters.C.getVal(idx) * x[col_idx];
            }
        }
        // obtain a new data point
        stoPoint xi_cur = generator_random(model_parameters, stochastic_map, generator, idx_scenario, e_DB, C_DB, idx_rep);
        // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
        for (int idx_C = 0; idx_C < xi_cur.C.size(); ++idx_C) {
            exprs_eq[stochastic_map.C_map[idx_C].first] += xi_cur.C[idx_C] * x[stochastic_map.C_map[idx_C].second];
        }
        // right hand side e
        for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
            exprs_eq[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
        }
        // right hand side (stochastic part) equality e_(i) equality
        for (int idx = 0; idx < xi_cur.e.size(); ++idx) {
            exprs_eq[stochastic_map.e_map[idx]] -= xi_cur.e[idx];
        }
        // add constraints
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            model.addConstr(exprs_eq[index_eq] == 0);
        }
        Xi.push_back(xi_cur);
    }
    // optimize
    model.optimize();
    // create outputs
    solverOutput res;
    // objective
    std::cout << "Optimal Value: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    res.obj = model.get(GRB_DoubleAttr_ObjVal);
    // solution
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        res.sol.push_back(x[idx_x].get(GRB_DoubleAttr_X));
    }
    res.Xi = Xi;
    //model.write("/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleLands/lands_model_saa.lp");
    return res;
}

solverOutput saa_sqqp(const std::string& folder_path, std::mt19937& generator, int sample_size) {
    // initialize uniform (0,1)
    std::uniform_real_distribution<double> dist_uniform(0.0,1.0);
    std::vector<stoPoint> Xi; // set of historic data
    // data file
    bool flag_data_e; // tell if e stochastic is generated
    bool flag_data_C; // tell if C stochastic is generated
    // create directory paths for database and model
    std::string e_DB_path = folder_path + "/e_DB.txt";
    std::string C_DB_path = folder_path + "/C_DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(rep_sqqp_SAAv1.0).txt";
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath, std::fstream::app); // append results to the end of the file
    // convert all the paths into constant chars
    const char* e_DB_path_const = e_DB_path.c_str();
    const char* C_DB_path_const = C_DB_path.c_str();
    // create stream object
    std::ifstream readFile_e(e_DB_path_const);
    std::ifstream readFile_C(C_DB_path_const);
    // create database
    std::vector<std::vector<dataPoint>> e_DB;
    std::vector<std::vector<dataPoint>> C_DB;
    // create model structure
    standardTwoStageParameters model_parameters;
    // create sto object
    stoMap stochastic_map;
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
    // read model file
    model_parameters = readStandardTwoStageParameters(model_path);
    // read sto file
    stochastic_map = readStochasticMap(sto_path);
    // --- Finish Reading ---
    // --- gurobi solver environment ---
    // Create the Gurobi environment
    GRBEnv env = GRBEnv();
    
    // Create an empty model
    GRBModel model = GRBModel(env);
    
    // create decision variables
    // first-stage decision variable
    std::vector<GRBVar> x;
    x.reserve(model_parameters.num_var_1stStage); // Pre-allocate memory
    
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        x.push_back(model.addVar(model_parameters.var_lb_1stStage, model_parameters.var_ub_1stStage, 0.0, GRB_CONTINUOUS));
    }
    // second-stage decision variable
    std::vector<std::vector<GRBVar>> y;
    // slack variable
    std::vector<std::vector<GRBVar>> y_slack;
    long total_num_var_2ndStage = model_parameters.num_var_2ndStage + model_parameters.num_slack_2ndStage;
    for (int idx_scenario = 0; idx_scenario < sample_size; ++idx_scenario) {
        std::vector<GRBVar> y_scenario;
        y.reserve(model_parameters.num_var_2ndStage);
        for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
            y_scenario.push_back(model.addVar(model_parameters.var_lb_2ndStage,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
        }
        y.push_back(y_scenario);
        if (model_parameters.num_slack_2ndStage > 0) {
            std::vector<GRBVar> y_slack_scenario;
            for (int idx_y_slack = 0; idx_y_slack < model_parameters.num_slack_2ndStage; ++idx_y_slack) {
                y_slack_scenario.push_back(model.addVar(0,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
            }
            y_slack.push_back(y_slack_scenario);
        }
    }
    
    // objective
    GRBQuadExpr expr_obj;
    // first stage objective
    for (int idx_x = 0; idx_x < model_parameters.c.getNzeroLen(); ++idx_x) {
        expr_obj += model_parameters.c.getVal(idx_x) * x[idx_x];
    }
    // quadratic terms
    for (int col_idx = 0; col_idx < model_parameters.Q.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.Q.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.Q.getClen(col_idx) + beg_idx; ++idx) {
            expr_obj += 0.5 * model_parameters.Q.getVal(idx) * x[model_parameters.Q.getRow(idx)] * x[col_idx];
        }
    }
    // second stage objective
    double one_over_N = 1.0 / ((double) sample_size);
    for (int idx_scenario = 0; idx_scenario < sample_size; ++idx_scenario) {
        for (int idx_y = 0; idx_y < model_parameters.d.getNzeroLen(); ++idx_y) {
            expr_obj += one_over_N * model_parameters.d.getVal(idx_y) * y[idx_scenario][model_parameters.d.getLoc(idx_y)];
        }
        // quadratic terms
        for (int col_idx = 0; col_idx < model_parameters.P.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.P.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.P.getClen(col_idx) + beg_idx; ++idx) {
                if (model_parameters.P.getRow(idx) < model_parameters.num_var_2ndStage) {
                    if (col_idx < model_parameters.num_var_2ndStage) {
                        // case 1: P_ij * y_i * y_j
                        expr_obj += 0.5 * one_over_N * model_parameters.P.getVal(idx) * y[idx_scenario][model_parameters.P.getRow(idx)] * y[idx_scenario][col_idx];
                    }
                    else {
                        // case 2: P_ij * y_i * yslack_j
                        expr_obj += 0.5 * one_over_N * model_parameters.P.getVal(idx) * y[idx_scenario][model_parameters.P.getRow(idx)] * y_slack[idx_scenario][col_idx - model_parameters.num_var_2ndStage];
                    }
                }
                else {
                    if (col_idx < model_parameters.num_var_2ndStage) {
                        // case 3: P_ij * yslack_i * y_j
                        expr_obj += 0.5 * one_over_N * model_parameters.P.getVal(idx) * y_slack[idx_scenario][model_parameters.P.getRow(idx) - model_parameters.num_var_2ndStage] * y[idx_scenario][col_idx];
                    }
                    else {
                        // case 4: P_ij * yslack_i * yslack_j
                        expr_obj += 0.5 * one_over_N * model_parameters.P.getVal(idx) * y_slack[idx_scenario][model_parameters.P.getRow(idx) - model_parameters.num_var_2ndStage] * y_slack[idx_scenario][col_idx - model_parameters.num_var_2ndStage];
                    }
                }
            }
        }
    }
    
    // add objective
    model.setObjective(expr_obj, GRB_MINIMIZE);
    
    // constrants
    // first stage constraints
    // <= constrants
    std::vector<GRBLinExpr> exprs_L;
    for (int index_cons = 0; index_cons < model_parameters.AL.getRowLength(); ++index_cons) {
        GRBLinExpr expr_L;
        exprs_L.push_back(expr_L);
    }
    for (int col_idx = 0; col_idx < model_parameters.AL.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AL.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AL.getClen(col_idx) + beg_idx; ++idx) {
            exprs_L[model_parameters.AL.getRow(idx)] += model_parameters.AL.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bL.getNzeroLen(); ++idx) {
        exprs_L[model_parameters.bL.getLoc(idx)] -= model_parameters.bL.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AL.getRowLength(); ++idx) {
        model.addConstr(exprs_L[idx] <= 0);
    }
    
    // == constraints
    std::vector<GRBLinExpr> exprs_E;
    for (int index_cons = 0; index_cons < model_parameters.AE.getRowLength(); ++index_cons) {
        GRBLinExpr expr_E;
        exprs_E.push_back(expr_E);
    }
    for (int col_idx = 0; col_idx < model_parameters.AE.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AE.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AE.getClen(col_idx) + beg_idx; ++idx) {
            exprs_E[model_parameters.AE.getRow(idx)] += model_parameters.AE.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bE.getNzeroLen(); ++idx) {
        exprs_E[model_parameters.bE.getLoc(idx)] -= model_parameters.bE.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AE.getRowLength(); ++idx) {
        model.addConstr(exprs_E[idx] == 0);
    }
    
    // >= constraints
    std::vector<GRBLinExpr> exprs_G;
    for (int index_cons = 0; index_cons < model_parameters.AG.getRowLength(); ++index_cons) {
        GRBLinExpr expr_G;
        exprs_G.push_back(expr_G);
    }
    for (int col_idx = 0; col_idx < model_parameters.AG.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AG.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AG.getClen(col_idx) + beg_idx; ++idx) {
            exprs_G[model_parameters.AG.getRow(idx)] += model_parameters.AG.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bG.getNzeroLen(); ++idx) {
        exprs_G[model_parameters.bG.getLoc(idx)] -= model_parameters.bG.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AG.getRowLength(); ++idx) {
        model.addConstr(exprs_G[idx] >= 0);
    }
    
    // second stage
    // equality constraints D y + C x + I_slack y_slack = e
    for (int idx_scenario = 0; idx_scenario < sample_size; ++idx_scenario) {
        std::vector<GRBLinExpr> exprs_eq;
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            GRBLinExpr expr;
            exprs_eq.push_back(expr);
        }
        // coefficients before y; Dy
        for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.D.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[idx_scenario][col_idx];
            }
        }
        // slack variables; I_slack y_slack
        if (model_parameters.num_slack_2ndStage > 0){
            for (int col_idx = 0; col_idx < model_parameters.ISlack.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.ISlack.getCbeg(col_idx);
                for (int idx = beg_idx; idx < model_parameters.ISlack.getClen(col_idx) + beg_idx; ++idx) {
                    exprs_eq[model_parameters.ISlack.getRow(idx)] += model_parameters.ISlack.getVal(idx) * y_slack[idx_scenario][col_idx];
                }
            }
        }
        // coefficients before x (deterministic part)
        for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.C.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.C.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.C.getRow(idx)] += model_parameters.C.getVal(idx) * x[col_idx];
            }
        }
        // obtain a new data point
        stoPoint xi_cur = generator_random(model_parameters, stochastic_map, generator, idx_scenario, e_DB, C_DB);
        // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
        for (int idx_C = 0; idx_C < xi_cur.C.size(); ++idx_C) {
            exprs_eq[stochastic_map.C_map[idx_C].first] += xi_cur.C[idx_C] * x[stochastic_map.C_map[idx_C].second];
        }
        // right hand side e
        for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
            exprs_eq[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
        }
        // right hand side (stochastic part) equality e_(i) equality
        for (int idx = 0; idx < xi_cur.e.size(); ++idx) {
            exprs_eq[stochastic_map.e_map[idx]] -= xi_cur.e[idx];
        }
        // add constraints
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            model.addConstr(exprs_eq[index_eq] == 0);
        }
        Xi.push_back(xi_cur);
    }
    // optimize
    model.optimize();
    // create outputs
    solverOutput res;
    // objective
    std::cout << "--- SAA of SQQP---\n";
    writeFile << "--- SAA of SQQP---\n";
    std::cout << "Number of replications: " << sample_size << std::endl;
    writeFile << "Number of replications: " << sample_size << std::endl;
    std::cout << "Optimal Value: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    writeFile << "Optimal Value: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    res.obj = model.get(GRB_DoubleAttr_ObjVal);
    // solution
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        res.sol.push_back(x[idx_x].get(GRB_DoubleAttr_X));
        writeFile << res.sol[idx_x] << ", ";
    }
    writeFile << std::endl;
    res.Xi = Xi;
    //model.write("/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleQP/saa_sqqp_model.lp");
    
    writeFile.close();
    return res;
    
}

// master problem

compromiseOutput saa_2slp_compromise_master(standardTwoStageParameters& model_parameters,
                                            const stoMap& stochastic_map,
                                            std::vector<solverOutput>& aggregate_outputs,
                                            int sample_size_per_replication,
                                            int num_replications,
                                            double regularizer_coefficient) {
    // --- gurobi solver environment ---
    // Create the Gurobi environment
    GRBEnv env = GRBEnv();
    
    // Create an empty model
    GRBModel model = GRBModel(env);
    
    // create decision variables
    // first-stage decision variable
    std::vector<GRBVar> x;
    x.reserve(model_parameters.num_var_1stStage); // Pre-allocate memory
    
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        x.push_back(model.addVar(model_parameters.var_lb_1stStage, model_parameters.var_ub_1stStage, 0.0, GRB_CONTINUOUS));
    }
    
    // second-stage decision variable
    std::vector<std::vector<std::vector<GRBVar>>> y;
    // slack variable
    std::vector<std::vector<std::vector<GRBVar>>> y_slack;
    long total_num_var_2ndStage = model_parameters.num_var_2ndStage + model_parameters.num_slack_2ndStage;
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        std::vector<std::vector<GRBVar>> y_rep;
        std::vector<std::vector<GRBVar>> y_slack_rep;
        for (int idx_scenario = 0; idx_scenario < sample_size_per_replication; ++idx_scenario) {
            std::vector<GRBVar> y_scenario;
            y.reserve(model_parameters.num_var_2ndStage);
            for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
                y_scenario.push_back(model.addVar(model_parameters.var_lb_2ndStage,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
            }
            y_rep.push_back(y_scenario);
            if (model_parameters.num_slack_2ndStage > 0) {
                std::vector<GRBVar> y_slack_scenario;
                for (int idx_y_slack = 0; idx_y_slack < model_parameters.num_slack_2ndStage; ++idx_y_slack) {
                    y_slack_scenario.push_back(model.addVar(0,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
                }
                y_slack_rep.push_back(y_slack_scenario);
            }
        } // end for (int idx_scenario = 0; idx_scenario < sample_size_per_replication; ++idx_scenario)
        y.push_back(y_rep);
        y_slack.push_back(y_slack_rep);
    } // for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep)
    
    // objective
    GRBQuadExpr expr_obj;
    
    double one_over_num_replications = 1.0 / ((double) num_replications);
    double one_over_N = 1.0 / ((double) (sample_size_per_replication * num_replications));
    // first stage objective
    for (int idx_x = 0; idx_x < model_parameters.c.getNzeroLen(); ++idx_x) {
        expr_obj += model_parameters.c.getVal(idx_x) * x[idx_x];
        // regularizer
        expr_obj += 0.5 * regularizer_coefficient * x[idx_x] * x[idx_x];
        double x_average = 0;
        // regularizer
        for (int idx_rep = 0; idx_rep < aggregate_outputs.size(); ++idx_rep) {
            x_average += aggregate_outputs[idx_rep].sol[idx_x];
        }
        x_average *= one_over_num_replications;
        expr_obj -= regularizer_coefficient * x_average * x[idx_x];
    }
    
    // second stage objective
    
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        for (int idx_scenario = 0; idx_scenario < sample_size_per_replication; ++idx_scenario) {
            for (int idx_y = 0; idx_y < model_parameters.d.getNzeroLen(); ++idx_y) {
                expr_obj += one_over_N * model_parameters.d.getVal(idx_y) * y[idx_rep][idx_scenario][model_parameters.d.getLoc(idx_y)];
            }
        }
    }
    
    // add objective
    model.setObjective(expr_obj, GRB_MINIMIZE);
    
    // constrants
    // first stage constraints
    // <= constrants
    std::vector<GRBLinExpr> exprs_L;
    for (int index_cons = 0; index_cons < model_parameters.AL.getRowLength(); ++index_cons) {
        GRBLinExpr expr_L;
        exprs_L.push_back(expr_L);
    }
    for (int col_idx = 0; col_idx < model_parameters.AL.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AL.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AL.getClen(col_idx) + beg_idx; ++idx) {
            exprs_L[model_parameters.AL.getRow(idx)] += model_parameters.AL.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bL.getNzeroLen(); ++idx) {
        exprs_L[model_parameters.bL.getLoc(idx)] -= model_parameters.bL.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AL.getRowLength(); ++idx) {
        model.addConstr(exprs_L[idx] <= 0);
    }
    
    // == constraints
    std::vector<GRBLinExpr> exprs_E;
    for (int index_cons = 0; index_cons < model_parameters.AE.getRowLength(); ++index_cons) {
        GRBLinExpr expr_E;
        exprs_E.push_back(expr_E);
    }
    for (int col_idx = 0; col_idx < model_parameters.AE.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AE.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AE.getClen(col_idx) + beg_idx; ++idx) {
            exprs_E[model_parameters.AE.getRow(idx)] += model_parameters.AE.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bE.getNzeroLen(); ++idx) {
        exprs_E[model_parameters.bE.getLoc(idx)] -= model_parameters.bE.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AE.getRowLength(); ++idx) {
        model.addConstr(exprs_E[idx] == 0);
    }
    
    // >= constraints
    std::vector<GRBLinExpr> exprs_G;
    for (int index_cons = 0; index_cons < model_parameters.AG.getRowLength(); ++index_cons) {
        GRBLinExpr expr_G;
        exprs_G.push_back(expr_G);
    }
    for (int col_idx = 0; col_idx < model_parameters.AG.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AG.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AG.getClen(col_idx) + beg_idx; ++idx) {
            exprs_G[model_parameters.AG.getRow(idx)] += model_parameters.AG.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bG.getNzeroLen(); ++idx) {
        exprs_G[model_parameters.bG.getLoc(idx)] -= model_parameters.bG.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AG.getRowLength(); ++idx) {
        model.addConstr(exprs_G[idx] >= 0);
    }
    
    // second-stage constraints
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        // equality constraints D y + C x + I_slack y_slack = e
        for (int idx_scenario = 0; idx_scenario < sample_size_per_replication; ++idx_scenario) {
            std::vector<GRBLinExpr> exprs_eq_rep;
            for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
                GRBLinExpr expr;
                exprs_eq_rep.push_back(expr);
            }
            // coefficients before y; Dy
            for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.D.getCbeg(col_idx);
                for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
                    exprs_eq_rep[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[idx_rep][idx_scenario][col_idx];
                }
            }
            // slack variables; I_slack y_slack
            if (model_parameters.num_slack_2ndStage > 0){
                for (int col_idx = 0; col_idx < model_parameters.ISlack.getColLength(); ++col_idx) {
                    int beg_idx = model_parameters.ISlack.getCbeg(col_idx);
                    for (int idx = beg_idx; idx < model_parameters.ISlack.getClen(col_idx) + beg_idx; ++idx) {
                        exprs_eq_rep[model_parameters.ISlack.getRow(idx)] += model_parameters.ISlack.getVal(idx) * y_slack[idx_rep][idx_scenario][col_idx];
                    }
                }
            }
            // coefficients before x (deterministic part)
            for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.C.getCbeg(col_idx);
                for (int idx = beg_idx; idx < model_parameters.C.getClen(col_idx) + beg_idx; ++idx) {
                    exprs_eq_rep[model_parameters.C.getRow(idx)] += model_parameters.C.getVal(idx) * x[col_idx];
                }
            }
            // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
            for (int idx_C = 0; idx_C < aggregate_outputs[idx_rep].Xi[idx_scenario].C.size(); ++idx_C) {
                exprs_eq_rep[stochastic_map.C_map[idx_C].first] += aggregate_outputs[idx_rep].Xi[idx_scenario].C[idx_C] * x[stochastic_map.C_map[idx_C].second];
            }
            // right hand side e
            for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
                exprs_eq_rep[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
            }
            // right hand side (stochastic part) equality e_(i) equality
            for (int idx = 0; idx < aggregate_outputs[idx_rep].Xi[idx_scenario].e.size(); ++idx) {
                exprs_eq_rep[stochastic_map.e_map[idx]] -= aggregate_outputs[idx_rep].Xi[idx_scenario].e[idx];
            }
            // add constraints
            for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
                model.addConstr(exprs_eq_rep[index_eq] == 0);
            }
        }
    } // end for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep)
    
    // optimize
    model.optimize();
    // create outputs
    compromiseOutput compromise_res;
    // solution
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        compromise_res.compromise_decision.push_back(x[idx_x].get(GRB_DoubleAttr_X));
    }
    compromise_res.replications = aggregate_outputs;
    //model.write("/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/lands/lands_compromise_master.lp");
    return compromise_res;
}

// SAA problem after aggragation, compromise decision
// master problem for optimization
compromiseOutput saa_sqqp_compromise_master(standardTwoStageParameters& model_parameters,
                                            const stoMap& stochastic_map,
                                            std::vector<solverOutput>& aggregate_outputs,
                                            int sample_size_per_replication,
                                            int num_replications,
                                            double regularizer_coefficient) {
    // --- gurobi solver environment ---
    // Create the Gurobi environment
    GRBEnv env = GRBEnv();
    
    // Create an empty model
    GRBModel model = GRBModel(env);
    
    // create decision variables
    // first-stage decision variable
    std::vector<GRBVar> x;
    x.reserve(model_parameters.num_var_1stStage); // Pre-allocate memory
    
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        x.push_back(model.addVar(model_parameters.var_lb_1stStage, model_parameters.var_ub_1stStage, 0.0, GRB_CONTINUOUS));
    }
    
    // second-stage decision variable
    std::vector<std::vector<std::vector<GRBVar>>> y;
    // slack variable
    std::vector<std::vector<std::vector<GRBVar>>> y_slack;
    long total_num_var_2ndStage = model_parameters.num_var_2ndStage + model_parameters.num_slack_2ndStage;
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        std::vector<std::vector<GRBVar>> y_rep;
        std::vector<std::vector<GRBVar>> y_slack_rep;
        for (int idx_scenario = 0; idx_scenario < sample_size_per_replication; ++idx_scenario) {
            std::vector<GRBVar> y_scenario;
            y.reserve(model_parameters.num_var_2ndStage);
            for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
                y_scenario.push_back(model.addVar(model_parameters.var_lb_2ndStage,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
            }
            y_rep.push_back(y_scenario);
            if (model_parameters.num_slack_2ndStage > 0) {
                std::vector<GRBVar> y_slack_scenario;
                for (int idx_y_slack = 0; idx_y_slack < model_parameters.num_slack_2ndStage; ++idx_y_slack) {
                    y_slack_scenario.push_back(model.addVar(0,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
                }
                y_slack_rep.push_back(y_slack_scenario);
            }
        } // end for (int idx_scenario = 0; idx_scenario < sample_size_per_replication; ++idx_scenario)
        y.push_back(y_rep);
        y_slack.push_back(y_slack_rep);
    } // for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep)
    
    // objective
    GRBQuadExpr expr_obj;
    
    double one_over_num_replications = 1.0 / ((double) num_replications);
    double one_over_N = 1.0 / ((double) (sample_size_per_replication * num_replications));
    // first-stage objective
    for (int idx_x = 0; idx_x < model_parameters.c.getNzeroLen(); ++idx_x) {
        expr_obj += model_parameters.c.getVal(idx_x) * x[idx_x];
        // regularizer
        expr_obj += 0.5 * regularizer_coefficient * x[idx_x] * x[idx_x];
        double x_average = 0;
        // regularizer
        for (int idx_rep = 0; idx_rep < aggregate_outputs.size(); ++idx_rep) {
            x_average += aggregate_outputs[idx_rep].sol[idx_x];
        }
        x_average *= one_over_num_replications;
        expr_obj -= regularizer_coefficient * x_average * x[idx_x];
    }
    // first-stage quadratic terms 0.5 x^\top Q x
    for (int col_idx = 0; col_idx < model_parameters.Q.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.Q.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.Q.getClen(col_idx) + beg_idx; ++idx) {
            expr_obj += 0.5 * model_parameters.Q.getVal(idx) * x[model_parameters.Q.getRow(idx)] * x[col_idx];
        }
    }
    
    // second stage objective
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        for (int idx_scenario = 0; idx_scenario < sample_size_per_replication; ++idx_scenario) {
            // d^\top y
            for (int idx_y = 0; idx_y < model_parameters.d.getNzeroLen(); ++idx_y) {
                expr_obj += one_over_N * model_parameters.d.getVal(idx_y) * y[idx_rep][idx_scenario][model_parameters.d.getLoc(idx_y)];
            }
            // quadratic term 0.5 y^\top P y
            for (int col_idx = 0; col_idx < model_parameters.P.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.P.getCbeg(col_idx);
                for (int idx = beg_idx; idx < model_parameters.P.getClen(col_idx) + beg_idx; ++idx) {
                    if (model_parameters.P.getRow(idx) < model_parameters.num_var_2ndStage) {
                        if (col_idx < model_parameters.num_var_2ndStage) {
                            // case 1: P_ij * y_i * y_j
                            expr_obj += 0.5 * one_over_N * model_parameters.P.getVal(idx) * y[idx_rep][idx_scenario][model_parameters.P.getRow(idx)] * y[idx_rep][idx_scenario][col_idx];
                        }
                        else {
                            // case 2: P_ij * y_i * yslack_j
                            expr_obj += 0.5 * one_over_N * model_parameters.P.getVal(idx) * y[idx_rep][idx_scenario][model_parameters.P.getRow(idx)] * y_slack[idx_rep][idx_scenario][col_idx - model_parameters.num_var_2ndStage];
                        }
                    }
                    else {
                        if (col_idx < model_parameters.num_var_2ndStage) {
                            // case 3: P_ij * yslack_i * y_j
                            expr_obj += 0.5 * one_over_N * model_parameters.P.getVal(idx) * y_slack[idx_rep][idx_scenario][model_parameters.P.getRow(idx) - model_parameters.num_var_2ndStage] * y[idx_rep][idx_scenario][col_idx];
                        }
                        else {
                            // case 4: P_ij * yslack_i * yslack_j
                            expr_obj += 0.5 * one_over_N * model_parameters.P.getVal(idx) * y_slack[idx_rep][idx_scenario][model_parameters.P.getRow(idx) - model_parameters.num_var_2ndStage] * y_slack[idx_rep][idx_scenario][col_idx - model_parameters.num_var_2ndStage];
                        }
                    }
                }
            }// end setting up quadratic terms
        }
    }
    // add objective
    model.setObjective(expr_obj, GRB_MINIMIZE);
    
    // constrants
    // first stage constraints
    // <= constrants
    std::vector<GRBLinExpr> exprs_L;
    for (int index_cons = 0; index_cons < model_parameters.AL.getRowLength(); ++index_cons) {
        GRBLinExpr expr_L;
        exprs_L.push_back(expr_L);
    }
    for (int col_idx = 0; col_idx < model_parameters.AL.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AL.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AL.getClen(col_idx) + beg_idx; ++idx) {
            exprs_L[model_parameters.AL.getRow(idx)] += model_parameters.AL.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bL.getNzeroLen(); ++idx) {
        exprs_L[model_parameters.bL.getLoc(idx)] -= model_parameters.bL.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AL.getRowLength(); ++idx) {
        model.addConstr(exprs_L[idx] <= 0);
    }
    
    // == constraints
    std::vector<GRBLinExpr> exprs_E;
    for (int index_cons = 0; index_cons < model_parameters.AE.getRowLength(); ++index_cons) {
        GRBLinExpr expr_E;
        exprs_E.push_back(expr_E);
    }
    for (int col_idx = 0; col_idx < model_parameters.AE.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AE.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AE.getClen(col_idx) + beg_idx; ++idx) {
            exprs_E[model_parameters.AE.getRow(idx)] += model_parameters.AE.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bE.getNzeroLen(); ++idx) {
        exprs_E[model_parameters.bE.getLoc(idx)] -= model_parameters.bE.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AE.getRowLength(); ++idx) {
        model.addConstr(exprs_E[idx] == 0);
    }
    
    // >= constraints
    std::vector<GRBLinExpr> exprs_G;
    for (int index_cons = 0; index_cons < model_parameters.AG.getRowLength(); ++index_cons) {
        GRBLinExpr expr_G;
        exprs_G.push_back(expr_G);
    }
    for (int col_idx = 0; col_idx < model_parameters.AG.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AG.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AG.getClen(col_idx) + beg_idx; ++idx) {
            exprs_G[model_parameters.AG.getRow(idx)] += model_parameters.AG.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bG.getNzeroLen(); ++idx) {
        exprs_G[model_parameters.bG.getLoc(idx)] -= model_parameters.bG.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AG.getRowLength(); ++idx) {
        model.addConstr(exprs_G[idx] >= 0);
    }
    
    // second-stage constraints
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        // equality constraints D y + C x + I_slack y_slack = e
        for (int idx_scenario = 0; idx_scenario < sample_size_per_replication; ++idx_scenario) {
            std::vector<GRBLinExpr> exprs_eq_rep;
            for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
                GRBLinExpr expr;
                exprs_eq_rep.push_back(expr);
            }
            // coefficients before y; Dy
            for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.D.getCbeg(col_idx);
                for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
                    exprs_eq_rep[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[idx_rep][idx_scenario][col_idx];
                }
            }
            // slack variables; I_slack y_slack
            if (model_parameters.num_slack_2ndStage > 0){
                for (int col_idx = 0; col_idx < model_parameters.ISlack.getColLength(); ++col_idx) {
                    int beg_idx = model_parameters.ISlack.getCbeg(col_idx);
                    for (int idx = beg_idx; idx < model_parameters.ISlack.getClen(col_idx) + beg_idx; ++idx) {
                        exprs_eq_rep[model_parameters.ISlack.getRow(idx)] += model_parameters.ISlack.getVal(idx) * y_slack[idx_rep][idx_scenario][col_idx];
                    }
                }
            }
            // coefficients before x (deterministic part)
            for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.C.getCbeg(col_idx);
                for (int idx = beg_idx; idx < model_parameters.C.getClen(col_idx) + beg_idx; ++idx) {
                    exprs_eq_rep[model_parameters.C.getRow(idx)] += model_parameters.C.getVal(idx) * x[col_idx];
                }
            }
            // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
            for (int idx_C = 0; idx_C < aggregate_outputs[idx_rep].Xi[idx_scenario].C.size(); ++idx_C) {
                exprs_eq_rep[stochastic_map.C_map[idx_C].first] += aggregate_outputs[idx_rep].Xi[idx_scenario].C[idx_C] * x[stochastic_map.C_map[idx_C].second];
            }
            // right hand side e
            for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
                exprs_eq_rep[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
            }
            // right hand side (stochastic part) equality e_(i) equality
            for (int idx = 0; idx < aggregate_outputs[idx_rep].Xi[idx_scenario].e.size(); ++idx) {
                exprs_eq_rep[stochastic_map.e_map[idx]] -= aggregate_outputs[idx_rep].Xi[idx_scenario].e[idx];
            }
            // add constraints
            for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
                model.addConstr(exprs_eq_rep[index_eq] == 0);
            }
        }
    } // end for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep)
    
    // optimize
    model.optimize();
    // create outputs
    compromiseOutput compromise_res;
    // solution
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        compromise_res.compromise_decision.push_back(x[idx_x].get(GRB_DoubleAttr_X));
    }
    compromise_res.replications = aggregate_outputs;
    //model.write("/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/lands/lands_model.lp");
    return compromise_res;
}


// compromise decision process to coordinate replication and aggregation
compromiseOutput saa_2slp_compromise(const std::string& folder_path,
                                     std::mt19937& generator,
                                     int sample_size_per_replication,
                                     int num_replications,
                                     double regularizer_coefficient) {
    // Initialization
    std::string model_path = folder_path + "/model.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(CD_SAA_2slpv1.0).txt";
    std::string sto_path = folder_path + "/sto.txt";
    // create sto object
    stoMap stochastic_map;
    // read sto file
    stochastic_map = readStochasticMap(sto_path);
    // read model file
    standardTwoStageParameters model_parameters = readStandardTwoStageParameters(model_path);
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    
    
    std::cout << "*******************************************\n";
    writeFile << "*******************************************\n";
    std::cout << "Compromise Decision for SAA for Solving Two-Stage Stochastic Linear Programs (v1.0)\n";
    writeFile << "Compromise Decision for SAA for Solving Two-Stage Stochastic Linear Programs (v1.0)\n";
    std::cout << "Parameters\n";
    writeFile << "Parameters\n";
    std::cout << "sample_size_per_replication: " << sample_size_per_replication << std::endl;
    writeFile << "sample_size_per_replication: " << sample_size_per_replication << std::endl;
    std::cout << "num_replications: " << num_replications << std::endl;
    writeFile << "num_replications: " << num_replications << std::endl;
    std::cout << "regularizer_coefficient: " << regularizer_coefficient << std::endl;
    writeFile << "regularizer_coefficient: " << regularizer_coefficient << std::endl;
    // Step 1: Replication
    std::vector<solverOutput> aggregate_outputs;
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        // *** algorithm in the replication step ***
        solverOutput replication_outputs = saa_2slp(folder_path, generator, sample_size_per_replication, idx_rep);
        // *** ***
        aggregate_outputs.push_back(replication_outputs);
    }
    
    // Step 2: Aggregation
    compromiseOutput compromise_res = saa_2slp_compromise_master(model_parameters, stochastic_map, aggregate_outputs, sample_size_per_replication, num_replications, regularizer_coefficient);
    std::cout << "Compromise Decision: \n";
    writeFile << "Compromise Decision: \n";
    for (int idx_x = 0; idx_x < compromise_res.compromise_decision.size() - 1; ++idx_x) {
        std::cout << compromise_res.compromise_decision[idx_x] << ", ";
        writeFile << compromise_res.compromise_decision[idx_x] << ", ";
    }
    std::cout << compromise_res.compromise_decision[compromise_res.compromise_decision.size()  - 1] << std::endl;
    writeFile << compromise_res.compromise_decision[compromise_res.compromise_decision.size()  - 1] << std::endl;
    std::cout << "*******************************************\n";
    writeFile << "*******************************************\n";
    writeFile.close();
    return compromise_res;
}

compromiseOutput saa_sqqp_compromise(const std::string& folder_path,
                                     std::mt19937& generator,
                                     int sample_size_per_replication,
                                     int num_replications,
                                     double regularizer_coefficient) {
    // Initialization
    std::string model_path = folder_path + "/model.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(CD_SAA_2sqqpv1.0).txt";
    std::string sto_path = folder_path + "/sto.txt";
    // create sto object
    stoMap stochastic_map;
    // read sto file
    stochastic_map = readStochasticMap(sto_path);
    // read model file
    standardTwoStageParameters model_parameters = readStandardTwoStageParameters(model_path);
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    
    
    std::cout << "*******************************************\n";
    writeFile << "*******************************************\n";
    std::cout << "Compromise Decision for Single-Cut Benders Decomposition (v1.0)\n";
    writeFile << "Compromise Decision for Single-Cut Benders Decomposition (v1.0)\n";
    std::cout << "Parameters\n";
    writeFile << "Parameters\n";
    std::cout << "sample_size_per_replication: " << sample_size_per_replication << std::endl;
    writeFile << "sample_size_per_replication: " << sample_size_per_replication << std::endl;
    std::cout << "num_replications: " << num_replications << std::endl;
    writeFile << "num_replications: " << num_replications << std::endl;
    std::cout << "regularizer_coefficient: " << regularizer_coefficient << std::endl;
    writeFile << "regularizer_coefficient: " << regularizer_coefficient << std::endl;
    // Step 1: Replication
    std::vector<solverOutput> aggregate_outputs;
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        // *** algorithm in the replication step ***
        solverOutput replication_outputs = saa_sqqp(folder_path, generator, sample_size_per_replication);
        // *** ***
        aggregate_outputs.push_back(replication_outputs);
    }
    
    // Step 2: Aggregation
    compromiseOutput compromise_res = saa_sqqp_compromise_master(model_parameters, stochastic_map, aggregate_outputs, sample_size_per_replication, num_replications, regularizer_coefficient);
    std::cout << "Compromise Decision: \n";
    writeFile << "Compromise Decision: \n";
    for (int idx_x = 0; idx_x < compromise_res.compromise_decision.size() - 1; ++idx_x) {
        std::cout << compromise_res.compromise_decision[idx_x] << ", ";
        writeFile << compromise_res.compromise_decision[idx_x] << ", ";
    }
    std::cout << compromise_res.compromise_decision[compromise_res.compromise_decision.size()  - 1] << std::endl;
    writeFile << compromise_res.compromise_decision[compromise_res.compromise_decision.size()  - 1] << std::endl;
    std::cout << "*******************************************\n";
    writeFile << "*******************************************\n";
    writeFile.close();
    return compromise_res;
}

// compute the ground truth (need to assume finite support)
solverOutput ground_truth_discrete(const std::string& folder_path) {
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string true_dist_path = folder_path + "/true_dist.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(groundTruth_v1.0).txt";
    // convert all the paths into constant chars
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    // read model file
    standardTwoStageParameters model_parameters = readStandardTwoStageParameters(model_path);
    // read sto file
    stoMap stochastic_map = readStochasticMap(sto_path);
    // read ground truth
    std::vector<stoPoint> true_dist = readGroundTruth_discreteRV(true_dist_path);
    // --- Finish Reading ---
    // --- gurobi solver environment ---
    // Create the Gurobi environment
    GRBEnv env = GRBEnv();
    
    // Create an empty model
    GRBModel model = GRBModel(env);
    
    // create decision variables
    // first-stage decision variable
    std::vector<GRBVar> x;
    x.reserve(model_parameters.num_var_1stStage); // Pre-allocate memory
    
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        x.push_back(model.addVar(model_parameters.var_lb_1stStage, model_parameters.var_ub_1stStage, 0.0, GRB_CONTINUOUS));
    }
    // second-stage decision variable
    std::vector<std::vector<GRBVar>> y;
    // slack variable
    std::vector<std::vector<GRBVar>> y_slack;
    long total_num_var_2ndStage = model_parameters.num_var_2ndStage + model_parameters.num_slack_2ndStage;
    for (int idx_scenario = 0; idx_scenario < true_dist.size(); ++idx_scenario) {
        std::vector<GRBVar> y_scenario;
        y.reserve(total_num_var_2ndStage);
        for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
            y_scenario.push_back(model.addVar(model_parameters.var_lb_2ndStage,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
        }
        y.push_back(y_scenario);
        if (model_parameters.num_slack_2ndStage > 0) {
            std::vector<GRBVar> y_slack_scenario;
            for (int idx_y_slack = 0; idx_y_slack < model_parameters.num_slack_2ndStage; ++idx_y_slack) {
                y_slack_scenario.push_back(model.addVar(0,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
            }
            y_slack.push_back(y_slack_scenario);
        }
    }
    // objective
    GRBLinExpr expr_obj;
    // first stage objective
    for (int idx_x = 0; idx_x < model_parameters.c.getNzeroLen(); ++idx_x) {
        expr_obj += model_parameters.c.getVal(idx_x) * x[idx_x];
    }
    // second stage objective
    for (int idx_scenario = 0; idx_scenario < true_dist.size(); ++idx_scenario) {
        for (int idx_y = 0; idx_y < model_parameters.d.getNzeroLen(); ++idx_y) {
            expr_obj += true_dist[idx_scenario].prob * model_parameters.d.getVal(idx_y) * y[idx_scenario][model_parameters.d.getLoc(idx_y)];
        }
    }
    
    // add objective
    model.setObjective(expr_obj, GRB_MINIMIZE);
    
    // constrants
    // first stage constraints
    // <= constrants
    std::vector<GRBLinExpr> exprs_L;
    for (int index_cons = 0; index_cons < model_parameters.AL.getRowLength(); ++index_cons) {
        GRBLinExpr expr_L;
        exprs_L.push_back(expr_L);
    }
    for (int col_idx = 0; col_idx < model_parameters.AL.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AL.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AL.getClen(col_idx) + beg_idx; ++idx) {
            exprs_L[model_parameters.AL.getRow(idx)] += model_parameters.AL.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bL.getNzeroLen(); ++idx) {
        exprs_L[model_parameters.bL.getLoc(idx)] -= model_parameters.bL.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AL.getRowLength(); ++idx) {
        model.addConstr(exprs_L[idx] <= 0);
    }
    
    // == constraints
    std::vector<GRBLinExpr> exprs_E;
    for (int index_cons = 0; index_cons < model_parameters.AE.getRowLength(); ++index_cons) {
        GRBLinExpr expr_E;
        exprs_E.push_back(expr_E);
    }
    for (int col_idx = 0; col_idx < model_parameters.AE.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AE.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AE.getClen(col_idx) + beg_idx; ++idx) {
            exprs_E[model_parameters.AE.getRow(idx)] += model_parameters.AE.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bE.getNzeroLen(); ++idx) {
        exprs_E[model_parameters.bE.getLoc(idx)] -= model_parameters.bE.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AE.getRowLength(); ++idx) {
        model.addConstr(exprs_E[idx] == 0);
    }
    
    // >= constraints
    std::vector<GRBLinExpr> exprs_G;
    for (int index_cons = 0; index_cons < model_parameters.AG.getRowLength(); ++index_cons) {
        GRBLinExpr expr_G;
        exprs_G.push_back(expr_G);
    }
    for (int col_idx = 0; col_idx < model_parameters.AG.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AG.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AG.getClen(col_idx) + beg_idx; ++idx) {
            exprs_G[model_parameters.AG.getRow(idx)] += model_parameters.AG.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bG.getNzeroLen(); ++idx) {
        exprs_G[model_parameters.bG.getLoc(idx)] -= model_parameters.bG.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AG.getRowLength(); ++idx) {
        model.addConstr(exprs_G[idx] >= 0);
    }
    
    // second stage
    // equality constraints D y + C x + I_slack y_slack = e
    for (int idx_scenario = 0; idx_scenario < true_dist.size(); ++idx_scenario) {
        std::vector<GRBLinExpr> exprs_eq;
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            GRBLinExpr expr;
            exprs_eq.push_back(expr);
        }
        // coefficients before y; Dy
        for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.D.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[idx_scenario][col_idx];
            }
        }
        // slack variables; I_slack y_slack
        if (model_parameters.num_slack_2ndStage > 0){
            for (int col_idx = 0; col_idx < model_parameters.ISlack.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.ISlack.getCbeg(col_idx);
                for (int idx = beg_idx; idx < model_parameters.ISlack.getClen(col_idx) + beg_idx; ++idx) {
                    exprs_eq[model_parameters.ISlack.getRow(idx)] += model_parameters.ISlack.getVal(idx) * y_slack[idx_scenario][col_idx];
                }
            }
        }
        // coefficients before x (deterministic part)
        for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.C.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.C.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.C.getRow(idx)] += model_parameters.C.getVal(idx) * x[col_idx];
            }
        }
        // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
        for (int idx_C = 0; idx_C < true_dist[idx_scenario].C.size(); ++idx_C) {
            exprs_eq[stochastic_map.C_map[idx_C].first] += true_dist[idx_scenario].C[idx_C] * x[stochastic_map.C_map[idx_C].second];
        }
        // right hand side e
        for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
            exprs_eq[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
        }
        // right hand side (stochastic part) equality e_(i) equality
        for (int idx = 0; idx < true_dist[idx_scenario].e.size(); ++idx) {
            exprs_eq[stochastic_map.e_map[idx]] -= true_dist[idx_scenario].e[idx];
        }
        // add constraints
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            model.addConstr(exprs_eq[index_eq] == 0);
        }
    }
    // optimize
    model.optimize();
    // create outputs
    solverOutput res;
    std::cout << "---Ground Truth (Finite Discrete Support)---\n";
    writeFile << "---Ground Truth (Finite Discrete Support)---\n";
    // objective
    std::cout << "Optimal Value (True Dist): " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    writeFile << "Optimal Value (True Dist): " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    res.obj = model.get(GRB_DoubleAttr_ObjVal);
    // solution
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        res.sol.push_back(x[idx_x].get(GRB_DoubleAttr_X));
    }
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage - 1; ++idx_x) {
        std::cout << res.sol[idx_x] << ", ";
        writeFile << res.sol[idx_x] << ", ";
    }
    std::cout << res.sol[res.sol.size() - 1] << ", " << std::endl;
    writeFile << res.sol[res.sol.size() - 1] << ", " << std::endl;
    res.Xi = true_dist;
    //model.write(folder_path + "/ground_truth_mod.lp");
    writeFile.close();
    return res;
}

solverOutput sqqp_ground_truth_discrete(const std::string& folder_path) {
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string true_dist_path = folder_path + "/true_dist.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(sqqp_groundTruth_v1.0).txt";
    // convert all the paths into constant chars
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    // read model file
    standardTwoStageParameters model_parameters = readStandardTwoStageParameters(model_path);
    // read sto file
    stoMap stochastic_map = readStochasticMap(sto_path);
    // read ground truth
    std::vector<stoPoint> true_dist = readGroundTruth_discreteRV(true_dist_path);
    // --- Finish Reading ---
    // --- gurobi solver environment ---
    // Create the Gurobi environment
    GRBEnv env = GRBEnv();
    
    // Create an empty model
    GRBModel model = GRBModel(env);
    
    // create decision variables
    // first-stage decision variable
    std::vector<GRBVar> x;
    x.reserve(model_parameters.num_var_1stStage); // Pre-allocate memory
    
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        x.push_back(model.addVar(model_parameters.var_lb_1stStage, model_parameters.var_ub_1stStage, 0.0, GRB_CONTINUOUS));
    }
    // second-stage decision variable
    std::vector<std::vector<GRBVar>> y;
    // slack variable
    std::vector<std::vector<GRBVar>> y_slack;
    long total_num_var_2ndStage = model_parameters.num_var_2ndStage + model_parameters.num_slack_2ndStage;
    for (int idx_scenario = 0; idx_scenario < true_dist.size(); ++idx_scenario) {
        std::vector<GRBVar> y_scenario;
        y.reserve(total_num_var_2ndStage);
        for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
            y_scenario.push_back(model.addVar(model_parameters.var_lb_2ndStage,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
        }
        y.push_back(y_scenario);
        if (model_parameters.num_slack_2ndStage > 0) {
            std::vector<GRBVar> y_slack_scenario;
            for (int idx_y_slack = 0; idx_y_slack < model_parameters.num_slack_2ndStage; ++idx_y_slack) {
                y_slack_scenario.push_back(model.addVar(0,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
            }
            y_slack.push_back(y_slack_scenario);
        }
    }
    // objective
    GRBQuadExpr expr_obj;
    
    // first stage objective
    for (int idx_x = 0; idx_x < model_parameters.c.getNzeroLen(); ++idx_x) {
        expr_obj += model_parameters.c.getVal(idx_x) * x[idx_x];
    }
    // quadratic terms
    for (int col_idx = 0; col_idx < model_parameters.Q.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.Q.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.Q.getClen(col_idx) + beg_idx; ++idx) {
            expr_obj += 0.5 * model_parameters.Q.getVal(idx) * x[model_parameters.Q.getRow(idx)] * x[col_idx];
        }
    }
    // second stage objective
    for (int idx_scenario = 0; idx_scenario < true_dist.size(); ++idx_scenario) {
        for (int idx_y = 0; idx_y < model_parameters.d.getNzeroLen(); ++idx_y) {
            expr_obj += true_dist[idx_scenario].prob * model_parameters.d.getVal(idx_y) * y[idx_scenario][model_parameters.d.getLoc(idx_y)];
        }
        // quadratic terms
        for (int col_idx = 0; col_idx < model_parameters.P.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.P.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.P.getClen(col_idx) + beg_idx; ++idx) {
                if (model_parameters.P.getRow(idx) < model_parameters.num_var_2ndStage) {
                    if (col_idx < model_parameters.num_var_2ndStage) {
                        // case 1: P_ij * y_i * y_j
                        expr_obj += 0.5 * true_dist[idx_scenario].prob * model_parameters.P.getVal(idx) * y[idx_scenario][model_parameters.P.getRow(idx)] * y[idx_scenario][col_idx];
                    }
                    else {
                        // case 2: P_ij * y_i * yslack_j
                        expr_obj += 0.5 * true_dist[idx_scenario].prob * model_parameters.P.getVal(idx) * y[idx_scenario][model_parameters.P.getRow(idx)] * y_slack[idx_scenario][col_idx - model_parameters.num_var_2ndStage];
                    }
                }
                else {
                    if (col_idx < model_parameters.num_var_2ndStage) {
                        // case 3: P_ij * yslack_i * y_j
                        expr_obj += 0.5 * true_dist[idx_scenario].prob * model_parameters.P.getVal(idx) * y_slack[idx_scenario][model_parameters.P.getRow(idx) - model_parameters.num_var_2ndStage] * y[idx_scenario][col_idx];
                    }
                    else {
                        // case 4: P_ij * yslack_i * yslack_j
                        expr_obj += 0.5 * true_dist[idx_scenario].prob * model_parameters.P.getVal(idx) * y_slack[idx_scenario][model_parameters.P.getRow(idx) - model_parameters.num_var_2ndStage] * y_slack[idx_scenario][col_idx - model_parameters.num_var_2ndStage];
                    }
                }
            }
        }
    }
    
    // add objective
    model.setObjective(expr_obj, GRB_MINIMIZE);
    
    // constrants
    // first stage constraints
    // <= constrants
    std::vector<GRBLinExpr> exprs_L;
    for (int index_cons = 0; index_cons < model_parameters.AL.getRowLength(); ++index_cons) {
        GRBLinExpr expr_L;
        exprs_L.push_back(expr_L);
    }
    for (int col_idx = 0; col_idx < model_parameters.AL.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AL.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AL.getClen(col_idx) + beg_idx; ++idx) {
            exprs_L[model_parameters.AL.getRow(idx)] += model_parameters.AL.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bL.getNzeroLen(); ++idx) {
        exprs_L[model_parameters.bL.getLoc(idx)] -= model_parameters.bL.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AL.getRowLength(); ++idx) {
        model.addConstr(exprs_L[idx] <= 0);
    }
    
    // == constraints
    std::vector<GRBLinExpr> exprs_E;
    for (int index_cons = 0; index_cons < model_parameters.AE.getRowLength(); ++index_cons) {
        GRBLinExpr expr_E;
        exprs_E.push_back(expr_E);
    }
    for (int col_idx = 0; col_idx < model_parameters.AE.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AE.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AE.getClen(col_idx) + beg_idx; ++idx) {
            exprs_E[model_parameters.AE.getRow(idx)] += model_parameters.AE.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bE.getNzeroLen(); ++idx) {
        exprs_E[model_parameters.bE.getLoc(idx)] -= model_parameters.bE.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AE.getRowLength(); ++idx) {
        model.addConstr(exprs_E[idx] == 0);
    }
    
    // >= constraints
    std::vector<GRBLinExpr> exprs_G;
    for (int index_cons = 0; index_cons < model_parameters.AG.getRowLength(); ++index_cons) {
        GRBLinExpr expr_G;
        exprs_G.push_back(expr_G);
    }
    for (int col_idx = 0; col_idx < model_parameters.AG.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.AG.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.AG.getClen(col_idx) + beg_idx; ++idx) {
            exprs_G[model_parameters.AG.getRow(idx)] += model_parameters.AG.getVal(idx) * x[col_idx];
        }
    }
    // right hand side
    for (int idx = 0; idx < model_parameters.bG.getNzeroLen(); ++idx) {
        exprs_G[model_parameters.bG.getLoc(idx)] -= model_parameters.bG.getVal(idx);
    }
    // add constraints
    for (int idx = 0; idx < model_parameters.AG.getRowLength(); ++idx) {
        model.addConstr(exprs_G[idx] >= 0);
    }
    
    // second stage
    // equality constraints D y + C x + I_slack y_slack = e
    for (int idx_scenario = 0; idx_scenario < true_dist.size(); ++idx_scenario) {
        std::vector<GRBLinExpr> exprs_eq;
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            GRBLinExpr expr;
            exprs_eq.push_back(expr);
        }
        // coefficients before y; Dy
        for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.D.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[idx_scenario][col_idx];
            }
        }
        // slack variables; I_slack y_slack
        if (model_parameters.num_slack_2ndStage > 0){
            for (int col_idx = 0; col_idx < model_parameters.ISlack.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.ISlack.getCbeg(col_idx);
                for (int idx = beg_idx; idx < model_parameters.ISlack.getClen(col_idx) + beg_idx; ++idx) {
                    exprs_eq[model_parameters.ISlack.getRow(idx)] += model_parameters.ISlack.getVal(idx) * y_slack[idx_scenario][col_idx];
                }
            }
        }
        // coefficients before x (deterministic part)
        for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.C.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.C.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.C.getRow(idx)] += model_parameters.C.getVal(idx) * x[col_idx];
            }
        }
        // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
        for (int idx_C = 0; idx_C < true_dist[idx_scenario].C.size(); ++idx_C) {
            exprs_eq[stochastic_map.C_map[idx_C].first] += true_dist[idx_scenario].C[idx_C] * x[stochastic_map.C_map[idx_C].second];
        }
        // right hand side e
        for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
            exprs_eq[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
        }
        // right hand side (stochastic part) equality e_(i) equality
        for (int idx = 0; idx < true_dist[idx_scenario].e.size(); ++idx) {
            exprs_eq[stochastic_map.e_map[idx]] -= true_dist[idx_scenario].e[idx];
        }
        // add constraints
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            model.addConstr(exprs_eq[index_eq] == 0);
        }
    }
    
    // optimize
    model.optimize();
    // create outputs
    solverOutput res;
    std::cout << "---Ground Truth (Finite Discrete Support)---\n";
    writeFile << "---Ground Truth (Finite Discrete Support)---\n";
    // objective
    std::cout << "Optimal Value (True Dist): " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    writeFile << "Optimal Value (True Dist): " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    res.obj = model.get(GRB_DoubleAttr_ObjVal);
    // solution
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        res.sol.push_back(x[idx_x].get(GRB_DoubleAttr_X));
    }
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage - 1; ++idx_x) {
        std::cout << res.sol[idx_x] << ", ";
        writeFile << res.sol[idx_x] << ", ";
    }
    std::cout << res.sol[res.sol.size() - 1] << ", " << std::endl;
    writeFile << res.sol[res.sol.size() - 1] << ", " << std::endl;
    res.Xi = true_dist;
    model.write("/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleQP/ground_truth_sqqp_model.lp");
    writeFile.close();
    return res;
}


// validation
double validation_ground_truth_discrete(const std::string& folder_path,
                                              std::vector<double>& x) {
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string true_dist_path = folder_path + "/true_dist.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(validationGroundTruth_v1.0).txt";
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    // convert all the paths into constant chars
    // read model file
    standardTwoStageParameters model_parameters = readStandardTwoStageParameters(model_path);
    // read sto file
    stoMap stochastic_map = readStochasticMap(sto_path);
    // read ground truth
    std::vector<stoPoint> true_dist = readGroundTruth_discreteRV(true_dist_path);
    // --- Finish Reading ---
    // --- gurobi solver environment ---
    // Create the Gurobi environment
    GRBEnv env = GRBEnv();
    
    // Create an empty model
    GRBModel model = GRBModel(env);
    
    // second-stage decision variable
    std::vector<std::vector<GRBVar>> y;
    // slack variable
    std::vector<std::vector<GRBVar>> y_slack;
    long total_num_var_2ndStage = model_parameters.num_var_2ndStage + model_parameters.num_slack_2ndStage;
    for (int idx_scenario = 0; idx_scenario < true_dist.size(); ++idx_scenario) {
        std::vector<GRBVar> y_scenario;
        y.reserve(total_num_var_2ndStage);
        for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
            y_scenario.push_back(model.addVar(model_parameters.var_lb_2ndStage,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
        }
        y.push_back(y_scenario);
        if (model_parameters.num_slack_2ndStage > 0) {
            std::vector<GRBVar> y_slack_scenario;
            for (int idx_y_slack = 0; idx_y_slack < model_parameters.num_slack_2ndStage; ++idx_y_slack) {
                y_slack_scenario.push_back(model.addVar(0,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
            }
            y_slack.push_back(y_slack_scenario);
        }
    }
    
    // objective
    GRBLinExpr expr_obj;
    // first stage objective
    for (int idx_x = 0; idx_x < model_parameters.c.getNzeroLen(); ++idx_x) {
        expr_obj += model_parameters.c.getVal(idx_x) * x[idx_x];
    }
    // second stage objective
    for (int idx_scenario = 0; idx_scenario < true_dist.size(); ++idx_scenario) {
        for (int idx_y = 0; idx_y < model_parameters.d.getNzeroLen(); ++idx_y) {
            expr_obj += true_dist[idx_scenario].prob * model_parameters.d.getVal(idx_y) * y[idx_scenario][model_parameters.d.getLoc(idx_y)];
        }
    }
    // add objective
    model.setObjective(expr_obj, GRB_MINIMIZE);
    
    // second stage
    // equality constraints D y + C x + I_slack y_slack = e
    for (int idx_scenario = 0; idx_scenario < true_dist.size(); ++idx_scenario) {
        std::vector<GRBLinExpr> exprs_eq;
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            GRBLinExpr expr;
            exprs_eq.push_back(expr);
        }
        // coefficients before y; Dy
        for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.D.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[idx_scenario][col_idx];
            }
        }
        // slack variables; I_slack y_slack
        if (model_parameters.num_slack_2ndStage > 0){
            for (int col_idx = 0; col_idx < model_parameters.ISlack.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.ISlack.getCbeg(col_idx);
                for (int idx = beg_idx; idx < model_parameters.ISlack.getClen(col_idx) + beg_idx; ++idx) {
                    exprs_eq[model_parameters.ISlack.getRow(idx)] += model_parameters.ISlack.getVal(idx) * y_slack[idx_scenario][col_idx];
                }
            }
        }
        // coefficients before x (deterministic part)
        for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.C.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.C.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.C.getRow(idx)] += model_parameters.C.getVal(idx) * x[col_idx];
            }
        }
        // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
        for (int idx_C = 0; idx_C < true_dist[idx_scenario].C.size(); ++idx_C) {
            exprs_eq[stochastic_map.C_map[idx_C].first] += true_dist[idx_scenario].C[idx_C] * x[stochastic_map.C_map[idx_C].second];
        }
        // right hand side e
        for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
            exprs_eq[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
        }
        // right hand side (stochastic part) equality e_(i) equality
        for (int idx = 0; idx < true_dist[idx_scenario].e.size(); ++idx) {
            exprs_eq[stochastic_map.e_map[idx]] -= true_dist[idx_scenario].e[idx];
        }
        // add constraints
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            model.addConstr(exprs_eq[index_eq] == 0);
        }
    }
    // optimize
    model.optimize();
    // create outputs
    double res = model.get(GRB_DoubleAttr_ObjVal);
    std::cout << "---Validation Results---\n";
    writeFile << "---Validation Results---\n";
    // objective
    std::cout << "Objective (True Dist): " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    writeFile << "Objective (True Dist): " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    // Decision
    std::cout << "Evaluated Decision: ";
    writeFile << "Evaluated Decision: ";
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage - 1; ++idx_x) {
        std::cout << x[idx_x] << ", ";
        writeFile << x[idx_x] << ", ";
    }
    std::cout << x[x.size() - 1] << std::endl;
    writeFile << x[x.size() - 1] << std::endl;
    std::cout << "--------------------------\n";
    writeFile << "--------------------------\n";
    //model.write("/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/lands/model_tureDist.lp");
    
    
    writeFile.close();
    return res;
}


// validation set is given
double validation_saa(const std::string& folder_path,
                          std::vector<double>& x) {
    // data file
    bool flag_data_e; // tell if e stochastic is generated
    bool flag_data_C; // tell if C stochastic is generated
    // create directory paths for database and model
    std::string e_DB_path = folder_path + "/e_DB.txt";
    std::string C_DB_path = folder_path + "/C_DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(validation_SAAv1.0).txt";
    // convert all the paths into constant chars
    const char* e_DB_path_const = e_DB_path.c_str();
    const char* C_DB_path_const = C_DB_path.c_str();
    // create stream object
    std::ifstream readFile_e(e_DB_path_const);
    std::ifstream readFile_C(C_DB_path_const);
    // create database
    std::vector<std::vector<dataPoint>> e_DB;
    std::vector<std::vector<dataPoint>> C_DB;
    // create model structure
    standardTwoStageParameters model_parameters;
    // create sto object
    stoMap stochastic_map;
    // read  be
    long sample_size = 0;
    if (readFile_e.is_open()) {
        std::cout << "e_DB (RHS in the 2nd stage) data file is found." << std::endl;
        readFile_e.close(); // close the file
        // read be database
        e_DB = readDB(e_DB_path);
        sample_size = e_DB[0].size();
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
        sample_size = C_DB[0].size();
        flag_data_C = true;
    }
    else {
        readFile_C.close(); // close the file
        flag_data_C = false;
        std::cout << "C_DB (2nd stage) data file is not found!" << std::endl;
    }
    // read model file
    model_parameters = readStandardTwoStageParameters(model_path);
    // read sto file
    stochastic_map = readStochasticMap(sto_path);
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    // --- Finish Reading ---
    // --- gurobi solver environment ---
    // Create the Gurobi environment
    GRBEnv env = GRBEnv();
    
    // Create an empty model
    GRBModel model = GRBModel(env);
    
    // second-stage decision variable
    std::vector<std::vector<GRBVar>> y;
    // slack variable
    std::vector<std::vector<GRBVar>> y_slack;
    long total_num_var_2ndStage = model_parameters.num_var_2ndStage + model_parameters.num_slack_2ndStage;
    for (int idx_scenario = 0; idx_scenario < sample_size; ++idx_scenario) {
        std::vector<GRBVar> y_scenario;
        y.reserve(total_num_var_2ndStage);
        for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
            y_scenario.push_back(model.addVar(model_parameters.var_lb_2ndStage,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
        }
        y.push_back(y_scenario);
        if (model_parameters.num_slack_2ndStage > 0) {
            std::vector<GRBVar> y_slack_scenario;
            for (int idx_y_slack = 0; idx_y_slack < model_parameters.num_slack_2ndStage; ++idx_y_slack) {
                y_slack_scenario.push_back(model.addVar(0,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
            }
            y_slack.push_back(y_slack_scenario);
        }
    }
    
    // objective
    double one_over_sample_size = 1.0 / ((double) sample_size);
    GRBLinExpr expr_obj;
    // first stage objective
    for (int idx_x = 0; idx_x < model_parameters.c.getNzeroLen(); ++idx_x) {
        expr_obj += model_parameters.c.getVal(idx_x) * x[idx_x];
    }
    // second stage objective
    for (int idx_scenario = 0; idx_scenario < sample_size; ++idx_scenario) {
        for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
            expr_obj += one_over_sample_size * model_parameters.d.getVal(idx_y) * y[idx_scenario][model_parameters.d.getLoc(idx_y)];
        }
    }
    // add objective
    model.setObjective(expr_obj, GRB_MINIMIZE);
    
    // second stage
    // equality constraints D y + C x + I_slack y_slack = e
    for (int idx_scenario = 0; idx_scenario < sample_size; ++idx_scenario) {
        std::vector<GRBLinExpr> exprs_eq;
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            GRBLinExpr expr;
            exprs_eq.push_back(expr);
        }
        // coefficients before y; Dy
        for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.D.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[idx_scenario][col_idx];
            }
        }
        // slack variables; I_slack y_slack
        if (model_parameters.num_slack_2ndStage > 0){
            for (int col_idx = 0; col_idx < model_parameters.ISlack.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.ISlack.getCbeg(col_idx);
                for (int idx = beg_idx; idx < model_parameters.ISlack.getClen(col_idx) + beg_idx; ++idx) {
                    exprs_eq[model_parameters.ISlack.getRow(idx)] += model_parameters.ISlack.getVal(idx) * y_slack[idx_scenario][col_idx];
                }
            }
        }
        // coefficients before x (deterministic part)
        for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.C.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.C.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.C.getRow(idx)] += model_parameters.C.getVal(idx) * x[col_idx];
            }
        }
        // *** stochastic part ***
        // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
        if (flag_data_C == true) {
            for (int idx_C = 0; idx_C < C_DB[0][idx_scenario].response.size(); ++idx_C) {
                exprs_eq[stochastic_map.C_map[idx_C].first] += C_DB[0][idx_scenario].response[idx_C] * x[stochastic_map.C_map[idx_C].second];
            }
        }
        
        // right hand side e
        for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
            exprs_eq[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
        }
        // right hand side (stochastic part) equality e_(i) equality
        if (flag_data_e == true) {
            for (int idx = 0; idx < e_DB[0][idx_scenario].response.size(); ++idx) {
                exprs_eq[stochastic_map.e_map[idx]] -= e_DB[0][idx_scenario].response[idx];
            }
        }
        
        // add constraints
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            model.addConstr(exprs_eq[index_eq] == 0);
        }
    }
    // optimize
    model.optimize();
    // create outputs
    double res = model.get(GRB_DoubleAttr_ObjVal);
    std::cout << "---Validation Results by SAA---\n";
    writeFile << "---Validation Results by SAA---\n";
    // objective
    std::cout << "Objective (Validation by SAA): " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    writeFile << "Objective (Validation by SAA): " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    // Decision
    std::cout << "Evaluated Decision: ";
    writeFile << "Evaluated Decision: ";
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage - 1; ++idx_x) {
        std::cout << x[idx_x] << ", ";
        writeFile << x[idx_x] << ", ";
    }
    std::cout << x[x.size() - 1] << std::endl;
    writeFile << x[x.size() - 1] << std::endl;
    std::cout << "--------------------------\n";
    writeFile << "--------------------------\n";
    //model.write("/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/lands/model_saa.lp");
    
    writeFile.close();
    return res;
}


double sqqp_validation_ground_truth_discrete(const std::string& folder_path, std::vector<double>& x) {
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string true_dist_path = folder_path + "/true_dist.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(sqqp_validationGroundTruth_v1.0).txt";
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    // convert all the paths into constant chars
    // read model file
    standardTwoStageParameters model_parameters = readStandardTwoStageParameters(model_path);
    // read sto file
    stoMap stochastic_map = readStochasticMap(sto_path);
    // read ground truth
    std::vector<stoPoint> true_dist = readGroundTruth_discreteRV(true_dist_path);
    // --- Finish Reading ---
    // --- gurobi solver environment ---
    // Create the Gurobi environment
    GRBEnv env = GRBEnv();
    
    // Create an empty model
    GRBModel model = GRBModel(env);
    
    // second-stage decision variable
    std::vector<std::vector<GRBVar>> y;
    // slack variable
    std::vector<std::vector<GRBVar>> y_slack;
    long total_num_var_2ndStage = model_parameters.num_var_2ndStage + model_parameters.num_slack_2ndStage;
    for (int idx_scenario = 0; idx_scenario < true_dist.size(); ++idx_scenario) {
        std::vector<GRBVar> y_scenario;
        y.reserve(total_num_var_2ndStage);
        for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
            y_scenario.push_back(model.addVar(model_parameters.var_lb_2ndStage,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
        }
        y.push_back(y_scenario);
        if (model_parameters.num_slack_2ndStage > 0) {
            std::vector<GRBVar> y_slack_scenario;
            for (int idx_y_slack = 0; idx_y_slack < model_parameters.num_slack_2ndStage; ++idx_y_slack) {
                y_slack_scenario.push_back(model.addVar(0,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
            }
            y_slack.push_back(y_slack_scenario);
        }
    }
    
    // objective
    GRBQuadExpr expr_obj;
    
    // first stage objective
    for (int idx_x = 0; idx_x < model_parameters.c.getNzeroLen(); ++idx_x) {
        expr_obj += model_parameters.c.getVal(idx_x) * x[idx_x];
    }
    // quadratic terms
    for (int col_idx = 0; col_idx < model_parameters.Q.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.Q.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.Q.getClen(col_idx) + beg_idx; ++idx) {
            expr_obj += 0.5 * model_parameters.Q.getVal(idx) * x[model_parameters.Q.getRow(idx)] * x[col_idx];
        }
    }
    // second stage objective
    for (int idx_scenario = 0; idx_scenario < true_dist.size(); ++idx_scenario) {
        for (int idx_y = 0; idx_y < model_parameters.d.getNzeroLen(); ++idx_y) {
            expr_obj += true_dist[idx_scenario].prob * model_parameters.d.getVal(idx_y) * y[idx_scenario][model_parameters.d.getLoc(idx_y)];
        }
        // quadratic terms
        for (int col_idx = 0; col_idx < model_parameters.P.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.P.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.P.getClen(col_idx) + beg_idx; ++idx) {
                if (model_parameters.P.getRow(idx) < model_parameters.num_var_2ndStage) {
                    if (col_idx < model_parameters.num_var_2ndStage) {
                        // case 1: P_ij * y_i * y_j
                        expr_obj += 0.5 * true_dist[idx_scenario].prob * model_parameters.P.getVal(idx) * y[idx_scenario][model_parameters.P.getRow(idx)] * y[idx_scenario][col_idx];
                    }
                    else {
                        // case 2: P_ij * y_i * yslack_j
                        expr_obj += 0.5 * true_dist[idx_scenario].prob * model_parameters.P.getVal(idx) * y[idx_scenario][model_parameters.P.getRow(idx)] * y_slack[idx_scenario][col_idx - model_parameters.num_var_2ndStage];
                    }
                }
                else {
                    if (col_idx < model_parameters.num_var_2ndStage) {
                        // case 3: P_ij * yslack_i * y_j
                        expr_obj += 0.5 * true_dist[idx_scenario].prob * model_parameters.P.getVal(idx) * y_slack[idx_scenario][model_parameters.P.getRow(idx) - model_parameters.num_var_2ndStage] * y[idx_scenario][col_idx];
                    }
                    else {
                        // case 4: P_ij * yslack_i * yslack_j
                        expr_obj += 0.5 * true_dist[idx_scenario].prob * model_parameters.P.getVal(idx) * y_slack[idx_scenario][model_parameters.P.getRow(idx) - model_parameters.num_var_2ndStage] * y_slack[idx_scenario][col_idx - model_parameters.num_var_2ndStage];
                    }
                }
            }
        }
        
    }
    // add objective
    model.setObjective(expr_obj, GRB_MINIMIZE);
    
    // second stage
    // equality constraints D y + C x + I_slack y_slack = e
    for (int idx_scenario = 0; idx_scenario < true_dist.size(); ++idx_scenario) {
        std::vector<GRBLinExpr> exprs_eq;
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            GRBLinExpr expr;
            exprs_eq.push_back(expr);
        }
        // coefficients before y; Dy
        for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.D.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[idx_scenario][col_idx];
            }
        }
        // slack variables; I_slack y_slack
        if (model_parameters.num_slack_2ndStage > 0){
            for (int col_idx = 0; col_idx < model_parameters.ISlack.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.ISlack.getCbeg(col_idx);
                for (int idx = beg_idx; idx < model_parameters.ISlack.getClen(col_idx) + beg_idx; ++idx) {
                    exprs_eq[model_parameters.ISlack.getRow(idx)] += model_parameters.ISlack.getVal(idx) * y_slack[idx_scenario][col_idx];
                }
            }
        }
        // coefficients before x (deterministic part)
        for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.C.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.C.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.C.getRow(idx)] += model_parameters.C.getVal(idx) * x[col_idx];
            }
        }
        // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
        for (int idx_C = 0; idx_C < true_dist[idx_scenario].C.size(); ++idx_C) {
            exprs_eq[stochastic_map.C_map[idx_C].first] += true_dist[idx_scenario].C[idx_C] * x[stochastic_map.C_map[idx_C].second];
        }
        // right hand side e
        for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
            exprs_eq[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
        }
        // right hand side (stochastic part) equality e_(i) equality
        for (int idx = 0; idx < true_dist[idx_scenario].e.size(); ++idx) {
            exprs_eq[stochastic_map.e_map[idx]] -= true_dist[idx_scenario].e[idx];
        }
        // add constraints
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            model.addConstr(exprs_eq[index_eq] == 0);
        }
    }
    
    // optimize
    model.optimize();
    // create outputs
    double res = model.get(GRB_DoubleAttr_ObjVal);
    std::cout << "---Validation Results---\n";
    writeFile << "---Validation Results---\n";
    // objective
    std::cout << "Objective (True Dist): " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    writeFile << "Objective (True Dist): " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    // Decision
    std::cout << "Evaluated Decision: ";
    writeFile << "Evaluated Decision: ";
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage - 1; ++idx_x) {
        std::cout << x[idx_x] << ", ";
        writeFile << x[idx_x] << ", ";
    }
    std::cout << x[x.size() - 1] << std::endl;
    writeFile << x[x.size() - 1] << std::endl;
    std::cout << "--------------------------\n";
    writeFile << "--------------------------\n";
    //model.write("/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleQP/validation_tureDist.lp");
    
    
    writeFile.close();
    return res;
}

double sqqp_validation_saa(const std::string& folder_path,
                           std::vector<double>& x) {
    // data file
    bool flag_data_e; // tell if e stochastic is generated
    bool flag_data_C; // tell if C stochastic is generated
    // create directory paths for database and model
    std::string e_DB_path = folder_path + "/e_DB.txt";
    std::string C_DB_path = folder_path + "/C_DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(sqqp_validation_SAAv1.0).txt";
    // convert all the paths into constant chars
    const char* e_DB_path_const = e_DB_path.c_str();
    const char* C_DB_path_const = C_DB_path.c_str();
    // create stream object
    std::ifstream readFile_e(e_DB_path_const);
    std::ifstream readFile_C(C_DB_path_const);
    // create database
    std::vector<std::vector<dataPoint>> e_DB;
    std::vector<std::vector<dataPoint>> C_DB;
    // create model structure
    standardTwoStageParameters model_parameters;
    // create sto object
    stoMap stochastic_map;
    // read  be
    long sample_size = 0;
    if (readFile_e.is_open()) {
        std::cout << "e_DB (RHS in the 2nd stage) data file is found." << std::endl;
        readFile_e.close(); // close the file
        // read be database
        e_DB = readDB(e_DB_path);
        sample_size = e_DB[0].size();
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
        sample_size = C_DB[0].size();
        flag_data_C = true;
    }
    else {
        readFile_C.close(); // close the file
        flag_data_C = false;
        std::cout << "C_DB (2nd stage) data file is not found!" << std::endl;
    }
    // read model file
    model_parameters = readStandardTwoStageParameters(model_path);
    // read sto file
    stochastic_map = readStochasticMap(sto_path);
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    // --- Finish Reading ---
    // --- gurobi solver environment ---
    // Create the Gurobi environment
    GRBEnv env = GRBEnv();
    
    // Create an empty model
    GRBModel model = GRBModel(env);
    
    // second-stage decision variable
    std::vector<std::vector<GRBVar>> y;
    // slack variable
    std::vector<std::vector<GRBVar>> y_slack;
    long total_num_var_2ndStage = model_parameters.num_var_2ndStage + model_parameters.num_slack_2ndStage;
    for (int idx_scenario = 0; idx_scenario < sample_size; ++idx_scenario) {
        std::vector<GRBVar> y_scenario;
        y.reserve(total_num_var_2ndStage);
        for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
            y_scenario.push_back(model.addVar(model_parameters.var_lb_2ndStage,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
        }
        y.push_back(y_scenario);
        if (model_parameters.num_slack_2ndStage > 0) {
            std::vector<GRBVar> y_slack_scenario;
            for (int idx_y_slack = 0; idx_y_slack < model_parameters.num_slack_2ndStage; ++idx_y_slack) {
                y_slack_scenario.push_back(model.addVar(0,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
            }
            y_slack.push_back(y_slack_scenario);
        }
    }
    
    // objective
    double one_over_sample_size = 1.0 / ((double) sample_size);
    GRBQuadExpr expr_obj;
    // first-stage objective c^\top x
    for (int idx_x = 0; idx_x < model_parameters.c.getNzeroLen(); ++idx_x) {
        expr_obj += model_parameters.c.getVal(idx_x) * x[idx_x];
    }
    // first-stage quadratic terms 0.5 * x^\top Q x
    for (int col_idx = 0; col_idx < model_parameters.Q.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.Q.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.Q.getClen(col_idx) + beg_idx; ++idx) {
            expr_obj += 0.5 * model_parameters.Q.getVal(idx) * x[model_parameters.Q.getRow(idx)] * x[col_idx];
        }
    }
    // second stage objective
    for (int idx_scenario = 0; idx_scenario < sample_size; ++idx_scenario) {
        for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
            expr_obj += one_over_sample_size * model_parameters.d.getVal(idx_y) * y[idx_scenario][model_parameters.d.getLoc(idx_y)];
        }
        // quadratic terms
        for (int col_idx = 0; col_idx < model_parameters.P.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.P.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.P.getClen(col_idx) + beg_idx; ++idx) {
                if (model_parameters.P.getRow(idx) < model_parameters.num_var_2ndStage) {
                    if (col_idx < model_parameters.num_var_2ndStage) {
                        // case 1: P_ij * y_i * y_j
                        expr_obj += 0.5 * one_over_sample_size * model_parameters.P.getVal(idx) * y[idx_scenario][model_parameters.P.getRow(idx)] * y[idx_scenario][col_idx];
                    }
                    else {
                        // case 2: P_ij * y_i * yslack_j
                        expr_obj += 0.5 * one_over_sample_size * model_parameters.P.getVal(idx) * y[idx_scenario][model_parameters.P.getRow(idx)] * y_slack[idx_scenario][col_idx - model_parameters.num_var_2ndStage];
                    }
                }
                else {
                    if (col_idx < model_parameters.num_var_2ndStage) {
                        // case 3: P_ij * yslack_i * y_j
                        expr_obj += 0.5 * one_over_sample_size * model_parameters.P.getVal(idx) * y_slack[idx_scenario][model_parameters.P.getRow(idx) - model_parameters.num_var_2ndStage] * y[idx_scenario][col_idx];
                    }
                    else {
                        // case 4: P_ij * yslack_i * yslack_j
                        expr_obj += 0.5 * one_over_sample_size * model_parameters.P.getVal(idx) * y_slack[idx_scenario][model_parameters.P.getRow(idx) - model_parameters.num_var_2ndStage] * y_slack[idx_scenario][col_idx - model_parameters.num_var_2ndStage];
                    }
                }
            }
        }
    }
    // add objective
    model.setObjective(expr_obj, GRB_MINIMIZE);
    
    // second stage
    // equality constraints D y + C x + I_slack y_slack = e
    for (int idx_scenario = 0; idx_scenario < sample_size; ++idx_scenario) {
        std::vector<GRBLinExpr> exprs_eq;
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            GRBLinExpr expr;
            exprs_eq.push_back(expr);
        }
        // coefficients before y; Dy
        for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.D.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[idx_scenario][col_idx];
            }
        }
        // slack variables; I_slack y_slack
        if (model_parameters.num_slack_2ndStage > 0){
            for (int col_idx = 0; col_idx < model_parameters.ISlack.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.ISlack.getCbeg(col_idx);
                for (int idx = beg_idx; idx < model_parameters.ISlack.getClen(col_idx) + beg_idx; ++idx) {
                    exprs_eq[model_parameters.ISlack.getRow(idx)] += model_parameters.ISlack.getVal(idx) * y_slack[idx_scenario][col_idx];
                }
            }
        }
        // coefficients before x (deterministic part)
        for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.C.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.C.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.C.getRow(idx)] += model_parameters.C.getVal(idx) * x[col_idx];
            }
        }
        // *** stochastic part ***
        // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
        if (flag_data_C == true) {
            for (int idx_C = 0; idx_C < C_DB[0][idx_scenario].response.size(); ++idx_C) {
                exprs_eq[stochastic_map.C_map[idx_C].first] += C_DB[0][idx_scenario].response[idx_C] * x[stochastic_map.C_map[idx_C].second];
            }
        }
        
        // right hand side e
        for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
            exprs_eq[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
        }
        // right hand side (stochastic part) equality e_(i) equality
        if (flag_data_e == true) {
            for (int idx = 0; idx < e_DB[0][idx_scenario].response.size(); ++idx) {
                exprs_eq[stochastic_map.e_map[idx]] -= e_DB[0][idx_scenario].response[idx];
            }
        }
        
        // add constraints
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            model.addConstr(exprs_eq[index_eq] == 0);
        }
    }
    // optimize
    model.optimize();
    // create outputs
    double res = model.get(GRB_DoubleAttr_ObjVal);
    std::cout << "---Validation Results by SAA---\n";
    writeFile << "---Validation Results by SAA---\n";
    // objective
    std::cout << "Objective (Validation by SAA): " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    writeFile << "Objective (Validation by SAA): " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    // Decision
    std::cout << "Evaluated Decision: ";
    writeFile << "Evaluated Decision: ";
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage - 1; ++idx_x) {
        std::cout << x[idx_x] << ", ";
        writeFile << x[idx_x] << ", ";
    }
    std::cout << x[x.size() - 1] << std::endl;
    writeFile << x[x.size() - 1] << std::endl;
    std::cout << "--------------------------\n";
    writeFile << "--------------------------\n";
    //model.write("/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/lands/model_saa.lp");
    
    writeFile.close();
    return res;
} // end sqqp_validation_saa

// random number generator
stoPoint generator_random(standardTwoStageParameters& model_parameters,
                          const stoMap& stochastic_map,
                          std::mt19937& generator,
                          int dataPoint_idx,
                          const std::vector<std::vector<dataPoint>>& e_DB,
                          const std::vector<std::vector<dataPoint>>& C_DB,
                          int idx_rep) {
    stoPoint est;
    // e
    if (e_DB.size() < 1) {
        if (stochastic_map.e_flag_dist.compare("DISCRETE") == 0) { // discrete random variable
            for (int idx_e = 0; idx_e < stochastic_map.e_discrete.size(); ++idx_e) {
                double e_component = 0;
                std::uniform_real_distribution<double> uniform(0,1);
                double temp_val = uniform(generator);
                e_component = discrete_random_number_generator(temp_val, stochastic_map.e_discrete[idx_e]);
                est.e.push_back(e_component);
            }
        } // end if (stochastic_map.e_flag_dist.compare("DISCRETE") == 0)
        else if (stochastic_map.e_flag_dist.compare("UNIFORM") == 0) { // uniform random variable
            for (int idx_e = 0; idx_e < stochastic_map.e_uniform_lower.size(); ++idx_e) {
                double e_component = uniform_random_number_generator(generator, stochastic_map.e_uniform_lower[idx_e], stochastic_map.e_uniform_upper[idx_e]);
                est.e.push_back(e_component);
            }
        } // end else if (stochastic_map.e_flag_dist.compare("UNIFORM") == 0)
        else if (stochastic_map.e_flag_dist.compare("NORMAL") == 0) { // normal random variable
            for (int idx_e = 0; idx_e < stochastic_map.e_normal_mean.size(); ++idx_e) {
                double e_component = normal_random_number_generator(generator, stochastic_map.e_normal_mean[idx_e], stochastic_map.e_normal_stddev[idx_e]);
                est.e.push_back(e_component);
            }
        }
        else if (stochastic_map.e_flag_dist.compare("TRUN_NORMAL") == 0) { // truncated normal random variable
            for (int idx_e = 0; idx_e < stochastic_map.e_normal_mean.size(); ++idx_e) {
                double e_component = truncated_normal_random_number_generator(generator, stochastic_map.e_normal_mean[idx_e], stochastic_map.e_normal_stddev[idx_e], stochastic_map.e_normal_lb[idx_e], stochastic_map.e_normal_ub[idx_e]);
                est.e.push_back(e_component);
            }
        }
    } // if (e_DB.size() < 1)
    else if (e_DB.size() > 0 && stochastic_map.e_flag_dist.compare("DATA") == 0){
        est.e = e_DB[idx_rep][dataPoint_idx].response;
    }
    // C
    if (C_DB.size() < 1) {
        if (stochastic_map.C_flag_dist.compare("DISCRETE") == 0) { // discrete random variable
            for (int idx_C = 0; idx_C < stochastic_map.C_discrete.size(); ++idx_C) {
                double C_component = 0;
                std::uniform_real_distribution<double> uniform(0,1);
                double temp_val = uniform(generator);
                C_component = discrete_random_number_generator(temp_val, stochastic_map.C_discrete[idx_C]);
                est.C.push_back(C_component);
            }
        } // end if (stochastic_map.C_flag_dist.compare("DISCRETE") == 0)
        else if (stochastic_map.C_flag_dist.compare("UNIFORM") == 0) { // uniform random variable
            for (int idx_C = 0; idx_C < stochastic_map.C_uniform_lower.size(); ++idx_C) {
                double C_component = uniform_random_number_generator(generator, stochastic_map.C_uniform_lower[idx_C], stochastic_map.C_uniform_upper[idx_C]);
                est.C.push_back(C_component);
            }
        } // end else if (stochastic_map.C_flag_dist.compare("UNIFORM") == 0)
        else if (stochastic_map.C_flag_dist.compare("NORMAL") == 0) { // normal random variable
            for (int idx_C = 0; idx_C < stochastic_map.C_normal_mean.size(); ++idx_C) {
                double C_component = normal_random_number_generator(generator, stochastic_map.C_normal_mean[idx_C], stochastic_map.C_normal_stddev[idx_C]);
                est.C.push_back(C_component);
            }
        }
        else if (stochastic_map.C_flag_dist.compare("TRUN_NORMAL") == 0) { // truncated normal random variable
            for (int idx_C = 0; idx_C < stochastic_map.e_normal_mean.size(); ++idx_C) {
                double C_component = truncated_normal_random_number_generator(generator, stochastic_map.C_normal_mean[idx_C], stochastic_map.C_normal_stddev[idx_C], stochastic_map.C_normal_lb[idx_C], stochastic_map.C_normal_ub[idx_C]);
                est.C.push_back(C_component);
            }
        }
    } // end if (C_DB.size() < 1 && stochastic_map.C_flag_dist.compare("DATA") != 0)
    else if (C_DB.size() > 0 && stochastic_map.C_flag_dist.compare("DATA") == 0){
        est.C = C_DB[idx_rep][dataPoint_idx].response;
    }
    return  est;
}

