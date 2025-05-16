//
//  CuttingPlane_solver.cpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/27/25.
//

#include "CuttingPlane_solver.hpp"

double SOLVER_INF = 1e9;

cp_solution cp_master(standardTwoStageParameters& model_parameters,
                      cp_subproblem_cuts& subproblem_constraints) {
    // --- gurobi solver environment ---
    // Create the Gurobi environment
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);  // Turn off all console logging
    // Create an empty model
    GRBModel model = GRBModel(env);
    
    // create decision variables
    // first-stage decision variable
    std::vector<GRBVar> x;
    x.reserve(model_parameters.num_var_1stStage); // Pre-allocate memory
    
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        x.push_back(model.addVar(model_parameters.var_lb_1stStage, model_parameters.var_ub_1stStage, 0.0, GRB_CONTINUOUS));
    }
    // epigraphical representation of the second stage
    GRBVar z = model.addVar(-1e9, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
    
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
    expr_obj += z;
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
    
    // subproblem cuts
    for (int idx_cut = 0; idx_cut < subproblem_constraints.alpha_array.size(); ++idx_cut) {
        GRBLinExpr expr_sub;
        // + beta * x
        for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
            expr_sub += subproblem_constraints.beta_array[idx_cut][idx_x] * x[idx_x];
        }
        // + alpha
        expr_sub += subproblem_constraints.alpha_array[idx_cut];
        // - z
        expr_sub -= z;
        //
        model.addConstr(expr_sub <= 0);
        
    }
    // optimize
    model.optimize();
    //
    cp_solution res;
    // solution
    // Check if optimal
    int status = model.get(GRB_IntAttr_Status);
    if (status == GRB_OPTIMAL) {
        for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
            res.x.push_back(x[idx_x].get(GRB_DoubleAttr_X));
        }
        res.obj = model.get(GRB_DoubleAttr_ObjVal);
        res.z = z.get(GRB_DoubleAttr_X);
        res.sol_flag = 0;
    }
    else if (status == GRB_INFEASIBLE) {
        for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
            res.x.push_back(0);
        }
        res.z = 0;
        res.sol_flag = 1;
        std::cout << "Warning: Model is infeasible. Return 0 vector." << std::endl;
    }
    else if (status == GRB_UNBOUNDED) {
        for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
            res.x.push_back(x[idx_x].get(GRB_DoubleAttr_X));
        }
        res.z = z.get(GRB_DoubleAttr_X);
        res.sol_flag = 1;
        std::cout << "Warning: Model is unbounded. " << std::endl;
    }
    else {
        for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
            res.x.push_back(0);
        }
        res.z = 0;
        res.sol_flag = 1;
        std::cout << "Warning: Other cases. " << std::endl;
    }
    //model.write("/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/lands/cp_master_large.lp");
    
    return res;
}


compromise_cp_output cp_compromise_master(standardTwoStageParameters& model_parameters,
                                         std::vector<cp_output>& aggregate_outputs,
                                         double regularizer_coefficient) {
    long num_replications = aggregate_outputs.size();
    double one_over_num_replications = 1.0 / ((double) num_replications);
    // --- gurobi solver environment ---
    // Create the Gurobi environment
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);  // Turn off all console logging
    // Create an empty model
    GRBModel model = GRBModel(env);
    
    // create decision variables
    // first-stage decision variable
    std::vector<GRBVar> x;
    x.reserve(model_parameters.num_var_1stStage); // Pre-allocate memory
    
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        x.push_back(model.addVar(model_parameters.var_lb_1stStage, model_parameters.var_ub_1stStage, 0.0, GRB_CONTINUOUS));
    }
    // epigraphical representation of the second stage
    std::vector<GRBVar> eta;
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        eta.push_back(model.addVar(-1e9, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
    }
    
    // quadratic objective
    GRBQuadExpr expr_obj;
    // first stage objective
    for (int idx_x = 0; idx_x < model_parameters.c.getNzeroLen(); ++idx_x) {
        expr_obj += model_parameters.c.getVal(idx_x) * x[idx_x];
        // regularizer
        expr_obj += 0.5 * regularizer_coefficient * x[idx_x] * x[idx_x];
        double x_average = 0;
        // regularizer
        for (int idx_rep = 0; idx_rep < aggregate_outputs.size(); ++idx_rep) {
            x_average += aggregate_outputs[idx_rep].x[idx_x];
        }
        x_average *= one_over_num_replications;
        expr_obj -= regularizer_coefficient * x_average * x[idx_x];
    }
    // quadratic terms
    for (int col_idx = 0; col_idx < model_parameters.Q.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.Q.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.Q.getClen(col_idx) + beg_idx; ++idx) {
            expr_obj += 0.5 * model_parameters.Q.getVal(idx) * x[model_parameters.Q.getRow(idx)] * x[col_idx];
        }
    }
    // epigraphical representation of each replication
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        expr_obj += one_over_num_replications * eta[idx_rep];
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
    
    // subproblem cuts
    for (int idx_rep = 0; idx_rep < aggregate_outputs.size(); ++idx_rep) {
        for (int idx_cut = 0; idx_cut < aggregate_outputs[idx_rep].subproblem_constraints.alpha_array.size(); ++idx_cut) {
            GRBLinExpr expr_sub;
            // + beta * x
            for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
                expr_sub += aggregate_outputs[idx_rep].subproblem_constraints.beta_array[idx_cut][idx_x] * x[idx_x];
            }
            // + alpha
            expr_sub += aggregate_outputs[idx_rep].subproblem_constraints.alpha_array[idx_cut];
            // - z
            expr_sub -= eta[idx_rep];
            //
            model.addConstr(expr_sub <= 0);
            
        }
    }
    // optimize
    model.optimize();
    
    //
    compromise_cp_output res;
    // solution
    // Check if optimal
    int status = model.get(GRB_IntAttr_Status);
    if (status == GRB_OPTIMAL) {
        for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
            res.compromise_decision.push_back(x[idx_x].get(GRB_DoubleAttr_X));
        }
    }
    else {
        std::cout << "Warning: Compromise Decision problem is not infeasible or unbounded." << std::endl;
    }
    //model.write("/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/lands/model_compromise_large.lp");
    res.replications = aggregate_outputs;
    return res;
}


// cutting plane algorithm (assume relative complete recourse)
cp_output cp_solver(const std::string& folder_path,
                    std::mt19937& generator,
                    int sample_size,
                    double error) {
    // timer
    std::clock_t time_start;
    time_start = std::clock();
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
    std::string resultsOutput_path = folder_path + "/computationalResults(rep_CPv1.0).txt";
    // convert all the paths into constant chars
    const char* e_DB_path_const = e_DB_path.c_str();
    const char* C_DB_path_const = C_DB_path.c_str();
    // create stream object
    std::ifstream readFile_e(e_DB_path_const);
    std::ifstream readFile_C(C_DB_path_const);
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
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
    // generate random variables
    for (int idx_scenario = 0; idx_scenario < sample_size; ++idx_scenario) {
        // obtain a new data point
        stoPoint xi_cur = generator_random(model_parameters, stochastic_map, generator, idx_scenario, e_DB, C_DB);
        Xi.push_back(xi_cur);
    } // end for (int idx_scenario = 0; idx_scenario < sample_size; ++idx_scenario)
    // --- Finish Reading ---
    // --- start settig up subproblem ---
    // --- gurobi solver environment for subproblem---
    // Create the Gurobi environment
    GRBEnv env_sub = GRBEnv();
    env_sub.set(GRB_IntParam_OutputFlag, 0);  // Turn off all console logging
    // Create an empty model
    GRBModel model_sub = GRBModel(env_sub);
    
    // second-stage decision variable
    std::vector<GRBVar> y;
    // slack variable
    std::vector<GRBVar> y_slack;
    // initialize decision variables
    y.reserve(model_parameters.num_var_2ndStage);
    for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
        y.push_back(model_sub.addVar(model_parameters.var_lb_2ndStage,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
    }
    if (model_parameters.num_slack_2ndStage > 0) {
        std::vector<GRBVar> y_slack_scenario;
        for (int idx_y_slack = 0; idx_y_slack < model_parameters.num_slack_2ndStage; ++idx_y_slack) {
            y_slack.push_back(model_sub.addVar(0,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
        }
    }
    // objective of the second-stage subproblem d^\top y
    GRBQuadExpr expr_obj_sub;
    for (int idx = 0; idx < model_parameters.d.getNzeroLen(); ++idx) {
        expr_obj_sub += model_parameters.d.getVal(idx) * y[model_parameters.d.getLoc(idx)];
    }
    // quadratic terms 0.5 y^\top P y
    for (int col_idx = 0; col_idx < model_parameters.P.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.P.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.P.getClen(col_idx) + beg_idx; ++idx) {
            if (model_parameters.P.getRow(idx) < model_parameters.num_var_2ndStage) {
                if (col_idx < model_parameters.num_var_2ndStage) {
                    // case 1: P_ij * y_i * y_j
                    expr_obj_sub += 0.5 * model_parameters.P.getVal(idx) * y[model_parameters.P.getRow(idx)] * y[col_idx];
                }
                else {
                    // case 2: P_ij * y_i * yslack_j
                    expr_obj_sub += 0.5 * model_parameters.P.getVal(idx) * y[model_parameters.P.getRow(idx)] * y_slack[col_idx - model_parameters.num_var_2ndStage];
                }
            }
            else {
                if (col_idx < model_parameters.num_var_2ndStage) {
                    // case 3: P_ij * yslack_i * y_j
                    expr_obj_sub += 0.5 * model_parameters.P.getVal(idx) * y_slack[model_parameters.P.getRow(idx) - model_parameters.num_var_2ndStage] * y[col_idx];
                }
                else {
                    // case 4: P_ij * yslack_i * yslack_j
                    expr_obj_sub += 0.5 * model_parameters.P.getVal(idx) * y_slack[model_parameters.P.getRow(idx) - model_parameters.num_var_2ndStage] * y_slack[col_idx - model_parameters.num_var_2ndStage];
                }
            }
        }
    }
    // add objective
    model_sub.setObjective(expr_obj_sub, GRB_MINIMIZE);
    
    // constraints of the second-stage subproblem
    // stndard form equality constraints D y + Islack y_slack + Cx = e
    std::vector<GRBLinExpr> exprs_eq_sub;
    for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
        GRBLinExpr expr;
        exprs_eq_sub.push_back(expr);
    }
    // coefficients before y; Dy
    for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.D.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
            exprs_eq_sub[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[col_idx];
        }
    }
    // slack variables; I_slack y_slack
    if (model_parameters.num_slack_2ndStage > 0){
        for (int col_idx = 0; col_idx < model_parameters.ISlack.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.ISlack.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.ISlack.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq_sub[model_parameters.ISlack.getRow(idx)] += model_parameters.ISlack.getVal(idx) * y_slack[col_idx];
            }
        }
    }
    // right hand side e
    for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
        exprs_eq_sub[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
    }
    // add the equality constraints
    std::vector<GRBConstr> constraintsEquality_sub;
    for (int idx = 0; idx< model_parameters.D.getRowLength(); ++idx) {
        constraintsEquality_sub.push_back(model_sub.addConstr(exprs_eq_sub[idx] == 0));
    }
    // intermediate values for rhs update
    // deterministic part of the rhs_bounds
    std::vector<double> e_det(model_parameters.D.getRowLength(), 0.0);
    for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
        e_det[model_parameters.e.getLoc(idx)] += model_parameters.e.getVal(idx);
    }
    std::vector<double> det_rhs_bounds = e_det;
    std::vector<double> rhs_bounds(model_parameters.D.getRowLength(), 0.0);
    // --- end settig up subproblem ---
    // STEP 3 MAIN LOOP
    std::vector<double> x_est;
    double z_est;
    cp_subproblem_cuts subproblem_cons;
    std::cout << "Start Main Loop\n";
    writeFile << "Start Main Loop\n";
    bool flag_init = true;
    bool flag_terminate = false;
    int it_num = 1;
    double res_max_gap = 1e10;
    double one_over_sample_size = 1.0 / ((double) sample_size);
    double opt_obj = SOLVER_INF;
    while (flag_terminate == false) {
        std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
        writeFile << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
        std::cout << "(Cutting Plane Algorithm) Iteration:  " << it_num << std::endl;
        writeFile << "(Cutting Plane Algorithm) Iteration:  " << it_num << std::endl;
        // solve the relaxed master problem
        // update estimated solution (step 1 of the cuttng plane algorithm)
        cp_solution est_sol = cp_master(model_parameters, subproblem_cons);
        opt_obj = est_sol.obj;
        if (est_sol.sol_flag == 0) {
            std::cout << "Master problem is feasible.\n";
            writeFile << "Master problem is feasible.\n";
            x_est = est_sol.x;
            z_est = est_sol.z;
        }
        else if (it_num > 1) {
            std::cout << "Warning: Master problem is infeasible. The program is terminated.\n";
            writeFile << "Warning: Master problem is infeasible. The program is terminated.\n";
            break;
        }
        else {
            std::cout << "First iteration, Master problem is infeasible or unbounded.\n";
            writeFile << "First iteration, Master problem is infeasible or unbounded.\n";
            std::cout << "Initial Solution is used as the starting point.\n";
            writeFile << "Initial Solution is used as the starting point.\n";
            z_est = -1e10;
            for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
                x_est = est_sol.x;
            }
        }
        flag_terminate = true;
        // update deterministic part of the rhs_bounds at current x_est
        // Dy + I_slack y_slack = [e - Cx]
        det_rhs_bounds = e_det;
        // det C x
        for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.C.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.C.getClen(col_idx) + beg_idx; ++idx) {
                det_rhs_bounds[model_parameters.C.getRow(idx)] -= model_parameters.C.getVal(idx) * x_est[col_idx];
            }
        }
        //
        double intercept_2nd_stage = 0;
        std::vector<double> slope_2nd_stage(model_parameters.num_var_1stStage, 0);
        double cost_2nd_stage = 0;
        for (int scenario_idx = 0; scenario_idx < sample_size; ++scenario_idx) {
            // update the subproblem
            for (int idx_row = 0; idx_row < model_parameters.D.getRowLength(); ++idx_row) {
                rhs_bounds[idx_row] = det_rhs_bounds[idx_row];
            } // update the deterministic part
            // update the stochastic parts of e
            for (int idx_e = 0; idx_e < Xi[scenario_idx].e.size(); ++idx_e) {
                rhs_bounds[stochastic_map.e_map[idx_e]] += Xi[scenario_idx].e[idx_e];
            }
            // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
            for (int idx_C = 0; idx_C < Xi[scenario_idx].C.size(); ++idx_C) {
                rhs_bounds[stochastic_map.C_map[idx_C].first] -= Xi[scenario_idx].C[idx_C] * x_est[stochastic_map.C_map[idx_C].second];
            }
            // update the RHS
            for (int idx_row = 0; idx_row < rhs_bounds.size(); ++idx_row) {
                constraintsEquality_sub[idx_row].set(GRB_DoubleAttr_RHS, rhs_bounds[idx_row]);
            }
            // solve the subproblem
            //model_sub.write("/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/lands/lands_model_sub.lp");
            model_sub.optimize();
            int status = model_sub.get(GRB_IntAttr_Status);
            if (status == GRB_OPTIMAL) {
                cost_2nd_stage += model_sub.get(GRB_DoubleAttr_ObjVal);
                // calculate the dual multipliers
                std::vector<double> pi;
                for (int idx_row = 0; idx_row < model_parameters.C.getRowLength(); ++idx_row) {
                    pi.push_back(constraintsEquality_sub[idx_row].get(GRB_DoubleAttr_Pi));
                }
                // compute the optimality cut
                /*
                // deterministic e
                double pi_e = model_parameters.e.fast_dotProduct(pi);
                // stochastic e
                // equality part
                for (int idx_e = 0; idx_e < Xi[scenario_idx].e.size(); ++idx_e) {
                    //std::cout << "Xi[scenario_idx].e[idx_e]: " << Xi[scenario_idx].e[idx_e] << std::endl;
                    pi_e += pi[stochastic_map.e_map[idx_e]] * Xi[scenario_idx].e[idx_e];
                }
                 */
                double intercept = model_sub.get(GRB_DoubleAttr_ObjVal);
                // slope
                std::vector<double> neg_pi_C = model_parameters.C.fast_rightMultiply(pi);
                // negate
                for (int idx2 = 0; idx2 < neg_pi_C.size(); ++idx2) {
                    neg_pi_C[idx2] = (-1.0) * neg_pi_C[idx2];
                }
                // stochastic C
                // equality
                for (int idx_C = 0; idx_C < Xi[scenario_idx].C.size(); ++idx_C) {
                    neg_pi_C[stochastic_map.C_map[idx_C].second] -= Xi[scenario_idx].C[idx_C] * pi[stochastic_map.C_map[idx_C].first];
                }
                // alpha = f(x') - \beta^\top x'
                intercept -= neg_pi_C * x_est;
                // update
                intercept_2nd_stage += intercept;
                for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
                    slope_2nd_stage[idx_x] += neg_pi_C[idx_x];
                }
            }
            else {
                writeFile << "Error: subproblem is infeasible" << std::endl;
                std::logic_error( "Error: subproblem is infeasible" );
            }
        } // end for (int scenario_idx = 0; scenario_idx < sample_size; ++scenario_idx)
        // multiply by sampling frequency
        intercept_2nd_stage = intercept_2nd_stage * one_over_sample_size;
        for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
            slope_2nd_stage[idx_x] = slope_2nd_stage[idx_x] * one_over_sample_size;
        }
        subproblem_cons.alpha_array.push_back(intercept_2nd_stage);
        subproblem_cons.beta_array.push_back(slope_2nd_stage);
        cost_2nd_stage = cost_2nd_stage * one_over_sample_size;
        // print cuts
        std::cout << "Cuts(Slope): \n";
        writeFile << "Cuts(Slope): \n";
        for (int idx_x = 0; idx_x < x_est.size() - 1; ++idx_x) {
            std::cout << slope_2nd_stage[idx_x] << ", ";
            writeFile << slope_2nd_stage[idx_x] << ", ";
        }
        std::cout << slope_2nd_stage[x_est.size() - 1] << std::endl;
        writeFile << slope_2nd_stage[x_est.size() - 1] << std::endl;
        
        std::cout << "Cuts(Intercept): " << intercept_2nd_stage << std::endl;
        writeFile << "Cuts(Intercept): " << intercept_2nd_stage << std::endl;
        double gap_2nd_stage = cost_2nd_stage - z_est;
        if (est_sol.sol_flag != 0 && flag_init == true) { // first iteration is infeasible or unbounded
            flag_terminate = false;
        }
        // termination criterion
        if (est_sol.sol_flag == 0 && cost_2nd_stage > z_est + error) {
            flag_terminate = false;
        }
        flag_init = false;
        std::cout << "Gap (2nd Stage Recourse): " << gap_2nd_stage << std::endl;
        writeFile << "Gap (2nd Stage Recourse): " << gap_2nd_stage << std::endl;
        std::cout << "Estimated First-Stage Solution: \n";
        writeFile << "Estimated First-Stage Solution: \n";
        for (int idx_x = 0; idx_x < x_est.size() - 1; ++idx_x) {
            std::cout << x_est[idx_x] << ", ";
            writeFile << x_est[idx_x] << ", ";
        }
        std::cout << x_est[x_est.size() - 1] << std::endl;
        writeFile << x_est[x_est.size() - 1] << std::endl;
        std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
        writeFile << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
        it_num += 1;
        // for debugging
        //if (it_num > 10) {
        //    break;
        //}
    } // end while (flag_terminate == false)
    // return estimated solutions
    std::cout << "Finish Solving Process\n";
    writeFile << "Finish Solving Process\n";
    std::cout << "Total Number of Iterations: " << it_num << std::endl;
    writeFile << "Total Number of Iterations: " << it_num << std::endl;
    std::cout << "Estimated First-Stage Solution: \n";
    writeFile << "Estimated First-Stage Solution: \n";
    for (int idx_x = 0; idx_x < x_est.size() - 1; ++idx_x) {
        std::cout << x_est[idx_x] << ", ";
        writeFile << x_est[idx_x] << ", ";
    }
    std::cout << x_est[x_est.size() - 1] << std::endl;
    writeFile << x_est[x_est.size() - 1] << std::endl;
    // write time elapsed
    double duration = (std::clock() - time_start ) / (double) CLOCKS_PER_SEC;
    //double opt_obj = model_parameters.c.fast_dotProduct(x_est) + z_est;
    std::cout << "Optimal Value: " << opt_obj << std::endl;
    writeFile << "Optimal Value: " << opt_obj << std::endl;
    writeFile << "Time elapsed(secs) : " << duration << "\n";
    cp_output res;
    res.x = x_est;
    res.max_gap = res_max_gap;
    if (flag_terminate == true) { // 0: Problem is solved; 1: Problem is not solved
        res.sol_flag = 0;
    }
    else {
        res.sol_flag = 1;
    }
    res.it_num = it_num;
    // store realizations
    res.Xi = Xi;
    // return cuts
    res.subproblem_constraints = subproblem_cons;
    // store solution time
    res.sol_time = duration;
    writeFile.close();
    return res;
}


// compromise decisioin
compromise_cp_output cp_compromise_solver(const std::string& folder_path,
                                          std::mt19937& generator,
                                          int sample_size_per_replication,
                                          int num_replications,
                                          double regularizer_coefficient,
                                          double error,
                                          double error2) {
    
    // Initialization
    std::string model_path = folder_path + "/model.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(CD_CPv1.0).txt";
    std::string sto_path = folder_path + "/sto.txt";
    // create sto object
    stoMap stochastic_map;
    // read sto file
    stochastic_map = readStochasticMap(sto_path);
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    // read model file
    standardTwoStageParameters model_parameters = readStandardTwoStageParameters(model_path);
    
    std::cout << "*******************************************\n";
    writeFile << "*******************************************\n";
    std::cout << "Compromise Decision for Cutting Plane Algorithm for Solving SQLP/SQQP(v1.0)\n";
    writeFile << "Compromise Decision for Cutting Plane Algorithm for Solving SQLP/SQQP(v1.0)\n";
    std::cout << "Parameters\n";
    writeFile << "Parameters\n";
    std::cout << "sample_size_per_replication: " << sample_size_per_replication << std::endl;
    writeFile << "sample_size_per_replication: " << sample_size_per_replication << std::endl;
    std::cout << "num_replications: " << num_replications << std::endl;
    writeFile << "num_replications: " << num_replications << std::endl;
    std::cout << "regularizer_coefficient: " << regularizer_coefficient << std::endl;
    writeFile << "regularizer_coefficient: " << regularizer_coefficient << std::endl;
    std::cout << "precision error in the replication step: " << error << std::endl;
    writeFile << "precision error in the replication step: " << error << std::endl;
    
    // Step 1: Replication
    std::vector<cp_output> aggregate_outputs;
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        cp_output replication_outputs = cp_solver(folder_path, generator, sample_size_per_replication, error);
        aggregate_outputs.push_back(replication_outputs);
    }
    
    // Step 2: Aggregation
    // set up subproblem
    // Create the Gurobi environment
    GRBEnv env_sub = GRBEnv();
    env_sub.set(GRB_IntParam_OutputFlag, 0);  // Turn off all console logging
    // Create an empty model
    GRBModel model_sub = GRBModel(env_sub);
    
    // second-stage decision variable
    std::vector<GRBVar> y;
    // slack variable
    std::vector<GRBVar> y_slack;
    // initialize decision variables
    y.reserve(model_parameters.num_var_2ndStage);
    for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
        y.push_back(model_sub.addVar(model_parameters.var_lb_2ndStage,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
    }
    if (model_parameters.num_slack_2ndStage > 0) {
        std::vector<GRBVar> y_slack_scenario;
        for (int idx_y_slack = 0; idx_y_slack < model_parameters.num_slack_2ndStage; ++idx_y_slack) {
            y_slack.push_back(model_sub.addVar(0,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
        }
    }
    // objective of the second-stage subproblem
    GRBQuadExpr expr_obj_sub;
    // linear terms d^\top y
    for (int idx = 0; idx < model_parameters.d.getNzeroLen(); ++idx) {
        expr_obj_sub += model_parameters.d.getVal(idx) * y[model_parameters.d.getLoc(idx)];
    }
    // quadratic terms 0.5 y^\top P y
    for (int col_idx = 0; col_idx < model_parameters.P.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.P.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.P.getClen(col_idx) + beg_idx; ++idx) {
            if (model_parameters.P.getRow(idx) < model_parameters.num_var_2ndStage) {
                if (col_idx < model_parameters.num_var_2ndStage) {
                    // case 1: P_ij * y_i * y_j
                    expr_obj_sub += 0.5 * model_parameters.P.getVal(idx) * y[model_parameters.P.getRow(idx)] * y[col_idx];
                }
                else {
                    // case 2: P_ij * y_i * yslack_j
                    expr_obj_sub += 0.5 * model_parameters.P.getVal(idx) * y[model_parameters.P.getRow(idx)] * y_slack[col_idx - model_parameters.num_var_2ndStage];
                }
            }
            else {
                if (col_idx < model_parameters.num_var_2ndStage) {
                    // case 3: P_ij * yslack_i * y_j
                    expr_obj_sub += 0.5 * model_parameters.P.getVal(idx) * y_slack[model_parameters.P.getRow(idx) - model_parameters.num_var_2ndStage] * y[col_idx];
                }
                else {
                    // case 4: P_ij * yslack_i * yslack_j
                    expr_obj_sub += 0.5 * model_parameters.P.getVal(idx) * y_slack[model_parameters.P.getRow(idx) - model_parameters.num_var_2ndStage] * y_slack[col_idx - model_parameters.num_var_2ndStage];
                }
            }
        }
    }
    // add objective
    model_sub.setObjective(expr_obj_sub, GRB_MINIMIZE);
    
    // constraints of the second-stage subproblem
    // stndard form equality constraints D y + Islack y_slack + Cx = e
    std::vector<GRBLinExpr> exprs_eq_sub;
    for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
        GRBLinExpr expr;
        exprs_eq_sub.push_back(expr);
    }
    // coefficients before y; Dy
    for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.D.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
            exprs_eq_sub[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[col_idx];
        }
    }
    // slack variables; I_slack y_slack
    if (model_parameters.num_slack_2ndStage > 0){
        for (int col_idx = 0; col_idx < model_parameters.ISlack.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.ISlack.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.ISlack.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq_sub[model_parameters.ISlack.getRow(idx)] += model_parameters.ISlack.getVal(idx) * y_slack[col_idx];
            }
        }
    }
    // right hand side e
    for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
        exprs_eq_sub[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
    }
    // add the equality constraints
    std::vector<GRBConstr> constraintsEquality_sub;
    for (int idx = 0; idx< model_parameters.D.getRowLength(); ++idx) {
        constraintsEquality_sub.push_back(model_sub.addConstr(exprs_eq_sub[idx] == 0));
    }
    // intermediate values for rhs update
    // deterministic part of the rhs_bounds
    std::vector<double> e_det(model_parameters.D.getRowLength(), 0.0);
    for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
        e_det[model_parameters.e.getLoc(idx)] += model_parameters.e.getVal(idx);
    }
    std::vector<double> det_rhs_bounds = e_det;
    std::vector<double> rhs_bounds(model_parameters.D.getRowLength(), 0.0);
    // compute cuts at each candidate solution
    for (int idx_rep1 = 0; idx_rep1 < num_replications; ++idx_rep1) {
        for (int idx_rep2 = 0; idx_rep2 < num_replications; ++idx_rep2) {
            // compute cuts
            // update deterministic part of the rhs_bounds at current x_est
            // Dy + I_slack y_slack = [e - Cx]
            det_rhs_bounds = e_det;
            // det C x
            for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.C.getCbeg(col_idx);
                for (int idx = beg_idx; idx < model_parameters.C.getClen(col_idx) + beg_idx; ++idx) {
                    det_rhs_bounds[model_parameters.C.getRow(idx)] -= model_parameters.C.getVal(idx) * aggregate_outputs[idx_rep2].x[col_idx];
                }
            }
            double intercept_2nd_stage = 0;
            std::vector<double> slope_2nd_stage(model_parameters.num_var_1stStage, 0);
            double cost_2nd_stage = 0;
            for (int scenario_idx = 0; scenario_idx < aggregate_outputs[idx_rep1].Xi.size(); ++scenario_idx) {
                // update the subproblem
                for (int idx_row = 0; idx_row < model_parameters.D.getRowLength(); ++idx_row) {
                    rhs_bounds[idx_row] = det_rhs_bounds[idx_row];
                } // update the deterministic part
                // update the stochastic parts of e
                for (int idx_e = 0; idx_e < aggregate_outputs[idx_rep1].Xi[scenario_idx].e.size(); ++idx_e) {
                    rhs_bounds[stochastic_map.e_map[idx_e]] += aggregate_outputs[idx_rep1].Xi[scenario_idx].e[idx_e];
                }
                // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
                for (int idx_C = 0; idx_C < aggregate_outputs[idx_rep1].Xi[scenario_idx].C.size(); ++idx_C) {
                    rhs_bounds[stochastic_map.C_map[idx_C].first] -= aggregate_outputs[idx_rep1].Xi[scenario_idx].C[idx_C] * aggregate_outputs[idx_rep2].x[stochastic_map.C_map[idx_C].second];
                }
                // update the RHS
                for (int idx_row = 0; idx_row < rhs_bounds.size(); ++idx_row) {
                    constraintsEquality_sub[idx_row].set(GRB_DoubleAttr_RHS, rhs_bounds[idx_row]);
                }
                model_sub.optimize();
                int status = model_sub.get(GRB_IntAttr_Status);
                if (status == GRB_OPTIMAL) {
                    cost_2nd_stage += model_sub.get(GRB_DoubleAttr_ObjVal);
                    // calculate the dual multipliers
                    std::vector<double> pi;
                    for (int idx_row = 0; idx_row < model_parameters.C.getRowLength(); ++idx_row) {
                        pi.push_back(constraintsEquality_sub[idx_row].get(GRB_DoubleAttr_Pi));
                    }
                    // compute the optimality cut
                    // deterministic e
                    double pi_e = model_parameters.e.fast_dotProduct(pi);
                    // stochastic e
                    // equality part
                    for (int idx_e = 0; idx_e < aggregate_outputs[idx_rep1].Xi[scenario_idx].e.size(); ++idx_e) {
                        //std::cout << "Xi[scenario_idx].e[idx_e]: " << Xi[scenario_idx].e[idx_e] << std::endl;
                        pi_e += pi[stochastic_map.e_map[idx_e]] * aggregate_outputs[idx_rep1].Xi[scenario_idx].e[idx_e];
                    }
                    // slope
                    std::vector<double> neg_pi_C = model_parameters.C.fast_rightMultiply(pi);
                    // negate
                    for (int idx2 = 0; idx2 < neg_pi_C.size(); ++idx2) {
                        neg_pi_C[idx2] = (-1.0) * neg_pi_C[idx2];
                    }
                    // stochastic C
                    // equality
                    for (int idx_C = 0; idx_C < aggregate_outputs[idx_rep1].Xi[scenario_idx].C.size(); ++idx_C) {
                        neg_pi_C[stochastic_map.C_map[idx_C].second] -= aggregate_outputs[idx_rep1].Xi[scenario_idx].C[idx_C] * pi[stochastic_map.C_map[idx_C].first];
                    }
                    // update
                    intercept_2nd_stage += pi_e;
                    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
                        slope_2nd_stage[idx_x] += neg_pi_C[idx_x];
                    }
                }
                else {
                    writeFile << "Error: subproblem is infeasible" << std::endl;
                    std::logic_error( "Error: subproblem is infeasible" );
                }
            } // for (int scenario_idx = 0; scenario_idx < sample_size; ++scenario_idx)
            double one_over_sample_size = 1.0 / ((double) aggregate_outputs[idx_rep1].Xi.size());
            // multiply by sampling frequency (epsilon-subgradient)
            intercept_2nd_stage = intercept_2nd_stage * one_over_sample_size - error2;
            for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
                slope_2nd_stage[idx_x] = slope_2nd_stage[idx_x] * one_over_sample_size;
            }
            aggregate_outputs[idx_rep1].subproblem_constraints.alpha_array.push_back(intercept_2nd_stage);
            aggregate_outputs[idx_rep1].subproblem_constraints.beta_array.push_back(slope_2nd_stage);
            cost_2nd_stage = cost_2nd_stage * one_over_sample_size;
            // print cuts
            std::cout << "Cuts(Slope): \n";
            writeFile << "Cuts(Slope): \n";
            for (int idx_x = 0; idx_x < aggregate_outputs[idx_rep2].x.size() - 1; ++idx_x) {
                std::cout << slope_2nd_stage[idx_x] << ", ";
                writeFile << slope_2nd_stage[idx_x] << ", ";
            }
            std::cout << slope_2nd_stage[aggregate_outputs[idx_rep2].x.size() - 1] << std::endl;
            writeFile << slope_2nd_stage[aggregate_outputs[idx_rep2].x.size() - 1] << std::endl;
            
            std::cout << "Cuts(Intercept): " << intercept_2nd_stage << std::endl;
            writeFile << "Cuts(Intercept): " << intercept_2nd_stage << std::endl;
        }
    }
    
    compromise_cp_output compromise_res = cp_compromise_master(model_parameters, aggregate_outputs, regularizer_coefficient);
    
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
