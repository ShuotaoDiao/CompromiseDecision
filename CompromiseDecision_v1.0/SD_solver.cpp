//
//  SD_solver.cpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/23/25.
//

#include "SD_solver.hpp"

// declare global variables (need to be defined in the source file)
const double SOLVER_PRECISION_LOWER = -1e-6;
const double SOLVER_PRECISION_UPPER = 1e-6;
const double SOLVER_PRECISION_LOWER2 = -1e-6; // alternative is -1e-6
const double SOLVER_PRECISION_UPPER2 = 1e-6; // alternative is 1e-6
const double SOLVER_INF = 1e10;

// compare duals
bool if_vec_equal(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    if (vec1.size() != vec2.size()) {
        return false;
    }
    // use L1 distance to measure the distance between two vectors
    double diff = 0;
    for (int index = 0; index < vec1.size(); ++index) {
        diff += std::abs(vec1[index] - vec2[index]);
        if (diff > SOLVER_PRECISION_UPPER) { // difference is significant
            //std::cout << "Debug check if two vector are equal.\n";
            //std::cout << "index: " << index;
            //std::cout << "  diff: " << diff << std::endl;
            return false;
        }
    }
    return true;
}

bool if_vec_equal(const std::vector<int>& vec1, const std::vector<int>& vec2) {
    if (vec1.size() != vec2.size()) {
        return false;
    }
    for (int index = 0; index < vec1.size(); ++index) {
        if (vec1[index] != vec2[index]) {
            return false;
        }
    }
    return true;
}

bool if_face_equal(const face& face1, const face& face2) {
    if (face1.dim == face2.dim) {
        return if_vec_equal(face1.axis, face2.axis);
    }
    else {
        return false;
    }
}

bool if_face_new(const std::vector<face>& face_collection, const face& face_candidate) {
    if (face_collection.size() < 1) { // vector "faces" have no face
        return true;
    }
    for (int idx = 0; idx < face_collection.size(); ++idx) {
        if (if_face_equal(face_collection[idx], face_candidate)) {
            return false;
        }
    }
    return true;
}

face compute_face(const std::vector<double> s) {
    face res;
    res.dim = s.size();
    for (int idx = 0; idx < s.size(); ++idx) {
        if (s[idx] >= SOLVER_PRECISION_LOWER2 && s[idx] <= SOLVER_PRECISION_UPPER2) {
            res.axis.push_back(idx);
        }
    }
    return res;
}

double compute_obj(std::vector<double> x, standardTwoStageParameters& model_parameters, std::vector<minorant>& minorant_collection) {
    double obj = 0;
    for (int idx = 0; idx < model_parameters.c.getNzeroLen(); ++idx) {
        obj += model_parameters.c.getVal(idx) * x[model_parameters.c.getLoc(idx)];
    }
    // quadratic terms
    for (int col_idx = 0; col_idx < model_parameters.Q.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.Q.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.Q.getClen(col_idx) + beg_idx; ++idx) {
            obj += 0.5 * model_parameters.Q.getVal(idx) * x[model_parameters.Q.getRow(idx)] * x[col_idx];
        }
    }
    // find the max value of minorants at x
    double max_val = 0;
    for (int idx_minorant = 0; idx_minorant < minorant_collection.size(); ++idx_minorant) {
        double cur_val = minorant_collection[idx_minorant].alpha;
        for (int idx_x = 0; idx_x < x.size(); ++idx_x) {
            cur_val += minorant_collection[idx_minorant].beta[idx_x] * x[idx_x];
        }
        if (idx_minorant > 0) {
            if (max_val < cur_val) {
                max_val = cur_val;
            }
        }
        else {
            max_val = cur_val;
        }
    }
    obj += max_val;
    return obj;
}

std::vector<double> sd_compromise_master(standardTwoStageParameters& model_parameters,
                                         const stoMap& stochastic_map,
                                         std::vector<sd_output> aggregate_outputs,
                                         double regularizer_coefficient,
                                         int num_replications,
                                         double lb_error) {
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
    // epigraphical representation of the objective
    std::vector<GRBVar> eta;
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        eta.push_back(model.addVar(aggregate_outputs[idx_rep].obj - lb_error, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
    }
    
    // quadratic objective
    GRBQuadExpr expr_obj;
    // regularizer
    for (int idx_x = 0; idx_x < model_parameters.c.getNzeroLen(); ++idx_x) {
        expr_obj += 0.5 * regularizer_coefficient * x[idx_x] * x[idx_x];
        double x_average = 0;
        // regularizer
        for (int idx_rep = 0; idx_rep < aggregate_outputs.size(); ++idx_rep) {
            x_average += aggregate_outputs[idx_rep].sol[idx_x];
        }
        x_average *= one_over_num_replications;
        expr_obj -= regularizer_coefficient * x_average * x[idx_x];
    }
    // epigraphical representation
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
        for (int idx_minorant = 0; idx_minorant < aggregate_outputs[idx_rep].minorant_collection.size(); ++idx_minorant) {
            GRBQuadExpr expr_sub;
            // + beta * x
            for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
                expr_sub += aggregate_outputs[idx_rep].minorant_collection[idx_minorant].beta[idx_x] * x[idx_x];
            }
            // + alpha
            expr_sub += aggregate_outputs[idx_rep].minorant_collection[idx_minorant].alpha;
            // c^\top x
            for (int idx_x = 0; idx_x < model_parameters.c.getNzeroLen(); ++idx_x) {
                expr_sub += model_parameters.c.getVal(idx_x) * x[idx_x];
            }
            // 0.5 x^\top Q x
            // quadratic terms
            for (int col_idx = 0; col_idx < model_parameters.Q.getColLength(); ++col_idx) {
                int beg_idx = model_parameters.Q.getCbeg(col_idx);
                for (int idx = beg_idx; idx < model_parameters.Q.getClen(col_idx) + beg_idx; ++idx) {
                    expr_sub += 0.5 * model_parameters.Q.getVal(idx) * x[model_parameters.Q.getRow(idx)] * x[col_idx];
                }
            }
            // - eta
            expr_sub -= eta[idx_rep];
            //
            model.addQConstr(expr_sub <= 0);
            
        }
    }
    
    // optimize
    model.optimize();
    
    //model.write("/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleQP/model_compromise_sd.lp");
    std::vector<double> res;
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        res.push_back(x[idx_x].get(GRB_DoubleAttr_X));
    }
    return res;
    
}

// SQLP
std::vector<double> sd_sqlp_presolve(standardTwoStageParameters& model_parameters,
                                 const stoMap& stochastic_map,
                                     stoPoint& xi) {
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
    
    // second-stage decision variable
    std::vector<GRBVar> y;
    // slack variable
    std::vector<GRBVar> y_slack;
    // initialize decision variables
    y.reserve(model_parameters.num_var_2ndStage);
    for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
        y.push_back(model.addVar(model_parameters.var_lb_2ndStage,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
    }
    if (model_parameters.num_slack_2ndStage > 0) {
        std::vector<GRBVar> y_slack_scenario;
        for (int idx_y_slack = 0; idx_y_slack < model_parameters.num_slack_2ndStage; ++idx_y_slack) {
            y_slack.push_back(model.addVar(0,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
        }
    }
    
    // objective
    GRBQuadExpr expr_obj;
    for (int idx = 0; idx < model_parameters.c.getNzeroLen(); ++idx) {
        expr_obj += model_parameters.c.getVal(idx) * x[model_parameters.c.getLoc(idx)];
    }
    for (int idx = 0; idx < model_parameters.d.getNzeroLen(); ++idx) {
        expr_obj += model_parameters.d.getVal(idx) * y[model_parameters.d.getLoc(idx)];
    }
    // quadratic terms
    for (int col_idx = 0; col_idx < model_parameters.Q.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.Q.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.Q.getClen(col_idx) + beg_idx; ++idx) {
            expr_obj += 0.5 * model_parameters.Q.getVal(idx) * x[model_parameters.Q.getRow(idx)] * x[col_idx];
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
    std::vector<GRBLinExpr> exprs_eq;
    for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
        GRBLinExpr expr;
        exprs_eq.push_back(expr);
    }
    // coefficients before y; Dy
    for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.D.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
            exprs_eq[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[col_idx];
        }
    }
    // slack variables; I_slack y_slack
    if (model_parameters.num_slack_2ndStage > 0){
        for (int col_idx = 0; col_idx < model_parameters.ISlack.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.ISlack.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.ISlack.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.ISlack.getRow(idx)] += model_parameters.ISlack.getVal(idx) * y_slack[col_idx];
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
    for (int idx_C = 0; idx_C < xi.C.size(); ++idx_C) {
        exprs_eq[stochastic_map.C_map[idx_C].first] += xi.C[idx_C] * x[stochastic_map.C_map[idx_C].second];
    }
    // right hand side e
    for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
        exprs_eq[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
    }
    // right hand side (stochastic part) equality e_(i) equality
    for (int idx = 0; idx < xi.e.size(); ++idx) {
        exprs_eq[stochastic_map.e_map[idx]] -= xi.e[idx];
    }
    // add constraints
    for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
        model.addConstr(exprs_eq[index_eq] == 0);
    }
    // optimize
    model.optimize();
    // create outputs
    std::vector<double> x_presolve;
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        x_presolve.push_back(x[idx_x].get(GRB_DoubleAttr_X));
    }
    //model.write("/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/lands_sqlp/sd_sqlp_presolve.lp");
    return x_presolve;
}

// master problem of SD for solving SQLP/SQQP problems
masterOutput sd_master(standardTwoStageParameters& model_parameters,
                                   std::vector<minorant>& minorant_collection,
                                   std::vector<double>& x_incumbent,
                                   double sigma) {
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
    GRBVar eta = model.addVar(-1e9, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
    
    // objective
    GRBQuadExpr expr_obj;
    // first stage objective c^\top x + 0.5 x^\top Q x + 0.5 * sigma \|x - x_incumbent\|^2 + eta
    for (int idx_x = 0; idx_x < model_parameters.c.getNzeroLen(); ++idx_x) {
        expr_obj += model_parameters.c.getVal(idx_x) * x[idx_x];
    }
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        // proximal mapping terms
        expr_obj += 0.5 * sigma * x[idx_x] * x[idx_x];
        expr_obj -= sigma * x_incumbent[idx_x] * x[idx_x];
    }
    // quadratic terms
    for (int col_idx = 0; col_idx < model_parameters.Q.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.Q.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.Q.getClen(col_idx) + beg_idx; ++idx) {
            expr_obj += 0.5 * model_parameters.Q.getVal(idx) * x[model_parameters.Q.getRow(idx)] * x[col_idx];
        }
    }
    expr_obj += eta;
    
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
    
    // minorants
    std::vector<GRBConstr> constraints_minorant;
    for (int idx_minorant = 0; idx_minorant < minorant_collection.size(); ++idx_minorant) {
        GRBLinExpr expr_minorant;
        expr_minorant += minorant_collection[idx_minorant].alpha - eta;
        for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
            expr_minorant += minorant_collection[idx_minorant].beta[idx_x] * x[idx_x];
        }
        constraints_minorant.push_back(model.addConstr(expr_minorant <= 0));
        //std::cout << "Debug: idx_minorant = " << idx_minorant << std::endl;
    }
    //model.write("/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/simpleQP/simpleQP_sd_master.lp");
    // optimize
    model.optimize();
    // initialize a masterOput object
    masterOutput res;
    // solution
    // Check if optimal
    int status = model.get(GRB_IntAttr_Status);
    if (status == GRB_OPTIMAL) {
        for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
            res.x.push_back(x[idx_x].get(GRB_DoubleAttr_X));
        }
        res.eta = eta.get(GRB_DoubleAttr_X);
        for (int idx_minorant = 0; idx_minorant < minorant_collection.size(); ++idx_minorant) {
            res.dual_minorant.push_back(constraints_minorant[idx_minorant].get(GRB_DoubleAttr_Pi));
        }
        res.sol_flag = 0;
    }
    else if (status == GRB_SUBOPTIMAL) {
        std::cout << "Warning: The problem is solved to sub-optimality by Gurobi." << std::endl;
        for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
            res.x.push_back(x[idx_x].get(GRB_DoubleAttr_X));
        }
        res.eta = eta.get(GRB_DoubleAttr_X);
        for (int idx_minorant = 0; idx_minorant < minorant_collection.size(); ++idx_minorant) {
            res.dual_minorant.push_back(constraints_minorant[idx_minorant].get(GRB_DoubleAttr_Pi));
        }
        res.sol_flag = 0;
    }
    else {
        std::cout << "Warning: The problem is not solved to optimality by Gurobi. Return x_incumbent and -SOLVER_INF." << std::endl;
        res.x = x_incumbent;
        res.eta = -SOLVER_INF;
        res.sol_flag = 1;
    }
    return res;
}


sd_output sd_sqlp_solver(const std::string& folder_path,
                            std::mt19937& generator,
                            stoPoint& point_est,
                            int max_iterations,
                            double f_lowerbound,
                            double sigma_init) {
    // for test
    //double sigma_lowerbound = 1;
    //double sigma_upperbound = 20;
    // timer
    std::clock_t time_start;
    time_start = std::clock();
    // current time
    std::time_t currTime = std::time(nullptr);
    // STEP 1: INITIALIZATION
    // algorithm parameters
    double sigma = sigma_init; // stepsize
    double q = 0.2;
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
    std::string resultsOutput_path = folder_path + "/computationalResults(rep_SD_SQLPv1.0).txt";
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
    // write initial setup
    std::cout << "*******************************************\n";
    writeFile << "*******************************************\n";
    std::cout << "SD-SQLP(v1.0) is initialized\n";
    writeFile << "SD-SQLP(v1.0) is initialized\n";
    std::cout << "Algorithmic Parameters\n";
    writeFile << "Algorithmic Parameters\n";
    std::cout << "sigma_init, q" << std::endl;
    writeFile << "sigma_init, q" << std::endl;
    std::cout << sigma << ", " << q << std::endl;
    writeFile << sigma << ", " << q << std::endl;
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
    
    // presolve problem to get x_incumbent
    std::vector<double> x_incumbent = sd_sqlp_presolve(model_parameters, stochastic_map, point_est);
    std::cout << "Incumbent solution after presolve:\n";
    writeFile << "Incumbent solution after presolve:\n";
    for (int idx_x = 0; idx_x < x_incumbent.size() - 1; ++idx_x) {
        std::cout << x_incumbent[idx_x] << ", ";
        writeFile << x_incumbent[idx_x] << ", ";
    }
    std::cout << x_incumbent[x_incumbent.size() - 1] << std::endl;
    writeFile << x_incumbent[x_incumbent.size() - 1] << std::endl;
    
    // initialize explored dual multipliers in the second stage
    std::vector<dualMultipliers> explored_duals;
    std::vector<double> pi_e_collection;
    std::vector<std::vector<double>> minus_pi_C_collection;
    
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
    // objective of the second-stage subproblem
    GRBLinExpr expr_obj_sub;
    for (int idx = 0; idx < model_parameters.d.getNzeroLen(); ++idx) {
        expr_obj_sub += model_parameters.d.getVal(idx) * y[model_parameters.d.getLoc(idx)];
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
    // *********************************
    // --- initialize the minorants ---
    // intermediate variables for computation
    std::vector<double> dual_workingVector(model_parameters.D.getRowLength(),0.0);
    double l1_norm = 0;
    bool flag_new = true;
    // initialize a collection of minorants
    std::vector<minorant> minorant_collection;
    // construct initial minorant
    std::cout << "Construct initial minorant.\n";
    writeFile << "Construct initial minorant.\n";
    minorant initial_minorant;
    initial_minorant.alpha = f_lowerbound; // should use lower bound for the intercept of the initial minorant
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        initial_minorant.beta.push_back(0);
    }
    minorant_collection.push_back(initial_minorant);
    // --- end initializing the minorants ---
    // --- main loop ---
    std::cout << "Maximum number of iterations: " << max_iterations << std::endl;
    writeFile << "Maximum number of iterations: " << max_iterations << std::endl;
    std::cout << "Start Solving Process\n";
    writeFile << "Start Solving Process\n";
    // intermediate values for the solver
    std::vector<double> old_candidate = x_incumbent;
    std::vector<double> old_incumbent = x_incumbent;
    minorant old_minorant_incumbent = initial_minorant;
    minorant old_minorant_candidate = initial_minorant;
    std::vector<double> x_candidate(model_parameters.num_var_1stStage, 0.0);
    std::vector<double> minus_pi_C_candidate;
    std::vector<double> minus_pi_C_incumbent;
    // sample size
    double N = 0;
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        std::cout << "***Iteration " << iteration << "***\n";
        writeFile << "***Iteration " << iteration << "***\n";
        std::cout << "sigma: " << sigma << std::endl;
        writeFile << "sigma: " << sigma << std::endl;
        std::cout << "Number of minorants used in the regularized master problem: " << minorant_collection.size() << std::endl;
        writeFile << "Number of minorants used in the regularized master problem: " << minorant_collection.size() << std::endl;
        N += 1; // increment sample size by one
        // solve the master problem
        masterOutput res_master = sd_master(model_parameters, minorant_collection, x_incumbent, sigma);
        // obtain the proximal point (condidate solution)
        for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
            x_candidate[idx_x] = res_master.x[idx_x];
            //std::cout << cplex.getValue(x_temp[idx_x]) << std::endl;
        }
        // --- compute the objective function value at the incumbent and candidate ---
        // objective value at the candidate
        double f_old_candidate = 0;
        // value at the incumbent
        double f_old_incumbent = 0;
        for (int idx = 0; idx < model_parameters.c.getNzeroLen(); ++idx) {
            f_old_candidate += model_parameters.c.getVal(idx) * x_candidate[model_parameters.c.getLoc(idx)];
            f_old_incumbent += model_parameters.c.getVal(idx) * x_incumbent[model_parameters.c.getLoc(idx)];
        }
        // quadratic terms
        for (int col_idx = 0; col_idx < model_parameters.Q.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.Q.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.Q.getClen(col_idx) + beg_idx; ++idx) {
                f_old_candidate += 0.5 * model_parameters.Q.getVal(idx) * x_candidate[model_parameters.Q.getRow(idx)] * x_candidate[col_idx];
                f_old_incumbent += 0.5 * model_parameters.Q.getVal(idx) * x_incumbent[model_parameters.Q.getRow(idx)] * x_incumbent[col_idx];
            }
        }
        double f_new_candidate = f_old_candidate;
        double f_new_incumbent = f_old_incumbent;
        f_old_candidate += res_master.eta;
        // second stage value at the incumbent
        double recourse_incumbent = f_lowerbound;
        for (int idx_minorant = 0; idx_minorant < minorant_collection.size(); ++idx_minorant) {
            double piece_val = minorant_collection[idx_minorant].alpha;
            piece_val += minorant_collection[idx_minorant].beta * x_incumbent;
            if (piece_val > recourse_incumbent) {
                recourse_incumbent = piece_val;
            }
        }
        f_old_incumbent += recourse_incumbent;
        // --- end computing the objective function value at the incumbent and candidate ---
        // --- update active minorants ---
        int dual_idx = -1;
        for (int idx_minorant = 0; idx_minorant < minorant_collection.size(); ++idx_minorant) {
            dual_idx++;
            if (res_master.dual_minorant[dual_idx] >= SOLVER_PRECISION_LOWER && res_master.dual_minorant[dual_idx] <= SOLVER_PRECISION_UPPER) { // remove inactive minorants
                minorant_collection.erase(minorant_collection.begin() + idx_minorant);
                idx_minorant--;
            }
        }
        // --- end updating active minorants ---
        // output candidate solution
        std::cout << "Candidate Solution: ";
        writeFile << "Candidate Solution: ";
        for (int x_index = 0; x_index < model_parameters.num_var_1stStage - 1; ++x_index) {
            std::cout << x_candidate[x_index] << ", ";
            writeFile << x_candidate[x_index] << ", ";
        }
        std::cout << x_candidate[model_parameters.num_var_1stStage - 1] << std::endl;
        writeFile << x_candidate[model_parameters.num_var_1stStage - 1] << std::endl;
        // --- end proximal mapping for updating candidate solutions ---
        // obtain a new data point
        int idx_scenario = iteration;
        stoPoint xi_cur = generator_random(model_parameters, stochastic_map, generator, idx_scenario, e_DB, C_DB);
        Xi.push_back(xi_cur);
        // --- second-stage subproblem ---
        // --- compute a new dual at x_candidate ---
        // update deterministic part of the rhs_bounds at current x_est
        // Dy + I_slack y_slack = [e - Cx]
        det_rhs_bounds = e_det;
        // det C x
        for (int col_idx = 0; col_idx < model_parameters.C.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.C.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.C.getClen(col_idx) + beg_idx; ++idx) {
                det_rhs_bounds[model_parameters.C.getRow(idx)] -= model_parameters.C.getVal(idx) * x_candidate[col_idx];
            }
        }
        // update the subproblem
        for (int idx_row = 0; idx_row < model_parameters.D.getRowLength(); ++idx_row) {
            rhs_bounds[idx_row] = det_rhs_bounds[idx_row];
        } // update the deterministic part
        // update the stochastic parts of e
        for (int idx_e = 0; idx_e < xi_cur.e.size(); ++idx_e) {
            rhs_bounds[stochastic_map.e_map[idx_e]] += xi_cur.e[idx_e];
        }
        // coefficients before x (stochastic part) equality (i.e., Cij * xj map: <i,j> )
        for (int idx_C = 0; idx_C < xi_cur.C.size(); ++idx_C) {
            rhs_bounds[stochastic_map.C_map[idx_C].first] -= xi_cur.C[idx_C] * x_candidate[stochastic_map.C_map[idx_C].second];
        }
        // update the RHS
        for (int idx_row = 0; idx_row < rhs_bounds.size(); ++idx_row) {
            constraintsEquality_sub[idx_row].set(GRB_DoubleAttr_RHS, rhs_bounds[idx_row]);
        }
        model_sub.optimize();
        int status = model_sub.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL) {
            // get duals
            std::vector<double> pi;
            l1_norm = 0; // reset l1_norm
            for (int idx_row = 0; idx_row < model_parameters.C.getRowLength(); ++idx_row) {
                dual_workingVector[idx_row] = constraintsEquality_sub[idx_row].get(GRB_DoubleAttr_Pi);
                l1_norm += abs(dual_workingVector[idx_row]);
            }
            // check if l1_norm is new
            flag_new = true;
            for (int idx_dual = 0; idx_dual < explored_duals.size(); ++idx_dual) {
                double diff = l1_norm - explored_duals[idx_dual].l1_norm;
                if (diff > SOLVER_PRECISION_LOWER && diff < SOLVER_PRECISION_UPPER) { // l1 norm is not new
                    flag_new = false;
                    break;
                }
            }
            if (flag_new == true) { // find a new dual extreme point
                dualMultipliers curr;
                curr.dual = dual_workingVector;
                curr.l1_norm = l1_norm;
                explored_duals.push_back(curr);
                // deterministic e
                double pi_e = model_parameters.e.fast_dotProduct(curr.dual);
                pi_e_collection.push_back(pi_e);
                // determinictic C
                std::vector<double> pi_C = model_parameters.C.fast_rightMultiply(curr.dual);
                // negate
                for (int idx = 0; idx < pi_C.size(); ++idx) {
                    pi_C[idx] = pi_C[idx] * (-1.0);
                }
                minus_pi_C_collection.push_back(pi_C);
            }
        }
        else {
            writeFile << "Error: subproblem is infeasible" << std::endl;
            std::logic_error( "Error: subproblem is infeasible" );
        }
        // erase working vector
        for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
            dual_workingVector[index_eq] = 0;
        }
        std::cout << "Number of unique duals: " << explored_duals.size() << std::endl;
        writeFile << "Number of unique duals: " << explored_duals.size() << std::endl;
        // --- end computing a new dual at x_candidate ---
        // --- construct minorants ---
        double one_over_N = 1.0 / N;
        minorant minorant_candidate;
        minorant minorant_incumbent;
        minorant_candidate.alpha = 0;
        minorant_incumbent.alpha = 0;
        for (int index_x = 0; index_x < model_parameters.num_var_1stStage; ++index_x) {
            minorant_candidate.beta.push_back(0.0);
            minorant_incumbent.beta.push_back(0.0);
        }
        minus_pi_C_candidate.clear();
        minus_pi_C_incumbent.clear();
        for (int idx_scenario = 0; idx_scenario < (int) N; ++idx_scenario) {
            // candidate
            double max_value = -SOLVER_INF;
            // incumbent
            double max_value_incumbent = -SOLVER_INF;
            int max_index = -1;
            int max_index_incumbent = -1;
            double alpha_candidate = 0;
            double alpha_incumbent = 0;
            std::vector<double> beta_candidate(model_parameters.num_var_1stStage, 0);
            std::vector<double> beta_incumbent(model_parameters.num_var_1stStage, 0);
            for (int dual_index = 0; dual_index < explored_duals.size(); ++dual_index) {
                // pi_e - pi_C
                // find optimal dual based on the given set of unique duals
                double current_value = 0;
                // deterministic e
                double pi_e = pi_e_collection[dual_index];
                // stochastic e
                if (flag_new == false) { // dual set is not changed
                    if (idx_scenario == iteration) { // new realization
                        double sto_pi_e = 0;
                        for (int idx_eq = 0; idx_eq < Xi[idx_scenario].e.size(); ++idx_eq) {
                            sto_pi_e += explored_duals[dual_index].dual[stochastic_map.e_map[idx_eq]] * Xi[idx_scenario].e[idx_eq];
                        }
                        explored_duals[dual_index].sto_pi_e.push_back(sto_pi_e);
                    }
                }
                else {
                    if (dual_index == explored_duals.size() - 1){ // new dual
                        double sto_pi_e = 0;
                        for (int idx_eq = 0; idx_eq < Xi[idx_scenario].e.size(); ++idx_eq) {
                            sto_pi_e += explored_duals[dual_index].dual[stochastic_map.e_map[idx_eq]] * Xi[idx_scenario].e[idx_eq];
                        }
                        explored_duals[dual_index].sto_pi_e.push_back(sto_pi_e);
                    }
                    else if (idx_scenario == iteration) {
                        double sto_pi_e = 0;
                        for (int idx_eq = 0; idx_eq < Xi[idx_scenario].e.size(); ++idx_eq) {
                            sto_pi_e += explored_duals[dual_index].dual[stochastic_map.e_map[idx_eq]] * Xi[idx_scenario].e[idx_eq];
                        }
                        explored_duals[dual_index].sto_pi_e.push_back(sto_pi_e);
                    }
                }
                pi_e += explored_duals[dual_index].sto_pi_e[idx_scenario];
                //std::cout << "Debug pi_e: " << pi_e << std::endl;
                current_value += pi_e;
                std::vector<double> sto_minus_pi_C(model_parameters.C.getColLength(), 0.0);
                // stochastic C
                // equality
                for (int idx_C = 0; idx_C < Xi[idx_scenario].C.size(); ++idx_C) {
                    sto_minus_pi_C[stochastic_map.C_map[idx_C].second] -= Xi[idx_scenario].C[idx_C] * explored_duals[dual_index].dual[stochastic_map.C_map[idx_C].first];
                }
                // incumbent
                double current_value_incumbent = pi_e;
                // deterministic part is only calculated at the first time
                /*
                if (idx_scenario < 1) {
                    minus_pi_C_candidate.push_back(minus_pi_C_collection[dual_index] * x_candidate);
                    if (iteration < 1) { // first iteration
                        minus_pi_C_incumbent.push_back(minus_pi_C_collection[dual_index] * x_incumbent);
                    }
                    else {
                        if (flag_new == true && dual_index == explored_duals.size() - 1) {
                            minus_pi_C_incumbent.push_back(minus_pi_C_collection[dual_index] * x_incumbent);
                        }
                    }
                }
                */
                if (idx_scenario < 1) {
                    minus_pi_C_candidate.push_back(minus_pi_C_collection[dual_index] * x_candidate);
                    minus_pi_C_incumbent.push_back(minus_pi_C_collection[dual_index] * x_incumbent);
                }
                current_value += minus_pi_C_candidate[dual_index];
                current_value_incumbent += minus_pi_C_incumbent[dual_index];
                if (Xi[idx_scenario].C.size() > 0) {
                    current_value += sto_minus_pi_C * x_candidate;
                    current_value_incumbent += sto_minus_pi_C * x_incumbent;
                }
                std::cout << "max_value: " << max_value << std::endl;
                std::cout << "current_value: " << current_value << std::endl;
                // update the maximum
                if (dual_index < 1) {
                    max_index = dual_index;
                    max_value = current_value;
                    max_index_incumbent = dual_index;
                    max_value_incumbent = current_value_incumbent;
                    // store the intercept and slope
                    alpha_candidate = pi_e; // \pi^\top e
                    alpha_incumbent = pi_e;
                    //beta_candidate = (-1.0) * pi_C; // -\pi^\top C
                    if (Xi[idx_scenario].C.size() > 0) {
                        beta_candidate = sto_minus_pi_C;
                        beta_incumbent = sto_minus_pi_C;
                    }
                } // if (dual_index < 1)
                else {
                    if (max_value < current_value) { // find the better dual for given candidate
                        max_index = dual_index;
                        max_value = current_value;
                        alpha_candidate = pi_e;
                        //beta_candidate = (-1.0) * pi_C;
                        if (Xi[idx_scenario].C.size() > 0) {
                            beta_candidate = sto_minus_pi_C;
                        }
                    }
                    if (max_value_incumbent < current_value_incumbent) { // find the better dual for given incumbent
                        max_index_incumbent = dual_index;
                        max_value_incumbent = current_value_incumbent;
                        alpha_incumbent = pi_e;
                        //beta_incumbent = (-1.0) * pi_C;
                        // only store the stochastic part, det part is stored as dual_index
                        if (Xi[idx_scenario].C.size() > 0) {
                            beta_incumbent = sto_minus_pi_C;
                        }
                    }
                }
            } // end for (int dual_index = 0; dual_index < explored_duals.size(); ++dual_index)
            // minorant on the candidate
            minorant_candidate.alpha += alpha_candidate;
            // minorant on the incumbent
            minorant_incumbent.alpha += alpha_incumbent;
            if (Xi[idx_scenario].C.size() > 0) {
                minorant_candidate.beta = minorant_candidate.beta + beta_candidate;
                minorant_candidate.beta = minorant_candidate.beta + minus_pi_C_collection[max_index];
                minorant_incumbent.beta = minorant_incumbent.beta + beta_incumbent;
                minorant_incumbent.beta = minorant_incumbent.beta + minus_pi_C_collection[max_index_incumbent];
            }
            else {
                minorant_candidate.beta = minorant_candidate.beta + minus_pi_C_collection[max_index];
                minorant_incumbent.beta = minorant_incumbent.beta + minus_pi_C_collection[max_index_incumbent];
            }
        }
        minorant_candidate.alpha *= one_over_N;
        minorant_candidate.beta = one_over_N * minorant_candidate.beta;
        minorant_incumbent.alpha *= one_over_N;
        minorant_incumbent.beta = one_over_N * minorant_incumbent.beta;
        // --- end constructing minorants ---
        for (int dual_idx = 0; dual_idx < explored_duals.size(); ++dual_idx) {
            std::cout << "dual[" << dual_idx << "](sto_pi_e): ";
            for (int idx_scenario = 0; idx_scenario < iteration + 1; ++idx_scenario) {
                std::cout << explored_duals[dual_idx].sto_pi_e[idx_scenario] << " ";
            }
            std::cout << std::endl;
        }
        // --- update minorant ---
        double N_over_Nplus1 = N / (N + 1.0);
        double one_over_Nplus1 = 1.0 / (N + 1.0);
        for (int idx_minorant = 0; idx_minorant < minorant_collection.size(); ++idx_minorant) {
            minorant_collection[idx_minorant].alpha = minorant_collection[idx_minorant].alpha * N_over_Nplus1 + one_over_Nplus1 * f_lowerbound;
        } // end for (int idx_minorant = 0; idx_minorant < minorant_collection.size(); ++idx_minorant)
        // --- end updating minorant ---
        // add new minorants
        minorant_collection.push_back(minorant_candidate);
        minorant_collection.push_back(minorant_incumbent);
        // report new minorants
        std::cout << "Minorant Candidate\n";
        writeFile << "Minorant Candidate\n";
        std::cout << "alpha: " << minorant_candidate.alpha << std::endl;
        writeFile << "alpha: " << minorant_candidate.alpha << std::endl;
        std::cout << "beta: ";
        writeFile << "beta: ";
        for (int x_index = 0; x_index < model_parameters.num_var_1stStage - 1; ++x_index) {
            std::cout << minorant_candidate.beta[x_index] << ", ";
            writeFile << minorant_candidate.beta[x_index] << ", ";
        }
        std::cout << minorant_candidate.beta[model_parameters.num_var_1stStage - 1] << std::endl;
        writeFile << minorant_candidate.beta[model_parameters.num_var_1stStage - 1] << std::endl;
        std::cout << "Minorant Incumbent\n";
        writeFile << "Minorant Incumbent\n";
        std::cout << "alpha: " << minorant_incumbent.alpha << std::endl;
        writeFile << "alpha: " << minorant_incumbent.alpha << std::endl;
        std::cout << "beta: ";
        writeFile << "beta: ";
        for (int x_index = 0; x_index < model_parameters.num_var_1stStage - 1; ++x_index) {
            std::cout << minorant_incumbent.beta[x_index] << ", ";
            writeFile << minorant_incumbent.beta[x_index] << ", ";
        }
        std::cout << minorant_incumbent.beta[model_parameters.num_var_1stStage - 1] << std::endl;
        writeFile << minorant_incumbent.beta[model_parameters.num_var_1stStage - 1] << std::endl;
        // --- incumbent selection ---
        // LHS
        // new function value at candidate solution
        // second stage value
        double new_recourse_candidate = f_lowerbound;
        double new_recourse_incumbent = f_lowerbound;
        for (int idx_minorant = 0; idx_minorant < minorant_collection.size(); ++idx_minorant) {
            double piece_val = minorant_collection[idx_minorant].alpha;
            piece_val += minorant_collection[idx_minorant].beta * x_candidate;
            if (piece_val > new_recourse_candidate) {
                new_recourse_candidate = piece_val;
            }
            double piece_val2 = minorant_collection[idx_minorant].alpha;
            piece_val2 += minorant_collection[idx_minorant].beta * x_incumbent;
            if (piece_val2 > new_recourse_incumbent) {
                new_recourse_incumbent = piece_val2;
            }
        }
        f_new_candidate += new_recourse_candidate;
        f_new_incumbent += new_recourse_incumbent;
        double LHS = f_new_candidate - f_new_incumbent;
        double RHS = q * (f_old_candidate - f_old_incumbent);
        std::cout << "LHS: " << LHS << std::endl;
        writeFile << "LHS: " << LHS << std::endl;
        std::cout << "RHS: " << RHS << std::endl;
        writeFile << "RHS: " << RHS << std::endl;
        if (LHS <= RHS) {
            std::cout << "Computation Log: Incumbent selection criterion is passed.\n";
            writeFile <<"Computation Log: Incumbent selection criterion is passed.\n";
            x_incumbent = x_candidate;
            // update stepsize
            sigma = ((double) iteration + 1.0) * sigma_init;
            //sigma = max(sigma * 0.5, sigma_lowerbound);
            // update minus_pi_C_inncumbent
            minus_pi_C_incumbent = minus_pi_C_candidate;
        }
        else {
            std::cout << "Computation Log: Incumbent selection criterion is not passed.\n";
            writeFile <<"Computation Log: Incumbent solution selection criterion is not passed.\n";
            //sigma = min(sigma * 2.0, sigma_upperbound);
        }
        // print out the incumbent solution
        std::cout << "Incumbent Solution: ";
        writeFile << "Incumbent Solution: ";
        for (int index = 0; index < model_parameters.num_var_1stStage -1; ++index) {
            std::cout << x_incumbent[index] << ", ";
            writeFile << x_incumbent[index] << ", ";
        }
        std::cout << x_incumbent[model_parameters.num_var_1stStage - 1] << std::endl;
        writeFile << x_incumbent[model_parameters.num_var_1stStage - 1] << std::endl;
        //
        std::cout << "sigma: " << sigma << std::endl;
        // update candidates and incuments
        old_candidate = x_candidate;
        old_incumbent = x_incumbent;
    } // end for (int iteration = 0; iteration < max_iterations; ++iteration)
    // --- end main loop ---
    std::cout << "*******************************************\n";
    writeFile << "*******************************************\n";
    std::cout << "Output Solution: ";
    writeFile << "Output Solution: ";
    for (int index = 0; index < model_parameters.num_var_1stStage-1; ++index) {
        std::cout << x_incumbent[index] << ", ";
        writeFile << x_incumbent[index] << ", ";
    }
    std::cout << x_incumbent[model_parameters.num_var_1stStage - 1] << std::endl;
    writeFile << x_incumbent[model_parameters.num_var_1stStage - 1] << std::endl;
    std::cout << "Computation Log: Finish Solving Process.\n";
    writeFile << "Computation Log: Finish Solving Process.\n";
    // write time elapsed
    double duration = (std::clock() - time_start ) / (double) CLOCKS_PER_SEC;
    writeFile << "Time elapsed(secs) : " << duration << "\n";
    writeFile << "*******************************************\n";
    writeFile.close();
    // return the results
    sd_output res;
    res.sol = x_incumbent;
    //res.Xi = Xi;
    res.it_num = max_iterations;
    res.sol_time = duration;
    res.obj = compute_obj(x_incumbent, model_parameters, minorant_collection);
    return res;
}


// compromise decision solver where SD is used in the replication step

compromise_sd_output compromise_sd_sqlp_solver(const std::string& folder_path,
                                               std::mt19937& generator,
                                               stoPoint& point_est,
                                               int max_iterations,
                                               double f_lowerbound,
                                               double sigma_init,
                                               double regularizer_coefficient,
                                               int num_replications,
                                               double lb_error) {
    // Initialization
    std::string model_path = folder_path + "/model.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(CD_SD_SQLPv1.0).txt";
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
    std::cout << "Compromise Decision for SD for Solving SQLP Problems (v1.0)\n";
    writeFile << "Compromise Decision for SD for Solving SQLP Problems (v1.0)\n";
    std::cout << "Parameters\n";
    writeFile << "Parameters\n";
    std::cout << "num_replications: " << num_replications << std::endl;
    writeFile << "num_replications: " << num_replications << std::endl;
    std::cout << "max_iterations: " << max_iterations << std::endl;
    writeFile << "max_iterations: " << max_iterations << std::endl;
    std::cout << "sigma_init: " << sigma_init << std::endl;
    writeFile << "sigma_init: " << sigma_init << std::endl;
    std::cout << "Lower bound error: " << lb_error << std::endl;
    writeFile << "Lower bound error: " << lb_error << std::endl;
    
    compromise_sd_output res;
    // Step 1: Replication
    std::vector<sd_output> aggregate_outputs;
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        sd_output replication_outputs = sd_sqlp_solver(folder_path, generator, point_est, max_iterations, f_lowerbound, sigma_init);
        aggregate_outputs.push_back(replication_outputs);
        res.replication_decision.push_back(replication_outputs.sol);
    }
    
    // Step 2: Aggregation
    std::vector<double> compromise_decision = sd_compromise_master(model_parameters, stochastic_map, aggregate_outputs, regularizer_coefficient, num_replications, lb_error);
    
    
    res.compromise_decision = compromise_decision;
    
    return res;
}


// --- SQQP ---
std::vector<double> sd_sqqp_presolve(standardTwoStageParameters& model_parameters,
                                     const stoMap& stochastic_map,
                                     stoPoint& xi) {
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
    
    // second-stage decision variable
    std::vector<GRBVar> y;
    // slack variable
    std::vector<GRBVar> y_slack;
    // initialize decision variables
    y.reserve(model_parameters.num_var_2ndStage);
    for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
        y.push_back(model.addVar(model_parameters.var_lb_2ndStage,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
    }
    if (model_parameters.num_slack_2ndStage > 0) {
        std::vector<GRBVar> y_slack_scenario;
        for (int idx_y_slack = 0; idx_y_slack < model_parameters.num_slack_2ndStage; ++idx_y_slack) {
            y_slack.push_back(model.addVar(0,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
        }
    }
    
    // objective
    GRBQuadExpr expr_obj;
    // linear terms
    for (int idx = 0; idx < model_parameters.c.getNzeroLen(); ++idx) {
        expr_obj += model_parameters.c.getVal(idx) * x[model_parameters.c.getLoc(idx)];
    }
    for (int idx = 0; idx < model_parameters.d.getNzeroLen(); ++idx) {
        expr_obj += model_parameters.d.getVal(idx) * y[model_parameters.d.getLoc(idx)];
    }
    // quadratic terms
    for (int col_idx = 0; col_idx < model_parameters.Q.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.Q.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.Q.getClen(col_idx) + beg_idx; ++idx) {
            expr_obj += 0.5 * model_parameters.Q.getVal(idx) * x[model_parameters.Q.getRow(idx)] * x[col_idx];
        }
    }
    
    for (int col_idx = 0; col_idx < model_parameters.P.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.P.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.P.getClen(col_idx) + beg_idx; ++idx) {
            if (model_parameters.P.getRow(idx) < model_parameters.num_var_2ndStage) {
                if (col_idx < model_parameters.num_var_2ndStage) {
                    // case 1: P_ij * y_i * y_j
                    expr_obj += 0.5 * model_parameters.P.getVal(idx) * y[model_parameters.P.getRow(idx)] * y[col_idx];
                }
                else {
                    // case 2: P_ij * y_i * yslack_j
                    expr_obj += 0.5 * model_parameters.P.getVal(idx) * y[model_parameters.P.getRow(idx)] * y_slack[col_idx - model_parameters.num_var_2ndStage];
                }
            }
            else {
                if (col_idx < model_parameters.num_var_2ndStage) {
                    // case 3: P_ij * yslack_i * y_j
                    expr_obj += 0.5 * model_parameters.P.getVal(idx) * y_slack[model_parameters.P.getRow(idx) - model_parameters.num_var_2ndStage] * y[col_idx];
                }
                else {
                    // case 4: P_ij * yslack_i * yslack_j
                    expr_obj += 0.5 * model_parameters.P.getVal(idx) * y_slack[model_parameters.P.getRow(idx) - model_parameters.num_var_2ndStage] * y_slack[col_idx - model_parameters.num_var_2ndStage];
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
    std::vector<GRBLinExpr> exprs_eq;
    for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
        GRBLinExpr expr;
        exprs_eq.push_back(expr);
    }
    // coefficients before y; Dy
    for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.D.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
            exprs_eq[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[col_idx];
        }
    }
    // slack variables; I_slack y_slack
    if (model_parameters.num_slack_2ndStage > 0){
        for (int col_idx = 0; col_idx < model_parameters.ISlack.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.ISlack.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.ISlack.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.ISlack.getRow(idx)] += model_parameters.ISlack.getVal(idx) * y_slack[col_idx];
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
    for (int idx_C = 0; idx_C < xi.C.size(); ++idx_C) {
        exprs_eq[stochastic_map.C_map[idx_C].first] += xi.C[idx_C] * x[stochastic_map.C_map[idx_C].second];
    }
    // right hand side e
    for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
        exprs_eq[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
    }
    // right hand side (stochastic part) equality e_(i) equality
    for (int idx = 0; idx < xi.e.size(); ++idx) {
        exprs_eq[stochastic_map.e_map[idx]] -= xi.e[idx];
    }
    // add constraints
    for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
        model.addConstr(exprs_eq[index_eq] == 0);
    }
    // optimize
    model.optimize();
    // create outputs
    std::vector<double> x_presolve;
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        x_presolve.push_back(x[idx_x].get(GRB_DoubleAttr_X));
    }
    //model.write("/Users/shuotaodiao/Documents/Research/numerical_experiments/CD/lands_sqlp/sd_sqqp_presolve.lp");
    return x_presolve;
}

dualMultipliersQP sd_sqqp_2ndStagePrimal(standardTwoStageParameters& model_parameters,
                                         const stoMap& stochastic_map,
                                         stoPoint& xi,
                                         const std::vector<double>& x) {
    // --- gurobi solver environment ---
    // Create the Gurobi environment
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);  // Turn off all console logging
    // Create an empty model
    GRBModel model = GRBModel(env);
    // second-stage decision variable
    std::vector<GRBVar> y;
    // slack variable
    std::vector<GRBVar> y_slack;
    // initialize decision variables
    y.reserve(model_parameters.num_var_2ndStage);
    for (int idx_y = 0; idx_y < model_parameters.num_var_2ndStage; ++idx_y) {
        y.push_back(model.addVar(model_parameters.var_lb_2ndStage,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
    }
    if (model_parameters.num_slack_2ndStage > 0) {
        std::vector<GRBVar> y_slack_scenario;
        for (int idx_y_slack = 0; idx_y_slack < model_parameters.num_slack_2ndStage; ++idx_y_slack) {
            y_slack.push_back(model.addVar(0,GRB_INFINITY, 0.0, GRB_CONTINUOUS));
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
    model.setObjective(expr_obj_sub, GRB_MINIMIZE);
    
    // constraints of the second-stage subproblem
    // stndard form equality constraints D y + Islack y_slack + Cx = e
    std::vector<GRBLinExpr> exprs_eq;
    for (int index_eq = 0; index_eq < model_parameters.D.getRowLength(); ++index_eq) {
        GRBLinExpr expr;
        exprs_eq.push_back(expr);
    }
    // coefficients before y; Dy
    for (int col_idx = 0; col_idx < model_parameters.D.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.D.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.D.getClen(col_idx) + beg_idx; ++idx) {
            exprs_eq[model_parameters.D.getRow(idx)] += model_parameters.D.getVal(idx) * y[col_idx];
        }
    }
    // slack variables; I_slack y_slack
    if (model_parameters.num_slack_2ndStage > 0){
        for (int col_idx = 0; col_idx < model_parameters.ISlack.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.ISlack.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.ISlack.getClen(col_idx) + beg_idx; ++idx) {
                exprs_eq[model_parameters.ISlack.getRow(idx)] += model_parameters.ISlack.getVal(idx) * y_slack[col_idx];
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
    for (int idx_C = 0; idx_C < xi.C.size(); ++idx_C) {
        exprs_eq[stochastic_map.C_map[idx_C].first] += xi.C[idx_C] * x[stochastic_map.C_map[idx_C].second];
    }
    // right hand side e
    for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
        exprs_eq[model_parameters.e.getLoc(idx)] -= model_parameters.e.getVal(idx);
    }
    // right hand side (stochastic part) equality e_(i) equality
    for (int idx = 0; idx < xi.e.size(); ++idx) {
        exprs_eq[stochastic_map.e_map[idx]] -= xi.e[idx];
    }
    // add the equality constraints
    std::vector<GRBConstr> constraintsEquality;
    for (int idx = 0; idx< model_parameters.D.getRowLength(); ++idx) {
        constraintsEquality.push_back(model.addConstr(exprs_eq[idx] == 0));
    }
    
    model.optimize();
    dualMultipliersQP res;
    for (int idx_row = 0; idx_row < model_parameters.C.getRowLength(); ++idx_row) {
        res.t.push_back(constraintsEquality[idx_row].get(GRB_DoubleAttr_Pi));
    }
    res.obj = model.get(GRB_DoubleAttr_ObjVal);
    return res;
}

// dual of the second stage problem, refer to (9) of the StQP paper
dualMultipliersQP sd_sqqp_2ndStageDual(standardTwoStageParameters& model_parameters,
                                       const stoMap& stochastic_map,
                                       stoPoint& xi,
                                       const std::vector<double>& x) {
    // --- gurobi solver environment ---
    // Create the Gurobi environment
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);  // Turn off all console logging
    // Create an empty model
    GRBModel model = GRBModel(env);
    
    // decision variables
    std::vector<GRBVar> s;
    long s_len = model_parameters.num_var_2ndStage + model_parameters.num_slack_2ndStage;
    s.reserve(s_len); // Pre-allocate memory
    for (int idx_s = 0; idx_s < s_len; ++idx_s) {
        s.push_back(model.addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
    }
    std::vector<GRBVar> t;
    t.reserve(model_parameters.D.getRowLength());
    for (int idx_t = 0; idx_t < model_parameters.D.getRowLength(); ++idx_t) {
        t.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
    }
    
    // objective function g_QQ
    GRBQuadExpr g_QQ;
    
    // - 0.5 d^\top Pinv d
    g_QQ -= 0.5 * model_parameters.dt_Pinv_d;
    
    // + [d^\top P^{-1} D^\top] t
    for (int idx = 0; idx < model_parameters.dt_Pinv_Dt.getNzeroLen(); ++idx) {
        g_QQ += model_parameters.dt_Pinv_Dt.getVal(idx) * t[model_parameters.dt_Pinv_Dt.getLoc(idx)];
    }
    // + [d^\top P^{-1}] s
    for (int idx = 0; idx < model_parameters.dt_Pinv.getNzeroLen(); ++idx) {
        g_QQ += model_parameters.dt_Pinv.getVal(idx) * s[model_parameters.dt_Pinv.getLoc(idx)];
    }
    // - 0.5 * t [D P^{-1} D^\top] t
    for (int col_idx = 0; col_idx < model_parameters.D_Pinv_Dt.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.D_Pinv_Dt.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.D_Pinv_Dt.getClen(col_idx) + beg_idx; ++idx) {
            g_QQ -= 0.5 * model_parameters.D_Pinv_Dt.getVal(idx) * t[model_parameters.D_Pinv_Dt.getRow(idx)] * t[col_idx];
        }
    }
    // - t^\top [D^\top P^{-1}] s
    for (int col_idx = 0; col_idx < model_parameters.D_Pinv.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.D_Pinv.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.D_Pinv.getClen(col_idx) + beg_idx; ++idx) {
            g_QQ -= model_parameters.D_Pinv.getVal(idx) * t[model_parameters.D_Pinv.getRow(idx)] * s[col_idx];
        }
    }
    // - 0.5 * s^\top [P^{-1}] s
    for (int col_idx = 0; col_idx < model_parameters.Pinv.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.Pinv.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.Pinv.getClen(col_idx) + beg_idx; ++idx) {
            g_QQ -= 0.5 * model_parameters.Pinv.getVal(idx) * s[model_parameters.Pinv.getRow(idx)] * s[col_idx];
        }
    }
    // --- (e - C x) t---
    // terms of x
    long num_row_2ndStage = model_parameters.num_E_2ndStage + model_parameters.num_G_2ndStage + model_parameters.num_L_2ndStage;
    std::vector<double> e(num_row_2ndStage, 0.0);
    for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
        e[model_parameters.e.getLoc(idx)] += model_parameters.e.getVal(idx);
    }
    // stochastic part
    for (int idx = 0; idx < xi.e.size(); ++idx) {
        e[stochastic_map.e_map[idx]] += xi.e[idx];
    }
    // deterministic part
    std::vector<double> Cx = model_parameters.C.fast_leftMultiply(x);
    // stochastic part
    for (int idx = 0; idx < xi.C.size(); ++idx) {
        Cx[stochastic_map.C_map[idx].first] = xi.C[idx] * x[stochastic_map.C_map[idx].second];
    }
    //
    for (int idx = 0; idx < num_row_2ndStage; ++idx) {
        g_QQ += e[idx] * t[idx];
        g_QQ -= Cx[idx] * t[idx];
    }
    // --- end (e - C x) t---
    // add objective
    model.setObjective(g_QQ, GRB_MAXIMIZE);
    // optimize
    model.optimize();
    dualMultipliersQP res;
    for (int idx = 0; idx < s_len; ++idx) {
        res.s.push_back(s[idx].get(GRB_DoubleAttr_X));
    }
    for (int idx = 0; idx < t.size(); ++idx) {
        res.t.push_back(t[idx].get(GRB_DoubleAttr_X));
    }
    res.obj = model.get(GRB_DoubleAttr_ObjVal);
    return res;
}


dualMultipliersQP sd_sqqp_2ndStageDual(standardTwoStageParameters& model_parameters,
                                       const stoMap& stochastic_map,
                                       stoPoint& xi,
                                       const std::vector<double>& x,
                                       face curFace) {
    // --- gurobi solver environment ---
    // Create the Gurobi environment
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);  // Turn off all console logging
    // Create an empty model
    GRBModel model = GRBModel(env);
    
    // set the precision
    //model.set(GRB_DoubleParam_OptimalityTol, 1e-3);
    // decision variables
    std::vector<GRBVar> s;
    long s_len = model_parameters.num_var_2ndStage + model_parameters.num_slack_2ndStage;
    s.reserve(s_len); // Pre-allocate memory
    int idx_axis = 0;
    for (int idx_s = 0; idx_s < s_len; ++idx_s) {
        
        if (idx_axis < curFace.axis.size()) {
            if (idx_s == curFace.axis[idx_axis]) {
                s.push_back(model.addVar(0, 0, 0.0, GRB_CONTINUOUS));
                idx_axis += 1;
            }
            else {
                s.push_back(model.addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
            }
        }
        else {
            s.push_back(model.addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
        }
         
        //s.push_back(model.addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
    }
    std::vector<GRBVar> t;
    t.reserve(model_parameters.D.getRowLength());
    for (int idx_t = 0; idx_t < model_parameters.D.getRowLength(); ++idx_t) {
        t.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
    }
    
    // objective function g_QQ
    GRBQuadExpr g_QQ;
    
    // - 0.5 d^\top Pinv d
    g_QQ -= 0.5 * model_parameters.dt_Pinv_d;
    
    // + [d^\top P^{-1} D^\top] t
    for (int idx = 0; idx < model_parameters.dt_Pinv_Dt.getNzeroLen(); ++idx) {
        g_QQ += model_parameters.dt_Pinv_Dt.getVal(idx) * t[model_parameters.dt_Pinv_Dt.getLoc(idx)];
    }
    // + [d^\top P^{-1}] s
    for (int idx = 0; idx < model_parameters.dt_Pinv.getNzeroLen(); ++idx) {
        g_QQ += model_parameters.dt_Pinv.getVal(idx) * s[model_parameters.dt_Pinv.getLoc(idx)];
    }
    // - 0.5 * t [D P^{-1} D^\top] t
    for (int col_idx = 0; col_idx < model_parameters.D_Pinv_Dt.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.D_Pinv_Dt.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.D_Pinv_Dt.getClen(col_idx) + beg_idx; ++idx) {
            g_QQ -= 0.5 * model_parameters.D_Pinv_Dt.getVal(idx) * t[model_parameters.D_Pinv_Dt.getRow(idx)] * t[col_idx];
        }
    }
    // - t^\top [D P^{-1}] s
    for (int col_idx = 0; col_idx < model_parameters.D_Pinv.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.D_Pinv.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.D_Pinv.getClen(col_idx) + beg_idx; ++idx) {
            g_QQ -= model_parameters.D_Pinv.getVal(idx) * t[model_parameters.D_Pinv.getRow(idx)] * s[col_idx];
        }
    }
    // - 0.5 * s^\top [P^{-1}] s
    for (int col_idx = 0; col_idx < model_parameters.Pinv.getColLength(); ++col_idx) {
        int beg_idx = model_parameters.Pinv.getCbeg(col_idx);
        for (int idx = beg_idx; idx < model_parameters.Pinv.getClen(col_idx) + beg_idx; ++idx) {
            g_QQ -= 0.5 * model_parameters.Pinv.getVal(idx) * s[model_parameters.Pinv.getRow(idx)] * s[col_idx];
        }
    }
    // --- (e - C x) t---
    // terms of x
    long num_row_2ndStage = model_parameters.num_E_2ndStage + model_parameters.num_G_2ndStage + model_parameters.num_L_2ndStage;
    std::vector<double> e(num_row_2ndStage, 0.0);
    for (int idx = 0; idx < model_parameters.e.getNzeroLen(); ++idx) {
        e[model_parameters.e.getLoc(idx)] += model_parameters.e.getVal(idx);
    }
    // stochastic part
    for (int idx = 0; idx < xi.e.size(); ++idx) {
        e[stochastic_map.e_map[idx]] += xi.e[idx];
    }
    // deterministic part
    std::vector<double> Cx = model_parameters.C.fast_leftMultiply(x);
    // stochastic part
    for (int idx = 0; idx < xi.C.size(); ++idx) {
        Cx[stochastic_map.C_map[idx].first] = xi.C[idx] * x[stochastic_map.C_map[idx].second];
    }
    //
    for (int idx = 0; idx < num_row_2ndStage; ++idx) {
        g_QQ += e[idx] * t[idx];
        g_QQ -= Cx[idx] * t[idx];
    }
    // --- end (e - C x) t---
    // add objective
    model.setObjective(g_QQ, GRB_MAXIMIZE);
    // fix faces
    /*
    for (int idx = 0; idx < curFace.axis.size(); ++idx) {
        model.addConstr(s[curFace.axis[idx]] == 0);
    }
     */
    // optimize
    model.optimize();
    dualMultipliersQP res;
    for (int idx = 0; idx < s_len; ++idx) {
        res.s.push_back(s[idx].get(GRB_DoubleAttr_X));
    }
    for (int idx = 0; idx < t.size(); ++idx) {
        res.t.push_back(t[idx].get(GRB_DoubleAttr_X));
    }
    res.obj = model.get(GRB_DoubleAttr_ObjVal);
    return res;
}


sd_output sd_sqqp_solver(const std::string& folder_path,
                            std::mt19937& generator,
                            stoPoint& point_est,
                            int max_iterations,
                            double f_lowerbound,
                            double sigma_init) {
    // timer
    std::clock_t time_start;
    time_start = std::clock();
    // current time
    std::time_t currTime = std::time(nullptr);
    // STEP 1: INITIALIZATION
    // algorithm parameters
    double sigma = sigma_init; // stepsize
    double q = 0.2;
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
    std::string resultsOutput_path = folder_path + "/computationalResults(rep_SD_SQQPv1.0).txt";
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
    // write initial setup
    std::cout << "*******************************************\n";
    writeFile << "*******************************************\n";
    std::cout << "SD-SQQP(v1.0) is initialized\n";
    writeFile << "SD-SQQP(v1.0) is initialized\n";
    std::cout << "Algorithmic Parameters\n";
    writeFile << "Algorithmic Parameters\n";
    std::cout << "sigma_init, q" << std::endl;
    writeFile << "sigma_init, q" << std::endl;
    std::cout << sigma << ", " << q << std::endl;
    writeFile << sigma << ", " << q << std::endl;
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
    // presolve problem to get x_incumbent
    std::vector<double> x_incumbent = sd_sqqp_presolve(model_parameters, stochastic_map, point_est);
    std::cout << "Incumbent solution after presolve:\n";
    writeFile << "Incumbent solution after presolve:\n";
    for (int idx_x = 0; idx_x < x_incumbent.size() - 1; ++idx_x) {
        std::cout << x_incumbent[idx_x] << ", ";
        writeFile << x_incumbent[idx_x] << ", ";
    }
    std::cout << x_incumbent[x_incumbent.size() - 1] << std::endl;
    writeFile << x_incumbent[x_incumbent.size() - 1] << std::endl;
    // initialize set of explored faces
    std::vector<face> explored_faces;
    // initialize a collection of minorants
    std::vector<minorant> minorant_collection;
    // construct initial minorant
    std::cout << "Construct initial minorant.\n";
    writeFile << "Construct initial minorant.\n";
    minorant initial_minorant;
    initial_minorant.alpha = f_lowerbound; // should use lower bound for the intercept of the initial minorant
    for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
        initial_minorant.beta.push_back(0);
    }
    minorant_collection.push_back(initial_minorant);
    // --- end initializing the minorants ---
    // --- main loop ---
    std::cout << "Maximum number of iterations: " << max_iterations << std::endl;
    writeFile << "Maximum number of iterations: " << max_iterations << std::endl;
    std::cout << "Start Solving Process\n";
    writeFile << "Start Solving Process\n";
    // intermediate values for the solver
    std::vector<double> old_candidate = x_incumbent;
    std::vector<double> old_incumbent = x_incumbent;
    minorant old_minorant_incumbent = initial_minorant;
    minorant old_minorant_candidate = initial_minorant;
    std::vector<double> x_candidate(model_parameters.num_var_1stStage, 0.0);
    // sample size
    double N = 0;
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        std::cout << "***Iteration " << iteration << "***\n";
        writeFile << "***Iteration " << iteration << "***\n";
        std::cout << "sigma: " << sigma << std::endl;
        writeFile << "sigma: " << sigma << std::endl;
        std::cout << "Number of minorants used in the regularized master problem: " << minorant_collection.size() << std::endl;
        writeFile << "Number of minorants used in the regularized master problem: " << minorant_collection.size() << std::endl;
        N += 1; // increment sample size by one
        // solve the master problem
        masterOutput res_master = sd_master(model_parameters, minorant_collection, x_incumbent, sigma);
        // obtain the proximal point (condidate solution)
        for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
            x_candidate[idx_x] = res_master.x[idx_x];
            //std::cout << cplex.getValue(x_temp[idx_x]) << std::endl;
        }
        // --- compute the objective function value at the incumbent and candidate ---
        // objective value at the candidate
        double f_old_candidate = 0;
        // value at the incumbent
        double f_old_incumbent = 0;
        for (int idx = 0; idx < model_parameters.c.getNzeroLen(); ++idx) {
            f_old_candidate += model_parameters.c.getVal(idx) * x_candidate[model_parameters.c.getLoc(idx)];
            f_old_incumbent += model_parameters.c.getVal(idx) * x_incumbent[model_parameters.c.getLoc(idx)];
        }
        // quadratic terms
        for (int col_idx = 0; col_idx < model_parameters.Q.getColLength(); ++col_idx) {
            int beg_idx = model_parameters.Q.getCbeg(col_idx);
            for (int idx = beg_idx; idx < model_parameters.Q.getClen(col_idx) + beg_idx; ++idx) {
                f_old_candidate += 0.5 * model_parameters.Q.getVal(idx) * x_candidate[model_parameters.Q.getRow(idx)] * x_candidate[col_idx];
                f_old_incumbent += 0.5 * model_parameters.Q.getVal(idx) * x_incumbent[model_parameters.Q.getRow(idx)] * x_incumbent[col_idx];
            }
        }
        double f_new_candidate = f_old_candidate;
        double f_new_incumbent = f_old_incumbent;
        f_old_candidate += res_master.eta;
        // second stage value at the incumbent
        double recourse_incumbent = f_lowerbound;
        for (int idx_minorant = 0; idx_minorant < minorant_collection.size(); ++idx_minorant) {
            double piece_val = minorant_collection[idx_minorant].alpha;
            piece_val += minorant_collection[idx_minorant].beta * x_incumbent;
            if (piece_val > recourse_incumbent) {
                recourse_incumbent = piece_val;
            }
        }
        f_old_incumbent += recourse_incumbent;
        // --- end computing the objective function value at the incumbent and candidate ---
        // --- update active minorants ---
        int dual_idx = -1;
        for (int idx_minorant = 0; idx_minorant < minorant_collection.size(); ++idx_minorant) {
            dual_idx++;
            if (res_master.dual_minorant[dual_idx] >= SOLVER_PRECISION_LOWER && res_master.dual_minorant[dual_idx] <= SOLVER_PRECISION_UPPER) { // remove inactive minorants
                minorant_collection.erase(minorant_collection.begin() + idx_minorant);
                idx_minorant--;
            }
        }
        // --- end updating active minorants ---
        // output candidate solution
        std::cout << "Candidate Solution: ";
        writeFile << "Candidate Solution: ";
        for (int x_index = 0; x_index < model_parameters.num_var_1stStage - 1; ++x_index) {
            std::cout << x_candidate[x_index] << ", ";
            writeFile << x_candidate[x_index] << ", ";
        }
        std::cout << x_candidate[model_parameters.num_var_1stStage - 1] << std::endl;
        writeFile << x_candidate[model_parameters.num_var_1stStage - 1] << std::endl;
        // --- end proximal mapping for updating candidate solutions ---
        // obtain a new data point
        int idx_scenario = iteration;
        stoPoint xi_cur = generator_random(model_parameters, stochastic_map, generator, idx_scenario, e_DB, C_DB);
        Xi.push_back(xi_cur);
        // --- second-stage subproblem ---
        // --- compute a new dual / face at x_candidate ---
        dualMultipliersQP newDual = sd_sqqp_2ndStageDual(model_parameters, stochastic_map, xi_cur, x_candidate);
        //dualMultipliersQP testDual = sd_sqqp_2ndStagePrimal(model_parameters, stochastic_map, xi_cur, x_candidate);
        
        face newFace = compute_face(newDual.s);
        if (if_face_new(explored_faces, newFace)) {
            explored_faces.push_back(newFace);
        }
        // --- end computing a new dual / face at x_candidate ---
        // --- construct new minorants
        // initialize new minorants
        minorant minorant_candidate;
        minorant minorant_incumbent;
        minorant_candidate.alpha = 0;
        minorant_incumbent.alpha = 0;
        for (int idx_x = 0; idx_x < model_parameters.num_var_1stStage; ++idx_x) {
            minorant_candidate.beta.push_back(0.0);
            minorant_incumbent.beta.push_back(0.0);
        }
        for (int idx_scenario = 0; idx_scenario < (int) N; ++idx_scenario) {
            double max_value_candiate = -SOLVER_INF;
            double max_value_incumbent = -SOLVER_INF;
            dualMultipliersQP max_dual_candidate;
            dualMultipliersQP max_dual_incumbent;
            int idx_face_max_candidate = 0;
            int idx_face_max_incumbent = 0;
            for (int idx_face = 0; idx_face < explored_faces.size(); ++idx_face) {
                // candidate
                dualMultipliersQP dual_candidate = sd_sqqp_2ndStageDual(model_parameters, stochastic_map, Xi[idx_scenario], x_candidate, explored_faces[idx_face]);
                // for test
                //dualMultipliersQP dual_candidate = sd_sqqp_2ndStagePrimal(model_parameters, stochastic_map, xi_cur, x_candidate);
                // for test 2
                //dualMultipliersQP dual_candidate = sd_sqqp_2ndStageDual(model_parameters, stochastic_map, xi_cur, x_candidate);
                // incumbent
                dualMultipliersQP dual_incumbent = sd_sqqp_2ndStageDual(model_parameters, stochastic_map, Xi[idx_scenario], x_incumbent, explored_faces[idx_face]);
                //dualMultipliersQP dual_incumbent = sd_sqqp_2ndStagePrimal(model_parameters, stochastic_map, xi_cur, x_incumbent);
                //dualMultipliersQP dual_incumbent = sd_sqqp_2ndStageDual(model_parameters, stochastic_map, xi_cur, x_incumbent);
                
                double current_value_candidate = dual_candidate.obj;
                double current_value_incumbent = dual_incumbent.obj;
                /*
                std::cout << "face: " << idx_face << ", max_value_candiate: " << max_value_candiate << std::endl;
                std::cout << "current_value_candidate: " << current_value_candidate << std::endl;
                std::cout << "max_value_incumbent: " << max_value_incumbent << std::endl;
                std::cout << "current_value_incumbent: " << current_value_incumbent << std::endl;
                */
                if (idx_face == 0) {
                    // candidate
                    max_value_candiate = current_value_candidate;
                    max_dual_candidate.obj = current_value_candidate;
                    max_dual_candidate.s = dual_candidate.s;
                    max_dual_candidate.t = dual_candidate.t;
                    idx_face_max_candidate = idx_face;
                    // incumbent
                    max_value_incumbent = current_value_incumbent;
                    max_dual_incumbent.obj = current_value_incumbent;
                    max_dual_incumbent.s = dual_incumbent.s;
                    max_dual_incumbent.t = dual_incumbent.t;
                    idx_face_max_incumbent = idx_face;
                }
                else { // if there is a tie, pick the latest one 
                    if (max_value_candiate <= current_value_candidate) {
                        max_value_candiate = current_value_candidate;
                        max_dual_candidate.obj = current_value_candidate;
                        max_dual_candidate.s = dual_candidate.s;
                        max_dual_candidate.t = dual_candidate.t;
                        idx_face_max_candidate = idx_face;
                    }
                    if (max_value_incumbent <= current_value_incumbent) {
                        max_value_incumbent = current_value_incumbent;
                        max_dual_incumbent.obj = current_value_incumbent;
                        max_dual_incumbent.s = dual_incumbent.s;
                        max_dual_incumbent.t = dual_incumbent.t;
                        idx_face_max_incumbent = idx_face;
                    }
                }
            } // end for (int idx_face = 0; idx_face < explored_faces.size(); ++idx_face)
            /*
            std::cout << "max index for candidate: " << idx_face_max_candidate << std::endl;
            std::cout << "max index for incumbent: " << idx_face_max_incumbent << std::endl;
            // for the test
            dualMultipliersQP dual_candidate2 = sd_sqqp_2ndStageDual(model_parameters, stochastic_map, Xi[idx_scenario], x_candidate);
            dualMultipliersQP dual_incumbent2 = sd_sqqp_2ndStageDual(model_parameters, stochastic_map, Xi[idx_scenario], x_incumbent);
            std::cout << "test(full faces) candidate: " << dual_candidate2.obj << std::endl;
            std::cout << "max_dual_candidate: " << max_dual_candidate.obj << std::endl;
            std::cout << "test(full faces) incumbent: " << dual_incumbent2.obj << std::endl;
            std::cout << "max_dual_incumbent: " << max_dual_incumbent.obj << std::endl;
            */
            // calculate alpha and beta
            // candidate
            // slope deterministic part
            std::vector<double> neg_t_C_candidate = model_parameters.C.fast_rightMultiply(max_dual_candidate.t);
            std::vector<double> neg_t_C_incumbent = model_parameters.C.fast_rightMultiply(max_dual_incumbent.t);
            // negate
            for (int idx = 0; idx < neg_t_C_candidate.size(); ++idx) {
                neg_t_C_candidate[idx] = (-1.0) * neg_t_C_candidate[idx];
                neg_t_C_incumbent[idx] = (-1.0) * neg_t_C_incumbent[idx];
            }
            // stochastic C
            // equality
            for (int idx_C = 0; idx_C < Xi[idx_scenario].C.size(); ++idx_C) {
                neg_t_C_candidate[stochastic_map.C_map[idx_C].second] -= Xi[idx_scenario].C[idx_C] * max_dual_candidate.t[stochastic_map.C_map[idx_C].first];
                neg_t_C_incumbent[stochastic_map.C_map[idx_C].second] -= Xi[idx_scenario].C[idx_C] * max_dual_incumbent.t[stochastic_map.C_map[idx_C].first];
            }
            // intercept
            double intercept_candidate = max_value_candiate - neg_t_C_candidate * x_candidate;
            double intercept_incumbent = max_value_incumbent - neg_t_C_incumbent * x_incumbent;
            // accumulate
            minorant_candidate.alpha += intercept_candidate;
            minorant_incumbent.alpha += intercept_incumbent;
            for (int idx = 0; idx < model_parameters.num_var_1stStage; ++idx) {
                minorant_candidate.beta[idx] += neg_t_C_candidate[idx];
                minorant_incumbent.beta[idx] += neg_t_C_incumbent[idx];
            }
        }
        // average
        double one_over_N = (1.0 / N);
        minorant_candidate.alpha *= one_over_N;
        minorant_incumbent.alpha *= one_over_N;
        for (int idx = 0; idx < model_parameters.num_var_1stStage; ++idx) {
            minorant_candidate.beta[idx] *= one_over_N;
            minorant_incumbent.beta[idx] *= one_over_N;
        }
        // --- end constructing new minorants
        // --- update old minorants ---
        double N_over_Nplus1 = N / (N + 1.0);
        double one_over_Nplus1 = 1.0 / (N + 1.0);
        for (int idx_minorant = 0; idx_minorant < minorant_collection.size(); ++idx_minorant) {
            minorant_collection[idx_minorant].alpha = minorant_collection[idx_minorant].alpha * N_over_Nplus1 + one_over_Nplus1 * f_lowerbound;
        } // end for (int idx_minorant = 0; idx_minorant < minorant_collection.size(); ++idx_minorant)
        // --- end updating minorant ---
        // add new minorants
        minorant_collection.push_back(minorant_candidate);
        minorant_collection.push_back(minorant_incumbent);
        // report new minorants
        std::cout << "Minorant Candidate\n";
        writeFile << "Minorant Candidate\n";
        std::cout << "alpha: " << minorant_candidate.alpha << std::endl;
        writeFile << "alpha: " << minorant_candidate.alpha << std::endl;
        std::cout << "beta: ";
        writeFile << "beta: ";
        for (int x_index = 0; x_index < model_parameters.num_var_1stStage - 1; ++x_index) {
            std::cout << minorant_candidate.beta[x_index] << ", ";
            writeFile << minorant_candidate.beta[x_index] << ", ";
        }
        std::cout << minorant_candidate.beta[model_parameters.num_var_1stStage - 1] << std::endl;
        writeFile << minorant_candidate.beta[model_parameters.num_var_1stStage - 1] << std::endl;
        std::cout << "Minorant Incumbent\n";
        writeFile << "Minorant Incumbent\n";
        std::cout << "alpha: " << minorant_incumbent.alpha << std::endl;
        writeFile << "alpha: " << minorant_incumbent.alpha << std::endl;
        std::cout << "beta: ";
        writeFile << "beta: ";
        for (int x_index = 0; x_index < model_parameters.num_var_1stStage - 1; ++x_index) {
            std::cout << minorant_incumbent.beta[x_index] << ", ";
            writeFile << minorant_incumbent.beta[x_index] << ", ";
        }
        std::cout << minorant_incumbent.beta[model_parameters.num_var_1stStage - 1] << std::endl;
        writeFile << minorant_incumbent.beta[model_parameters.num_var_1stStage - 1] << std::endl;
        // --- end updating old minorants ---
        // --- incumbent selection ---
        // LHS
        // new function value at candidate solution
        // second stage value
        double new_recourse_candidate = f_lowerbound;
        double new_recourse_incumbent = f_lowerbound;
        for (int idx_minorant = 0; idx_minorant < minorant_collection.size(); ++idx_minorant) {
            double piece_val = minorant_collection[idx_minorant].alpha;
            piece_val += minorant_collection[idx_minorant].beta * x_candidate;
            if (piece_val > new_recourse_candidate) {
                new_recourse_candidate = piece_val;
            }
            double piece_val2 = minorant_collection[idx_minorant].alpha;
            piece_val2 += minorant_collection[idx_minorant].beta * x_incumbent;
            if (piece_val2 > new_recourse_incumbent) {
                new_recourse_incumbent = piece_val2;
            }
        }
        f_new_candidate += new_recourse_candidate;
        f_new_incumbent += new_recourse_incumbent;
        double LHS = f_new_candidate - f_new_incumbent;
        double RHS = q * (f_old_candidate - f_old_incumbent);
        std::cout << "LHS: " << LHS << std::endl;
        writeFile << "LHS: " << LHS << std::endl;
        std::cout << "RHS: " << RHS << std::endl;
        writeFile << "RHS: " << RHS << std::endl;
        if (LHS <= RHS) {
            std::cout << "Computation Log: Incumbent selection criterion is passed.\n";
            writeFile <<"Computation Log: Incumbent selection criterion is passed.\n";
            x_incumbent = x_candidate;
            // update stepsize
            sigma = ((double) iteration + 1.0) * sigma_init;
            //sigma = max(sigma * 0.5, sigma_lowerbound);
        }
        else {
            std::cout << "Computation Log: Incumbent selection criterion is not passed.\n";
            writeFile <<"Computation Log: Incumbent solution selection criterion is not passed.\n";
            //sigma = min(sigma * 2.0, sigma_upperbound);
        }
        // print out the incumbent solution
        std::cout << "Incumbent Solution: ";
        writeFile << "Incumbent Solution: ";
        for (int index = 0; index < model_parameters.num_var_1stStage -1; ++index) {
            std::cout << x_incumbent[index] << ", ";
            writeFile << x_incumbent[index] << ", ";
        }
        std::cout << x_incumbent[model_parameters.num_var_1stStage - 1] << std::endl;
        writeFile << x_incumbent[model_parameters.num_var_1stStage - 1] << std::endl;
        //
        std::cout << "sigma: " << sigma << std::endl;
        writeFile << "sigma: " << sigma << std::endl;
        std::cout << "Number of faces: " << explored_faces.size() << std::endl;
        writeFile << "Number of faces: " << explored_faces.size() << std::endl;
        // update candidates and incuments
        old_candidate = x_candidate;
        old_incumbent = x_incumbent;
    } // end for (int iteration = 0; iteration < max_iterations; ++iteration)
    // --- end main loop ---
    std::cout << "*******************************************\n";
    writeFile << "*******************************************\n";
    std::cout << "Output Solution: ";
    writeFile << "Output Solution: ";
    for (int index = 0; index < model_parameters.num_var_1stStage-1; ++index) {
        std::cout << x_incumbent[index] << ", ";
        writeFile << x_incumbent[index] << ", ";
    }
    std::cout << x_incumbent[model_parameters.num_var_1stStage - 1] << std::endl;
    writeFile << x_incumbent[model_parameters.num_var_1stStage - 1] << std::endl;
    std::cout << "Computation Log: Finish Solving Process.\n";
    writeFile << "Computation Log: Finish Solving Process.\n";
    // write time elapsed
    double duration = (std::clock() - time_start ) / (double) CLOCKS_PER_SEC;
    writeFile << "Time elapsed(secs) : " << duration << "\n";
    writeFile << "*******************************************\n";
    writeFile.close();
    // return the results
    sd_output res;
    res.sol = x_incumbent;
    res.it_num = max_iterations;
    res.sol_time = duration;
    res.minorant_collection = minorant_collection;
    res.obj = compute_obj(x_incumbent, model_parameters, minorant_collection);
    // end cplex environment
    return res;
}

compromise_sd_output compromise_sd_sqqp_solver(const std::string& folder_path,
                                               std::mt19937& generator,
                                               stoPoint& point_est,
                                               int max_iterations,
                                               double f_lowerbound,
                                               double sigma_init,
                                               double regularizer_coefficient,
                                               int num_replications,
                                               double lb_error) {
    // Initialization
    std::string model_path = folder_path + "/model.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(CD_SD_SQQPv1.0).txt";
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
    std::cout << "Compromise Decision for SD for Solving SQQP Problems (v1.0)\n";
    writeFile << "Compromise Decision for SD for Solving SQQP Problems (v1.0)\n";
    std::cout << "Parameters\n";
    writeFile << "Parameters\n";
    std::cout << "num_replications: " << num_replications << std::endl;
    writeFile << "num_replications: " << num_replications << std::endl;
    std::cout << "max_iterations: " << max_iterations << std::endl;
    writeFile << "max_iterations: " << max_iterations << std::endl;
    std::cout << "sigma_init: " << sigma_init << std::endl;
    writeFile << "sigma_init: " << sigma_init << std::endl;
    std::cout << "Lower bound error: " << lb_error << std::endl;
    writeFile << "Lower bound error: " << lb_error << std::endl;
    
    compromise_sd_output res;
    // Step 1: Replication
    std::vector<sd_output> aggregate_outputs;
    for (int idx_rep = 0; idx_rep < num_replications; ++idx_rep) {
        sd_output replication_outputs = sd_sqqp_solver(folder_path, generator, point_est, max_iterations, f_lowerbound, sigma_init);
        aggregate_outputs.push_back(replication_outputs);
        res.replication_decision.push_back(replication_outputs.sol);
    }
    
    // Step 2: Aggregation
    std::vector<double> compromise_decision = sd_compromise_master(model_parameters, stochastic_map, aggregate_outputs, regularizer_coefficient, num_replications, lb_error);
    
    
    res.compromise_decision = compromise_decision;
    
    return res;
}
