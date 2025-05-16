//
//  ioModel.cpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/17/25.
//

#include "ioModel.hpp"

// declare global variables (need to be defined in the source file)
double MODEL_PRECISION_LOWER = -1e-6;
double MODEL_PRECISION_UPPER = 1e-6;

standardTwoStageParameters readStandardTwoStageParameters(const std::string& parameterPath) {
    // initialize variables
    standardTwoStageParameters parameters;
    // constant strings
    const std::string nameBeginModel("<model:twoStageSLP>");
    const std::string nameEndModel("</model:twoStageSLP>");
    // first stage
    const std::string nameBeginParameter_c("<c>");
    const std::string nameEndParameter_c("</c>");
    const std::string nameBeginParameter_Q("<Q>"); // quadratic terms in the objective
    const std::string nameEndParameter_Q("</Q>"); // quadratic terms in the objective
    const std::string nameBeginParameter_AL("<AL>");
    const std::string nameEndParameter_AL("</AL>");
    const std::string nameBeginParameter_AE("<AE>");
    const std::string nameEndParameter_AE("</AE>");
    const std::string nameBeginParameter_AG("<AG>");
    const std::string nameEndParameter_AG("</AL=G>");
    const std::string nameBeginParameter_bL("<bL>");
    const std::string nameEndParameter_bL("</bL>");
    const std::string nameBeginParameter_bE("<bE>");
    const std::string nameEndParameter_bE("</bE>");
    const std::string nameBeginParameter_bG("<bG>");
    const std::string nameEndParameter_bG("</bG>");
    // second stage
    const std::string nameBeginParameter_d("<d>");
    const std::string nameEndParameter_d("</d>");
    const std::string nameBeginParameter_P("<P>"); // quadratic terms in the objective
    const std::string nameEndParameter_P("</P>"); // quadratic terms in the objective
    const std::string nameBeginParameter_D("<D>");
    const std::string nameEndParameter_D("</D>");
    const std::string nameBeginParameter_C("<C>");
    const std::string nameEndParameter_C("</C>");
    const std::string nameBeginParameter_e("<e>");
    const std::string nameEndParameter_e("</e>");
    // intermediate matrices
    //Pinv Pinv_Dt; // here, D' = [D Islack] D_Pinv_Dt; D_Pinv;
    const std::string nameBeginParameter_Pinv("<Pinv>");
    const std::string nameEndParameter_Pinv("</Pinv>");
    const std::string nameBeginParameter_Pinv_Dt("<Pinv_Dt>");
    const std::string nameEndParameter_Pinv_Dt("</Pinv_Dt>");
    const std::string nameBeginParameter_D_Pinv_Dt("<D_Pinv_Dt>");
    const std::string nameEndParameter_D_Pinv_Dt("</D_Pinv_Dt>");
    const std::string nameBeginParameter_D_Pinv("<D_Pinv>");
    const std::string nameEndParameter_D_Pinv("</D_Pinv>");
    // intermediate value
    const std::string nameBeginParameter_dt_Pinv_d("<dt_Pinv_d>");
    const std::string nameEndParameter_dt_Pinv_d("</dt_Pinv_d>");
    // intermediate vectors
    const std::string nameBeginParameter_dt_Pinv_Dt("<dt_Pinv_Dt>");
    const std::string nameEndParameter_dt_Pinv_Dt("</dt_Pinv_Dt>");
    const std::string nameBeginParameter_dt_Pinv("<dt_Pinv>");
    const std::string nameEndParameter_dt_Pinv("</dt_Pinv>");
    // slack
    const std::string nameBeginParameter_slack("<Slack>");
    const std::string nameEndParameter_slack("</Slack>");
    // specs
    const std::string nameBeginParameter_specs("<specs>");
    const std::string nameEndParameter_specs("</specs>");
    std::string readCondition("null");
    const char* parameterPathConst = parameterPath.c_str(); // convert a string path into a constant path
    std::ifstream readFile(parameterPathConst);
    int count_line = 0;
    int count_line_specs = 0;
    if (readFile.is_open()) {
        std::string line1;
        while (getline(readFile,line1)) {// get the whole line
            std::stringstream ss1(line1); // convert a string into stream
            if (nameBeginModel.compare(line1) != 0 && nameEndModel.compare(line1) != 0) {
                // main content
                if (nameBeginParameter_c.compare(line1) == 0) { // beign reading parameter c
                    readCondition = "c";
                    }
                else if (nameEndParameter_c.compare(line1) == 0) {
                    readCondition = "null";
                }
                else if (nameBeginParameter_Q.compare(line1) == 0) { // begin reading parameter Q
                    readCondition = "Q";
                }
                else if (nameEndParameter_Q.compare(line1) == 0) {
                    readCondition = "null";
                }
                else if (nameBeginParameter_AL.compare(line1) == 0) { // begin reading parameter AL
                    readCondition = "AL";
                }
                else if (nameEndParameter_AL.compare(line1) == 0) {
                    readCondition = "null";
                }
                else if (nameBeginParameter_AG.compare(line1) == 0) { // begin reading AG
                    readCondition = "AG";
                }
                else if (nameEndParameter_AG.compare(line1) == 0) { // end reading AG
                    readCondition = "null";
                }
                else if (nameBeginParameter_AE.compare(line1) == 0) { // begin reading AE
                    readCondition = "AE";
                }
                else if (nameEndParameter_AE.compare(line1) == 0) { // end reading AE
                    readCondition = "null";
                }
                else if (nameBeginParameter_bL.compare(line1) == 0) { // begin reading bL
                    readCondition = "bL";
                }
                else if (nameEndParameter_bL.compare(line1) == 0) { // end reading bL
                    readCondition = "null";
                }
                else if (nameBeginParameter_bE.compare(line1) == 0) { // begin reading bE
                    readCondition = "bE";
                }
                else if (nameEndParameter_bE.compare(line1) == 0) { // end reading bE
                    readCondition = "null";
                }
                else if (nameBeginParameter_bG.compare(line1) == 0) { // begin reading bG
                    readCondition = "bG";
                }
                else if (nameEndParameter_bG.compare(line1) == 0) { // end reading bG
                    readCondition = "null";
                }
                else if (nameBeginParameter_d.compare(line1) == 0) { // begin reading d
                    readCondition = "d";
                }
                else if (nameEndParameter_d.compare(line1) == 0) { // begin reading d
                    readCondition = "null";
                }
                else if (nameBeginParameter_P.compare(line1) == 0) { // end reading P
                    readCondition = "P";
                }
                else if (nameEndParameter_P.compare(line1) == 0) {
                    readCondition = "null"; 
                }
                else if (nameBeginParameter_D.compare(line1) == 0) { // begin reading D
                    readCondition = "D";
                }
                else if (nameEndParameter_D.compare(line1) == 0) { // end reading D
                    readCondition = "null";
                }
                else if (nameBeginParameter_C.compare(line1) == 0) { // begin reading C
                    readCondition = "C";
                }
                else if (nameEndParameter_C.compare(line1) == 0) { // emd reading C
                    readCondition = "null";
                }
                else if (nameBeginParameter_e.compare(line1) == 0) { // begin reading d
                    readCondition = "e";
                }
                else if (nameEndParameter_e.compare(line1) == 0) { // end reading d
                    readCondition = "null";
                }
                else if (nameBeginParameter_Pinv.compare(line1) == 0) { // begin reading Pinv
                    readCondition = "Pinv";
                }
                else if (nameEndParameter_Pinv.compare(line1) == 0) {
                    readCondition = "null";
                }
                else if (nameBeginParameter_D_Pinv.compare(line1) == 0) {
                    readCondition = "D_Pinv";
                }
                else if (nameEndParameter_D_Pinv.compare(line1) == 0) {
                    readCondition = "null";
                }
                else if (nameBeginParameter_Pinv_Dt.compare(line1) == 0) {
                    readCondition = "Pinv_Dt";
                }
                else if (nameEndParameter_Pinv_Dt.compare(line1) == 0) {
                    readCondition = "null";
                }
                else if (nameBeginParameter_D_Pinv_Dt.compare(line1) == 0) {
                    readCondition = "D_Pinv_Dt";
                }
                else if (nameEndParameter_D_Pinv_Dt.compare(line1) == 0) {
                    readCondition = "null";
                }
                else if (nameBeginParameter_dt_Pinv_d.compare(line1) == 0) {
                    readCondition = "dt_Pinv_d";
                }
                else if (nameEndParameter_dt_Pinv_d.compare(line1) == 0) {
                    readCondition = "null";
                }
                else if (nameBeginParameter_dt_Pinv_Dt.compare(line1) == 0) {
                    readCondition = "dt_Pinv_Dt";
                }
                else if (nameEndParameter_dt_Pinv_Dt.compare(line1) == 0) {
                    readCondition = "null";
                }
                else if (nameBeginParameter_dt_Pinv.compare(line1) == 0) {
                    readCondition = "dt_Pinv";
                }
                else if (nameEndParameter_dt_Pinv.compare(line1) == 0) {
                    readCondition = "null";
                }
                else if (nameBeginParameter_slack.compare(line1) == 0) { // begin reading slack
                    readCondition = "slack";
                }
                else if (nameEndParameter_slack.compare(line1) == 0) { // end reading slack
                    readCondition = "null";
                }
                else if (nameBeginParameter_specs.compare(line1) == 0) { // begin reading slack
                    readCondition = "specs";
                }
                else if (nameEndParameter_specs.compare(line1) == 0) { // begin reading slack
                    readCondition = "null";
                }
                else {
                    if (readCondition.compare("specs") == 0) { // input sytax: name:val;
                        if (count_line_specs == 0) { // number of decision variables in the 1st stage
                            while (getline(ss1, line1, ';')) {
                                std::stringstream ss2(line1);
                                int count_item = 0;
                                long num = 0;
                                while (getline(ss2, line1, ':')) {
                                    std::stringstream ss3(line1);
                                    if (count_item == 1) {
                                        ss3 >> num;
                                        parameters.num_var_1stStage = num;
                                        }
                                    count_item += 1;
                                }
                            }
                        } // end if (count_line_specs == 0)
                        else if (count_line_specs == 1) { // number of decision variables in the 2nd stage
                            while (getline(ss1, line1, ';')) {
                                std::stringstream ss2(line1);
                                int count_item = 0;
                                long num = 0;
                                while (getline(ss2, line1, ':')) {
                                    std::stringstream ss3(line1);
                                    if (count_item == 1) {
                                        ss3 >> num;
                                        parameters.num_var_2ndStage = num;
                                        }
                                    count_item += 1;
                                }
                            }
                        } // end else if (count_line_specs == 1)
                        else if (count_line_specs == 2) { // number of slack variables in the 2nd stage
                            while (getline(ss1, line1, ';')) {
                                std::stringstream ss2(line1);
                                int count_item = 0;
                                long num = 0;
                                while (getline(ss2, line1, ':')) {
                                    std::stringstream ss3(line1);
                                    if (count_item == 1) {
                                        ss3 >> num;
                                        parameters.num_slack_2ndStage = num;
                                        }
                                    count_item += 1;
                                }
                            }
                        } // end else if (count_line_specs == 2)
                        else if (count_line_specs == 3) { // number of equality constraints in the 1st stage
                            while (getline(ss1, line1, ';')) {
                                std::stringstream ss2(line1);
                                int count_item = 0;
                                long num = 0;
                                while (getline(ss2, line1, ':')) {
                                    std::stringstream ss3(line1);
                                    if (count_item == 1) {
                                        ss3 >> num;
                                        parameters.num_E_1stStage = num;
                                        }
                                    count_item += 1;
                                }
                            }
                        } // end else if (count_line_specs == 3)
                        else if (count_line_specs == 4) { // number of L constraints in the 1st stage
                            while (getline(ss1, line1, ';')) {
                                std::stringstream ss2(line1);
                                int count_item = 0;
                                long num = 0;
                                while (getline(ss2, line1, ':')) {
                                    std::stringstream ss3(line1);
                                    if (count_item == 1) {
                                        ss3 >> num;
                                        parameters.num_L_1stStage = num;
                                        }
                                    count_item += 1;
                                }
                            }
                        } // end else if (count_line_specs == 4)
                        else if (count_line_specs == 5) { // number of G constraints in the 1st stage
                            while (getline(ss1, line1, ';')) {
                                std::stringstream ss2(line1);
                                int count_item = 0;
                                long num = 0;
                                while (getline(ss2, line1, ':')) {
                                    std::stringstream ss3(line1);
                                    if (count_item == 1) {
                                        ss3 >> num;
                                        parameters.num_G_1stStage = num;
                                        }
                                    count_item += 1;
                                }
                            }
                        } // end else if (count_line_specs == 5)
                        else if (count_line_specs == 6) { // number of E constraints in the 2nd stage
                            while (getline(ss1, line1, ';')) {
                                std::stringstream ss2(line1);
                                int count_item = 0;
                                long num = 0;
                                while (getline(ss2, line1, ':')) {
                                    std::stringstream ss3(line1);
                                    if (count_item == 1) {
                                        ss3 >> num;
                                        parameters.num_E_2ndStage = num;
                                        }
                                    count_item += 1;
                                }
                            }
                        } // end else if (count_line_specs == 6)
                        else if (count_line_specs == 7) { // number of L constraints in the 2nd stage
                            while (getline(ss1, line1, ';')) {
                                std::stringstream ss2(line1);
                                int count_item = 0;
                                long num = 0;
                                while (getline(ss2, line1, ':')) {
                                    std::stringstream ss3(line1);
                                    if (count_item == 1) {
                                        ss3 >> num;
                                        parameters.num_L_2ndStage = num;
                                        }
                                    count_item += 1;
                                }
                            }
                        } // end else if (count_line_specs == 7)
                        else if (count_line_specs == 8) { // number of G constraints in the 2nd stage
                            while (getline(ss1, line1, ';')) {
                                std::stringstream ss2(line1);
                                int count_item = 0;
                                long num = 0;
                                while (getline(ss2, line1, ':')) {
                                    std::stringstream ss3(line1);
                                    if (count_item == 1) {
                                        ss3 >> num;
                                        parameters.num_G_2ndStage = num;
                                        }
                                    count_item += 1;
                                }
                            }
                        } // end else if (count_line_specs == 8)
                        else if (count_line_specs == 9) { // upper bound of decision variable in the 1st stage
                            while (getline(ss1, line1, ';')) {
                                std::stringstream ss2(line1);
                                int count_item = 0;
                                double ub = 0;
                                while (getline(ss2, line1, ':')) {
                                    std::stringstream ss3(line1);
                                    if (count_item == 1) {
                                        ss3 >> ub;
                                        parameters.var_ub_1stStage = ub;
                                        }
                                    count_item += 1;
                                }
                            }
                        } // end else if (count_line_specs == 8)
                        count_line_specs += 1;
                    } // end if (readCondition.compare("specs") == 0)
                    else if (readCondition.compare("c") == 0) { // input syntax: pos,val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
                            int pos = 0;
                            while (getline(ss2, line1, ',')) {
                                std::stringstream ss3(line1);
                                if (count_item == 0) {
                                    ss3 >> pos;
                                }
                                else if (count_item == 1) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) { // update c
                                parameters.c.insert(pos, val);
                            }
                        }
                    } // end if (readCondition.compare("c") == 0)
                    else if (readCondition.compare("Q") == 0) { // input syntax: row,col,val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
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
                                else if (count_item == 2) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) {
                                parameters.Q.insert(row, col, val);
                            }
                        }
                    }
                    else if (readCondition.compare("AL") == 0) { // input syntax: row,col,val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
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
                                else if (count_item == 2) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) {
                                //std::pair<int, int> loc(row,col);
                                parameters.AL.insert(row, col, val);
                            }
                        }
                    } // end else if (readCondition.compare("AL") == 0)
                    else if (readCondition.compare("AE") == 0) { // input syntax: row,col,val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
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
                                else if (count_item == 2) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) {
                                parameters.AE.insert(row, col, val);
                            }
                        }
                    } // end else if (readCondition.compare("AE") == 0)
                    else if (readCondition.compare("AG") == 0) { // input syntax: row,col,val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
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
                                else if (count_item == 2) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) {
                                parameters.AG.insert(row, col, val);
                            }
                        }
                    } // end else if (readCondition.compare("AG") == 0)
                    else if (readCondition.compare("bL") == 0) { // input syntax: pos,val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
                            int pos = 0;
                            while (getline(ss2, line1, ',')) {
                                std::stringstream ss3(line1);
                                if (count_item == 0) {
                                    ss3 >> pos;
                                }
                                else if (count_item == 1) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) { // update c
                                parameters.bL.insert(pos, val);
                            }
                        }
                    } // end  else if (readCondition.compare("bL") == 0)
                    else if (readCondition.compare("bE") == 0) { // input syntax: pos,val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
                            int pos = 0;
                            while (getline(ss2, line1, ',')) {
                                std::stringstream ss3(line1);
                                if (count_item == 0) {
                                    ss3 >> pos;
                                }
                                else if (count_item == 1) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) { // update c
                                parameters.bE.insert(pos, val);
                            }
                        }
                    } // end else if (readCondition.compare("bE") == 0)
                    else if (readCondition.compare("bG") == 0) { // input syntax: pos,val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
                            int pos = 0;
                            while (getline(ss2, line1, ',')) {
                                std::stringstream ss3(line1);
                                if (count_item == 0) {
                                    ss3 >> pos;
                                }
                                else if (count_item == 1) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) { // update c
                                parameters.bG.insert(pos, val);
                            }
                        }
                    } // end else if (readCondition.compare("bG") == 0)
                    else if (readCondition.compare("d") == 0) {// input syntax: pos,val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
                            int pos = 0;
                            while (getline(ss2, line1, ',')) {
                                std::stringstream ss3(line1);
                                if (count_item == 0) {
                                    ss3 >> pos;
                                }
                                else if (count_item == 1) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) { // update c
                                parameters.d.insert(pos, val);
                            }
                        }
                    } // end else if (readCondition.compare("d") == 0)
                    else if (readCondition.compare("P") == 0) { // input syntax: row,col,val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
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
                                else if (count_item == 2) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) {
                                parameters.P.insert(row, col, val);
                            }
                        }
                    }
                    else if (readCondition.compare("D") == 0) { // input syntax: row,col,val;
                        while (getline(ss1, line1, ';')) {
                            //std::cout << line1 << std::endl;
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
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
                                else if (count_item == 2) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) {
                                std::pair<int, int> loc(row,col);
                                parameters.D.insert(row, col, val);
                            }
                        }
                    } // end else if (readCondition.compare("D") == 0)
                    else if (readCondition.compare("C") == 0) { // input syntax: row,col,val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
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
                                else if (count_item == 2) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) {
                                std::pair<int, int> loc(row,col);
                                parameters.C.insert(row, col, val);
                            }
                        }
                    } // end else if (readCondition.compare("C") == 0)
                    else if (readCondition.compare("slack") == 0) { // input syntax: row,col,val; notes: revisit after setting up the solver env
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
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
                                else if (count_item == 2) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) {
                                std::pair<int, int> loc(row,col);
                                parameters.ISlack.insert(row, col, val);
                            }
                        }
                    } // end else if (readCondition.compare("slack") == 0)
                    else if (readCondition.compare("e") == 0) { // input syntax: pos,val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
                            int pos = 0;
                            while (getline(ss2, line1, ',')) {
                                std::stringstream ss3(line1);
                                if (count_item == 0) {
                                    ss3 >> pos;
                                }
                                else if (count_item == 1) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) { // update e
                                parameters.e.insert(pos, val);
                            }
                        }
                    } // end if (readCondition.compare("e") == 0)
                    else if (readCondition.compare("Pinv") == 0) { // input syntax: row,col,val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
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
                                else if (count_item == 2) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) {
                                parameters.Pinv.insert(row, col, val);
                            }
                        }
                    }
                    else if (readCondition.compare("Pinv_Dt") == 0) { // input syntax: row,col,val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
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
                                else if (count_item == 2) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) {
                                parameters.Pinv_Dt.insert(row, col, val);
                            }
                        }
                    }
                    else if (readCondition.compare("D_Pinv") == 0) { // input syntax: row,col,val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
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
                                else if (count_item == 2) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) {
                                parameters.D_Pinv.insert(row, col, val);
                            }
                        }
                    }
                    else if (readCondition.compare("D_Pinv_Dt") == 0) { // input syntax: row,col,val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
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
                                else if (count_item == 2) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) {
                                parameters.D_Pinv_Dt.insert(row, col, val);
                            }
                        }
                    }
                    else if (readCondition.compare("dt_Pinv_d") == 0) { // input syntax: val;
                        double val = 0;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            ss2 >> val;
                        }
                        parameters.dt_Pinv_d = val;
                    }
                    else if (readCondition.compare("dt_Pinv_Dt") == 0) { // input syntax: pos, val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
                            int pos = 0;
                            while (getline(ss2, line1, ',')) {
                                std::stringstream ss3(line1);
                                if (count_item == 0) {
                                    ss3 >> pos;
                                }
                                else if (count_item == 1) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) { // update e
                                parameters.dt_Pinv_Dt.insert(pos, val);
                            }
                        }
                    }
                    else if (readCondition.compare("dt_Pinv") == 0) { // input syntax: pos, val;
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            int count_item = 0;
                            double val = 0;
                            int pos = 0;
                            while (getline(ss2, line1, ',')) {
                                std::stringstream ss3(line1);
                                if (count_item == 0) {
                                    ss3 >> pos;
                                }
                                else if (count_item == 1) {
                                    ss3 >> val;
                                }
                                count_item += 1;
                            }
                            if (val > MODEL_PRECISION_UPPER || val < MODEL_PRECISION_LOWER) { // update e
                                parameters.dt_Pinv.insert(pos, val);
                            }
                        }
                    }
                } // end else
            }
        } // end while (getline(readFile,line1))
    }
    // finalize inserting
    // first stage
    parameters.c.setLen(parameters.num_var_1stStage);
    if (parameters.num_E_1stStage > 0) {
        parameters.AE.insert_end(parameters.num_E_1stStage, parameters.num_var_1stStage);
    }
    if (parameters.num_L_1stStage > 0) {
        parameters.AL.insert_end(parameters.num_L_1stStage, parameters.num_var_1stStage);
    }
    if (parameters.num_G_1stStage > 0) {
        parameters.AG.insert_end(parameters.num_G_1stStage, parameters.num_var_1stStage);
    }
    if (parameters.Q.getNonZeroEntries() > 0) {
        parameters.Q.insert_end(parameters.num_var_1stStage, parameters.num_var_1stStage);
    }
    // second stage subproblem
    parameters.d.setLen(parameters.num_var_2ndStage);
    long num_total_var_2ndStage = parameters.num_var_2ndStage + parameters.num_slack_2ndStage;
    if (parameters.P.getNonZeroEntries() > 0) {
        parameters.P.insert_end(num_total_var_2ndStage, num_total_var_2ndStage);
    }
    if (parameters.Pinv.getNonZeroEntries() > 0) {
        parameters.Pinv.insert_end(num_total_var_2ndStage, num_total_var_2ndStage);
    }
    if (parameters.dt_Pinv.getNzeroLen() > 0) {
        parameters.dt_Pinv.setLen(num_total_var_2ndStage);
    }
    long num_con_2ndStage = parameters.num_E_2ndStage + parameters.num_L_2ndStage + parameters.num_G_2ndStage;
    if (parameters.D_Pinv_Dt.getNonZeroEntries() > 0) {
        parameters.D_Pinv_Dt.insert_end(num_con_2ndStage, num_con_2ndStage);
    }
    if (parameters.D_Pinv.getNonZeroEntries() > 0) {
        parameters.D_Pinv.insert_end(num_con_2ndStage, num_total_var_2ndStage);
    }
    if (parameters.dt_Pinv_Dt.getNzeroLen()) {
        parameters.dt_Pinv_Dt.setLen(num_con_2ndStage);
    }
    parameters.D.insert_end(num_con_2ndStage, parameters.num_var_2ndStage);
    parameters.C.insert_end(num_con_2ndStage, parameters.num_var_1stStage);
    if (parameters.num_slack_2ndStage > 0) {
        parameters.ISlack.insert_end(num_con_2ndStage, parameters.num_slack_2ndStage);
    }
    parameters.e.setLen(num_con_2ndStage);
    return parameters;
}
