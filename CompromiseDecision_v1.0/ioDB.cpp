//
//  ioDB.cpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/18/25.
//

#include "ioDB.hpp"

std::vector<std::vector<dataPoint>> readDB(std::string readPath){
    std::cout << "Read Database" << std::endl;
    std::vector<std::vector<dataPoint>> dataPointDB; // database
    std::vector<dataPoint> dataPointDS; // dataset
    std::string nameBeginDB("<database>");
    std::string nameEndDB("</database>");
    std::string nameBeginDS("<dataset>");
    std::string nameEndDS("</dataset>");
    bool flag_ifDS = false; // flag of whether the current dataset has been reading
    const char* readPathConst = readPath.c_str(); // convert the string type path to constant
    std::ifstream readFile(readPathConst); // create a readFile object
    if (readFile.is_open()) {
        std::string line1;
        while (getline(readFile, line1)) { // get the whole line
            std::stringstream ss1(line1); // convert a string into stream
            dataPoint dataPointTemp;
            unsigned int index_position = 0; // 1 for response, 3 for weight
            if (nameBeginDB.compare(line1) != 0 && nameEndDB.compare(line1) != 0) { // main
                // contents of DB
                //std::cout << line1 << std::endl;
                //std::cout << "Main Content" << std::endl;
                if (nameBeginDS.compare(line1) == 0 && !flag_ifDS) {// if it is a new dataset
                    //std::cout << "Begin of Dataset" << std::endl;
                    flag_ifDS = true; // start reading new dataset and set the flag to be true
                }
                else if (flag_ifDS && nameEndDS.compare(line1) != 0){ // if it is the current dataset
                    //std::cout << "Data Point" << std::endl;
                    // syntax:  response:(numbers);weight:(number);
                    while (getline(ss1, line1, ';')) {
                        std::stringstream ss2(line1);
                        while (getline(ss2, line1, ':')) {
                            if (index_position == 1){ // read vector
                                std::stringstream ss3(line1);
                                while (getline(ss3, line1, ',')) {
                                    double value;
                                    std::stringstream ss4(line1);
                                    ss4 >> value;
                                    dataPointTemp.response.push_back(value);
                                }
                            }
                            else if (index_position == 3){
                                double value;
                                std::stringstream ss3(line1);
                                ss3 >> value;
                                dataPointTemp.weight = value;
                            }
                            index_position++;
                        }
                    }
                    dataPointDS.push_back(dataPointTemp);
                } // end else if (flag_ifDS && nameEndDS.compare(line1) != 0)
                else if (nameEndDS.compare(line1) == 0){ // if it is the end of a dataset
                    //std::cout << "End of Dataset" << std::endl;
                    flag_ifDS = false; // finishing reading current dataset and set the flag to be false
                    dataPointDB.push_back(dataPointDS); // add the dataset the database
                    dataPointTemp = {}; // clear a struct
                    dataPointDS.clear(); // clear the temporary dataset
                }
            } // end if (nameBeginDB.compare(line1) != 0 && nameEndDB.compare(line1) != 0)
        }
    }// end if (readFile.is_open())
    readFile.close(); // close the file
    return dataPointDB;
}


// print dataPoint
void printDB(const std::vector<std::vector<dataPoint>>& dataPointDB) {
    long DS_number = dataPointDB.size();
    std::cout << DS_number << std::endl;
    for (int DS_index = 0; DS_index < DS_number; DS_index++) {
        long dataPoint_number = dataPointDB[DS_index].size();
        std::cout << "Dataset: " << DS_index + 1 << std::endl;
        for (int dataPoint_index = 0; dataPoint_index < dataPoint_number; dataPoint_index++) {
            long response_size = dataPointDB[DS_index][dataPoint_index].response.size(); // size of predictor
            std::cout << "Response " << dataPoint_index + 1 << ": ";
            for (int response_index = 0; response_index < response_size; response_index++) {
                std::cout << dataPointDB[DS_index][dataPoint_index].response[response_index] << ", ";
            }
            std::cout << "Weight: " << dataPointDB[DS_index][dataPoint_index].weight << std::endl;
        }
    }
}


// print dataPoint
void printDataPoint(const dataPoint& dataPoint01) {
    // size
    long response_size = dataPoint01.response.size();
    std::cout << "Response: ";
    for (int index = 0; index < response_size; ++index) {
        std::cout << dataPoint01.response[index] << " ";
    }
    std::cout << "Weight: " << dataPoint01.weight << std::endl;
}
