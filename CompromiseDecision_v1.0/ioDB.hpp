//
//  ioDB.hpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/18/25.
//

#ifndef ioDB_hpp
#define ioDB_hpp

#include <stdio.h>
#include <fstream>
#include <sstream>

#include "DataStructure.hpp"

// read database from a text file
std::vector<std::vector<dataPoint>> readDB(std::string readPath);

// print database
void printDB(const std::vector<std::vector<dataPoint>>& dataPointDB);
// print dataPoint
void printDataPoint(const dataPoint& dataPoint01);
#endif /* ioDB_hpp */
