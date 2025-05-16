//
//  ioModel.hpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/17/25.
//

#ifndef ioModel_hpp
#define ioModel_hpp

#include <stdio.h>

#include "DataStructure.hpp"

// standard form in the second stage LP
// |DE|     |CE|     |0  0 |     |eE|
// |DL| y + |CL| x + |IL 0 | s = |eL|
// |DG|     |CG|     |0 -IG|     |eG|
// y and s share the indices
standardTwoStageParameters readStandardTwoStageParameters(const std::string& parameterPath);


#endif /* ioModel_hpp */
