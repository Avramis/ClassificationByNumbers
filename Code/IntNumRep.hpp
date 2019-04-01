//
//  IntNumRep.hpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#ifndef IntNumRep_hpp
#define IntNumRep_hpp

#include <stdio.h>
#include <vector>
#include <string>
class IntNumRep{
private:
public:
    IntNumRep(std::string read , std::vector<std::vector<double>> *rep);
    ~IntNumRep();
};
#endif /* IntNumRep_hpp */
