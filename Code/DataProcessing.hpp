//
//  DataProcessing.hpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#ifndef DataProcessing_hpp
#define DataProcessing_hpp

#include <stdio.h>
#include <vector>
#include <math.h>

class DataProcessing{
private:
public:
    DataProcessing();
    void repAccumulation(std::vector<std::vector<double>> *rep);
    
    void repAccumulation(std::vector<std::vector<double>> *rep, int s, int l);
    
    void repZnormalisation(std::vector<std::vector<double>> *rep);
    
    void repZnormalisation(std::vector<std::vector<double>> *rep, int s, int l);
    
    ~DataProcessing();
};
#endif /* DataProcessing_hpp */
