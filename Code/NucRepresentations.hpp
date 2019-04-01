//
//  NucRepresentations.hpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#ifndef NucRepresentations_hpp
#define NucRepresentations_hpp

#include <stdio.h>
#include <string>
#include <vector>

#include "AtomNumRep.hpp"
#include "CompNumRep.hpp"
#include "DnaWNumRep.hpp"
#include "EIIPNumRep.hpp"
#include "IntNumRep.hpp"
#include "PairNumRep.hpp"
#include "RealNumRep.hpp"
#include "TetrNumRep.hpp"
#include "VossNumRep.hpp"
#include "ZCurNumRep.hpp"
#include "DataProcessing.hpp"
class NucRepresentations{
private:
public:
    NucRepresentations();
    NucRepresentations(std::string method, std::string read, std::vector<std::vector<double>> *rep);
    NucRepresentations(std::string method, std::string read, std::vector<std::vector<double>> *rep, bool accu, bool norm);
    
    void createRepresentation(std::string method, std::string read, std::vector<std::vector<double>> *rep);
    
    void createRepresentation(std::string method, std::string read, std::vector<std::vector<double>> *rep, bool accu, bool norm);
    ~NucRepresentations();
};

#endif /* NucRepresentations_hpp */
