//
//  FastaParser.hpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#ifndef FastaParser_hpp
#define FastaParser_hpp

#include <stdio.h>
#include <iostream>
#include<vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include "NGS.hpp"
#include "NucRepresentations.hpp"
#include "DataTransformations.hpp"
#include "DataProcessing.hpp"

class FastaParser{
private:
    NucRepresentations rep;
    DataTransformations tra;
    DataProcessing DaPr;
public:
    FastaParser(std::string filepath, std::vector<std::pair<std::string, std::string>> *ReferenceList);
    
    ~FastaParser();
};
#endif /* FastaParser_hpp */
