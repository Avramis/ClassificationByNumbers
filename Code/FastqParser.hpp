//
//  FastqParser.hpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#ifndef FastqParser_hpp
#define FastqParser_hpp

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

class FastqParser{
private:
    NucRepresentations rep;
    DataTransformations tra;
    DataProcessing DaPr;
public:
    FastqParser(std::string filepath, std::string repmeth, int refkmerlen, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences, int sa, int ss, bool conper, int rcount);
    
    ~FastqParser();
};
#endif /* FastqParser_hpp */
