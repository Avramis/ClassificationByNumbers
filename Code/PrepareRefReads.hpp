//
//  PrepareRefReads.hpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#ifndef PrepareRefReads_hpp
#define PrepareRefReads_hpp

#include <stdio.h>
#include <vector>
#include <map>
#include "NGS.hpp"

class PrepareRefReads{
private:
public:
    PrepareRefReads();
    PrepareRefReads(std::vector<std::pair<std::string, std::string>> *RefList, std::vector<NGS> *RefNGS, std::vector<NGS> *TreeData, int kmer);
    
    ~PrepareRefReads();
};
#endif /* PrepareRefReads_hpp */
