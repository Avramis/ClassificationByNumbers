//
//  ReverseCompliment.hpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#ifndef ReverseCompliment_hpp
#define ReverseCompliment_hpp

#include <stdio.h>
#include <stdio.h>
#include <string>
#include <algorithm>

class ReverseCompliment{
private:
    
    char NucleotideComparison(char x);
    
public:
    ~ReverseCompliment();
    
    ReverseCompliment();
    
    std::string returnReversCompliment(std::string read);
    
};
#endif /* ReverseCompliment_hpp */
