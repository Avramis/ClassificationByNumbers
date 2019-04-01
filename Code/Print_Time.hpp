//
//  Print_Time.hpp
//  ClassificationTool
//
//  Created by Avraam Tapinos on 24/04/2018.
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#ifndef Print_Time_hpp
#define Print_Time_hpp

#include <stdio.h>
#include <iostream>
class Print_Time{
private:
    
public:
    Print_Time();
    
    std::string returnTimestr(clock_t start, clock_t end);
};

#endif /* Print_Time_hpp */
