//
//  Print_Time.cpp
//  ClassificationTool
//
//  Created by Avraam Tapinos on 24/04/2018.
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#include "Print_Time.hpp"

Print_Time::Print_Time(){};

std::string Print_Time::returnTimestr(clock_t start, clock_t end){
    std::string timestr = "";
    double finaltime = (double)(end - start)/CLOCKS_PER_SEC;
    int h = (int) (finaltime/3600.00);
    finaltime -= ((double)h * 3600.00);
    int m = (int) (finaltime/60.00);
    finaltime -= ((double)m * 60.00);
    
    timestr = "";
    if(h < 10){
        timestr += "0";
    }
    timestr += std::to_string(h) +":";
    
    if(m < 10){
        timestr += "0";
    }
    timestr += std::to_string(m) + ":" ;
    
    if(finaltime < 10){
        timestr += "0";
    }
    timestr += std::to_string(finaltime);
    
    return timestr;
};
