//
//  TestFileExistance.cpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#include "TestFileExistance.hpp"
TestFileExistance::TestFileExistance(){};

bool TestFileExistance::exists_test(std::string name){
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
};
