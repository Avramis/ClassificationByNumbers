//
//  TestFileExistance.hpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#ifndef TestFileExistance_hpp
#define TestFileExistance_hpp

#include <stdio.h>
#include <sys/stat.h>
#include <string>
class TestFileExistance{
private:
public:
    TestFileExistance();
    //bool exists_test3 (const std::string& name) {
    bool exists_test(std::string name);
};
#endif /* TestFileExistance_hpp */
