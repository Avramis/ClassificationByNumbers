//
//  SVDClass.hpp
//  ClassificationTool
//
//  Created by Samaneh Kouchaki
//  Copyright © 2018 Samaneh Kouchaki. All rights reserved.
//

#ifndef SVDClass_hpp
#define SVDClass_hpp

#include <stdio.h>
#include <vector>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <algorithm>


#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX(x,y) ((x)>(y)?(x):(y))

class SVDClass{
private:
    static bool doubleGreater(double a, double b);
public:
    SVDClass();
    std::vector<double> SVD(std::vector<std::vector<double>> *m);
    double PYTHAG(double a, double b);
    int dsvd(std::vector<std::vector<double>> *a, int m, int n, std::vector<double> *w, std::vector<std::vector<double>> *v);
    
    
};

#endif /* SVDClass_hpp */
