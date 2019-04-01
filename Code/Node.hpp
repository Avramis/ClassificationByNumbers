//
//  Node.hpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#ifndef Node_hpp
#define Node_hpp

#include <stdio.h>
#include <random>
#include <stdlib.h>
#include <memory>
#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

#include "NGS.hpp"
#include "DistanceMethods.hpp"

class Node{
private:
    Node *lNode;
    Node *rNode;
    double tau = 0.00, thresh, maxdis = 0.00;
    bool Ntype = false;
    int pp = 1;
    
    NGS *vantagepoint;
    DistanceMethods DM;
    
    void setNtype(bool b);
    
    void setINode(double t, double c, double m);
    
    void QuickSort(std::vector<NGS> *dataPoints, int leftmost, int rightmost);
    
public:
    Node();
    int randominitiation(std::vector<NGS> *dataPoints, int m, int n);
    
    void generateNode(std::vector<NGS> *dataPoints, int m, int n, bool r);
    
    void generateNode(std::vector<NGS> *dataPoints, int m, int n, bool r, int l);
    
    void generateNode(std::vector<NGS> *dataPoints, int m, int n, bool r, int l, std::map<int, std::pair<int, double>> *levDetails);
    
    
    void initiateSearch(NGS *q, int k, std::map<double, std::vector<NGS>, std::greater<double>> *KNN, int *sc, int *absc, bool fs);
    
    void KnnSearch(NGS *q, int k, std::map<double, std::vector<NGS>, std::greater<double>> *KNN, int *sc, int *absc, bool fs);
    
    void initiateRangeSearch(NGS *q, double t, std::map<double, std::vector<NGS>, std::greater<double>> *NN, int *sc, int *absc);
    
    void RangeSearch(NGS *q, double t, std::map<double, std::vector<NGS>, std::greater<double>> *NN, int *sc, int *absc);
    
    void saveNode(std::ofstream *sout);
    
    int loadNode(std::ifstream *sin, std::vector<NGS> *dataPoints);
    
    void averageTau(double *x, double y);
    
    void totalTau(long double *x);
    
};
#endif /* Node_hpp */
