//
//  VPTree.hpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//
#ifndef VPTree_hpp
#define VPTree_hpp

#include <stdio.h>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

#include "Node.hpp"
#include "NGS.hpp"
class VPTree{
private:
    Node vp;
    std::map<int, std::pair<int, double>, std::less<int>> levDetails;

public:
    
    VPTree();
    void TreeGeneration(std::vector<NGS> *dataPoints);
    
    void buildTree(Node *vp, int n, std::vector<NGS> *dataPoints);
    
    void buildTree(Node *vp, int n, std::vector<NGS> *dataPoints, std::map<int, std::pair<int, double>> *levdetails);
    
    void KnnSearch(NGS *q, int k, std::map<double, std::vector<NGS>, std::greater<double>> *KNN, int *sc, int *absc, bool sen);
    
    void FullKnnSearch(NGS *q, int k, std::map<double, std::vector<NGS>, std::greater<double>> *KNN, int *sc, int *absc, bool sen);
    
    
    void RangeSearch(NGS *q, double t, std::map<double, std::vector<NGS>, std::greater<double>> *NN, int *sc, int *absc);
    
    void FullRangeSearch(NGS *q, double t, std::map<double, std::vector<NGS>, std::greater<double>> *NN, int *sc, int *absc);
    
    void saveIndex(std::ofstream *sout);
    
    int loadVPTreee(std::ifstream *sin, std::vector<NGS> *dataPoints);
    
    int loadLevDetails(std::string lineContents);
    
    
    long double totalTau();
    
    double averageTau();
    
    double levAvTau(int y);
    
    void printAllTau();
    
    int NodesNum();
};
#endif /* VPTree_hpp */
