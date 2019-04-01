//
//  ClassificationTree.hpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#ifndef ClassificationTree_hpp
#define ClassificationTree_hpp

#include <stdio.h>
#include <string>
#include <map>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "DataTransformations.hpp"
#include "NucRepresentations.hpp"
#include "VPTree.hpp"
#include "FastaParser.hpp"
#include "NucRepresentations.hpp"
#include "DataTransformations.hpp"
#include "DataProcessing.hpp"
#include "DistanceMethods.hpp"
#include "ReverseCompliment.hpp"
#include "ssw_cpp.h"

class ClassificationTree{
private:
    int Kmer, Comlvl;
    std::string Repmeth, Tranmet;
    bool Accum, Norm;
    ReverseCompliment RVC;
    VPTree VP;
    std::vector<std::pair<std::string, std::string>> RefList;
    std::vector<NGS> uniqueRefNGS;
    NucRepresentations rep;
    DataTransformations tra;
    DataProcessing DaPr;
    DistanceMethods DM;
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    
    void genUnique();
    
    int loadRef(std::ifstream *sin);
    
    int loadReads(std::ifstream *sin);
    
public:
    
    ClassificationTree();
    ClassificationTree(std::string fastadir);
    ClassificationTree(std::string fastadir, int kmer);
    ClassificationTree(std::string fastadir, int kmer, std::string rep, std::string tra, int lvl, bool accum, bool norm);
    
    void BuildTree(std::string fastadir, int kmer, std::string rep, std::string tra, int lvl, bool accum, bool norm);

    
    void setKmer(int kmer);
    
    void setRepresentation(std::string rep);
    
    void setTansformationApproach(std::string tra);
    
    void setCompressionLevel(int lvl);
    
    void setAccummulation(bool accum);
    
    void setNormalisation(bool norm);
    
    std::pair<bool, bool> returnAccumNorm();
    
    void KNNSearch(std::string query, int knn, std::vector<NGS> *KNNList, long long *ovsc, long long *ovas, bool sen);
    
    void RangeSearch(std::string query, double tau, std::vector<NGS> *NNList, long long *ovsc, long long *ovas);
    
    std::vector<std::pair<std::string, std::string>> *returnRefList();
    
    std::string returnRepmeth();
    
    void saveTree(std::string savedir);
    
    int loadTree(std::string loaddir);
    
    double TreeaverageTau();
    
    double TreelevAvTau(int y);
    
    void TreeprintAllTau();
    
    int TreeNodesNum();
    
    long double totalTau();
    
    ~ClassificationTree();
};
#endif /* ClassificationTree_hpp */
