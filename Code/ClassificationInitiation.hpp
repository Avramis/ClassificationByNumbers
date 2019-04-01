//
//  ClassificationInitiation.hpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#ifndef ClassificationInitiation_hpp
#define ClassificationInitiation_hpp

#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <iomanip>

#include "ClassificationTree.hpp"
#include "DistanceMethods.hpp"
#include "ReverseCompliment.hpp"
#include "ssw_cpp.h"

class ClassificationInitiation{
private:
    ClassificationTree *CVPTree;
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    NucRepresentations rep;
    DataTransformations tra;
    DataProcessing DaPr;
    DistanceMethods DM;
    ReverseCompliment RVC;
    std::vector<std::pair<std::string, std::string>> *RefList;
    std::string Repmeth;
    double maxswscore = 1.0;
    bool Accum, Norm;
    static bool sortNGSDiscore(NGS n, NGS m);
    
    static bool sortNGSSWscore(NGS n, NGS m);
    
    void eliminateDuplicateMatching(std::vector<NGS> *NNList);
    
    void EDEvaluation(std::vector<std::vector<double>> *Fwrep, std::vector<std::vector<double>> *Rvrep, std::vector<NGS> *KNNList, bool r);
    
    void DWTEvaluation(std::vector<std::vector<double>> *Fwrep, std::vector<std::vector<double>> *Rvrep, std::vector<NGS> *KNNList, bool r);
    
    void SWEvaluation(std::string Fwread, std::string Rvread, std::vector<NGS> *KNNList, bool r);
    
    void SWFinalEvaluation(std::string Fwread, std::vector<NGS> *KNNList);
    
    void indexRead(std::string searchmeth, std::string read, std::vector<NGS> *KNNList, bool r);

    int testfile(std::string filepath);

    void fastaParserfile(ClassificationTree *Tree, std::string shorreadsdir, std::string outdir, std::string searchmeth, int knn, bool r, bool repall, bool sen);
    
    void fastaParserfile(ClassificationTree *Tree, std::string shorreadsdir, std::string outdir, std::string searchmeth, double tau, bool r, bool repall);
    
    void fastqParserfile(ClassificationTree *Tree, std::string shorreadsdir, std::string outdir, std::string searchmeth, int knn, bool r, bool repall, bool sen);
    
    void fastqParserfile(ClassificationTree *Tree, std::string shorreadsdir, std::string outdir, std::string searchmeth, double tau, bool r, bool repall);
    
    std::string returnAligstring(std::string read, std::string ref, std::string cigar, int refstart, bool dir);
    
    std::string DecomposeCigar(std::string cigar);
    
    std::string adjustString(int i, std::string title);
    
    void setAccumNorm(std::pair<bool, bool>);
    
public:
    ClassificationInitiation();
    
    ClassificationInitiation(ClassificationTree *Tree, std::string shorreadsdir, std::string outdir, std::string searchmeth, int knn, bool r, bool repall, bool sen);
    
    ClassificationInitiation(ClassificationTree *Tree, std::string shorreadsdir, std::string outdir, std::string searchmeth, double tau, bool r, bool repall);
    
    void setTree(ClassificationTree *Tree);
    
    ~ClassificationInitiation();
};
#endif /* ClassificationInitiation_hpp */
