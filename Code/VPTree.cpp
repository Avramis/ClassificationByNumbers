//
//  VPTree.cpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#include "VPTree.hpp"

VPTree::VPTree(){};

void VPTree::TreeGeneration(std::vector<NGS> *dataPoints){
    //buildTree(&vp, (int)dataPoints->size(), dataPoints);
    buildTree(&vp, (int)dataPoints->size(), dataPoints, &levDetails);
    
    
    for(std::map<int, std::pair<int, double>, std::less<int>>::iterator it = levDetails.begin(); it != levDetails.end(); it++){
        levDetails[it->first].second = levDetails[it->first].second / ((double)levDetails[it->first].first);
    }
};

void VPTree::buildTree(Node *vp, int n, std::vector<NGS> *dataPoints){
    vp->generateNode(dataPoints, 0, n, true);
};

void VPTree::buildTree(Node *vp, int n, std::vector<NGS> *dataPoints, std::map<int, std::pair<int, double>> *levdetails){
    vp->generateNode(dataPoints, 0, n, true, 0 ,levdetails);
};

void VPTree::KnnSearch(NGS *q, int k, std::map<double, std::vector<NGS>, std::greater<double>> *KNN, int *sc, int *absc, bool sen){
    vp.KnnSearch(q, k, KNN, sc, absc, sen);
};

void VPTree::FullKnnSearch(NGS *q, int k,std::map<double, std::vector<NGS>, std::greater<double>> *KNN, int *sc, int *absc, bool sen){
    vp.initiateSearch(q, k, KNN, sc, absc, sen);
};



void VPTree::RangeSearch(NGS *q, double t, std::map<double, std::vector<NGS>, std::greater<double>> *NN, int *sc, int *absc){
    vp.RangeSearch(q, t, NN, sc, absc);
    
};

void VPTree::FullRangeSearch(NGS *q, double t, std::map<double, std::vector<NGS>, std::greater<double>> *NN, int *sc, int *absc){
    vp.initiateRangeSearch(q, t, NN, sc, absc);
};


void VPTree::saveIndex(std::ofstream *sout){
    int sp = 100;
    (*sout) << "levDet_{";
    for (std::map<int, std::pair<int, double>>::iterator it = levDetails.begin(); it != levDetails.end(); it++){
        (*sout) << std::setprecision(sp)  << it->first << "|" << std::setprecision(sp)  << levDetails[it->first].first << "|" << std::setprecision(sp)  <<  levDetails[it->first].second << "?|?";
    }
    (*sout) << "}_levDet\n";
    vp.saveNode(sout);
};

int VPTree::loadVPTreee(std::ifstream *sin, std::vector<NGS> *dataPoints){
    int t = 0;
    std::string lineContents = "";
    size_t didx;
    while (!(*sin).eof()){
        getline((*sin), lineContents);
        didx = lineContents.find_first_of("{");
        if(lineContents == "}_VPTree"){
            t = 0;
             break;
        }
        else{
            if(lineContents.substr(0,didx+1) == "levDet_{"){
                loadLevDetails(lineContents);
            }
            else if(lineContents.substr(0,didx+1) == "Node_S{"){
                t = vp.loadNode(sin, dataPoints);
            }
    }
    }
    return t;
};


int VPTree::loadLevDetails(std::string lineContents){
    int t = 0;
    size_t didx = lineContents.find_first_of("{");;
    lineContents = lineContents.substr(didx+1);
    didx = lineContents.find_last_of("}_levDet");
    lineContents = lineContents.substr(0, didx-7);
    didx = lineContents.find_first_of("?");
    while (didx != std::string::npos){
        std::string subseq = lineContents.substr(0, didx) + "|";
        lineContents = lineContents.substr(didx + 3);
        int x, y;
        double t;
        didx = subseq.find_first_of("|");
        x = stoi(subseq.substr(0, didx));
        subseq = subseq.substr(didx+1);
        didx = subseq.find_first_of("|");
        y= stoi(subseq.substr(0, didx));
        subseq = subseq.substr(didx+1);
        didx = subseq.find_first_of("|");
        t = stod(subseq.substr(0, didx));
        levDetails[x] = std::make_pair(y, t);
        didx = lineContents.find_first_of("?");
    }
    return t;
};



long double VPTree::totalTau(){
    long double x = 0.00;
    vp.totalTau(&x);
    return x;
};


double VPTree::averageTau(){
    double x = 0.0;
    vp.averageTau(&x, NodesNum());
    
    return x;
};

double  VPTree::levAvTau(int y){
    return levDetails[y].second;
};

void  VPTree::printAllTau(){
    for (std::map<int, std::pair<int, double>>::iterator it = levDetails.begin(); it != levDetails.end(); it++){
        std::cout << "> Level: " << it->first << " Nodes num: " << std::setprecision(100) << it->second.first << " | Average Tau: " << std::setprecision(100) <<  it->second.second << "\n";
    }
};

int VPTree::NodesNum(){
    int nn = 0;
    for (std::map<int, std::pair<int, double>>::iterator it = levDetails.begin(); it != levDetails.end(); it++){
        nn += it->second.first;
    }
    return nn;
};
