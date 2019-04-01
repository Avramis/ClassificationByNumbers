//
//  DataProcessing.cpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#include "DataProcessing.hpp"
DataProcessing::DataProcessing(){};

void DataProcessing::repAccumulation(std::vector<std::vector<double>> *rep){
    for (int i = 0; i < (int)rep->size(); i++){
        for (int j = 1; j < (int)rep->at(i).size(); j++){
            (*rep)[i][j] += (*rep)[i][j-1];
        }
    }
};

void DataProcessing::repAccumulation(std::vector<std::vector<double>> *rep, int s, int l){
    for (int i = 0; i < (int)rep->size(); i++){
        for (int j = s+1; j < s + l; j++){
            (*rep)[i][j] += (*rep)[i][j-1];
        }
    }
};


void DataProcessing::repZnormalisation(std::vector<std::vector<double>> *rep){
    double meanv, stvar;
    int rsize = (int)rep->at(0).size();
    double dds = (double) rsize;
    
    for (int i = 0; i < (int)rep->size(); i++){
        meanv = 0.0;
        stvar = 0.0;
        for (int j = 0; j < rsize; j++){
            meanv += (*rep)[i][j]/dds;
        }
        if(meanv != meanv){
            meanv = 0;
        }
        
        for (int j = 0; j < rsize; j++){
            stvar += ((pow((*rep)[i][j] - meanv , 2.0))/dds);
        }
        if(stvar != stvar){
            stvar = 1.0;
        }
        for (int j = 0; j < rsize; j++){
            (*rep)[i][j] = ((*rep)[i][j]-meanv)/stvar;
        }
    }
};

void DataProcessing::repZnormalisation(std::vector<std::vector<double>> *rep, int s, int l){
    double meanv, stvar;
    double dds = (double) l;
    
    for (int i = 0; i < (int)rep->size(); i++){
        meanv = 0.0;
        stvar = 0.0;
        for (int j = s + 1; j < s + l; j++){
            meanv += (*rep)[i][j]/dds;
        }
        if(meanv != meanv){
            meanv = 0;
        }
        
        //for (int j = 0; j < rsize; j++){
        for (int j = s + 1; j < s + l; j++){
            stvar += ((pow((*rep)[i][j] - meanv , 2.0))/dds);
        }
        if(stvar != stvar){
            stvar = 1.0;
        }
        for (int j = s; j < s + l; j++){
            (*rep)[i][j] = ((*rep)[i][j]-meanv)/stvar;
        }
    }
};

DataProcessing::~DataProcessing(){};
