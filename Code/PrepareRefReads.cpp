//
//  PrepareRefReads.cpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#include "PrepareRefReads.hpp"
PrepareRefReads::PrepareRefReads(){};
PrepareRefReads::PrepareRefReads(std::vector<std::pair<std::string, std::string>> (*RefList), std::vector<NGS> *RefNGS, std::vector<NGS> *TreeData, int kmer){
    std::map<std::string, std::vector<std::pair<int, int>>> uniquereads;
    for(int i = 0; i < (int)(*RefList).size(); i++){
        std::map<std::string, std::vector<std::pair<int, int>>>::iterator it;
        NGS n((*RefList)[i].first, (*RefList)[i].second);
        n.setIdx(i);
        (*RefNGS).push_back(n);
        //for (int j = 0; j < (x - k) + 1; j++){
        if((int)(*RefList)[i].second.size() >= kmer){
            for (int j = 0; j < ((int)(*RefList)[i].second.size() - kmer) + 1; j++){
                uniquereads[(*RefList)[i].second.substr(j, kmer)].push_back(std::make_pair(i, j));
            }
        }
        else{
            uniquereads[(*RefList)[i].second].push_back(std::make_pair(i, 0));
        }
    }
};
