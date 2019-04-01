//
//  ClassificationTree.cpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#include "ClassificationTree.hpp"

void ClassificationTree::genUnique(){
    std::map <std::string, std::vector<std::pair<std::pair<int, int>, bool>>> uniquekmers;
    for(int i = 0; i < (int)RefList.size(); i++){
        std::string reverscompliment = RVC.returnReversCompliment(RefList[i].second);
        long long refsize = (long long)reverscompliment.size();
        if((int)RefList[i].second.size() > Kmer){
            for (long long j = 0; j < ((long long)RefList[i].second.size() - Kmer) + 1; j++){
                uniquekmers[RefList[i].second.substr(j, Kmer)].push_back(std::make_pair(std::make_pair(i, j), false));
                uniquekmers[reverscompliment.substr((refsize - j)-Kmer, Kmer)].push_back(std::make_pair(std::make_pair(i, j), true));
            }
        }
        else{
            uniquekmers[RefList[i].second].push_back(std::make_pair(std::make_pair(i, 0), false));
            uniquekmers[reverscompliment].push_back(std::make_pair(std::make_pair(i, 0), true));
        }
    }
    std::map <std::string, std::vector<std::pair<std::pair<int, int>, bool>>>::iterator it;
    int count = 0;
    for (it = uniquekmers.begin(); it != uniquekmers.end(); it++){
        NGS n(uniquekmers[it->first][0].first.first, uniquekmers[it->first][0].first.second);
        n.setAlDire(uniquekmers[it->first][0].second);
        for (long long i = 1; i < (long long)uniquekmers[it->first].size(); i++){
            n.addCoordinates(uniquekmers[it->first][i].first.first, uniquekmers[it->first][i].first.second);
            n.setAlDire(uniquekmers[it->first][i].second);
        }
        rep.createRepresentation(Repmeth,it->first, n.returnRep());
        int k = Kmer;
        if(k > (int)(*n.returnRep())[0].size()){
            k = (int)(*n.returnRep())[0].size();
        }
        if(Accum == true){
            DaPr.repAccumulation(n.returnRep(), 0, k);
        }
        if(Norm == true){
            DaPr.repZnormalisation(n.returnRep(), 0, k);
        }
        tra.createTransformation(Tranmet, n.returnRep(), n.returnTra(), 0, Kmer,  Comlvl);
        
        n.setIdx(count);
        n.clearRep();
        uniqueRefNGS.push_back(n);
        count++;
    }
    std::map <std::string, std::vector<std::pair<std::pair<int, int>, bool>>>().swap(uniquekmers);
};

ClassificationTree::ClassificationTree(){};

ClassificationTree::ClassificationTree(std::string fastadir){
    
};

ClassificationTree::ClassificationTree(std::string fastadir, int kmer){
    setKmer(kmer);
};

ClassificationTree::ClassificationTree(std::string fastadir, int kmer, std::string rep, std::string tra, int lvl, bool accum, bool norm){
    setKmer(kmer);
    setRepresentation(rep);
    setTansformationApproach(tra);
    setCompressionLevel(lvl);
    setAccummulation(accum);
    setNormalisation(norm);
};

void ClassificationTree::BuildTree(std::string fastadir, int kmer, std::string rep, std::string tra, int lvl, bool accum, bool norm){
    setKmer(kmer);
    setRepresentation(rep);
    setTansformationApproach(tra);
    setCompressionLevel(lvl);
    setAccummulation(accum);
    setNormalisation(norm);
    FastaParser(fastadir, &RefList);
    genUnique();
    VP.TreeGeneration(&uniqueRefNGS);
};

void ClassificationTree::setKmer(int kmer){
    Kmer = kmer;
};

void ClassificationTree::setRepresentation(std::string rep){
    Repmeth = rep;
};

void ClassificationTree::setTansformationApproach(std::string tra){
    Tranmet = tra;
    
};

void ClassificationTree::setCompressionLevel(int lvl){
    Comlvl = lvl;
};

void ClassificationTree::setAccummulation(bool accum){
    Accum = accum;
};

void ClassificationTree::setNormalisation(bool norm){
    Norm = norm;
};

std::pair<bool, bool> ClassificationTree::returnAccumNorm(){
    return std::make_pair(Accum, Norm);
};

void ClassificationTree::KNNSearch(std::string query, int knn, std::vector<NGS> *KNNList, long long *ovsc, long long *ovas, bool sen){
    int k = Kmer;
    if (k > (int)query.size()){
        k = (int)query.size();
    }
    
    NGS q(query, k);
    rep.createRepresentation(Repmeth, q.returnRead(), q.returnRep());
    if(Accum == true){
        DaPr.repAccumulation(q.returnRep(), 0, k);
    }
    if(Norm == true){
        DaPr.repZnormalisation(q.returnRep(), 0, k);
    }
    tra.createTransformation(Tranmet, q.returnRep(), q.returnTra(), 0, k,  Comlvl);
    int sc = 0;
    int absc = 0;
    std::map<double, std::vector<NGS>, std::greater<double>> Neighbours;
    VP.FullKnnSearch(&q, knn, &Neighbours, &sc, &absc, sen);
    (*ovsc) += (long long) sc;
    (*ovas) += (long long) absc;
    
    std::map<double, std::vector<NGS>, std::greater<double>>::iterator it;
    for (it = Neighbours.begin(); it != Neighbours.end(); it++){
        for (int i = 0; i < (int)Neighbours[it->first].size(); i++){
            for (int j = 0; j < (int)Neighbours[it->first][i].returnCoordinates().size(); j++){
                NGS n;
                n.setAlDire(Neighbours[it->first][i].returnAldirectionList()[j]);
                
                if(Neighbours[it->first][i].returnAldirectionList()[j] == false){
                    n.addCoordinates(Neighbours[it->first][i].returnCoordinates()[j].first, Neighbours[it->first][i].returnCoordinates()[j].second);
                    if(Neighbours[it->first][i].returnCoordinates()[j].second < 0){
                        std::cout << "";
                    }
                }
                else{
                    int alp = Neighbours[it->first][i].returnCoordinates()[j].second - ((int)q.returnRead().size() - Kmer);
                    n.addCoordinates(Neighbours[it->first][i].returnCoordinates()[j].first,  alp);
                    if(alp < 0){
                        std::cout << "";
                    }
                }
                (*KNNList).push_back(n);
            }
        }
        
    }
};

void ClassificationTree::RangeSearch(std::string query, double tau, std::vector<NGS> *NNList, long long *ovsc, long long *ovas){
    int k = Kmer;
    if (k > (int)query.size()){
        k = (int)query.size();
    }
    
    NGS q(query, k);
    rep.createRepresentation(Repmeth, q.returnRead(), q.returnRep());
    if(Accum == true){
        DaPr.repAccumulation(q.returnRep(), 0 ,k);
    }
    if(Norm == true){
        DaPr.repZnormalisation(q.returnRep(), 0 ,k);
    }
    tra.createTransformation(Tranmet, q.returnRep(), q.returnTra(), 0, k,  Comlvl);
    int sc = 0;
    int absc = 0;

    std::map<double, std::vector<NGS>, std::greater<double>> Neighbours;
    VP.FullRangeSearch(&q, tau, &Neighbours, &sc, &absc);
    (*ovsc) += (long long) sc;
    (*ovas) += (long long) absc;
    std::map<double, std::vector<NGS>, std::greater<double>>::iterator it;
    for (it = Neighbours.begin(); it != Neighbours.end(); it++){
        for (int i = 0; i < (int)Neighbours[it->first].size(); i++){
            for (int j = 0; j < (int)Neighbours[it->first][i].returnCoordinates().size(); j++){
                NGS n;
                n.setAlDire(Neighbours[it->first][i].returnAldirectionList()[j]);
                
                if(Neighbours[it->first][i].returnAldirectionList()[j] == false){
                    n.addCoordinates(Neighbours[it->first][i].returnCoordinates()[j].first, Neighbours[it->first][i].returnCoordinates()[j].second);
                    if(Neighbours[it->first][i].returnCoordinates()[j].second < 0){
                        std::cout << "";
                    }
                }
                else{
                    int alp = Neighbours[it->first][i].returnCoordinates()[j].second - ((int)q.returnRead().size() - Kmer);
                    n.addCoordinates(Neighbours[it->first][i].returnCoordinates()[j].first,  alp);
                    if(alp < 0){
                        std::cout << "";
                    }
                }
                (*NNList).push_back(n);
            }
        }
    }
};




std::vector<std::pair<std::string, std::string>> * ClassificationTree::returnRefList(){
    return &RefList;
};

std::string ClassificationTree::returnRepmeth(){
    return Repmeth;
};

void ClassificationTree::saveTree(std::string savedir){
    savedir+=".vptstruct";
    std::ofstream sout(savedir);
    sout << "InitiateSaveTree{\n";
    sout << "RefList_{\n";
    for(int i = 0; i < (int)RefList.size(); i++){
        sout << "{" << RefList[i].first << "^*^" << RefList[i].second << "}\n";
    }
    sout << "}_RefList\n";
    sout << "Repmeth_{" << Repmeth << "}_Repmeth\n";
    sout << "Tranmet_{" << Tranmet << "}_Tranmet\n";
    sout << "Kmer_{"<< Kmer << "}_Kmer\n";
    sout << "Comlvl_{"<< Comlvl << "}_Comlvl\n";
    if(Accum == false){
        sout << "Accum_{False}_Accum\n";
    }
    else{
        sout << "Accum_{True}_Accum\n";
    }
    if(Norm == false){
        sout << "Norm_{False}_Norm\n";
    }
    else{
        sout << "Norm_{True}_Norm\n";
    }
    
    for(int i = 0; i < (int)uniqueRefNGS.size(); i++){
        uniqueRefNGS[i].setIdx(i);
    }
    
    sout << "uniqueKmersSize_{" <<std::setprecision(100) << (long long)uniqueRefNGS.size() << "}_uniqueKmersSize\n";
    sout << "uniqueRefNGS_{\n";
    for(long long i = 0; i < (long long)uniqueRefNGS.size(); i++){
        uniqueRefNGS[i].savePoint(&sout);
    }
    sout << "}_uniqueRefNGS\n";
    sout << "VPTree_{\n";
    VP.saveIndex(&sout);
    sout << "}_VPTree\n";
    sout << "}SaveTreeCompleted\n";
    sout.close();
};


int ClassificationTree::loadTree(std::string loaddir){
    int tl = 0;
    size_t didx = loaddir.find_last_of(".");
    if(didx == std::string::npos){
        std::cout << loaddir << " has an invalide format for loading VPTree.\n";
        tl = 1;
        return tl;
    }
    else{
        if(loaddir.substr(didx) != ".vptstruct"){
            std::cout << loaddir << " has an invalide format for loading VPTree.\n";
            tl = 1;
            return tl;
        }
    }
    std::ifstream sin(loaddir);
    std::string lineContents;
    while(!sin.eof()){
        getline(sin, lineContents);
        didx = lineContents.find_first_of("{");
        if(lineContents.substr(0,didx+1) == "RefList_{"){
            std::cout << "> Load reference(s).\n";
            tl = loadRef(&sin);
        }
        else if(lineContents.substr(0,didx+1) == "Repmeth_{"){
            lineContents = lineContents.substr(didx+1);
            didx = lineContents.find_last_of("}_Repmeth");
            Repmeth = lineContents.substr(0, didx-8);
        }
        else if(lineContents.substr(0,didx+1) == "Tranmet_{"){
            lineContents = lineContents.substr(didx+1);
            didx = lineContents.find_last_of("}_Tranmet");
            Tranmet = lineContents.substr(0, didx-8);
        }
        else if(lineContents.substr(0,didx+1) == "Kmer_{"){
            lineContents = lineContents.substr(didx+1);
            didx = lineContents.find_last_of("}_Kmer");
            Kmer = stoi(lineContents.substr(0, didx-5));
        }
        else if(lineContents.substr(0,didx+1) == "Comlvl_{"){
            lineContents = lineContents.substr(didx+1);
            didx = lineContents.find_last_of("}_Comlvl");
            Comlvl = stoi(lineContents.substr(0,didx-7));
        }
        else if(lineContents.substr(0,didx+1) == "Accum_{"){
            lineContents = lineContents.substr(didx+1);
            didx = lineContents.find_last_of("}_Accum");
            if(lineContents.substr(0,didx-6) == "False"){
                Accum = false;
            }
            else{
                Accum = true;
            }
        }
        else if(lineContents.substr(0,didx+1) == "Norm_{"){
            lineContents = lineContents.substr(didx+1);
            didx = lineContents.find_last_of("}_Norm");
            if(lineContents.substr(0,didx-5) == "False"){
                Norm = false;
            }
            else{
                Norm = true;
            }
        }
        else if(lineContents.substr(0,didx+1) == "uniqueKmersSize_{"){
            lineContents = lineContents.substr(didx+1);
            didx = lineContents.find_last_of("}_uniqueKmersSize");
            std::vector<NGS>(stoll(lineContents.substr(0,didx-16)), NGS()).swap(uniqueRefNGS);
        }
        else if(lineContents.substr(0,didx+1) == "uniqueRefNGS_{"){
            std::cout << "> Load unique kmers.\n";
            loadReads(&sin);
        }
        else if(lineContents.substr(0,didx+1) == "VPTree_{"){
            std::cout << "> Load tree.\n";
            VP.loadVPTreee(&sin, &uniqueRefNGS);
        }
    };
    return tl;
};


int ClassificationTree::loadRef(std::ifstream *sin){
    int tl = 0;
    std::string lineContents;
    size_t sidx;
    while(!(*sin).eof()){
        getline((*sin), lineContents);
        sidx = lineContents.find_last_of("}");
        if(sidx != std::string::npos){
            if(lineContents.substr(sidx) == "}_RefList"){
                tl = 0;
                break;
            }
            else{
                lineContents = lineContents.substr(0, sidx);
                sidx = lineContents.find_first_of("{");
                lineContents = lineContents.substr(sidx+1);
                sidx = lineContents.find("^*^");
                RefList.push_back(std::make_pair(lineContents.substr(0, sidx), lineContents.substr(sidx+3)));
            }
        }
        else{
            tl = 1;
            break;
        }
    }
    return tl;
};


int ClassificationTree::loadReads(std::ifstream *sin){
    int tl = 0;
    std::string lineContents;
    size_t sidx;
    while(!(*sin).eof()){
        getline((*sin), lineContents);
        sidx = lineContents.find_last_of("}");
        if(sidx != std::string::npos){
            if(lineContents.substr(sidx) == "}_uniqueRefNGS"){
                tl = 0;
                break;
            }
            else{
                lineContents = lineContents.substr(0, sidx);
                sidx = lineContents.find_first_of("{");
                lineContents = lineContents.substr(sidx+1);
                uniqueRefNGS[stoll(lineContents)].loadPoint(sin);
            }
        }
        else{
            tl = 1;
            break;
        }
    }
    
    return tl;
};


double ClassificationTree::TreeaverageTau(){
    return VP.averageTau();
};

double ClassificationTree::TreelevAvTau(int y){
    return VP.levAvTau(y);
};

void ClassificationTree::TreeprintAllTau(){
    VP.printAllTau();
};

int ClassificationTree::TreeNodesNum(){
    return VP.NodesNum();
};

long double ClassificationTree::totalTau(){
    return VP.totalTau();
};

ClassificationTree::~ClassificationTree(){};
