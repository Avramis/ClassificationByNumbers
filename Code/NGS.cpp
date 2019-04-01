//
//  NGS.cpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#include "NGS.hpp"
NGS::NGS(){};
//For short reads

NGS::NGS(std::string u, std::string r){
    setSRid(u);
    setSread(r);
};

NGS::NGS(std::string u, std::string r, std::string q){
    setSRid(u);
    setSread(r);
    setQual(q);
};

NGS::NGS(std::string u, std::string r, std::string q, int k){
    setSRid(u);
    setSread(r);
    setQual(q);
    setKmer(k);
};

NGS::NGS(std::string r, int k){
    setSread(r);
    setKmer(k);
};


//For reference dubliacted reads

NGS::NGS(int g, int s){
    addCoordinates(g, s);
};

NGS::NGS( int g, int s, int k){
    
};

NGS::NGS(std::string *u, int g, int s){
    addGnp(u);
    addCoordinates(g, s);
    
};

NGS::NGS(std::string *u, int g, int s, int k){
    addGnp(u);
    addCoordinates(g, s);
    setKmer(k);
};


NGS::NGS(int g, int s, std::string *r){
    setRefreafp(r);
    addCoordinates(g, s);
};

NGS::NGS(int g, int s, int k, std::string *r){
    setRefreafp(r);
    addCoordinates(g, s);
    setKmer(k);
};

NGS::NGS(std::string *u, std::string *r, int g, int s){
    addGnp(u);
    setRefreafp(r);
    addCoordinates(g, s);
};

NGS::NGS(std::string *u, std::string *r, int g, int s, int k){
    addGnp(u);
    setRefreafp(r);
    addCoordinates(g, s);
    setKmer(k);
};

void NGS::setSRid(std::string s){
    rid = s;
};

void NGS::setSread(std::string s){
    read = s;
};

void NGS::setQual(std::string s){
    qual = s;
};

void NGS::setKmer(int i){
    kmer = i;
};

void NGS::setIdx(int i){
    idx = i;
};

void NGS::setAldetails(int i, int j, std::string s, bool b){
    addCoordinates(i, j);
    cigars.push_back(s);
    setAlDire(b);
};

void NGS::setAlDire(bool d){
    alDir.push_back(d);
};

void NGS::setAlDire(int i, bool d){
    alDir[i] = d;
};

void NGS::setAldetails(int i, int j, std::string s, std::string r){
    addCoordinates(i, j);
    cigars.push_back(s);
    alRef.push_back(r);
};


void NGS::setAldetails(int i, int j, std::string s, std::string r, bool d){
    addCoordinates(i, j);
    cigars.push_back(s);
    alRef.push_back(r);
    setAlDire(d);
}


void NGS::addGnp(std::string *s){
    gid.push_back(s);
};

void NGS::addCoordinates(int i, int p){
    coordinates.push_back(std::make_pair(i, p));
};

void NGS::clearRep(){
    std::vector<std::vector<double>>().swap(nFwrep);
    std::vector<std::vector<double>>().swap(nRvrep);
};

void NGS::clearTra(){
    std::vector<std::vector<double>>().swap(nFwtra);
    std::vector<std::vector<double>>().swap(nRvtra);
};


void NGS::clearEndpos(){
    std::vector<int>().swap(endpos);
};

void NGS::setRefreafp(std::string *r){
    refread = r;
};

void NGS::deleteRead(){
    read = "";
};

void NGS::deleteUid(){
    rid = "";
};


void NGS::setAlscore(double i){
    alSc = i;
};



void NGS::setAlVecscore(double i){
    alScVec.push_back(i);
};


void NGS::setEnpos(int i){
    endpos.push_back(i);
};

std::vector<int>* NGS::returnEndpos(){
    return &endpos;
};

void NGS::clearAldetails(){
    alSc = 0;
    std::vector <double>().swap(alScVec);
    std::vector <std::string>().swap(alRef);
    std::vector<bool>().swap(alDir);
    std::vector<std::pair<int, int>>().swap(coordinates);
    std::vector<std::string>().swap(cigars);
    std::vector<bool>().swap(alDir);
    clearEndpos();
};

void NGS::setCigar(std::string c){
    cigars.push_back(c);
};

void NGS::setAlRef(std::string s){
    alRef.push_back(s);
};

void NGS::editCigar(int i, std::string c){
    cigars.at(i) = c;
};

void NGS::editCoordinates(int i, int c){
    coordinates.at(i).second = c;
};


void NGS::reverseCigar(int i){
    
    editCigar(i, returnReverseCigar(i));

    
};


std::string NGS::returnReverseCigar(int i){
    std::string c, ec;
    
    c = cigars.at(i);
    size_t spos;
    std::map<size_t, std::string, std::greater<size_t>> cidec;
    
    spos = c.size()-2;
    while (c.size() > 0  && spos > 0){
        if(!isdigit(*c.substr(spos,1).c_str())){
            ec = ec +  c.substr(spos+1);
            c = c.substr(0, spos+1);
            spos = c.size()-2;
        }
        else{
            
            spos--;
            
        }
    }
    ec = ec +  c;
    return ec;
};

std::string NGS::returnRead(){
    return read;
};

std::string NGS::returnSRid(){
    return rid ;
};

std::string NGS::returnQual(){
    return qual;
};

int NGS::returnKmer(){
    return kmer;
};

std::vector<std::vector<double>> * NGS::returnRep(){
    return &nFwrep;
};

std::vector<std::vector<double>> * NGS::returnRvRep(){
    return &nRvrep;
};


std::vector<std::vector<double>> * NGS::returnTra(){
    return &nFwtra;
};

std::vector<std::vector<double>> * NGS::returnRvTra(){
    return &nRvtra;
};


std::string* NGS::returnrefReads(){
    return refread;
};

std::vector<std::string*> NGS::returnrefUID(){
    return gid;
};

std::string *NGS::returnSRidP(){
    return &read;
};

std::string *NGS::returnUIP(){
    return &rid;
};

double NGS::returnAlscore(){
    return alSc;
};

std::vector<double> NGS::returnAlVec(){
    return alScVec;
};

std::vector<std::pair<int, int>> NGS::returnCoordinates(){
    return  coordinates;
};


std::vector <std::string*> NGS::returnGID(){
    return gid;
};


std::vector <std::string> NGS::returnAlref(){
    return alRef;
};

std::vector <std::string> NGS::returnCigar(){
    return cigars;
};

int NGS::returnIdx(){
    return idx;
};



std::string NGS::returnRvRead(){
    std::string rvread = read;
    std::reverse(rvread.begin(),rvread.end());
    for(int i = 0; i < (int) rvread.size(); i++){
        rvread[i]=NucleotideComparison(rvread[i]);
    }
    return rvread;
};


char NGS::NucleotideComparison(char x){
    if(toupper(x) == 'A' )
    {
        x='T';
    }
    else if(toupper(x) == 'T' || toupper(x) == 'U' )
    {
        x='A';
    }
    else if(toupper(x) == 'C' )
    {
        x='G';
    }
    else if(toupper(x) == 'G' )
    {
        x='C';
    }
    else if(toupper(x) == 'N' )
    {
        x='N';
    }
    else if(toupper(x) == 'M' )
    {
        x='K';
    }
    else if(toupper(x) == 'R' )
    {
        x='Y';
    }
    else if(toupper(x) == 'W' )
    {
        x='W';
    }
    else if(toupper(x) == 'S' )
    {
        x='S';
    }
    else if(toupper(x) == 'Y' )
    {
        x='R';
    }
    else if(toupper(x) == 'K' )
    {
        x='M';
    }
    else if(toupper(x) == 'V' )
    {
        x='B';
    }
    else if(toupper(x) == 'H' )
    {
        x='D';
    }
    else if(toupper(x) == 'D' )
    {
        x='H';
    }
    else if(toupper(x) == 'B' )
    {
        x='V';
    }
    return x;
};



void NGS::setDirection(bool i){
    rdir = i;
};

bool NGS::returnDirection(){
    return rdir;
};

bool NGS::returnAldirection(int i){
    return alDir.at(i);
};


std::vector<bool> NGS::returnAldirectionList(){
    return alDir;
};



void NGS::savePoint(std::ofstream *sout){
    int sp = 100;
    (*sout) << "NGSPoint_{" << idx <<"}\n";
    (*sout) << "idx_{" << idx << "}_idx\n";
    (*sout) << "nFwtra_{";
    (*sout) << "[" << std::setprecision(sp) << (int)nFwtra.size() << "|" << std::setprecision(sp) << (int)nFwtra[0].size() << "]";
    for (int i = 0; i < (int)nFwtra.size(); i++){
        for (int j = 0; j <(int)nFwtra[i].size()-1; j++){
            (*sout) << std::setprecision(sp) << nFwtra[i][j] << "|";
        }
        (*sout) << std::setprecision(sp) << nFwtra[i][(int)nFwtra[i].size()-1] << "?|?";
    }
    (*sout) << "}_nFwtra\n";
    (*sout) << "coordinates_{";
    for (int i = 0; i < (int) coordinates.size(); i++){
        (*sout) << coordinates[i].first << "|" << coordinates[i].second <<"?|?";
    }
    (*sout) << "}_coordinates\n";
    (*sout) << "alDir_{";
    for (int i = 0; i < (int)alDir.size(); i++){
        if(alDir[i] == false){
            (*sout) << "false";
        }
        else{
            (*sout) << "true";
        }
        (*sout) << "?|?";
    }
    (*sout) << "}_alDir\n";
    (*sout) << "}_NGSPoint\n";
};



int NGS::loadPoint(std::ifstream *sin){
    int t = 0;
    std::string lineContents = "";
    size_t didx;
    while (!(*sin).eof()){
        getline((*sin), lineContents);
        if(lineContents == "}_NGSPoint"){
            t = 0;
            break;
        }
        else{
            didx = lineContents.find_first_of("{");
            if(lineContents.substr(0,didx+1) == "idx_{"){
                lineContents = lineContents.substr(didx+1);
                didx = lineContents.find_last_of("}_idx");
                idx = stoi(lineContents.substr(0,didx-4));
            }
            else if(lineContents.substr(0,didx+1) == "nFwtra_{"){
                loadFwTran(lineContents);
            }
            else if(lineContents.substr(0,didx+1) == "coordinates_{"){
                loadCoordiantes(lineContents);
            }
            else if(lineContents.substr(0,didx+1) == "alDir_{"){
                loadDiretrions(lineContents);
            }
        }
    }
    return t;
};


void NGS::loadFwTran(std::string lineContents){
    size_t didx;
    int xe, ye;
    int xs = 0, ys = 0;
    didx = lineContents.find_first_of("{");
    lineContents = lineContents.substr(didx+1);
    didx = lineContents.find_last_of("}_nFwtra");
    lineContents = lineContents.substr(0, didx-7);
    didx = lineContents.find_first_of("]");
    std::string vecdim = lineContents.substr(1,didx-1);
    lineContents = lineContents.substr(didx + 1);
    didx = vecdim.find("|");
    xe = stoi(vecdim.substr(0,didx));
    ye = stoi(vecdim.substr(didx+1));
    std::vector<std::vector<double>>(xe, std::vector<double>(ye, 0.0)).swap(nFwtra);
    didx = lineContents.find_first_of("?");
    while (didx != std::string::npos){
        std::string subseq = lineContents.substr(0, didx) + "|";
        
        lineContents = lineContents.substr(didx + 3);
        didx = subseq.find_first_of("|");
        ys = 0;
        while (didx != std::string::npos){
            nFwtra[xs][ys] = stod(subseq.substr(0, didx));
            ys++;
            if(didx < subseq.size()-1){
                subseq = subseq.substr(didx+1);
            }
            else{
                subseq = "";
            }
            didx = subseq.find_first_of("|");
        }
        xs++;
        
        didx = lineContents.find_first_of("?");
    }
};

void NGS::loadCoordiantes(std::string lineContents){
    size_t didx;
    didx = lineContents.find_first_of("{");
    lineContents = lineContents.substr(didx+1);
    didx = lineContents.find_last_of("}_coordinates");
    lineContents = lineContents.substr(0, didx-12);
    didx = lineContents.find_first_of("?");
    while (didx != std::string::npos){
        std::string subseq = lineContents.substr(0, didx);
        lineContents = lineContents.substr(didx + 3);
        didx = subseq.find_first_of("|");
        coordinates.push_back(std::make_pair(stoi(subseq.substr(0, didx)), stoi(subseq.substr(didx+1))));
        didx = lineContents.find_first_of("?");
    }
};

void NGS::loadDiretrions(std::string lineContents){
    size_t didx;
    didx = lineContents.find_first_of("{");
    lineContents = lineContents.substr(didx+1);
    didx = lineContents.find_last_of("}_alDir");
    lineContents = lineContents.substr(0, didx-6);
    didx = lineContents.find_first_of("?");
    while (didx != std::string::npos){
        if(lineContents.substr(0, didx) == "false"){
            alDir.push_back(false);
        }
        else{
            alDir.push_back(true);
        }
        lineContents = lineContents.substr(didx + 3);
        didx = lineContents.find_first_of("?");
    }
};

