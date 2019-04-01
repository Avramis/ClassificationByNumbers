//
//  NGS.hpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#ifndef NGS_hpp
#define NGS_hpp

#include <stdio.h>
#include <vector>
#include<string>
#include <map>
#include <ctype.h>
#include <algorithm> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

class NGS{
private:
    std::vector<std::vector<double>> nFwrep, nFwtra,  nRvrep, nRvtra;
    std::vector <std::string*> gid;
    std::string *refread;
    bool rdir = false;
    std::vector<std::pair<int, int>> coordinates;
    std::vector<bool> alDir;
    std::vector<std::string> cigars, alRef;
    std::string rid, read, qual;
    int kmer, idx;
    double alSc = 0;
    std::vector<double> alScVec;
    std::vector<int> endpos;
    
public:
    NGS();
    //For short reads
    NGS(std::string u, std::string r);
    NGS(std::string u, std::string r, std::string q);
    NGS(std::string r, int k);
    NGS(std::string u, std::string r, std::string q, int k);
    
    //For reference dubliacted reads
    NGS(int g, int s);
    NGS(int g, int s, int k);
    
    NGS(std::string *u, int g, int s);
    NGS(std::string *u, int g, int s, int k);
    
    NGS(int g, int s, std::string *r);
    NGS(int g, int s, int k, std::string *r);
    
    NGS(std::string *u, std::string *r, int g, int s);
    NGS(std::string *u, std::string *r, int g, int s, int k);
    
    
    void setSRid(std::string s);
    void setSread(std::string s);
    void setQual(std::string s);
    void setKmer(int i);
    void setIdx(int i);
    
    
    void setAlDire(bool d);
    void setAlDire(int i, bool d);

    
    void setAldetails(int i, int j, std::string s, bool b);
    void setAldetails(int i, int j, std::string s, std::string r);
    void setAldetails(int i, int j, std::string s, std::string r, bool d);
    
    void setRefreafp(std::string *r);
    
    void addGnp(std::string *s);
    
    void addCoordinates(int i, int p);
    
    void clearRep();
    void clearTra();
    void clearEndpos();
    
    
    void deleteRead();
    void deleteUid();
    
    void setAlscore(double i);
    
    
    void setAlVecscore(double i);
    void setEnpos(int i);
    
    void clearAldetails();
    
    void setCigar(std::string c);
    void setAlRef(std::string s);
    void editCigar(int i, std::string c);
    void editCoordinates(int i, int c);
    
    
    void reverseCigar(int i);
    
    std::string returnReverseCigar(int i);
    
    std::string returnRead();
    
    std::string returnQual();
    std::string returnSRid();
    
    int returnKmer();
    
    std::vector<std::vector<double>> *returnRep();
    std::vector<std::vector<double>> *returnTra();
    
    std::vector<std::vector<double>> *returnRvRep();
    std::vector<std::vector<double>> *returnRvTra();
    
    std::string* returnrefReads();
    std::vector<std::string*> returnrefUID();
    
    std::string *returnSRidP();
    std::string *returnUIP();
    
    double returnAlscore();
    
    std::vector<double> returnAlVec();
    
    std::vector<std::pair<int, int>> returnCoordinates();
    std::vector <std::string*> returnGID();
    std::vector <std::string> returnCigar();
    std::vector <std::string> returnAlref();
    
    int returnIdx();
    
    std::string returnRvRead();
    
    char NucleotideComparison(char x);
    
    void setDirection(bool i);
    bool returnDirection();
    
    bool returnAldirection(int i);
    std::vector<bool>returnAldirectionList();
    std::vector<int>* returnEndpos();
    void savePoint(std::ofstream *sout);
    
    int loadPoint(std::ifstream *sin);
    
    void loadFwTran(std::string lineContents);
    
    void loadCoordiantes(std::string lineContents);
    
    void loadDiretrions(std::string lineContents);
};

#endif /* NGS_hpp */

