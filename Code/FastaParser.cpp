//
//  FastaParser.cpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#include "FastaParser.hpp"
// Called
FastaParser::FastaParser(std::string filepath, std::vector<std::pair<std::string, std::string>> *ReferenceList){
    std::ifstream testpecase(filepath);
    std::string lineContents;
    std::string name1 = "";
    
    std::ifstream fastqFile(filepath);
    lineContents = "";
    std::string read, uid, qualstr;
    int count = 0, count1 = 0, readlength = 0;
    std::vector<std::vector<double>> saxrep ,saxtra;
    bool ident = false;
    while(!fastqFile.eof()){
        
        getline(fastqFile, lineContents);
        std::stringstream readstream(lineContents);
        if(count == 0 && lineContents[0] == '>'){
            if(ident == true){
                (*ReferenceList).push_back(std::make_pair(uid, read));
            }
            uid = *new std::string;
            std::getline (readstream, uid, '\n');
            uid.erase(0,1);
            readlength = 0;
            count++;
            read = *new std::string;;
            ident = false;
            count1++;
            if((count1%100000) == 0){
                std::cout << "> Fasta read: " << count1 << "\n";
            }
        }
        else{
            ident = true;
            std::string curstr;
            readstream >> curstr;
            curstr.erase(std::remove(curstr.begin(), curstr.end(), '\n'), curstr.end());
            readlength += (int)curstr.size();
            read+=curstr;
            count=0;
            
        }
    }
    fastqFile.close();
    
    if(ident == true){
        (*ReferenceList).push_back(std::make_pair(uid, read));
    }
};


// Called
FastaParser::~FastaParser(){
    rep.~NucRepresentations();
    tra.~DataTransformations();
};

