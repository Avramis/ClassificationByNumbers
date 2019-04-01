//
//  FastqParser.cpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#include "FastqParser.hpp"
// Called
FastqParser::FastqParser(std::string filepath, std::string repmeth, int refkmerlen, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences, int sa, int ss, bool conper, int rcount){
    std::ifstream testpecase(filepath);
    std::string lineContents;
    std::string name1 = "", name2 = "";
    
    std::ifstream fastqFile(filepath);
    lineContents = "";
    std::string read, uid, qualstr;
    int count = 0, count1 = rcount, readlength = 0;
    std::vector<std::vector<double>> saxrep ,saxtra;
    int  c = 0;
    while(!fastqFile.eof()){
        getline(fastqFile, lineContents);
        std::stringstream readstream(lineContents);
        if(count == 0 && lineContents[0] == '@'){
            uid = *new std::string;
            std::getline (readstream, uid, '\n');
            uid.erase(0,1);
            count++;
        }
        else{
            if(count == 1)
            {
                readlength = *new int;
                read = *new std::string;
                readstream >> read;
                readlength = (int)read.size();
                count++;
            }
            else if(count == 2){
                count++;
            }
            else if(count == 3){
                qualstr = *new std::string;
                readstream >> qualstr;
                if(conper == true){
                    int rsize = (int)read.size()/2;
                    
                    int reqlen = rsize;
                    if(reqlen > 0)
                    {
                        if(reqlen >= refkmerlen){
                            reqlen = refkmerlen;
                        }
                    }
                    NGS n1(uid+"/1", read.substr(0, rsize), qualstr.substr(0, rsize));
                    NGS n2(uid+"/2", read.substr(rsize), qualstr.substr(rsize, rsize));
                    n1.setIdx(c);
                    n2.setIdx(c);
                    rep.createRepresentation(repmeth, n1.returnRead(), n1.returnRep());
                    rep.createRepresentation(repmeth, n2.returnRead(), n2.returnRep());
                    tra.createTransformation(tranmeth, n1.returnRep(), n1.returnTra(), 0, reqlen,  tranlev);
                    tra.createTransformation(tranmeth, n2.returnRep(), n2.returnTra(), 0, reqlen,  tranlev);

                    
                    rep.createRepresentation(repmeth,n1.returnRvRead(), n1.returnRvRep());
                    tra.createTransformation(tranmeth, n1.returnRvRep(), n1.returnRvTra(), 0, reqlen,  tranlev);
                    rep.createRepresentation(repmeth,n2.returnRvRead(), n2.returnRvRep());
                    tra.createTransformation(tranmeth, n2.returnRvRep(), n2.returnRvTra(), 0, reqlen,  tranlev);
                    
                    (*NuSequences).push_back(n1);
                    (*NuSequences).push_back(n2);
                    c++;
                    
                }
                else{
                    int reqlen = readlength;
                    if(reqlen > 0)
                    {
                        if(reqlen >= refkmerlen){
                            reqlen = refkmerlen;
                        }
                    }
                    
                    NGS n(uid, read, qualstr);
                    n.setIdx(c);
                    rep.createRepresentation(repmeth, n.returnRead(), n.returnRep());
                    tra.createTransformation(tranmeth, n.returnRep(), n.returnTra(), 0, reqlen,  tranlev);
                    
                    rep.createRepresentation(repmeth, n.returnRvRead(), n.returnRvRep());
                    tra.createTransformation(tranmeth, n.returnRvRep(), n.returnRvTra(), 0, reqlen,  tranlev);
                    
                    (*NuSequences).push_back(n);
                    c++;
                }
                
                count = 0;
                count1++;
                if((count1%100000) == 0){
                    std::cout << "> Fastq read: " << count1 << "\n";
                }
            }
        }
    }
    fastqFile.close();
};

// Called
FastqParser::~FastqParser(){
    rep.~NucRepresentations();
    tra.~DataTransformations();
};
