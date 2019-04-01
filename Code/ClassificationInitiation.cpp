//
//  ClassificationInitiation.cpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#include "ClassificationInitiation.hpp"

bool ClassificationInitiation::sortNGSDiscore(NGS n, NGS m){
    return n.returnAlscore() < m.returnAlscore();
};

bool ClassificationInitiation::sortNGSSWscore(NGS n, NGS m){
    return n.returnAlscore() > m.returnAlscore();
};

void ClassificationInitiation::eliminateDuplicateMatching(std::vector<NGS> *NNList){
    std::map<std::pair<std::pair<int, int>, bool>, NGS> eliminatedublicate;
    for (int i = 0; i < (int)(*NNList).size(); i++){
        eliminatedublicate[std::make_pair(std::make_pair((*NNList)[i].returnCoordinates()[0].first , (*NNList)[i].returnCoordinates()[0].second), (*NNList)[i].returnAldirection(0))] = (*NNList)[i];
    }
    std::vector<NGS>().swap((*NNList));
    
    for(std::map<std::pair<std::pair<int, int>, bool>, NGS>::iterator it = eliminatedublicate.begin(); it != eliminatedublicate.end(); it++){
        (*NNList).push_back(eliminatedublicate[it->first]);
    }
    
    std::map<std::pair<std::pair<int, int>, bool>, NGS>().swap(eliminatedublicate);
};


void ClassificationInitiation::EDEvaluation(std::vector<std::vector<double>> *Fwrep, std::vector<std::vector<double>> *Rvrep, std::vector<NGS> *KNNList, bool r){
    std::vector<NGS> tempKNN;
    tempKNN.swap((*KNNList));
    std::vector<NGS> ().swap((*KNNList));
    int qlen = (int)(*Fwrep)[0].size();
    
    std::vector<std::pair<std::string, std::string>> *RefList = (*CVPTree).returnRefList();
    setAccumNorm((*CVPTree).returnAccumNorm());
    
    std::vector<std::vector<double>> refrep;
    for (int i = 0; i < (int) tempKNN.size(); i++){
        int ridx = tempKNN[i].returnCoordinates()[0].first;
        int didx = tempKNN[i].returnCoordinates()[0].second;
        if(didx < 0){
            didx = 0;
        }
        std::string refread;
        if((*RefList)[ridx].second.substr(didx).size() > qlen)
        {
            refread = (*RefList)[ridx].second.substr(didx, didx+qlen);
        }
        else{
            refread = (*RefList)[ridx].second.substr(didx);
        }
        rep.createRepresentation(Repmeth, refread, &refrep);
        if(Accum == true){
            DaPr.repAccumulation(&refrep);
        }
        if(Norm == true){
            DaPr.repZnormalisation(&refrep);
        }
        
        double dist = 0;
        if(tempKNN[i].returnAldirectionList()[0] == false){
            
            dist = DM.returnDistance( &refrep, Fwrep);
        }
        else{
            dist = DM.returnDistance( &refrep, Rvrep);
        }
        tempKNN[i].setAlscore(dist);
        tempKNN[i].setAlVecscore(dist);
        std::vector<std::vector<double>>().swap(refrep);
    }
    
    sort(tempKNN.begin(), tempKNN.end(), sortNGSDiscore);
    
    std::vector<double> allscores;
    double aveScore = 0;;
    for (int i = 0; i < (int) tempKNN.size(); i++){
        aveScore += tempKNN[i].returnAlscore();
        allscores.push_back(tempKNN[i].returnAlscore());
    }
    aveScore = aveScore/ (double)tempKNN.size();
    
    if(r == true){
        for (int i = 0; i < (int) tempKNN.size(); i++){
            if(tempKNN[i].returnAlscore() <= aveScore){
                (*KNNList).push_back(tempKNN[i]);
            }
            else{
                break;
            }
        }
    }
    else{
        (*KNNList).swap(tempKNN);
    }
    std::vector<NGS> ().swap(tempKNN);
};

void ClassificationInitiation::DWTEvaluation(std::vector<std::vector<double>> *Fwrep, std::vector<std::vector<double>> *Rvrep, std::vector<NGS> *KNNList, bool r){
    std::vector<NGS> tempKNN;
    tempKNN.swap((*KNNList));
    std::vector<NGS> ().swap((*KNNList));
    int qlen = (int)(*Fwrep)[0].size();
    std::vector<std::vector<double>> refrep;
    int wi = 6;
    std::vector<double> allscores;
    double aveScore = 0;
    for (int i = 0; i < (int) tempKNN.size(); i++){
        int ridx = tempKNN[i].returnCoordinates()[0].first;
        int didx = tempKNN[i].returnCoordinates()[0].second;
        if(didx < 0){
            didx = 0;
        }
        std::string refread;
        if((*RefList)[ridx].second.substr(didx).size() > qlen)
        {
            refread = (*RefList)[ridx].second.substr(didx, didx + qlen);
        }
        else{
            refread = (*RefList)[ridx].second.substr(didx);
        }
        rep.createRepresentation(Repmeth, refread, &refrep);
        if(Accum == true){
            DaPr.repAccumulation(&refrep);
        }
        if(Norm == true){
            DaPr.repZnormalisation(&refrep);
        }
        
        double dist = 0;
        if(tempKNN[i].returnAldirectionList()[0] == false){
            
            dist = DM.returnDTW( &refrep, Fwrep, wi);
        }
        else{
            dist = DM.returnDTW( &refrep, Rvrep, wi);
        }
        tempKNN[i].setAlscore(dist);
        tempKNN[i].setAlVecscore(dist);
        aveScore += dist;
        std::vector<std::vector<double>>().swap(refrep);
    }
    sort(tempKNN.begin(), tempKNN.end(), sortNGSDiscore);
    
    aveScore = aveScore/ (double)tempKNN.size();
    
    if(r == true){
        for (int i = 0; i < (int) tempKNN.size(); i++){
            if(tempKNN[i].returnAlscore() <= aveScore){
                (*KNNList).push_back(tempKNN[i]);
            }
            else{
                break;
            }
        }
    }
    else{
        (*KNNList).swap(tempKNN);
    }
    std::vector<NGS> ().swap(tempKNN);
    
};

void ClassificationInitiation::SWEvaluation(std::string Fwread, std::string Rvread, std::vector<NGS> *KNNList, bool r){
    int qlen = (int) Fwread.size();
    std::vector<NGS> tempKNN;
    tempKNN.swap((*KNNList));
    std::vector<NGS> ().swap((*KNNList));
    
    for (int i = 0; i < (int) tempKNN.size(); i++){
        int ridx = tempKNN[i].returnCoordinates()[0].first;
        int didx = tempKNN[i].returnCoordinates()[0].second;
        if(didx < 0){
            didx = 0;
        }
        int subsize = 0;
        std::string refread;
        
        subsize = (int)(*RefList)[ridx].second.substr(didx).size();
        if(subsize > (qlen+10)){
            subsize = qlen+10;
        }
        
        if(didx - 10 > 0 ){
            didx = didx - 10;
            subsize = subsize + 10;
        }
        else{
            subsize = subsize + didx;
            didx = 0;
        }
        
        refread = (*RefList)[ridx].second.substr(didx, subsize);
        
        if(tempKNN[i].returnAldirectionList()[0] == false){
            aligner.Align(Fwread.c_str(),refread.c_str(), (int)refread.size(),filter, &alignment);
        }
        else{
            aligner.Align(Rvread.c_str(),refread.c_str(), (int)refread.size(),filter, &alignment);
        }
        tempKNN[i].setAlscore(alignment.sw_score / (double) qlen);
        tempKNN[i].setAlVecscore(alignment.sw_score / (double) qlen);
    }
    sort(tempKNN.begin(), tempKNN.end(), sortNGSSWscore);
    
    std::vector<double> allscores;
    double aveScore = 0;;
    for (int i = 0; i < (int) tempKNN.size(); i++){
        aveScore += tempKNN[i].returnAlscore();
        allscores.push_back(tempKNN[i].returnAlscore());
    }
    aveScore = aveScore/ (double)tempKNN.size();
    
    if(r == true){
        for (int i = 0; i < (int) tempKNN.size(); i++){
            if(tempKNN[i].returnAlscore() >= aveScore){
                (*KNNList).push_back(tempKNN[i]);
            }
            else{
                break;
            }
        }
    }
    else{
        (*KNNList).swap(tempKNN);
    }
    std::vector<NGS> ().swap(tempKNN);
    
};

void ClassificationInitiation::SWFinalEvaluation(std::string Fwread, std::vector<NGS> *KNNList){
    int qlen = (int) Fwread.size();
    std::vector<NGS> tempKNN;
    tempKNN.swap((*KNNList));
    std::vector<NGS> ().swap((*KNNList));
    int cadd = 0;
    for (int i = 0; i < (int) tempKNN.size(); i++){
        int ridx = tempKNN[i].returnCoordinates()[0].first;
        int didx = tempKNN[i].returnCoordinates()[0].second;
        if(didx < 0){
            didx = 0;
        }
        int subsize = 0;
        std::string refread;
        
        subsize = (int)(*RefList)[ridx].second.substr(didx).size();
        if(subsize > (qlen+10)){
            subsize = qlen+10;
        }
        
        if((didx - 10) > 0 ){
            didx = (didx - 10);
            subsize = (subsize + 10);
        }
        else{
            subsize = subsize + (didx+1);
            didx = 0;
        }
        
        refread = (*RefList)[ridx].second.substr(didx, subsize);
        if(tempKNN[i].returnAldirectionList()[0] == true){
            refread = RVC.returnReversCompliment(refread);
        }
        
        aligner.Align(Fwread.c_str(),refread.c_str(), (int)refread.size(),filter, &alignment);
        
        if((alignment.sw_score/(double)qlen) >= maxswscore){
            int endpos;
            if(tempKNN[i].returnAldirectionList()[0] == false){
                didx = didx + alignment.ref_begin;
            }
            else{
                didx = didx + (subsize  - alignment.ref_end);
            }
            endpos = didx + ((int)alignment.ref_end - (int)alignment.ref_begin);
            
            cadd++;
            (*KNNList).push_back(tempKNN[i]);
            (*KNNList)[cadd-1].setCigar(alignment.cigar_string);
            (*KNNList)[cadd-1].setAlscore(alignment.sw_score);
            (*KNNList)[cadd-1].editCoordinates(0, didx);
            (*KNNList)[cadd-1].setEnpos(endpos);
        }
        
    }
    sort((*KNNList).begin(), (*KNNList).end(), sortNGSSWscore);
    std::vector<NGS> ().swap(tempKNN);
    
};

int ClassificationInitiation::testfile(std::string filepath){
    int t = 3;
    std::string lineContents;
    std::ifstream fastqFile(filepath);
    lineContents = "";
    while(!fastqFile.eof()){
        getline(fastqFile, lineContents);
        std::stringstream readstream(lineContents);
        if(lineContents[0] == '@'){
            t = 0;
            break;
        }
        else if(lineContents[0] == '>'){
            t = 1;
            break;
        }
        else{
            t = 2;
        }
    }
    return t;
};


void ClassificationInitiation::indexRead(std::string searchmeth, std::string read, std::vector<NGS> *KNNList, bool r){
    std::string rvread = RVC.returnReversCompliment(read);
    std::vector<std::vector<double>> Fwrep, Rvrep;
    rep.createRepresentation(Repmeth, read, &Fwrep);
    rep.createRepresentation(Repmeth, rvread, &Rvrep);
    
    if(Accum == true){
        DaPr.repAccumulation(&Fwrep);
        DaPr.repZnormalisation(&Fwrep);
    }
    if(Norm == true){
        DaPr.repAccumulation(&Rvrep);
        DaPr.repZnormalisation(&Rvrep);
    }
    
    if(searchmeth == "ED"){
        EDEvaluation(&Fwrep, &Rvrep, KNNList, r);
    }
    else if(searchmeth == "DTW"){
        DWTEvaluation(&Fwrep, &Rvrep, KNNList, r);
    }
    else{
        SWEvaluation(read, rvread, KNNList, r);
    }
    rvread = "";
    std::vector<std::vector<double>>().swap(Fwrep);
    std::vector<std::vector<double>>().swap(Rvrep);
    
    SWFinalEvaluation(read, KNNList);
};


void ClassificationInitiation::fastaParserfile(ClassificationTree *Tree, std::string shorreadsdir, std::string outdir, std::string searchmeth, int knn, bool r, bool repall, bool sen){
    
    long long ovsc = 0;
    long long ovas = 0;
    int uncreads = 0;
    int creads = 0;
    
    std::string lineContents;
    std::string name1 = "";
    std::ifstream fastaFile(shorreadsdir);
    std::string details = outdir;
    size_t didx = outdir.find_last_of(".");
    std::string csvoutdir ="";
    if(didx != std::string::npos){
        csvoutdir = outdir.substr(0,didx) + ".tsv";
    }
    else{
        csvoutdir = outdir + ".tsv";
        outdir+=".txt";
    }
    
    std::ofstream sout(outdir);
    std::ofstream sout2(csvoutdir);
    sout2 << "Read_ID\tNum_of_Matches\tHighest_Score\tBest_Match\tDirection\tStart_Pos\n";
    
    lineContents = "";
    std::string read, uid;
    int count = 0, count1 = 0, readlength = 0;
    bool ident = false;
    
    while(!fastaFile.eof()){
        getline(fastaFile, lineContents);
        std::stringstream readstream(lineContents);
        if(count == 0 && lineContents[0] == '>'){
            if(ident == true){
                std::vector<NGS> KNNList;
                (*CVPTree).KNNSearch(read, knn, &KNNList, &ovsc, &ovas, sen);
                indexRead(searchmeth, read, &KNNList, r);
                
                eliminateDuplicateMatching(&KNNList);
                
                
                if((int)KNNList.size() > 0){
                    sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
                    sout << uid << "\n";
                    sout << "-------------------------------------------------------------------------------------------\n";
                    sout << "-------------------------------------------------------------------------------------------\n\n";
                    int bestscor = std::numeric_limits<int>::lowest();
                    int bestidx = 0;
                    for (int i = 0; i < (int)KNNList.size(); i++){
                        if(bestscor <  (int)KNNList[i].returnAlscore()){
                            bestscor =  (int)KNNList[i].returnAlscore();
                            bestidx = i;
                        }
                        int stapos = KNNList[i].returnCoordinates()[0].second;
                        int enpos = (*KNNList[i].returnEndpos())[0];
                        
                        sout << "-------------------------------------------------------------------------------------------\n";
                        sout << (*RefList)[KNNList[i].returnCoordinates()[0].first].first << "\n";
                        if(KNNList[i].returnAldirection(0) == false){
                            sout << "Direction: Forward\n";
                        }
                        else{
                            sout << "Direction: Reverse\n";
                        }
                        sout << "Range from " << std::to_string(stapos+1) << " to " << std::to_string (enpos+1) << "\n";
                        sout << "Alignment Score: " << std::to_string((int)KNNList[i].returnAlscore()) << "/" << std::to_string((int)read.size()*2)<< "\n";
                        sout << "Cigar: " << KNNList[i].returnCigar()[0] << "\n\n";
                        
                        
                        
                        if(KNNList[i].returnAldirection(0) == false){
                            std::string refread = (*RefList)[KNNList[i].returnCoordinates()[0].first].second.substr(stapos, (enpos - stapos)+1);
                            sout << returnAligstring(read, refread, KNNList[i].returnCigar()[0], stapos+1, false) << "\n";
                        }
                        else{
                            std::string refread = (*RefList)[KNNList[i].returnCoordinates()[0].first].second.substr(stapos-1, (enpos - stapos)+1);
                            
                            sout << returnAligstring(read, RVC.returnReversCompliment(refread), KNNList[i].returnCigar()[0], enpos + 1, true) << "\n";
                        }
                        sout << "-------------------------------------------------------------------------------------------\n\n";
                    }
                    sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
                    std::string bdir = "Forward";
                    if(KNNList[bestidx].returnAldirection(0) == true){
                        bdir = "Reverse";
                    }
                    sout2 << uid << "\t" << std::to_string((int) KNNList.size()) << "\t" << std::to_string(bestscor) <<"\t"<< (*RefList)[KNNList[bestidx].returnCoordinates()[0].first].first <<"\t" << bdir << "\t" << std::to_string(KNNList[bestidx].returnCoordinates()[0].second) <<"\n";
                    creads++;
                }
                else{
                    if(repall == true){
                        sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
                        sout << uid << "\n";
                        sout << "-------------------------------------------------------------------------------------------\n";
                        sout << "-------------------------------------------------------------------------------------------\n\n";
                        sout << "Unclassified read.\n";
                        sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
                        sout2 << uid << "\t" << 0 << "\t"  << "\t" << "\t" << "\t" <<"\n";
                    }
                    uncreads++;
                }
                
                
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
    fastaFile.close();
    
    if(ident == true){
        std::vector<NGS> KNNList;
        (*CVPTree).KNNSearch(read, knn, &KNNList, &ovsc, &ovas, sen);
        indexRead(searchmeth, read, &KNNList, r);
        
        eliminateDuplicateMatching(&KNNList);
        
        
        if((int)KNNList.size() > 0){
            sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
            sout << uid << "\n";
            sout << "-------------------------------------------------------------------------------------------\n";
            sout << "-------------------------------------------------------------------------------------------\n\n";
            int bestscor = std::numeric_limits<int>::lowest();
            int bestidx = 0;
            for (int i = 0; i < (int)KNNList.size(); i++){
                if(bestscor <  (int)KNNList[i].returnAlscore()){
                    bestscor =  (int)KNNList[i].returnAlscore();
                    bestidx = i;
                }
                int stapos = KNNList[i].returnCoordinates()[0].second;
                int enpos = (*KNNList[i].returnEndpos())[0];
                
                sout << "-------------------------------------------------------------------------------------------\n";
                sout << (*RefList)[KNNList[i].returnCoordinates()[0].first].first << "\n";
                if(KNNList[i].returnAldirection(0) == false){
                    sout << "Direction: Forward\n";
                }
                else{
                    sout << "Direction: Reverse\n";
                }
                sout << "Range from " << std::to_string(stapos+1) << " to " << std::to_string (enpos+1) << "\n";
                sout << "Alignment Score: " << std::to_string((int)KNNList[i].returnAlscore()) << "/" << std::to_string((int)read.size()*2)<< "\n";
                sout << "Cigar: " << KNNList[i].returnCigar()[0] << "\n\n";
                
                
                
                if(KNNList[i].returnAldirection(0) == false){
                    std::string refread = (*RefList)[KNNList[i].returnCoordinates()[0].first].second.substr(stapos, (enpos - stapos)+1);
                    sout << returnAligstring(read, refread, KNNList[i].returnCigar()[0], stapos+1, false) << "\n";
                }
                else{
                    std::string refread = (*RefList)[KNNList[i].returnCoordinates()[0].first].second.substr(stapos-1, (enpos - stapos)+1);
                    
                    sout << returnAligstring(read, RVC.returnReversCompliment(refread), KNNList[i].returnCigar()[0], enpos + 1, true) << "\n";
                }
                sout << "-------------------------------------------------------------------------------------------\n\n";
            }
            sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
            std::string bdir = "Forward";
            if(KNNList[bestidx].returnAldirection(0) == true){
                bdir = "Reverse";
            }
            sout2 << uid << "\t" << std::to_string((int) KNNList.size()) << "\t" << std::to_string(bestscor) <<"\t"<< (*RefList)[KNNList[bestidx].returnCoordinates()[0].first].first <<"\t" << bdir << "\t" << std::to_string(KNNList[bestidx].returnCoordinates()[0].second) <<"\n";
            creads++;
        }
        else{
            
            if(repall == true){
                sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
                sout << uid << "\n";
                sout << "-------------------------------------------------------------------------------------------\n";
                sout << "-------------------------------------------------------------------------------------------\n\n";
                sout << "Unclassified read.\n";
                sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
                sout2 << uid << "\t" << 0 << "\t"  << "\t" << "\t" << "\t" <<"\n";
            }
            uncreads++;
        }
        
    }
    sout.close();
    sout2 << "Total number of reads:\t" << std::to_string(creads) <<"\tMatched reads:\t" << std::to_string(creads - uncreads) << "\tUnmatched reads:\t" << std::to_string(uncreads)<<"\n";
    sout2.close();
    std::cout << "> Num of classified reads: " << creads << "\n";
    std::cout << "> Num of unclassified reads: " << uncreads << "\n";
};


void ClassificationInitiation::fastaParserfile(ClassificationTree *Tree, std::string shorreadsdir, std::string outdir, std::string searchmeth, double tau, bool r, bool repall){
    
    long long ovsc = 0;
    long long ovas = 0;
    int uncreads = 0;
    int creads = 0;
    
    std::string lineContents;
    std::string name1 = "";
    std::ifstream fastaFile(shorreadsdir);
    size_t didx = outdir.find_last_of(".");
    std::string csvoutdir ="";
    if(didx != std::string::npos){
        csvoutdir = outdir.substr(0,didx) + ".tsv";
    }
    else{
        csvoutdir = outdir + ".tsv";
        outdir+=".txt";
    }
    std::ofstream sout(outdir);
    std::ofstream sout2(csvoutdir);
    sout2 << "Read_ID\tNum_of_Matches\tHighest_Score\tBest_Match\tDirection\tStart_Pos\n";
    lineContents = "";
    std::string read, uid;
    int count = 0, count1 = 0, readlength = 0;
    bool ident = false;
    
    while(!fastaFile.eof()){
        getline(fastaFile, lineContents);
        std::stringstream readstream(lineContents);
        if(count == 0 && lineContents[0] == '>'){
            if(ident == true){
                std::vector<NGS> NNList;
                (*CVPTree).RangeSearch(read, tau, &NNList, &ovsc, &ovas);
                indexRead(searchmeth, read, &NNList, r);
                
                eliminateDuplicateMatching(&NNList);
                
                
                if((int)NNList.size() > 0){
                    sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
                    sout << uid << "\n";
                    sout << "-------------------------------------------------------------------------------------------\n";
                    sout << "-------------------------------------------------------------------------------------------\n\n";
                    int bestscor = std::numeric_limits<int>::lowest();
                    int bestidx = 0;
                    for (int i = 0; i < (int)NNList.size(); i++){
                        if(bestscor <  (int)NNList[i].returnAlscore()){
                            bestscor =  (int)NNList[i].returnAlscore();
                            bestidx = i;
                        }
                        int stapos = NNList[i].returnCoordinates()[0].second;
                        int enpos = (*NNList[i].returnEndpos())[0];
                        
                        sout << "-------------------------------------------------------------------------------------------\n";
                        sout << (*RefList)[NNList[i].returnCoordinates()[0].first].first << "\n";
                        if(NNList[i].returnAldirection(0) == false){
                            sout << "Direction: Forward\n";
                        }
                        else{
                            sout << "Direction: Reverse\n";
                        }
                        sout << "Range from " << std::to_string(stapos+1) << " to " << std::to_string (enpos+1) << "\n";
                        sout << "Alignment Score: " << std::to_string((int)NNList[i].returnAlscore()) << "/" << std::to_string((int)read.size()*2)<< "\n";
                        sout << "Cigar: " << NNList[i].returnCigar()[0] << "\n\n";
                        
                        
                        
                        if(NNList[i].returnAldirection(0) == false){
                            std::string refread = (*RefList)[NNList[i].returnCoordinates()[0].first].second.substr(stapos, (enpos - stapos)+1);
                            sout << returnAligstring(read, refread, NNList[i].returnCigar()[0], stapos+1, false) << "\n";
                        }
                        else{
                            std::string refread = (*RefList)[NNList[i].returnCoordinates()[0].first].second.substr(stapos-1, (enpos - stapos)+1);
                            
                            sout << returnAligstring(read, RVC.returnReversCompliment(refread), NNList[i].returnCigar()[0], enpos + 1, true) << "\n";
                        }
                        sout << "-------------------------------------------------------------------------------------------\n\n";
                    }
                    sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
                    std::string bdir = "Forward";
                    if(NNList[bestidx].returnAldirection(0) == true){
                        bdir = "Reverse";
                    }
                    sout2 << uid << "\t" << std::to_string((int) NNList.size()) << "\t" << std::to_string(bestscor) <<"\t"<< (*RefList)[NNList[bestidx].returnCoordinates()[0].first].first <<"\t" << bdir << "\t" << std::to_string(NNList[bestidx].returnCoordinates()[0].second) <<"\n";
                    creads++;
                }
                else{
                    if(repall == true){
                        sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
                        sout << uid << "\n";
                        sout << "-------------------------------------------------------------------------------------------\n";
                        sout << "-------------------------------------------------------------------------------------------\n\n";
                        sout << "Unclassified read.\n";
                        sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
                        sout2 << uid << "\t" << 0 << "\t"  << "\t" << "\t" << "\t" <<"\n";
                    }
                    uncreads++;
                }
                
                
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
    fastaFile.close();
    
    if(ident == true){
        
        std::vector<NGS> NNList;
        (*CVPTree).RangeSearch(read, tau, &NNList, &ovsc, &ovas);
        
        indexRead(searchmeth, read, &NNList, r);
        eliminateDuplicateMatching(&NNList);
        
        if((int)NNList.size() > 0){
            sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
            sout << uid << "\n";
            sout << "-------------------------------------------------------------------------------------------\n";
            sout << "-------------------------------------------------------------------------------------------\n\n";
            int bestscor = std::numeric_limits<int>::lowest();
            int bestidx = 0;
            for (int i = 0; i < (int)NNList.size(); i++){
                if(bestscor <  (int)NNList[i].returnAlscore()){
                    bestscor =  (int)NNList[i].returnAlscore();
                    bestidx = i;
                }
                int stapos = NNList[i].returnCoordinates()[0].second;
                int enpos = (*NNList[i].returnEndpos())[0];
                
                sout << "-------------------------------------------------------------------------------------------\n";
                sout << (*RefList)[NNList[i].returnCoordinates()[0].first].first << "\n";
                if(NNList[i].returnAldirection(0) == false){
                    sout << "Direction: Forward\n";
                }
                else{
                    sout << "Direction: Reverse\n";
                }
                sout << "Range from " << std::to_string(stapos+1) << " to " << std::to_string (enpos+1) << "\n";
                sout << "Alignment Score: " << std::to_string((int)NNList[i].returnAlscore()) << "/" << std::to_string((int)read.size()*2)<< "\n";
                sout << "Cigar: " << NNList[i].returnCigar()[0] << "\n\n";
                
                
                
                if(NNList[i].returnAldirection(0) == false){
                    std::string refread = (*RefList)[NNList[i].returnCoordinates()[0].first].second.substr(stapos, (enpos - stapos)+1);
                    sout << returnAligstring(read, refread, NNList[i].returnCigar()[0], stapos+1, false) << "\n";
                }
                else{
                    std::string refread = (*RefList)[NNList[i].returnCoordinates()[0].first].second.substr(stapos-1, (enpos - stapos)+1);
                    sout << returnAligstring(read, RVC.returnReversCompliment(refread), NNList[i].returnCigar()[0], enpos + 1, true) << "\n";
                }
                sout << "-------------------------------------------------------------------------------------------\n\n";
            }
            sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
            std::string bdir = "Forward";
            if(NNList[bestidx].returnAldirection(0) == true){
                bdir = "Reverse";
            }
            sout2 << uid << "\t" << std::to_string((int) NNList.size()) << "\t" << std::to_string(bestscor) <<"\t"<< (*RefList)[NNList[bestidx].returnCoordinates()[0].first].first <<"\t" << bdir << "\t" << std::to_string(NNList[bestidx].returnCoordinates()[0].second) <<"\n";
            creads++;
        }
        else{
            if(repall == true){
                sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
                sout << uid << "\n";
                sout << "-------------------------------------------------------------------------------------------\n";
                sout << "-------------------------------------------------------------------------------------------\n\n";
                sout << "Unclassified read.\n";
                sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
                sout2 << uid << "\t" << 0 << "\t"  << "\t" << "\t" << "\t" <<"\n";
            }
            uncreads++;
        }
        
    }
    sout.close();
    sout2 << "Total number of reads:\t" << std::to_string(creads) <<"\tMatched reads:\t" << std::to_string(creads - uncreads) << "\tUnmatched reads:\t" << std::to_string(uncreads)<<"\n";
    
    sout2.close();
    std::cout << "> Num of classified reads: " << creads << "\n";
    std::cout << "> Num of unclassified reads: " << uncreads << "\n";
};

void ClassificationInitiation::fastqParserfile(ClassificationTree *Tree, std::string shorreadsdir, std::string outdir, std::string searchmeth, int knn, bool r, bool repall, bool sen){
    
    long long ovsc = 0;
    long long ovas = 0;
    int uncreads = 0;
    int creads = 0;
    
    std::string lineContents;
    std::string name1 = "";
    std::ifstream fastqFile(shorreadsdir);
    size_t didx = outdir.find_last_of(".");
    std::string csvoutdir ="";
    if(didx != std::string::npos){
        csvoutdir = outdir.substr(0,didx) + ".tsv";
    }
    else{
        csvoutdir = outdir + ".tsv";
        outdir+=".txt";
    }
    
    std::ofstream sout(outdir);
    std::ofstream sout2(csvoutdir);
    
    sout2 << "Read_ID\tNum_of_Matches\tHighest_Score\tBest_Match\tDirection\tStart_Pos\n";
    
    
    lineContents = "";
    std::string read, uid;
    int count = 0, count1 = 0, readlength = 0;
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
                
                std::vector<NGS> KNNList;
                (*CVPTree).KNNSearch(read, knn, &KNNList, &ovsc, &ovas, sen);
                indexRead(searchmeth, read, &KNNList, r);
                eliminateDuplicateMatching(&KNNList);
                
                if((int)KNNList.size() > 0){
                    sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
                    sout << uid << "\n";
                    sout << "-------------------------------------------------------------------------------------------\n";
                    sout << "-------------------------------------------------------------------------------------------\n\n";
                    int bestscor = std::numeric_limits<int>::lowest();
                    int bestidx = 0;
                    for (int i = 0; i < (int)KNNList.size(); i++){
                        if(bestscor <  (int)KNNList[i].returnAlscore()){
                            bestscor =  (int)KNNList[i].returnAlscore();
                            bestidx = i;
                        }
                        int stapos = KNNList[i].returnCoordinates()[0].second;
                        int enpos = (*KNNList[i].returnEndpos())[0];
                        
                        sout << "-------------------------------------------------------------------------------------------\n";
                        sout << (*RefList)[KNNList[i].returnCoordinates()[0].first].first << "\n";
                        if(KNNList[i].returnAldirection(0) == false){
                            sout << "Direction: Forward\n";
                        }
                        else{
                            sout << "Direction: Reverse\n";
                        }
                        sout << "Range from " << std::to_string(stapos+1) << " to " << std::to_string (enpos+1) << "\n";
                        sout << "Alignment Score: " << std::to_string((int)KNNList[i].returnAlscore()) << "/" << std::to_string((int)read.size()*2)<< "\n";
                        sout << "Cigar: " << KNNList[i].returnCigar()[0] << "\n\n";
                        
                        
                        
                        if(KNNList[i].returnAldirection(0) == false){
                            std::string refread = (*RefList)[KNNList[i].returnCoordinates()[0].first].second.substr(stapos, (enpos - stapos)+1);
                            sout << returnAligstring(read, refread, KNNList[i].returnCigar()[0], stapos+1, false) << "\n";
                        }
                        else{
                            std::string refread = (*RefList)[KNNList[i].returnCoordinates()[0].first].second.substr(stapos-1, (enpos - stapos)+1);
                            
                            sout << returnAligstring(read,  RVC.returnReversCompliment(refread), KNNList[i].returnCigar()[0], enpos + 1, true) << "\n";
                        }
                        sout << "-------------------------------------------------------------------------------------------\n\n";
                    }
                    sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
                    std::string bdir = "Forward";
                    if(KNNList[bestidx].returnAldirection(0) == true){
                        bdir = "Reverse";
                    }
                    sout2 << uid << "\t" << std::to_string((int) KNNList.size()) << "\t" << std::to_string(bestscor) <<"\t"<< (*RefList)[KNNList[bestidx].returnCoordinates()[0].first].first <<"\t" << bdir << "\t" << std::to_string(KNNList[bestidx].returnCoordinates()[0].second) <<"\n";
                    creads++;
                }
                else{
                    if(repall == true){
                        sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
                        sout << uid << "\n";
                        sout << "-------------------------------------------------------------------------------------------\n";
                        sout << "-------------------------------------------------------------------------------------------\n\n";
                        sout << "Unclassified read.\n";
                        sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
                        sout2 << uid << "\t" << 0 << "\t"  << "\t" << "\t" << "\t" <<"\n";
                    }
                    uncreads++;
                }
                
                
            }
            else if(count == 2){
                count++;
            }
            else if(count == 3){
                count = 0;
                count1++;
                if((count1%100000) == 0){
                    std::cout << "> Fastq read: " << count1 << "\n";
                }
            }
        }
    }
    fastqFile.close();
    sout.close();
    sout2 << "Total number of reads:\t" << std::to_string(count1) <<"\tMatched reads:\t" << std::to_string(count1 - uncreads) << "\tUnmatched reads:\t" << std::to_string(uncreads)<<"\n";
    sout2.close();
    std::cout << "> Num of classified reads: " << creads << "\n";
    std::cout << "> Num of unclassified reads: " << uncreads << "\n";
};

void ClassificationInitiation::fastqParserfile(ClassificationTree *Tree, std::string shorreadsdir, std::string outdir, std::string searchmeth, double tau, bool r, bool repall){
    
    long long ovsc = 0;
    long long ovas = 0;
    int uncreads = 0;
    int creads = 0;
    
    std::string lineContents;
    std::string name1 = "";
    std::ifstream fastqFile(shorreadsdir);
    size_t didx = outdir.find_last_of(".");
    std::string csvoutdir ="";
    if(didx != std::string::npos){
        csvoutdir = outdir.substr(0,didx) + ".tsv";
    }
    else{
        csvoutdir = outdir + ".tsv";
        outdir+=".txt";
    }
    
    std::ofstream sout(outdir);
    std::ofstream sout2(csvoutdir);
    sout2 << "Read_ID\tNum_of_Matches\tHighest_Score\tBest_Match\tDirection\tStart_Pos\n";
    
    lineContents = "";
    std::string read, uid;
    int count = 0, count1 = 0, readlength = 0;
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
                
                std::vector<NGS> NNList;
                
                (*CVPTree).RangeSearch(read, tau, &NNList, &ovsc, &ovas);
                indexRead(searchmeth, read, &NNList, r);
                eliminateDuplicateMatching(&NNList);
                
                if((int)NNList.size() > 0){
                    sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
                    sout << uid << "\n";
                    sout << "-------------------------------------------------------------------------------------------\n";
                    sout << "-------------------------------------------------------------------------------------------\n\n";
                    int bestscor = std::numeric_limits<int>::lowest();
                    int bestidx = 0;
                    for (int i = 0; i < (int)NNList.size(); i++){
                        if(bestscor <  (int)NNList[i].returnAlscore()){
                            bestscor =  (int)NNList[i].returnAlscore();
                            bestidx = i;
                        }
                        int stapos = NNList[i].returnCoordinates()[0].second;
                        int enpos = (*NNList[i].returnEndpos())[0];
                        sout << "-------------------------------------------------------------------------------------------\n";
                        sout << (*RefList)[NNList[i].returnCoordinates()[0].first].first << "\n";
                        if(NNList[i].returnAldirection(0) == false){
                            sout << "Direction: Forward\n";
                        }
                        else{
                            sout << "Direction: Reverse\n";
                        }
                        sout << "Range from " << std::to_string(stapos+1) << " to " << std::to_string (enpos+1) << "\n";
                        sout << "Alignment Score: " << std::to_string((int)NNList[i].returnAlscore()) << "/" << std::to_string((int)read.size()*2)<< "\n";
                        sout << "Cigar: " << NNList[i].returnCigar()[0] << "\n\n";
                        
                        if(NNList[i].returnAldirection(0) == false){
                            std::string refread = (*RefList)[NNList[i].returnCoordinates()[0].first].second.substr(stapos, (enpos - stapos)+1);
                            sout << returnAligstring(read, refread, NNList[i].returnCigar()[0], stapos+1, false) << "\n";
                        }
                        else{
                            std::string refread = (*RefList)[NNList[i].returnCoordinates()[0].first].second.substr(stapos-1, (enpos - stapos)+1);
                            sout << returnAligstring(read, RVC.returnReversCompliment(refread), NNList[i].returnCigar()[0], enpos + 1, true) << "\n";
                        }
                        sout << "-------------------------------------------------------------------------------------------\n\n";
                    }
                    sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
                    std::string bdir = "Forward";
                    if(NNList[bestidx].returnAldirection(0) == true){
                        bdir = "Reverse";
                    }
                    sout2 << uid << "\t" << std::to_string((int) NNList.size()) << "\t" << std::to_string(bestscor) <<"\t"<< (*RefList)[NNList[bestidx].returnCoordinates()[0].first].first <<"\t" << bdir << "\t" << std::to_string(NNList[bestidx].returnCoordinates()[0].second) <<"\n";
                    
                    creads++;
                }
                else{
                    if(repall == true){
                        sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
                        sout << uid << "\n";
                        sout << "-------------------------------------------------------------------------------------------\n";
                        sout << "-------------------------------------------------------------------------------------------\n\n";
                        sout << "Unclassified read.\n";
                        sout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
                        sout2 << uid << "\t" << 0 << "\t"  << "\t" << "\t" << "\t" <<"\n";
                    }
                    
                    uncreads++;
                }
                
                
            }
            else if(count == 2){
                count++;
            }
            else if(count == 3){
                count = 0;
                count1++;
                if((count1%100000) == 0){
                    std::cout << "> Fastq read: " << count1 << "\n";
                }
            }
        }
    }
    fastqFile.close();
    sout.close();
    sout2 << "Total number of reads:\t" << std::to_string(count1) <<"\tMatched reads:\t" << std::to_string(count1 - uncreads) << "\tUnmatched reads:\t" << std::to_string(uncreads)<<"\n";
    sout2.close();
    std::cout << "> Num of classified reads: " << creads << "\n";
    std::cout << "> Num of unclassified reads: " << uncreads << "\n";
}

std::string ClassificationInitiation::returnAligstring(std::string read, std::string ref, std::string cigar, int refstart, bool dir){
    std::string alstr;
    
    std::string refsubread = "";
    std::string readsubread = "";
    std::string pipesubread = "";
    
    std::string decig = DecomposeCigar(cigar);
    int refcount = 0;
    int readcount = 0;
    for (int i = 0; i < (int)decig.size(); i++){
        
        if(decig.substr(i,1) == "M"){
            refsubread += ref.substr(refcount,1);
            if(ref.substr(refcount,1) == read.substr(readcount,1)){
                pipesubread += "|";
            }
            else{
                pipesubread += "x";
            }
            readsubread += read.substr(readcount,1);
            refcount++;
            readcount++;
        }
        else if(decig.substr(i,1) == "D"){
            refsubread += ref.substr(refcount,1);
            pipesubread += " ";
            readsubread += "-";
            refcount++;
        }
        else if(decig.substr(i,1) == "I"){
            refsubread += "-";
            pipesubread += " ";
            readsubread += read.substr(readcount,1);
            readcount++;
        }
        else if(decig.substr(i,1) == "S"){
            readcount++;
        }
        if((int)pipesubread.size() == 60){
            std::string Refr = "Ref: ";
            std::string Quer = "Que: ";
            std::string Pipr = "               ";
            
            int inreadcount = readcount - 60;
            if(inreadcount < 1){
                inreadcount = 1;
            }
            
            int refs2use, refe2use;
            if(dir == false){
                refs2use = refstart + (refcount - 60) + 1;
                if(refs2use < 1){
                    refs2use = 1;
                }
                refe2use = refstart+ refcount;
            }
            else{
                refs2use = refstart - (refcount - (60 - 1)) + 1;
                refe2use = (refstart + 1) - (refcount - 1);
                if(refe2use < 1){
                    refe2use = 1;
                }
            }
            alstr += adjustString(refs2use, "Ref") +  refsubread + "   " + std::to_string(refe2use-1)  + "   " + "\n" + Pipr + pipesubread + "\n" + adjustString(inreadcount, "Que") + readsubread + "   " + std::to_string(readcount-1) + "\n\n";
            refsubread = "";
            pipesubread = "";
            readsubread = "";
        }
    }
    
    if((int)refsubread.size() > 0){
        int inreadcount = readcount - ((int)readsubread.size());
        if(inreadcount < 1){
            inreadcount = 1;
        }
        
        int refs2use, refe2use;
        if(dir == false){
            refs2use = refstart + (refcount - (int)refsubread.size());
            if(refs2use < 1){
                refs2use = 1;
            }
            refe2use = refstart+ refcount;
        }
        else{
            refs2use = refstart - (refcount - ((int)refsubread.size()-1)) + 1;
            refe2use = (refstart + 1) - (refcount - 1);
            if(refe2use < 1){
                refe2use = 1;
            }
        }
        std::string Pipr = "               ";
        alstr += adjustString(refs2use, "Ref") +  refsubread + "   " + std::to_string(refe2use)  + "   " + "\n" + Pipr + pipesubread + "\n" + adjustString(inreadcount, "Que") + readsubread + "   " + std::to_string(readcount) + "\n\n";
    }
    else{
        alstr +="\n";
    }
    return alstr;
};


std::string ClassificationInitiation::adjustString(int i, std::string title){
    std::string r = title + ": ";
    int l = 10;
    for(int j = 0; j < l - (int) std::to_string(i).size() ; j++){
        title +=" ";
    }
    title += std::to_string(i) + "  ";
    return title;
};

void ClassificationInitiation::setAccumNorm(std::pair<bool, bool> AccumNorm){
    Accum = AccumNorm.first;
    Norm  = AccumNorm.second;
};


std::string ClassificationInitiation::DecomposeCigar(std::string cigar){
    std::string decompose = "";
    std::map<size_t, std::string,  std::less<int>> mapkeys;
    size_t lpos = 0;
    size_t plpos = 0;
    while (lpos !=  std::string::npos){
        lpos = cigar.substr(plpos+1).find_first_of("M");
        if(lpos != std::string::npos){
            mapkeys[plpos + lpos+1] = "M";
            plpos += lpos + 1;
        }
    }
    
    lpos = 0;
    plpos = 0;
    while (lpos !=  std::string::npos){
        lpos = cigar.substr(plpos+1).find_first_of("D");
        if(lpos != std::string::npos){
            mapkeys[plpos + lpos+1] = "D";
            plpos += lpos + 1;
        }
    }
    
    lpos = 0;
    plpos = 0;
    while (lpos !=  std::string::npos){
        lpos = cigar.substr(plpos+1).find_first_of("I");
        if(lpos != std::string::npos){
            mapkeys[plpos + lpos+1] = "I";
            plpos += lpos + 1;
        }
    }
    
    lpos = 0;
    plpos = 0;
    while (lpos !=  std::string::npos){
        lpos = cigar.substr(plpos+1).find_first_of("S");
        if(lpos != std::string::npos){
            mapkeys[plpos + lpos+1] = "S";
            plpos += lpos + 1;
        }
    }
    
    lpos = 0;
    for (std::map<size_t, std::string,  std::greater<int>>::iterator it = mapkeys.begin(); it != mapkeys.end(); it++){
        for (int i = 0; i < std::stoi(cigar.substr(lpos, (it->first - lpos))); i++){
            decompose += it->second;
        }
        lpos = it->first+1;
    }
    
    return decompose;
};

ClassificationInitiation::ClassificationInitiation(){};

ClassificationInitiation::ClassificationInitiation(ClassificationTree *Tree, std::string shorreadsdir, std::string outdir, std::string searchmeth, int knn, bool r, bool repall, bool sen){
    setTree(Tree);
    int filetype = testfile(shorreadsdir);
    if(searchmeth != "ED" && searchmeth != "DTW" && searchmeth != "SW"){
        searchmeth = "DTW";
    }
    
    if (filetype == 0){
        //Fastq
        fastqParserfile(Tree, shorreadsdir, outdir, searchmeth, knn, r, repall, sen);
    }
    else if (filetype == 1){
        //Fasta
        fastaParserfile(Tree, shorreadsdir, outdir, searchmeth, knn, r, repall, sen);
    }
    else{
        //Wrong File type
        std::cout << "Wrong file was provided.\n";
    }
    
};


ClassificationInitiation::ClassificationInitiation(ClassificationTree *Tree, std::string shorreadsdir, std::string outdir, std::string searchmeth, double tau, bool r, bool repall){
    setTree(Tree);
    int filetype = testfile(shorreadsdir);
    
    if(searchmeth != "ED" && searchmeth != "DTW" && searchmeth != "SW"){
        searchmeth = "DTW";
    }
    
    if (filetype == 0){
        fastqParserfile(Tree, shorreadsdir, outdir, searchmeth, tau, r, repall);
    }
    else if (filetype == 1){
        
        fastaParserfile(Tree, shorreadsdir, outdir, searchmeth, tau, r, repall);
    }
    else{
        std::cout << "Wrong file was provided.\n";
    }
    
};

void ClassificationInitiation::setTree(ClassificationTree *Tree){
    CVPTree = Tree;
    RefList = (*CVPTree).returnRefList();
    Repmeth = (*CVPTree).returnRepmeth();
    setAccumNorm((*CVPTree).returnAccumNorm());
};


ClassificationInitiation::~ClassificationInitiation(){};


