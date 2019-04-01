//
//  main.cpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>

#include "ClassificationTree.hpp"
#include "ClassificationInitiation.hpp"
#include "InstructionsClass.hpp"
#include "Print_Time.hpp"
#include "TestFileExistance.hpp"
int main(int argc, const char * argv[]) {
    ClassificationTree CVPT;
    Print_Time PT;
    TestFileExistance TE;
    
    bool vptprov = false, faprov = false, rprov = false , outprov = false, pKNN = true, pRange = false, vpts = false, inddirprov = false, preporc = true, accu = false, znorm = false, ptd = false, repall = true, isen = false;
    std::string fadir = "", vptdir = "", readsdir = "", outdir = "", rep = "Tetrahedron", tra = "DWT", indexdir = "", sapp = "DTW";
    int kmer = 64, clvl = 4, knn = 5;
    double range = 0.15;
    clock_t start, stop;
    size_t sidx = readsdir.find_last_of(".");
    std::cout << readsdir.substr(0,sidx+1) << "\n";
    
    std::cout << "#############################################################\n";
    std::cout << "#                                                           #\n";
    std::cout << "#                            CBN                            #\n";
    std::cout << "#                  ClassificationByNumbers                  #\n";
    std::cout << "#                                                           #\n";
    std::cout << "#############################################################\n";
    
    int civ = (argc - 1) % 2;
    if (civ || argc - 1 == 0) {
        for (int i = 1; i < argc; i++) {
            if ((std::string) (argv[i]) == "-h"
                || (std::string) (argv[i]) == "-help"
                || (std::string) (argv[i]) == "-H"
                || (std::string) (argv[i]) == "-Help"
                || (std::string) (argv[i]) == "-HELP") {
                InstructionsClass instructions(1);
                return 1;
            }
        }
        InstructionsClass instreuctions(0);
        return 1;
    }
    else {
        for (int i = 1; i < argc; i += 2) {
            if ((std::string) (argv[i]) == "-h"
                || (std::string) (argv[i]) == "-help"
                || (std::string) (argv[i]) == "-H"
                || (std::string) (argv[i]) == "-Help"
                || (std::string) (argv[i]) == "-HELP") {
                InstructionsClass instructions(0);
            }
            else if ((std::string) (argv[i]) == "-o") {
                outdir = (std::string) (argv[i + 1]);
                outprov = true;
            }
            else if ((std::string) (argv[i]) == "-r") {
                readsdir = (std::string) (argv[i + 1]);
                rprov = TE.exists_test(readsdir);
                if (rprov == false){
                    std::cout << "> No such file or directory exist for the provided fastq file\n";
                    InstructionsClass instructions(0);
                    return 1;
                }
            }
            else if ((std::string) (argv[i]) == "-fa"){
                fadir = (std::string) (argv[i + 1]);
                faprov = TE.exists_test(fadir);
                if(faprov == false){
                    std::cout << "> No such file or directory exist for the provided fasta file\n";
                    InstructionsClass instructions(0);
                    return 1;
                }
            }
            else if ((std::string) (argv[i]) == "-vpt_load"){
                vptdir = (std::string) (argv[i + 1]);
                size_t didx = vptdir.find(".vptstruct");
                if(!(didx != std::string::npos)){
                    vptprov = TE.exists_test(vptdir+".vptstruct");
                    if(vptprov == true){
                        vptdir+=".vptstruct";
                    }
                    else{
                        std::cout << "> No such file or directory exist for the provided VPTree index file\n";
                        InstructionsClass instructions(1);
                        return 1;
                    }
                }
                else{
                    vptprov = TE.exists_test(vptdir);
                    if(vptprov == false){
                        std::cout << "> No such file or directory exist for the provided VPTree index file\n";
                        InstructionsClass instructions(1);
                        return 1;
                    }
                }
            }
            else if ((std::string) (argv[i]) == "-kmer") {
                kmer = std::stoi(argv[i + 1]);
            }
            else if ((std::string) (argv[i]) == "-rep") {
                rep = (std::string) (argv[i + 1]);
            }
            else if ((std::string) (argv[i]) == "-tra") {
                tra = (std::string) (argv[i + 1]);
            }
            else if ((std::string) (argv[i]) == "-clvl") {
                clvl = std::stoi(argv[i + 1]);
            }
            else if ((std::string) (argv[i]) == "-knn") {
                knn = std::stoi(argv[i + 1]);
                pKNN = true;
                pRange = false;
            }
            else if ((std::string) (argv[i]) == "-ra") {
                range = std::stod(argv[i + 1]);
                pKNN = false;
                pRange = true;
            }
            else if ((std::string) (argv[i]) == "-s") {
                if((std::string) (argv[i + 1]) == "y" || (std::string) (argv[i + 1]) == "Y" || (std::string) (argv[i + 1]) == "Yes"){
                    vpts = true;
                }
                else{
                    vpts = false;
                }
            }
            else if ((std::string) (argv[i]) == "-vpt_save") {
                indexdir = (std::string) (argv[i + 1]);
                inddirprov = true;
                vpts = true;
            }
            else if((std::string) (argv[i]) == "-sen"){
                if((std::string) (argv[i + 1]) == "y" || (std::string) (argv[i + 1]) == "Y" || (std::string) (argv[i + 1]) == "Yes"){
                    isen = true;
                }
                else{
                    isen = false;
                }
            }
            else if ((std::string) (argv[i]) == "-evs") {
                if((std::string) (argv[i + 1]) == "s" ||(std::string) (argv[i + 1]) == "S" || (std::string) (argv[i + 1]) == "Slow" || (std::string) (argv[i + 1]) == "SLOW" || (std::string) (argv[i + 1]) == "slow"){
                    sapp = "SW";
                }
                else if((std::string) (argv[i + 1]) == "m" ||(std::string) (argv[i + 1]) == "M" || (std::string) (argv[i + 1]) == "Medium" || (std::string) (argv[i + 1]) == "MEDIUM" || (std::string) (argv[i + 1]) == "medium"){
                    sapp = "DTW";
                }
                else{
                    sapp = "ED";
                }
            }
            else if ((std::string) (argv[i]) == "-proc") {
                if((std::string) (argv[i + 1]) == "y" || (std::string) (argv[i + 1]) == "Y" || (std::string) (argv[i + 1]) == "Yes"){
                    preporc = true;
                }
                else{
                    preporc = false;
                }
            }
            else if ((std::string) (argv[i]) == "-accu") {
                if((std::string) (argv[i + 1]) == "y" || (std::string) (argv[i + 1]) == "Y" || (std::string) (argv[i + 1]) == "Yes"){
                    accu = true;
                }
                else{
                    accu = false;
                }
            }
            else if ((std::string) (argv[i]) == "-znorm") {
                if((std::string) (argv[i + 1]) == "y" || (std::string) (argv[i + 1]) == "Y" || (std::string) (argv[i + 1]) == "Yes"){
                    znorm = true;
                }
                else{
                    znorm = false;
                }
            }
            else if ((std::string) (argv[i]) == "-ptd") {
                if((std::string) (argv[i + 1]) == "y" || (std::string) (argv[i + 1]) == "Y" || (std::string) (argv[i + 1]) == "Yes"){
                    ptd = true;
                }
                else{
                    ptd = false;
                }
            }
            
            else if ((std::string) (argv[i]) == "-repall") {
                if((std::string) (argv[i + 1]) == "y" || (std::string) (argv[i + 1]) == "Y" || (std::string) (argv[i + 1]) == "Yes"){
                    repall = true;
                }
                else{
                    repall = false;
                }
            }   
        }
    }
    
    if( vptprov == false && faprov == false){
        std::cout << "> No index file or fasta file reference file has been provided.\n";
        std::cout << "> Process is termianted.\n";
        InstructionsClass instructions(1);
    }
    
    
    if(rprov == false){
        std::cout << "> No reads file has been provided.\n";
        std::cout << "> Process is termianted.\n";
        InstructionsClass instructions(1);
    }
    
    if(outprov == false){
        size_t sidx = readsdir.find_last_of(".");
        outdir =  readsdir.substr(0,sidx);
        if(pKNN == true){
            outdir += "> KNN_Classification.txt";
        }
        else{
            outdir += "> Range_Search_Classification.txt";
        }
        std::cout << "> No output directory has been provided for storing the classification results.\n";
        std::cout << "> Results will be stored in:\n";
        std::cout << "> " << outdir << "\n";
        
    }
    
    if (inddirprov == true){
        vpts = true;
        
    };
    
    
    if(vpts == true){
        if(vptprov == true){
            std::cout << "> Indexing tree already exist in loaction:\n";
            std::cout << "> " << vptdir << "n";
            std::cout << "> Tree will not be saved.\n";
            vpts = false;
        }
        else{
            if(inddirprov == false){
            size_t sidx = fadir.find_last_of(".");
            indexdir =  fadir.substr(0,sidx);
            std::cout << "> No directory has been provided for storing the indexing tree.\n";
            std::cout << "> Indexing tree  will be storied in:\n";
            std::cout <<  "> " << indexdir << "\n";
            }
        }
    }
    
    if(vptprov == true){
        std::cout << "> Initiate loading indexing structure.\n";
        start = clock();
        CVPT.loadTree(vptdir);
        stop = clock();
        std::cout << "> Loading indexing structure completed.\n";
        std::cout << "> Loading tree elapsed time: " << PT.returnTimestr(start, stop) << "\n";
    }
    else{
        start = clock();
        std::cout << "> Initiate building indexing structure.\n";
        start = clock();
        CVPT.BuildTree(fadir, kmer, rep, tra, clvl, accu, znorm);
        stop = clock();
        std::cout << "> Building indexing structure completed.\n";
        std::cout << "> Building tree elapsed time: " << PT.returnTimestr(start, stop) << "\n";
    }
    if(ptd == true){
        std::cout << "> Index tree contains " << CVPT.TreeNodesNum() << " nodes\n";
        CVPT.TreeprintAllTau();
    }
    if(pKNN == true){
        std::cout << "> Initiate KNN search.\n";
        start = clock();
        ClassificationInitiation CIKNN(&CVPT, readsdir, outdir, sapp, knn, preporc, repall, isen);
        stop = clock();
        std::cout << "> KNN search completed.\n";
        std::cout << "> KNN search elapsed time: " << PT.returnTimestr(start, stop) << "\n";
    }
    else if(pRange == true){
        std::cout << "> Initiate range search.\n";
        start = clock();
        ClassificationInitiation CIRange(&CVPT, readsdir, outdir ,sapp, range, preporc, repall);
        stop = clock();
        std::cout << "> Range search completed.\n";
        std::cout << "> Range search elapsed time: " << PT.returnTimestr(start, stop) << "\n";
    }
    
    
    if(vpts == true){
        std::cout << "> Initiate saving indexing structure.\n";
        start = clock();
        CVPT.saveTree(indexdir);
        stop = clock();
        std::cout << "> Saving indexing structure completed.\n";
        std::cout << "> Saving indexing structure elapsed time: " << PT.returnTimestr(start, stop) << "\n";
    }
    std::cout << "#############################################################\n\n";
    return 0;
}
