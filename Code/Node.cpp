//
//  Node.cpp
//  ClassificationTool
//
//  Created by Avraam Tapinos
//  Copyright © 2018 Avraam Tapinos. All rights reserved.
//

#include "Node.hpp"
void Node::setNtype(bool b){
    
    Ntype = b;
};

void Node::setINode(double t, double c, double m){
    tau = t;
    thresh = c;
    maxdis = m;
};

void Node::QuickSort(std::vector<NGS> *dataPoints, int leftmost, int rightmost){
    int i = leftmost, j = rightmost;
    double  pivot = (*dataPoints)[(leftmost + rightmost)/2].returnAlscore();
    
    while (i <= j) {
        
        while ((*dataPoints)[i].returnAlscore() < pivot){
            
            i++;
        };
        
        while ((*dataPoints)[j].returnAlscore() > pivot){
            
            j--;
        };
        
        if (i <= j) {
            std::swap((*dataPoints)[i], (*dataPoints)[j]);
            i++;
            j--;
        };
    };
    
    // recursion //
    if (leftmost < j){
        QuickSort(dataPoints, leftmost, j);
    };
    if (i < rightmost){
        QuickSort(dataPoints, i, rightmost);
    };
    
};

Node::Node(){
    tau = 0;
    thresh = 0;
    Ntype = false;
    pp=1;
    
};

int Node::randominitiation(std::vector<NGS> *dataPoints, int m, int n){
    srand ((unsigned)time(NULL));
    int v = rand() % n;
    vantagepoint  = &dataPoints->at(v);
    for(int i = 0;  i < n; i++){
        (*dataPoints)[i].setAlscore(DM.returnDistance(vantagepoint->returnTra(), (*dataPoints)[i].returnTra()));
    };
    
    
    QuickSort(dataPoints, m, n-1);
    return n-1;
};

void Node::generateNode(std::vector<NGS> *dataPoints, int m, int n, bool r){
    setNtype(true);
    if(r == true){
        n  = randominitiation(dataPoints, m, n);
    };
    vantagepoint  = &(*dataPoints)[n];
    
    if(m < n){
        
        for(int i = m;  i < n; i++){
            (*dataPoints)[i].setAlscore(DM.returnDistance(vantagepoint->returnTra(), (*dataPoints)[i].returnTra()));
        };
        QuickSort(dataPoints, m, n-1);
        
        int med = m +(((n-m))/2);
        for (int i = med+1; i < n; i++){
            if((*dataPoints)[i].returnAlscore() <= (*dataPoints)[med].returnAlscore()){
                med++;
            }
            else{
                break;
            }
        }
        n--;
        double tau = (*dataPoints)[med].returnAlscore();
        double thre =0;
        
        int q2 = med+((n-med)/2);
        q2 = ((n-m)*3/4);
        if (q2 > n){
            q2 = n;
        }
        
        
        
        if((*dataPoints)[m+q2].returnAlscore() >= tau){
            thre =  (*dataPoints)[m+q2].returnAlscore() - tau;
        }
        else{
            thre = tau - (*dataPoints)[m+q2].returnAlscore();
        }
        setINode(tau, thre, (*dataPoints)[n].returnAlscore());
        if(m < med){
            lNode = new Node;
            lNode->generateNode(dataPoints, m, med-1, false);
        }
        else{
            lNode = NULL;
        };
        
        if(med <= n){
            rNode = new Node;
            rNode->generateNode(dataPoints, med, n, false);
        }
        else{
            rNode = NULL;
        };
    }
    else{
        lNode = NULL;
        rNode = NULL;
    }
};


void Node::generateNode(std::vector<NGS> *dataPoints, int m, int n, bool r, int l){
    pp = l;
    setNtype(true);
    if(r == true){
        n  = randominitiation(dataPoints, m, n);
    };
    vantagepoint  = &(*dataPoints)[n];
    
    if(m < n){
        
        for(int i = m;  i < n; i++){
            (*dataPoints)[i].setAlscore(DM.returnDistance(vantagepoint->returnTra(), (*dataPoints)[i].returnTra()));
        };
        QuickSort(dataPoints, m, n-1);
        
        int med = m +(((n-m))/2);
        for (int i = med+1; i < n; i++){
            if((*dataPoints)[i].returnAlscore() <= (*dataPoints)[med].returnAlscore()){
                med++;
            }
            else{
                break;
            }
        }
        n--;
        double tau = (*dataPoints)[med].returnAlscore();
        double thre =0;
        
        int q2 = med+((n-med)/2);
        q2 = ((n-m)*3/4);
        if (q2 > n){
            q2 = n;
        }
        if((*dataPoints)[m+q2].returnAlscore() >= tau){
            thre =  (*dataPoints)[m+q2].returnAlscore() - tau;
        }
        else{
            thre = tau - (*dataPoints)[m+q2].returnAlscore();
        }
        setINode(tau, thre, (*dataPoints)[n].returnAlscore());
        if(m < med){
            lNode = new Node;
            lNode->generateNode(dataPoints, m, med-1, false, pp+1);
        }
        else{
            lNode = NULL;
        };
        
        if(med <= n){
            rNode = new Node;
            rNode->generateNode(dataPoints, med, n, false, pp+1);
        }
        else{
            rNode = NULL;
        };
    }
    else{
        lNode = NULL;
        rNode = NULL;
    }
};

void Node::generateNode(std::vector<NGS> *dataPoints, int m, int n, bool r, int l, std::map<int, std::pair<int, double>> *levDetails){
    pp = l;
    
    setNtype(true);
    if(r == true){
        n  = randominitiation(dataPoints, m, n);
    };
    vantagepoint  = &(*dataPoints)[n];
    
    if(m < n){
        
        for(int i = m;  i < n; i++){
            (*dataPoints)[i].setAlscore(DM.returnDistance(vantagepoint->returnTra(), (*dataPoints)[i].returnTra()));
        };
        QuickSort(dataPoints, m, n-1);
        
        int med = m +(((n-m))/2);
        for (int i = med+1; i < n; i++){
            if((*dataPoints)[i].returnAlscore() <= (*dataPoints)[med].returnAlscore()){
                med++;
            }
            else{
                break;
            }
        }
        n--;
        double tau = (*dataPoints)[med].returnAlscore();
        
        
        std::map<int, std::pair<int, double>> ::iterator lDit = (*levDetails).find(pp);
        if(lDit != (*levDetails).end()){
            (*levDetails)[pp].first++;
            (*levDetails)[pp].second += tau;
        }
        else{
            (*levDetails)[pp].first = 1;
            (*levDetails)[pp].second = tau;
        }
        
        
        double thre = 0;
        
        int q2 = med+((n-med)/2);
        q2 = ((n-m)*3/4);
        if (q2 > n){
            q2 = n;
        }
        if((*dataPoints)[m+q2].returnAlscore() >= tau){
            thre =  (*dataPoints)[m+q2].returnAlscore() - tau;
        }
        else{
            thre = tau - (*dataPoints)[m+q2].returnAlscore();
        }
        setINode(tau, thre, (*dataPoints)[n].returnAlscore());
        if(m < med){
            lNode = new Node;
            lNode->generateNode(dataPoints, m, med-1, false, pp+1, levDetails);
        }
        else{
            lNode = NULL;
        };
        
        if(med <= n){
            rNode = new Node;
            rNode->generateNode(dataPoints, med, n, false, pp+1, levDetails);
        }
        else{
            rNode = NULL;
        };
    }
    else{
        
        std::map<int, std::pair<int, double>> ::iterator lDit = (*levDetails).find(pp);
        if(lDit != (*levDetails).end()){
            (*levDetails)[pp].first++;
            (*levDetails)[pp].second += tau;
        }
        else{
            (*levDetails)[pp].first = 1;
            (*levDetails)[pp].second = tau;
        }
        
        lNode = NULL;
        rNode = NULL;
    }
};

void Node::KnnSearch(NGS *q, int k, std::map<double, std::vector<NGS>, std::greater<double>> *KNN, int *sc, int *absc, bool sen){
    double dist;
    dist = DM.returnDistance(q->returnTra(), vantagepoint->returnTra());
    double ftau = thresh;
    std::map<double, std::vector<NGS>, std::greater<double>>::iterator it;
    
    if ((int)vantagepoint->returnGID().size()  ==  1 && q->returnSRid() == *vantagepoint->returnGID().at(0)){
        
    }
    else{
        if((int)KNN->size() < k){
            (*KNN)[dist].push_back(*vantagepoint);
        }
        else{
            it = KNN->begin();
            if(it->first >= dist){
                (*KNN)[dist].push_back(*vantagepoint);
                if((int)KNN->size() > k){
                    it = KNN->begin();
                    (*KNN).erase (it);
                }
            }
        }
    }

    
    it = KNN->begin();
    for (int i = 0; i < ((int)KNN->size())-1; i++){
        it++;
    }
    
    if(sen == true){
        if(ftau < it->first){
            ftau = it->first;
        }
    }
    else{
        if(ftau > it->first){
            ftau = it->first;
        }
    }
    
    ftau = ftau*(2.0/4.0);

    if(dist-ftau <= maxdis){
        (*sc)++;
        if(dist-ftau <= tau){
            if(lNode == NULL){}
            else{
                lNode->KnnSearch(q, k, KNN, sc, absc, sen);
            }
        }
        if(dist + ftau >= tau){
            if(rNode == NULL){}
            else{
                rNode->KnnSearch(q, k, KNN, sc, absc, sen);
            }
        }
    }
    else{
        (*absc)++;
    }
};



void Node::RangeSearch(NGS *q, double t, std::map<double, std::vector<NGS>, std::greater<double>> *NN, int *sc, int *absc){
    double dist;
    dist = DM.returnDistance(q->returnTra(), vantagepoint->returnTra());
    
    if(dist <= t){
        (*NN)[dist].push_back(*vantagepoint);
    }
    
    double ftau = thresh;
    
    
    if(ftau < dist*(2.0/4.0)){
        ftau = dist*(2.0/4.0);
    }
    
    if(dist-ftau <= maxdis){
        (*sc)++;
        if(dist - ftau <= tau){
            if(lNode == NULL){}
            else{
                lNode->RangeSearch(q, t, NN, sc, absc);
            }
        }
        if(dist + ftau >= tau){
            if(rNode == NULL){}
            else{
                rNode->RangeSearch(q, t, NN, sc, absc);
            }
        }
        
    }
    else{
         (*absc)++;
    }
};

void Node::initiateSearch(NGS *q, int k, std::map<double, std::vector<NGS>, std::greater<double>> *KNN, int *sc, int *absc, bool sen){
    KnnSearch(q, k, KNN, sc, absc, sen);
};

void Node::initiateRangeSearch(NGS *q, double t, std::map<double, std::vector<NGS>, std::greater<double>> *NN, int *sc, int *absc){
    RangeSearch(q, t, NN, sc, absc);
};


void Node::saveNode(std::ofstream *sout){
    std::string str2use;
    int sp = 100;
    (*sout) << "Node_S{\n";
    (*sout) << "Tau_{" << std::setprecision(sp) << tau << "}_Tau\n";
    (*sout) << "Thresh_{" << std::setprecision(sp) << thresh << "}_Thresh\n";
    (*sout) << "Maxdis_{" <<std::setprecision(sp) << maxdis << "}_Maxdis\n";
    if(Ntype == false){
        str2use = "false";
    }
    else{
        str2use = "true";
    }
    (*sout) << "Ntype_{" << str2use<< "}_Ntype\n";
    str2use = "";
    (*sout) << "pp_{" << std::setprecision(sp) << pp << "}_pp\n";
    (*sout) << "VantageIdx_{" << std::setprecision(sp) << (*vantagepoint).returnIdx() << "}_VantageIdx\n";
    
    if(lNode == NULL){
        (*sout) << "lNode_{null|\n";
        (*sout) << "}_lNode\n";
    }
    else{
        (*sout) << "lNode_{exist|\n";
        (*lNode).saveNode(sout);
        (*sout) << "}_lNode\n";
    }
    if(rNode == NULL){
        (*sout) << "rNode_{null|\n";
        (*sout) << "}_rNode\n";
    }
    else{
        (*sout) << "rNode_{exist|\n";
        (*rNode).saveNode(sout);
        (*sout) << "}_rNode\n";
    }
    (*sout) << "}Node_S\n";
};


int Node::loadNode(std::ifstream *sin, std::vector<NGS> *dataPoints){
    int t = 0;
    std::string lineContents = "";
    size_t didx;
    while (!(*sin).eof()){
        getline((*sin), lineContents);
        if(lineContents == "}Node_S"){
            t = 0;
            break;
        }
        else{
            didx = lineContents.find_first_of("{");
            if(lineContents.substr(0,didx+1) == "Tau_{"){
                lineContents = lineContents.substr(didx+1);
                didx = lineContents.find_last_of("}_Tau");
                tau = stod(lineContents.substr(0,didx-4));
            }
            else if(lineContents.substr(0,didx+1) == "Thresh_{"){
                lineContents = lineContents.substr(didx+1);
                didx = lineContents.find_last_of("}_Thresh");
                thresh = stod(lineContents.substr(0,didx-7));
            }
            else if(lineContents.substr(0,didx+1) == "Maxdis_{"){
                lineContents = lineContents.substr(didx+1);
                didx = lineContents.find_last_of("}_Maxdis");
                maxdis = stod(lineContents.substr(0,didx-7));
            }
            else if(lineContents.substr(0,didx+1) == "Ntype_{"){
                lineContents = lineContents.substr(didx+1);
                didx = lineContents.find_last_of("}_Ntype");
                if(lineContents.substr(0,didx-6) == "false"){
                    Ntype = false;
                }
                else{
                    Ntype = true;
                }
            }
            else if(lineContents.substr(0,didx+1) == "pp_{"){
                lineContents = lineContents.substr(didx+1);
                didx = lineContents.find_last_of("}_pp");
                pp = stoi(lineContents.substr(0,didx-3));
            }
            else if(lineContents.substr(0,didx+1) == "VantageIdx_{"){
                lineContents = lineContents.substr(didx+1);
                didx = lineContents.find_last_of("}_VantageIdx");
                int vidx = stoi(lineContents.substr(0,didx-11));
                vantagepoint = &(*dataPoints)[vidx];
            }
            else if(lineContents.substr(0,didx+1) == "lNode_{"){
                lineContents = lineContents.substr(didx+1);
                didx = lineContents.find_last_of("|");
                if(lineContents.substr(0,didx) == "exist"){
                    lNode = new Node;
                    t  = (*lNode).loadNode(sin, dataPoints);
                }
                else{
                    lNode = NULL;
                }
                
            }
            else if(lineContents.substr(0,didx+1) == "rNode_{"){
                lineContents = lineContents.substr(didx+1);
                didx = lineContents.find_last_of("|");
                if(lineContents.substr(0,didx) == "exist"){
                    rNode = new Node;
                    t  = (*rNode).loadNode(sin, dataPoints);
                }
                else{
                    rNode = NULL;
                }
            }
        }
    }
    return t;
};


void Node::averageTau(double *x, double y){
    (*x) += (tau/y);
    if(lNode == NULL){}
    else{
        (*lNode).averageTau(x, y);
    }
    if(rNode == NULL){}
    else{
        (*rNode).averageTau(x, y);
    }
};


void Node::totalTau(long double *x){
    (*x) += (long double) tau;
    if(lNode == NULL){}
    else{
        (*lNode).totalTau(x);
    }
    if(rNode == NULL){}
    else{
        (*rNode).totalTau(x);
    }
};
