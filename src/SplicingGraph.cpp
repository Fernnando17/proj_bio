#include "SplicingGraph.hpp"

std::string getExonID(int s, int e) {
    return std::to_string(s) + ":" + std::to_string(e);
}

SplicingGraph::SplicingGraph(const char *genomic,
                             const std::string &gtf) {
    refLen = strlen(genomic);

    std::ifstream gtfFile;

/*   edges.push_back({0});
    edges.push_back({0});

    parents.push_back({});
    sons.push_back({});
    
    T = "|";

    int exID = 1;
    std::map<std::string, int> id2index;
    int curr_i = 1;
    int last_i = -1;
 */  
    int index_vector = -1;
    std::string line;
    gtfFile.open(gtf);
    if(gtfFile.is_open()) {
        while(getline(gtfFile,line)) {
            Feature feat (line);
            
            if(feat.type.compare("gene") == 0) {
                index_vector++;
                Gene g;
                vector_gene.push_back(g);
                vector_gene[index_vector].reference_g = feat.seqid;
                vector_gene[index_vector].strand_g = feat.strand;
            }
            else if(feat.type.compare("transcript") == 0 || feat.type.compare("mRNA") == 0) {
                vector_gene[index_vector].last_i_g = -1;
            }
            else if(feat.type.compare("exon") == 0) {
                std::string posID = getExonID(feat.start, feat.end);

                if (vector_gene[index_vector].id2index_g.find(posID) == vector_gene[index_vector].id2index_g.end()) {
                    vector_gene[index_vector].id2index_g[posID] = vector_gene[index_vector].exID_g;
                    std::string currExString(genomic + feat.start-1, feat.end-feat.start+1);
                    vector_gene[index_vector].T_g += currExString + "|";
                    vector_gene[index_vector].ExonsPos_g.push_back(std::make_pair(feat.start, feat.end));
                    vector_gene[index_vector].parents_g.push_back({});
                    vector_gene[index_vector].sons_g.push_back({});
                    vector_gene[index_vector].exID_g++;
                }

                vector_gene[index_vector].curr_i_g = vector_gene[index_vector].id2index_g[posID];
                if(vector_gene[index_vector].last_i_g != -1) {
                    if(vector_gene[index_vector].last_i_g != vector_gene[index_vector].curr_i_g) {
                        if(vector_gene[index_vector].last_i_g >= (int)vector_gene[index_vector].edges_g.size()) {
                            int i = vector_gene[index_vector].edges_g.size();
                            while(i<=vector_gene[index_vector].last_i_g + 1) {
                                vector_gene[index_vector].edges_g.push_back({0});
                                ++i;
                            }
                        }

                        if(vector_gene[index_vector].curr_i_g >= (int)vector_gene[index_vector].edges_g[vector_gene[index_vector].last_i_g].size()) {
                            int i = vector_gene[index_vector].edges_g[vector_gene[index_vector].last_i_g].size();
                            while(i<=vector_gene[index_vector].curr_i_g+1) {
                                vector_gene[index_vector].edges_g[vector_gene[index_vector].last_i_g].push_back(0);
                                ++i;
                            }
                        }
                        //edges[last_i] vector<vector<int>> rightValidVariants(1, {i});
                        vector_gene[index_vector].edges_g[vector_gene[index_vector].last_i_g][vector_gene[index_vector].curr_i_g] = 1;
                        if(vector_gene[index_vector].curr_i_g >= (int)vector_gene[index_vector].parents_g.size()) {
                            int i = 0;
                            while(i<=vector_gene[index_vector].curr_i_g+1) {
                                vector_gene[index_vector].parents_g.push_back({});
                                ++i;
                            }
                        }
                        vector_gene[index_vector].parents_g[vector_gene[index_vector].curr_i_g].push_back(vector_gene[index_vector].last_i_g);

                        if(vector_gene[index_vector].last_i_g >= (int)vector_gene[index_vector].sons_g.size()) {
                            int i = 0;
                            while(i<=vector_gene[index_vector].last_i_g + 1) {
                                vector_gene[index_vector].sons_g.push_back({});
                                ++i;
                            }
                        }
                        vector_gene[index_vector].sons_g[vector_gene[index_vector].last_i_g].push_back(vector_gene[index_vector].curr_i_g);
                    }
                }
                vector_gene[index_vector].last_i_g = vector_gene[index_vector].curr_i_g;
            }

        }


    /*print gei geni*/
    /*    for(unsigned int i = 0; i < vector_gene.size(); i++){
        
            std::cout << vector_gene[i].T_g << std::endl;

           std::cout << i << std::endl;
            int j = 0;
            int range = vector_gene[i].bitVector_g.size();
            while(j < range){
                std::cout << vector_gene[i].bitVector_g[j];
                j++;
            }

            std::cout << std::endl;
        
        }

        std::cout << vector_gene.size()<< std::endl;    */
    //---------------------------

    }


    // Transitive closure on the graph
    for(unsigned int k = 0; k < vector_gene.size(); k++){
        int i = 1;
        for(const std::pair<int,int>& p1 : vector_gene[k].ExonsPos_g) {
            if(i>=(int)vector_gene[k].edges_g.size()) {
                int i_ = vector_gene[k].edges_g.size();
                while(i_<=i+1) {
                    vector_gene[k].edges_g.push_back({0});
                    ++i_;
                }
            }
            int j = 1;
            for(const std::pair<int,int>& p2 : vector_gene[k].ExonsPos_g) {
                if(p1.second <= p2.first) {
                    if(j>=(int)vector_gene[k].edges_g[i].size()) {
                        int j_ = vector_gene[k].edges_g[i].size();
                        while(j_<=j+1) {
                            vector_gene[k].edges_g[i].push_back(0);
                            ++j_;
                        }
                    }

                    if(vector_gene[k].edges_g[i][j] == 0) {
                        vector_gene[k].edges_g[i][j] = 2;
                        vector_gene[k].parents_g[j].push_back(i);
                        vector_gene[k].sons_g[i].push_back(j);
                    }
                }
                ++j;
            }
            ++i;
        }
    }

    gtfFile.close();
    setupBitVector();
    /*save(gtf);*/
}

void SplicingGraph::setupBitVector() {
   
    /*costruzione del bitvector di ogni gene del vector_gene*/

    for(unsigned int j = 0; j < vector_gene.size(); j++){

        sdsl::bit_vector BV (vector_gene[j].T_g.length(), 0);
        unsigned int i = 0;
        while(i<vector_gene[j].T_g.length()) {
            if(vector_gene[j].T_g[i] == '|') {
                BV[i] = 1;
            }
            i++;
        }

        vector_gene[j].bitVector_g = sdsl::rrr_vector<>(BV);
        vector_gene[j].selectBV_g = sdsl::rrr_vector<>::select_1_type(&vector_gene[j].bitVector_g);
        vector_gene[j].rankBV_g = sdsl::rrr_vector<>::rank_1_type(&vector_gene[j].bitVector_g);
    
    }

    /*concate delle stringhe T*/
    concatSqe(T);

/*--------------------------------*/

    /*creazione dei bitvector su T*/
    sdsl::bit_vector BV (T.size(), 0);
    unsigned int i = 0;
    while(i < T.length()){
        if(T[i] == '|'){
            BV[i] =1;
        }
        i++;
    }

    bitVector = sdsl::rrr_vector<>(BV);
    selectBV = sdsl::rrr_vector<>::select_1_type(&bitVector);
    rankBV = sdsl::rrr_vector<>::rank_1_type(&bitVector);

/*--------------------------*/

    /*bitvector delle concatenazioni dei bitvector*/
    sdsl::bit_vector BV_2 (T.size(), 0);

    BV_2[0] = 1;
    int pos_concat = 0;

   for(unsigned int i = 0; i  < vector_gene.size(); i++){    
        pos_concat += (vector_gene[i].T_g.length() - 1 );
      //  std::cout << T[pos_concat]<< std::endl;
        if(T[pos_concat]== '|'){   
            BV_2[pos_concat] = 1; 
        }
    }

    bitVector_sg = sdsl::rrr_vector<>(BV_2);
    selectBV_sg = sdsl::rrr_vector<>::select_1_type(&bitVector_sg);
    rankBV_sg = sdsl::rrr_vector<>::rank_1_type(&bitVector_sg);

}

void SplicingGraph::concatSqe(std::string &T){
    
    T = "";
    for(unsigned int i = 0; i < vector_gene.size(); i++){
        T += vector_gene[i].T_g.substr(0, vector_gene[i].T_g.length()-1); 
    }
    T += "|";

}

std::string SplicingGraph::getText() const {
    return T;
}

/*std::list<int> SplicingGraph::getParents(const int& i) const {
    return parents[i];
}*/

/*std::list<int> SplicingGraph::getSons(const int& i) const {
    return sons[i];
}*/

/*std::string SplicingGraph::getExon(const int& i) const {
    int s = selectBV(i)+1;
    int e = selectBV(i+1);
    std::string exonText = T.substr(s, e-s);
    return exonText;
}*/

int SplicingGraph::rank(const int& i) const {
    return rankBV(i);
}

int SplicingGraph::select(const int& i) const {
    return selectBV(i);
}

/*bool SplicingGraph::contain(const int& x, const int& y) const {
    if(edges[x][y] >= 1) {
        return true;
    } else {
        return false;
    }
}*/

/*bool SplicingGraph::isNew(const int& x, const int& y) const {
    if(edges[x][y] > 1) {
        return true;
    } else {
        return false;
    }
}*/

void SplicingGraph::print() const {
    
    std::cout << T << std::endl;
    std::cout << std::endl;
    
    /*stampa bitvector*/
    std::cout << "bitvector" << std::endl; 
    unsigned int i = 0;
    while(i<bitVector_sg.size()) {
        std::cout << bitVector[i];
        i++;
    }
    std::cout << std::endl;
    std::cout << std::endl;


    std::cout << std::endl;
    
    /*stampa bitvector_sg*/
     std::cout<< "bitvector_sg"<< std::endl;
    unsigned int i_2 = 0;
    while(i_2 <bitVector_sg.size()) {
        std::cout << bitVector_sg[i_2];
        i_2++;
    }
    std::cout << std::endl;
    std::cout << std::endl;


    /*stamp di (parents_g, edges_g, sons_g) dei geni*/
    for(unsigned int i = 0; i < vector_gene.size(); i++){

        std::cout << "gene " << i << std::endl;
        std::cout << "edges_g"<< i <<  std::endl;
        
        std::cout << std::endl;
        
        for(std::vector<int> v : vector_gene[i].edges_g) {
            for(const int& e : v) {
                std::cout << e << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        std::cout << "parents_g " << i << std::endl;
        std::cout << std::endl;

        for(auto e : vector_gene[i].parents_g) {
            for(auto p : e) {
                std::cout << p << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

        std::cout << "sons_g" << i << std::endl;
        std::cout << std::endl;
        for(auto e : vector_gene[i].sons_g) {
            for(auto s : e) {
                std::cout << s << " ";
            }
            std::cout << std::endl;
        }

        std::cout << std::endl;
    }

}

/*void SplicingGraph::save(const std::string path) {
    std::ofstream ofile;
    ofile.open(path + ".sg");
    ofile << reference << " " << refLen << " ";
    if(strand)
      ofile << "+";
    else
      ofile << "-";
    ofile << "\n";
    ofile << T << "\n";
    ofile << exsN << "\n";
    for(const std::vector<int>& v : edges) {
        for(const int& e : v) {
            ofile << e << " ";
        }
        ofile << "; ";
    }
    ofile << "\n";
    for(const std::pair<int,int>& p : ExonsPos) {
        ofile << p.first << "," << p.second << " ";
    }
    ofile << "\n";
    // for(const std::string& name : ExonsName) {
    //     ofile << name << " ";
    // }
    ofile << "\n";
    ofile.close();
}

int SplicingGraph::getExonsNumber() const {
    return exsN;
}*/
