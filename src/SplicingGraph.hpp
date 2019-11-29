#ifndef _SPLICING_GRAPH_HPP_
#define _SPLICING_GRAPH_HPP_

#include <iostream>
#include <fstream>
#include <list>
#include <utility>
#include <string>
#include <map>
#include <unordered_map>


#include "sdsl/bit_vectors.hpp"
#include "utils.hpp"


struct Feature {
    int i;
    int index;
    std::string seqid;
    std::string type;
    int start;
    int end;
    bool strand; //+: true; -:false
    std::string id;

    Feature() {
        i = 1;
    }

    Feature(std::string line) {
        i = 1;
        std::string token;
        std::string delimiter = "\t";
        std::size_t pos;
        while ((pos = line.find(delimiter)) != std::string::npos) {
            token = line.substr(0, pos);
            line.erase(0, pos + delimiter.length());
            add(token);
        }
        add(line);
    }

    void add(std::string elem) {
        switch(i) {
        case 1:
            seqid = elem;
            break;
        case 3:
            type = elem;
            break;
        case 4:
            start = std::stoi(elem);
            break;
        case 5:
            end = std::stoi(elem);
            break;
        case 7:
            if(elem.compare("+") == 0) {
                strand = true;
            } else {
                strand = false;
            }
            break;
        case 9:
            bool flag = true;
            //GFF
            //std::string delimiter = ";";
            //std::string string_to_search = "ID=";

            //GTF
            std::string delimiter = "; ";
            std::string string_to_search;
            if(type.compare("gene") == 0) {
                string_to_search = "gene_id \"";
            } else if(type.compare("transcript") == 0) {
                string_to_search = "transcript_id \"";
            } else if(type.compare("exon") == 0) {
                string_to_search = "exon_id \"";
            }
            std::size_t pos;
            std::string token;
            std::string subtoken;
            while((pos = elem.find(delimiter)) != std::string::npos) {
                token = elem.substr(0, pos);
                if(token.substr(0,string_to_search.size()).compare(string_to_search) == 0) {
                    //GTF
                    id = token.substr(string_to_search.size(),token.size()-string_to_search.size()-1);
                    //GFF
                    //id = token.substr(string_to_search.size(),token.size()-string_to_search.size());
                    flag = false;
                    break;
                }
                elem.erase(0, pos + delimiter.length());
            }
            if(flag && elem.substr(0,string_to_search.size()).compare(string_to_search) == 0) {
                //GTF
                id = elem.substr(string_to_search.size(),elem.size()-string_to_search.size()-2);
                //GFF
                //id = elem.substr(string_to_search.size(),elem.size()-string_to_search.size()-1);
                flag = false;
            }
            if(flag) {
                id = ".";
            }
        }
        i++;
    }

    std::string toStr() {
        return seqid + " " +
            type + " " +
            std::to_string(start) + " " +
            std::to_string(end) + " " +
            std::to_string(strand) + " " +
            id;
    }
};

class SplicingGraph {

    struct Gene
    {
        /* data */
        std::string T_g;
        std::string reference_g;
        bool strand_g;
        std::vector<std::pair<int, int> > ExonsPos_g;
        std::vector<std::list<int> > parents_g;
        std::vector<std::list<int> > sons_g;
        int exsN_g;

        int exID_g ;
        std::map<std::string, int> id2index_g;
        int curr_i_g;
        int last_i_g;

        std::vector<std::vector<int> > edges_g;
        sdsl::rrr_vector<> bitVector_g;
        sdsl::rrr_vector<>::select_1_type selectBV_g;
        sdsl::rrr_vector<>::rank_1_type rankBV_g;

        Gene(){
            T_g = "|";
            edges_g.push_back({0});
            edges_g.push_back({0});

            parents_g.push_back({});
            sons_g.push_back({});

            exID_g = 1;
            curr_i_g = 1;
            last_i_g = -1;
        }

    };

   /*sarà il risultato della concatenazione della lisa Gene e Tn(concatazione delle T_g)*/
    std::string T;
    int exsN;


   /*bitvector che contiene il bitvector della T */
    sdsl::rrr_vector<> bitVector;
    sdsl::rrr_vector<>::select_1_type selectBV;
    sdsl::rrr_vector<>::rank_1_type rankBV;

   /*bitvector_sg che concatenazione dei bitvector dei geni*/
    sdsl::rrr_vector<> bitVector_sg;
    sdsl::rrr_vector<>::select_1_type selectBV_sg;
    sdsl::rrr_vector<>::rank_1_type rankBV_sg;

    
    int refLen;  

    /*colezione di geni*/
    std::vector<Gene > vector_gene;

    void setupBitVector();
   // void save(const std::string);
    //void load(const std::string);
    
    public:
        //SplicingGraph(const std::string&);
        SplicingGraph(const char*, const std::string&);
        void concatSqe(std::string &tmp);
        std::string getText() const;
        std::list<int> getParents(const int&,const int&) const;
        std::list<int> getSons(const int&,const int&) const;
        std::string getExon( const int&) const;
        int rank(const int&) const;
        int select(const int&) const;
        int rank_sg(const int&) const;
        int select_sg(const int&) const;
        bool contain(const int&, const int&, const int&) const;
        bool isNew(const int&, const int&, const int&) const;
    
        void print() const;
        int getExonsNumber() const;
        int getGeneNumber() const;
};

#endif
