#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <utility>
#include <list>
#include <stdio.h>
#include <zlib.h>
#include <tuple>

#include "kseq.h"
#include "SplicingGraph.hpp"
#include "bMEM.hpp"
#include "utils.hpp"
#include "MEMsGraph.hpp"

KSEQ_INIT(gzFile, gzread)

std::pair<char, std::list<std::pair<int, std::list<Mem> > > >  analyzeRead(std::list<Mem>& l_memb,
                                                        const SplicingGraph sg,
                                                        const std::string& read,
                                                        const int& L,
                                                        const int& eps,
                                                        const int& exsN,
                                                        const bool& verbose){
    
    std::list<std::pair<int, std::list<Mem> > > paths;

   if(l_memb.size() > 0){
        MemsGraph mg (read, L, eps, exsN,verbose);
        mg.build(sg, l_memb);
        paths = mg.visit(sg);
    }

    char strand  = '+';
    return std::make_pair(strand, paths);

    }

int main(int argc, char* argv[]){

    std::string genomic;
    std::string annotation;
    std::string rna_seqs;

    int L = 0;
    int eps = -1;
    std::string out;
    bool verbose = false;
 

    int c;
    while (1) {
        static struct option long_options[] =
            {
                {"genomic", required_argument, 0, 'g'},
                {"annotation", required_argument, 0, 'a'},
                {"sample",  required_argument, 0, 's'},
                {"L",  required_argument, 0, 'l'},
                {"erate",    required_argument, 0, 'e'},
                {"output", required_argument, 0, 'o'},
                {"help", no_argument, 0, 'h'},
                //{"verbose", no_argument, 0, 'v'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "g:a:s:l:e:o:h", long_options, &option_index);

        if (c == -1) {
            break;
        }

        switch(c) {
        case 'g':
            genomic = optarg;
            break;
        case 'a':
            annotation = optarg;
            break;
        case 's':
            rna_seqs = optarg;
            break;
        case 'l':
            L = std::stoi(optarg);
            break;
        case 'e':
            eps = std::stoi(optarg);
            break;
        case 'o':
            out = std::string(optarg);
            break;
        // case 'v':
        //     verbose = true;
        //     break;
        case 'h':
          //  printHelp();
            exit(EXIT_SUCCESS);
        default:
           // printHelp();
            exit(EXIT_FAILURE);
        }
    }
    if(L == 0) {
        L = 7;
    }
    if(eps < 0 || eps > 100) {
        eps = 3;
    }


    gzFile fastain = gzopen(genomic.c_str(), "r");
    kseq_t *reference = kseq_init(fastain);
    kseq_read(reference);

    /*------------*/

    SplicingGraph sg (reference->seq.s, annotation);

    //sg.print();

   /*std::cout <<"prova rank select" <<std::endl;
   std::cout << sg.rank(210) << std::endl;
   std::cout << sg.select(2) << std::endl;*/

   /*il rank mi da una poiszione in piu rispetto al select*/

    kseq_destroy(reference);
    gzclose(fastain);

    /*setting up MEMs index*/
    BackwardMEM bm (sg.getText(), genomic);

    std::ofstream outFile;
    outFile.open(out);

    fastain = gzopen(rna_seqs.c_str(), "r");

    kseq_t *seqs = kseq_init(fastain);
    
    int l;
    std::string head;
    std::string read;
  

    /**main loop*/

    std::list<Mem> mems;
    std::vector<std::list<Mem> > list_memb;

    int n_gene = sg.getGeneNumber();
    
    while(n_gene > 0){
        list_memb.push_back({mems});
        n_gene--;
    }

    L = 7;

    std::pair<char, std::list<std::pair<int , std::list<Mem> > > >  paths;
    std::vector<std::tuple< std::string, std::string,std::pair<char,std::list <std::pair<int, std::list<Mem> > > > > > vector_path;
    std::list<int>pesi_paths;
    
    while((l = kseq_read(seqs)) >= 0){

        head = seqs -> name.s;
        read = seqs -> seq.s;

        mems = bm.getMEMs(read,L);

        
        for(Mem m : mems) {          
            list_memb[sg.rank_sg(m.t) - 1].push_back(m);
        }

        for(unsigned int i = 0; i < list_memb.size(); i++){
            
            std::pair<char, std::list<std::pair<int , std::list<Mem> > > >  tmp;
            paths = analyzeRead(list_memb[i], sg, read, L, eps, sg.getExonsNumber(), verbose);

            if(paths.second.size() > 0){
                
                vector_path.push_back(std::make_tuple(head, read, paths));
                /*lista dei pesi per poi trovare il peso minimo complessivo*/
                for(auto& x : paths.second){
                    pesi_paths.push_back(x.first);
                }

            }
            
        }

        /*pulisco le liste del dei vettori list_memb per evitare di ciclare sugli elementi gia elaborati*/
        for(unsigned int i = 0; i < list_memb.size(); i++){
            list_memb[i].clear();
        }

    }


    std::vector<std::tuple< std::string, std::string, std::pair<char,std::list <std::pair<int, std::list<Mem> > > > > > optimized_vector_path;
    int min = *min_element(pesi_paths.begin(), pesi_paths.end());
    //std::cout << min << std::endl;

/*scelgo quello col cammino minimo*/
    for(unsigned int i = 0; i < vector_path.size();i++){
        std::pair<char, std::list<std::pair<int , std::list<Mem> > > > tmp = std::get<2>(vector_path[i]);
        for(auto& x : tmp.second){
            if(x.first == min){
                optimized_vector_path.push_back(vector_path[i]);
            }
        }
    }

/*scrivo su file*/
    for(unsigned int i = 0; i < optimized_vector_path.size(); i++){
        std::pair<char, std::list<std::pair<int , std::list<Mem> > > > tmp = std::get<2>(optimized_vector_path[i]);
        outFile << tmp.first << " ";
        outFile << std::get<0>(optimized_vector_path[i]) << " ";
        for(auto& x : tmp.second){
            outFile << x.first << " ";
            for(auto& y :x.second){
                outFile << y.toStr() << " ";
            }
        }

        outFile <<std::endl;
        outFile<<std::get<1>(optimized_vector_path[i]) << std::endl;
    }

    
/*da fare ancora lo strand*/
    kseq_destroy(seqs);
    gzclose (fastain);
   
    return 0;
}
