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

int peso_min(std::vector<std::pair<char,std::list <std::pair<int, std::list<Mem> > > > > & vt){
    std::list<int> n_pesi;
    bool t = true;
    for(auto& x : vt){

        for(auto& y : x.second){
            /*cerco il minimo*/
   //         std::cout << y.first << " ";
            n_pesi.push_back(y.first);
        }
  //          std::cout << std::endl;
    }
        
    int min_peso = *min_element(n_pesi.begin(), n_pesi.end());
   // std::cout << "min : " << min_peso << "." << std::endl;

    return min_peso;

}

std::vector<std::list<Mem>> get_list_memb(std::list<Mem>& m, const SplicingGraph& sg, int& n_gene){
    
    std::vector<std::list<Mem> > tmp;

    while(n_gene > 0){
        tmp.push_back({});
        n_gene--;
    }

    for(Mem x : m){
        tmp[sg.rank_sg(x.t) -1].push_back(x);
    }

    return tmp;
} 

std::vector<std::pair<char, std::list<std::pair<int, std::list<Mem> > > > >  analyzeRead(BackwardMEM& bm,
                                                        const SplicingGraph sg,
                                                        const std::string& read,
                                                        const int& L,
                                                        const int& eps,
                                                        const int& exsN,
                                                        const bool& verbose){
        
        std::list<Mem> mems=bm.getMEMs(read, L);
        std::vector<std::list<Mem> > list_mem;
        std::list<std::pair<int, std::list<Mem> > > paths;
        std::vector<std::list<std::pair<int, std::list<Mem> > > >   vector_paths;
        std::vector<std::pair<char, std::list<std::pair<int, std::list<Mem> > > > > vector_result;
        char strand = '/';

        int n_gene = sg.getGeneNumber();

        if(!mems.empty()){    
            list_mem = get_list_memb(mems, sg, n_gene);

            for(auto& x : list_mem){
                MemsGraph mg(read, L, eps, sg.getExonsNumber(), verbose);
                mg.build(sg,x);
                paths = mg.visit(sg);
                vector_paths.push_back(paths);
            }
        }

        /*reversed and complementd read*/
        n_gene = sg.getGeneNumber();
        std::string readRC = reverseAndComplement(read);
        std::list<Mem> memsRC = bm.getMEMs(readRC, L);
        std::list<std::pair<int, std::list<Mem> > > pathsRC;
        std::vector<std::list<std::pair<int, std::list<Mem> > > > vector_paths_RC;
        if(!memsRC.empty()){
            list_mem = get_list_memb(memsRC, sg, n_gene);
            for(auto& x : list_mem){
                MemsGraph mgRC(readRC, L, eps, sg.getExonsNumber(), verbose);
                mgRC.build(sg, x);
                pathsRC = mgRC.visit(sg);
                vector_paths_RC.push_back(pathsRC);
            }
        }

      //  std::cout << "grandezza vector_paths " << vector_paths.size() << " gramdezza vector_paths_RC " << vector_paths_RC.size() << std::endl;
            for(unsigned int i = 0; i < sg.getGeneNumber(); i++){
                if(!vector_paths[i].empty() || !vector_paths_RC[i].empty()){
                    if(!vector_paths[i].empty() && vector_paths_RC[i].empty()){
                        strand = '+';
                        vector_result.push_back(std::make_pair(strand, vector_paths[i]));
                    }
                    else{
                        if(vector_paths[i].empty() && !vector_paths_RC[i].empty()){
                            strand = '-';
                            vector_result.push_back(std::make_pair(strand, vector_paths_RC[i]));
                        }
                        else{
                            std::list<std::pair<int, std::list<Mem> > > pt = vector_paths[i];
                            std::list<std::pair<int, std::list<Mem> > > pt_RC = vector_paths_RC[i];
                            if(pt.front().first <= pt_RC.front().first){
                                strand = '+';
                                vector_result.push_back(std::make_pair(strand, vector_paths[i]));
                            }
                            else{
                                strand = '-';
                                vector_result.push_back(std::make_pair(strand, vector_paths_RC[i]));
                            }
                        }
                    }
                }
            }
    return vector_result;
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

    std::pair<char, std::list<std::pair<int , std::list<Mem> > > >  paths;
    std::pair<char, std::list<std::pair<int , std::list<Mem> > > >  paths_rc;
    std::vector<std::pair<char,std::list <std::pair<int, std::list<Mem> > > > >  vector_path;
    std::list<int>pesi_paths;
    
    std::map<std::string, std::vector<std::list<Mem>>> control_double;
    L = 7;

    while((l = kseq_read(seqs)) >= 0){

        head = seqs -> name.s;
        read = seqs -> seq.s;

        vector_path = analyzeRead(bm, sg, read, L, eps, sg.getExonsNumber(), verbose);

        int min = peso_min(vector_path);

        /********manca la scelta del minore tra i vettori di path****************/
 
        /*scrittura su file*/
        int index_vector = 0;
        for(auto&x : vector_path){
            bool t = true;
            if(!x.second.empty()){
                //outFile << x.first << " " << head << " ";                
                for(auto& y: x.second){
                    if(y.first == min){
                        if(t){
                            outFile << x.first << " " << head << " ";
                            t = false;
                        }
                        /*peso*/
                        outFile << y.first << " ";
    
                        for(auto& z : y.second){
                            index_vector = sg.rank_sg(z.t) -1;
                            outFile << z.toStr() << " ";
                        } 
                        outFile << index_vector;
                        outFile << std::endl;
                        if(x.first == '+'){
                            outFile << read;
                            }
                        else{
                            outFile << reverseAndComplement(read);
                        }
                    }                        
                outFile<<std::endl;
                }
            }
        }   
    }

    kseq_destroy(seqs);
    gzclose (fastain);
   
    return 0;
}