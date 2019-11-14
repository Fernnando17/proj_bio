#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <utility>
#include <list>
#include <stdio.h>



#include "SplicingGraph.hpp"

int main(){

    std::string genomic;
    std::string annotation;
    std::string rna_seqs;

    std::ifstream input_file;
    input_file.open("./input/Danio_rerio.GRCz11.dna.chromosome.23.fa");
    std::stringstream strStream;
    strStream << input_file.rdbuf();
    genomic = strStream.str();

    annotation = "./input/chr23.small.gtf";
    SplicingGraph sg (genomic.c_str(), annotation);

    sg.print();
/*prova
    std::cout << sg.select(2) << std::endl;
    std::cout << sg.rank(210) << std::endl;
*/
    input_file.close();
   
    return 0;
}