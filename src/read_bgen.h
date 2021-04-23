#ifndef READBGEN_H
#define READBGEN_H

#include <Rcpp.h>
#include <fstream>
#include <string>
#include <vector>

typedef unsigned int   uint;

using namespace Rcpp;
class BGEN {
  
  
  public:
    
    // For file
    std::string bgenFile;
    // For BGEN offset
    uint offset;
    
    // For BGEN header block
    uint Mbgen;
    uint Nbgen;
    uint Layout;
    uint Compression;
    uint SampleIdentifiers;
    Rcpp::CharacterVector sampleID;

    uint Counter;
};

#endif



