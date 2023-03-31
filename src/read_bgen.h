#ifndef READBGEN_H
#define READBGEN_H

#include <Rcpp.h>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <zlib.h>

#include "libdeflate-1.6/libdeflate.h"
#include "zstd-1.4.5/lib/zstd.h"

typedef unsigned int  uint;
typedef unsigned char uchar;
typedef long long unsigned int llui;
typedef unsigned long long int ulli;

using namespace Rcpp;
class BGEN {
  public:
    std::string bgenFile;
    
    uint offset;
    uint Mbgen;
    uint Nbgen;
    uint Layout;
    uint Compression;
    uint SampleIdentifiers;
    Rcpp::CharacterVector sampleID;

    uint Counter;
};

Rcpp::List query_bgen11(SEXP bgenR_in, SEXP seek_in);
Rcpp::List query_bgen13(SEXP bgenR_in, SEXP seek_in);
Rcpp::List query_bgen13_zlib(SEXP bgenR_in, SEXP seek_in);
std::vector<llui>* get_bytes(bool getIndices_in, FILE* & bStream_in, uint offset_in, uint layout_in, uint M_in, uint N_in, uint compression_in);
#endif



