#include "read_bgen.h"
#include "libdeflate-1.6/libdeflate.h"
#include "zstd-1.4.5/lib/zstd.h"
#include <Rcpp.h>
#include <stdio.h>
#include <zlib.h>




typedef unsigned int   uint;
typedef unsigned char  uchar;
bool bgenIsOpen = false;

using namespace std;
using namespace Rcpp;



Rcpp::List query_bgen11(SEXP bgenR_in, SEXP seek_in);
Rcpp::List query_bgen13(SEXP bgenR, SEXP seek_in);

vector<long long unsigned int>* get_bytes(FILE* & bStream_in, uint offset_in, uint Layout_in, uint M_in, uint Compression_in);

Rcpp::String close_bgen(SEXP bgenR_in){

    Rcpp::List bgenR_list(bgenR_in);
    FILE* bStream = (FILE*)R_ExternalPtrAddr(bgenR_list["fin"]);
    vector<long long unsigned int>* bytes = (vector<long long unsigned int>*)R_ExternalPtrAddr(bgenR_list["fbytes"]);
    
    if(bStream == NULL) {
       Rcpp::stop("BGEN file is already closed.");
    } else {
      fseek(bStream, 0, SEEK_END);
      fclose(bStream);
      
      R_ClearExternalPtr(bgenR_list["fin"]);
    }

    
    R_ClearExternalPtr(bgenR_list["fcounter"]);
    if (bytes != NULL) {
      R_ClearExternalPtr(bgenR_list["fbytes"]);
    }
    
    bgenIsOpen = false;
    return("BGEN file closed successfully.");
}



Rcpp::List open_bgen(SEXP bgenfile_in, SEXP bytes_in){
  
  
    if (bgenIsOpen) {
      Rcpp::stop("A BGEN file is already open. Use close_bgen() to close the open BGEN file.");
    }
    
    uint maxLA = 65536;
    string bgenfile = Rcpp::as<string>(bgenfile_in);
    bool getBytes   = Rcpp::as<bool>(bytes_in);
    
    FILE* bStream = fopen(bgenfile.c_str(), "rb");
    if (!bStream) { 
        Rcpp::stop("Cannot not open BGEN file: " + bgenfile + "."); 
    }
    bgenIsOpen = true;
    
    // Header Block
    uint offset;
    if (!fread(&offset, 4, 1, bStream)) {
        Rcpp::stop("Cannot read BGEN header block (offset).");
    }
  
    uint L_H; 
    if (!fread(&L_H, 4, 1, bStream)) {
        Rcpp::stop("Cannot read BGEN header block (LH).");
    }
  
    uint Mbgen;
    if (!fread(&Mbgen, 4, 1, bStream)) {
        Rcpp::stop("Cannot read BGEN header block (M).");
    }
  
    uint Nbgen;
    if (!fread(&Nbgen, 4, 1, bStream)) {
        Rcpp::stop("Cannot read BGEN header block (N).");
    }
  
    char magic[5]; 
    if (!fread(magic, 1, 4, bStream)) {
        Rcpp::stop("Cannot read BGEN header block (magic numbers).");
    }
    magic[4] = '\0';
    if (!(magic[0] == 'b' && magic[1] == 'g' && magic[2] == 'e' && magic[3] == 'n')) {
        Rcpp::stop("BGEN file's magic number does not match 'b', 'g', 'e', 'n'.");
    }
 
    fseek(bStream, L_H - 20, SEEK_CUR);
  
    uint flags; 
    if (!fread(&flags, 4, 1, bStream)) {
        Rcpp::stop("Cannot read BGEN header block (flags).");
    }
    
    
    // Header Block Flags
    uint Compression = flags & 3;
    if (Compression != 0 && Compression != 1 && Compression != 2) {
        Rcpp::stop("BGEN compression flag should be 0, 1, or 2.");
    }
  
    uint Layout = (flags >> 2) & 0xf;
    if (Layout != 1 && Layout != 2) {
        Rcpp::stop("BGEN layout flag should be 1 or 2.");
    }
  
    uint SampleIdentifiers = flags >> 31;
    if (SampleIdentifiers != 0 && SampleIdentifiers != 1) {
        Rcpp::stop("BGEN sample identifier flag should be 0 or 1.");
    }
    
    Rcpp::CharacterVector sampleID;
    if (SampleIdentifiers == 1){
        sampleID = CharacterVector(Nbgen);
        uint LS1;
        if (!fread(&LS1, 4, 1, bStream)) {
            Rcpp::stop("Cannot read BGEN sample block (LS).");
        }
       
        uint Nsamples;
        if (!fread(&Nsamples, 4, 1, bStream)) {
            Rcpp::stop("Cannot read BGEN sample block (N).");
        }
        if (Nsamples != Nbgen) {
            Rcpp::stop("Number of sample identifiers (" + std::to_string(Nsamples) + ") does not match the number of BGEN samples (" + std::to_string(Nbgen) + ") in the header block.");
        }
      
      
        char* samID = new char[maxLA + 1];
        int ret;
        for (uint n = 0; n < Nsamples; n++) {
             uint16_t LSID;
             ret = fread(&LSID, 2, 1, bStream);
             ret = fread(samID, 1, LSID, bStream);
             samID[LSID] = '\0';
             Rcpp::String samID_2 = string(samID);
             sampleID[n] = samID_2;
        }
        (void)ret;
        delete[] samID;
    } else {
        sampleID = CharacterVector(1);
    }

    vector<long long unsigned int>* bytesPtr = nullptr;
    if (getBytes) {
      bytesPtr = get_bytes(bStream, offset, Layout, Mbgen, Compression);
    }    

    std::vector<uint>* Counter = new std::vector<uint>(1);
    SEXP fin      = R_MakeExternalPtr(bStream, R_NilValue, R_NilValue);
    SEXP fbytes   = R_MakeExternalPtr(bytesPtr, R_NilValue, R_NilValue);
    SEXP fcounter = R_MakeExternalPtr(Counter, R_NilValue, R_NilValue);
    
    fseek(bStream, offset + 4, SEEK_SET);
    return(Rcpp::List::create(Named("offset")           = offset,
                              Named("M")                = Mbgen,
                              Named("N")                = Nbgen,
                              Named("Compression")      = Compression,
                              Named("Layout")           = Layout,
                              Named("SampleIdentifier") = SampleIdentifiers,
                              Named("SampleID")         = sampleID,
                              Named("fin")              = fin,
                              Named("fbytes")           = fbytes,
                              Named("fcounter")         = fcounter
                              ));
  
}




void Bgen13GetTwoVals(const unsigned char* prob_start, uint32_t bit_precision, uintptr_t offset, uintptr_t* first_val_ptr, uintptr_t* second_val_ptr) {
  
    switch(bit_precision) {
    case 8:
      *first_val_ptr  = prob_start[0];
      prob_start += offset;
      *second_val_ptr = prob_start[0];
      break;
    case 16:
      *first_val_ptr  = prob_start[0]|(prob_start[1]<<8);
      prob_start += offset;
      *second_val_ptr = prob_start[0]|(prob_start[1]<<8);
      break;
    case 24:
      *first_val_ptr  = prob_start[0]|(prob_start[1]<<8)|(prob_start[2]<<16);
      prob_start += offset;
      *second_val_ptr = prob_start[0]|(prob_start[1]<<8)|(prob_start[2]<<16);
      break;
    case 32:
      *first_val_ptr  = prob_start[0]|(prob_start[1]<<8)|(prob_start[2]<<16)|(prob_start[3]<<24);
      prob_start += offset;
      *second_val_ptr = prob_start[0]|(prob_start[1]<<8)|(prob_start[2]<<16)|(prob_start[3]<<24);
    default:
      break;
    }

}


Rcpp::List query_bgen(SEXP bgenR_in, SEXP seek_in){
    Rcpp::List bgenR_list(bgenR_in);
    uint Layout = Rcpp::as<uint>(bgenR_list["Layout"]);
    if(Layout == 2){
      return(query_bgen13(bgenR_in, seek_in));
    }
    else if(Layout == 1){
      return(query_bgen11(bgenR_in, seek_in));
    }

    return(NA_REAL);
} 



Rcpp::List query_bgen13(SEXP bgenR_in, SEXP seek_in){

       Rcpp::List bgenR_list(bgenR_in);
       FILE* bStream = (FILE*)R_ExternalPtrAddr(bgenR_list["fin"]);
       if (bStream == NULL) {
         Rcpp::stop("BGEN file is not open. Please reopen the BGEN file with open_bgen().");
       }
       vector<uint>* counter = (vector<uint>*)R_ExternalPtrAddr(bgenR_list["fcounter"]);
       vector<unsigned long long int>* bytes= (vector<unsigned long long int>*)R_ExternalPtrAddr(bgenR_list["fbytes"]);
       vector<uint>&vr  = *counter;
       vector<unsigned long long int>&br  = *bytes;
       uint nvars  = Rcpp::as<uint>(bgenR_list["M"]);
       uint Nsamples    = Rcpp::as<uint>(bgenR_list["N"]);
       uint Compression = Rcpp::as<uint>(bgenR_list["Compression"]);
       std::vector<uchar> zBuf12;
       std::vector<uchar> shortBuf12;
       uint maxLA = 65536;
       
       vr[0]++;
       if (vr[0] >= (nvars + 1)){
           Rcpp::stop("End of BGEN file has already been reached. Please close the file with close_bgen().");
       }

       NumericMatrix probs(Nsamples, 2);
       NumericVector dosVec(Nsamples);
       int ret;
       
       uint seek = Rcpp::as<uint>(seek_in);
       if (seek != 0) {
         if (bytes == NULL) {
           Rcpp::stop("\nERROR: Reopen the BGEN file with `bytes=TRUE` in open_bgen() to use the `seek` argument in this function.");
         }
         fseek(bStream, br[seek-1], SEEK_SET);
       }
       
       struct libdeflate_decompressor* decompressor = libdeflate_alloc_decompressor();
       
       uint16_t LS;
       ret = fread(&LS, 2, 1, bStream);
       char* snpID = new char[LS + 1];
       ret = fread(snpID, 1, LS, bStream);
       snpID[LS] = '\0';
       Rcpp::String snpID_2 = string(snpID);

       uint16_t LR;
       ret = fread(&LR, 2, 1, bStream);
       char* rsID = new char[LR + 1];
       ret = fread(rsID, 1, LR, bStream);
       rsID[LR] = '\0';
       Rcpp::String rsID_2 = string(rsID);

       uint16_t LC;
       ret = fread(&LC, 2, 1, bStream);
       char* chrStr = new char[LC + 1];
       ret = fread(chrStr, 1, LC, bStream);
       chrStr[LC] = '\0';
       Rcpp::String chrStr_2 = string(chrStr);

       uint32_t physpos;
       ret = fread(&physpos, 4, 1, bStream);

       uint16_t LKnum;
       ret = fread(&LKnum,   2, 1, bStream);
       if (LKnum != 2U) {
           char* allele1 = new char[maxLA + 1];
           for (uint16_t a = 0; a < LKnum; a++) {
                uint32_t LA;
                ret = fread(&LA, 4, 1, bStream);
                ret = fread(allele1, 1, LA, bStream);
           }

           if (Compression > 0) {
               uint zLen;
               ret = fread(&zLen, 4, 1, bStream);
               fseek(bStream, 4 + zLen - 4, SEEK_CUR);

           }
           else {
               uint zLen;
               ret = fread(&zLen, 4, 1, bStream);
               fseek(bStream, zLen, SEEK_CUR);
           }
           delete[] allele1;
           return(Rcpp::List::create(Named("SNPID") = snpID_2,
                                     Named("RSID") = rsID_2,
                                     Named("Chromosome") = chrStr_2,
                                     Named("Position") = physpos,
                                     Named("Alleles") = LKnum,
                                     Named("Allele1") = NA_REAL,
                                     Named("Allele2") = NA_REAL,
                                     Named("AF") = NA_REAL,
                                     Named("Probabilities") = NA_REAL,
                                     Named("Dosages") = NA_REAL));
       }

       uint32_t LA;
       ret = fread(&LA, 4, 1, bStream);
       char* allele1 = new char[LA + 1];
       ret = fread(allele1, 1, LA, bStream);
       allele1[LA] = '\0';
       Rcpp::String allele1_2 = string(allele1);

       uint32_t LB;
       ret = fread(&LB, 4, 1, bStream);
       char* allele2 = new char[LB + 1];
       ret = fread(allele2, 1, LB, bStream);
       allele2[LB] = '\0';
       Rcpp::String allele2_2 = string(allele2);

       uchar* prob_start;
       uint cLen;
       ret = fread(&cLen, 4, 1, bStream);
       if (Compression == 1) {
           zBuf12.resize(cLen - 4);
           uint dLen;
           ret = fread(&dLen, 4, 1, bStream);
           ret = fread(&zBuf12[0], 1, cLen - 4, bStream);
           shortBuf12.resize(dLen);

           uLongf destLen = dLen;
           if (libdeflate_zlib_decompress(decompressor, &zBuf12[0], cLen - 4, &shortBuf12[0], destLen, NULL) != LIBDEFLATE_SUCCESS) {
               Rcpp::stop("Decompressing " + string(rsID) + "genotype block failed with libdeflate.");
           }
           prob_start = &shortBuf12[0];
       }
       else if (Compression == 2) {
           zBuf12.resize(cLen - 4);
           uint dLen;
           ret = fread(&dLen, 4, 1, bStream);
           ret = fread(&zBuf12[0], 1, cLen - 4, bStream);
           shortBuf12.resize(dLen);

           uLongf destLen = dLen;
           size_t ret = ZSTD_decompress(&shortBuf12[0], destLen, &zBuf12[0], cLen - 4);
           if (ret > destLen) {
               if (ZSTD_isError(ret)) {
                   Rcpp::stop("Decompressing " + string(rsID) + "genotype block failed with zstd.");
               }
           }
           prob_start = &shortBuf12[0];
       }
       else {
           zBuf12.resize(cLen);
           ret = fread(&zBuf12[0], 1, cLen, bStream);
           prob_start = &zBuf12[0];
       }
       (void)ret;

       uint32_t N;
       memcpy(&N, prob_start, sizeof(int32_t));
       uint16_t K;
       memcpy(&K, &(prob_start[4]), sizeof(int16_t));

       const uint32_t min_ploidy = prob_start[6];
       if (min_ploidy != 2) {
           Rcpp::stop("Variants with ploidy != 2 is currently not supported.");
       }
       const uint32_t max_ploidy = prob_start[7];
       if (max_ploidy != 2) {
           Rcpp::stop("Variants with ploidy != 2 is currently not supported.");
       }

       const unsigned char* missing_and_ploidy_info = &(prob_start[8]);
       const unsigned char* probs_start = &(prob_start[10 + N]);
       const uint32_t is_phased = probs_start[-2];
       if (is_phased != 0 && is_phased != 1) {
           Rcpp::stop("phased value must be 0 or 1.");
       } 

       const uint32_t B = probs_start[-1];
       if (B != 8 && B != 16 && B != 24 && B != 32) {
           Rcpp::stop("Bits to store probabilities must be 8, 16, 24, or 32.");
       }
       const uintptr_t numer_mask = (1U << B) - 1;
       const uintptr_t probs_offset = B / 8;


       double gmean = 0.0;
       double nmiss = 0.0;
       if(!is_phased){

          for (uint32_t i = 0; i < N; i++) {
               const uint32_t missing_and_ploidy = missing_and_ploidy_info[i];


               if(missing_and_ploidy == 2){
                  uintptr_t numer_aa;
                  uintptr_t numer_ab;
 
                  Bgen13GetTwoVals(probs_start, B, probs_offset, &numer_aa, &numer_ab);
                  probs_start += (probs_offset * 2);
 
                  double p11 = numer_aa / double(1.0 * (numer_mask));
                  double p10 = numer_ab / double(1.0 * (numer_mask));
                  double dosage = 2 * (1 - p11 - p10) + p10;

                  probs(i, 0) = p11;
                  probs(i, 1) = p10;
                  dosVec[i] = dosage;
                  gmean += dosage;

              }
              else if (missing_and_ploidy == 130) {
                  probs(i, 0) = NA_REAL;
                  probs(i, 1) = NA_REAL;
                  dosVec[i]   = NA_REAL;
                  nmiss+=1.0;
              }
              else {
                  Rcpp::stop("Variants with ploidy != 2 is currently not supported.");
              }
          }

       } else if(is_phased){
            for (uint32_t i = 0; i < N; i++) {
 
                 const uint32_t missing_and_ploidy = missing_and_ploidy_info[i];
 
                 if (missing_and_ploidy == 2) {
                     uintptr_t numer_aa;
                     uintptr_t numer_ab;

                     Bgen13GetTwoVals(probs_start, B, probs_offset, &numer_aa, &numer_ab);
                     probs_start += (probs_offset * 2);

                     double p11 = numer_aa / double(1.0 * (numer_mask));
                     double p10 = numer_ab / double(1.0 * (numer_mask));
                     double dosage = 2 - (p11 + p10);

                     probs(i, 0) = p11;
                     probs(i, 1) = p10;
                     dosVec[i] = dosage;
                     gmean += dosage;

                 }
                 else if (missing_and_ploidy == 130) {
                     probs(i, 0) = NA_REAL;
                     probs(i, 1) = NA_REAL;
                     dosVec[i] = NA_REAL;
                     nmiss+=1;
                 }
                 else {
                     Rcpp::stop("Variants with ploidy != 2 is currently not supported.");
                 }
            }

       }

       if (vr[0] == (nvars)){
           Rcout << "End of BGEN file has been reached. Close the file with close_bgen().\n";
       }

       double AF = gmean / (Nsamples-nmiss) / 2.0;

       delete[] snpID;
       delete[] rsID;
       delete[] chrStr;
       delete[] allele1;
       delete[] allele2;
       libdeflate_free_decompressor(decompressor);
       return(Rcpp::List::create(Named("SNPID") = snpID_2,
                                 Named("RSID") = rsID_2,
                                 Named("Chromosome") = chrStr_2,
                                 Named("Position") = physpos,
                                 Named("Alleles") = LKnum,
                                 Named("Allele1")  = allele1_2,
                                 Named("Allele2")  = allele2_2,
                                 Named("AF") = AF,
                                 Named("Probabilities") = probs,
                                 Named("Missing") = nmiss,
                                 Named("Dosages") = dosVec));
}





Rcpp::List query_bgen11(SEXP bgenR_in, SEXP seek_in){
  
  Rcpp::List bgenR_list(bgenR_in);
  FILE* bStream = (FILE*)R_ExternalPtrAddr(bgenR_list["fin"]);
  if (bStream == NULL) {
    Rcpp::stop("BGEN file is not open. Please reopen the BGEN file with open_bgen().");
  }
  uint nvars = Rcpp::as<uint>(bgenR_list["M"]);
  vector<uint>* counter = (vector<uint>*)R_ExternalPtrAddr(bgenR_list["fcounter"]);
  vector<unsigned long long int>* bytes= (vector<unsigned long long int>*)R_ExternalPtrAddr(bgenR_list["fbytes"]);
  vector<uint>&vr  = *counter;
  vector<unsigned long long int>&br  = *bytes;
  vr[0]++;
  if (vr[0] >= (nvars + 1)){
    Rcpp::stop("End of BGEN file has already been reached. Please close the file with close_bgen().");
  }


  uint Nsamples    = Rcpp::as<uint>(bgenR_list["N"]);
  uint Compression = Rcpp::as<uint>(bgenR_list["Compression"]);
  NumericMatrix probs(Nsamples, 3);
  NumericVector dosVec(Nsamples);
  int ret;
  
  uint seek = Rcpp::as<uint>(seek_in);
  if (seek != 0) {
    if (bytes == NULL) {
      Rcpp::stop("\nERROR: Reopen the BGEN file with `bytes=TRUE` in open_bgen() to use the `seek` argument in this function.");
    }
    fseek(bStream, br[seek-1], SEEK_SET);
  }
  
  uLongf destLen1 = 6 * Nsamples;
  std::vector<uchar> zBuf11;
  std::vector<uint16_t> shortBuf11;
  if (Compression == 0) {
    zBuf11.resize(destLen1);
  } else {
    shortBuf11.resize(destLen1);
  }

  struct libdeflate_decompressor* decompressor = libdeflate_alloc_decompressor();
  
  uint Nrow2; 
  ret = fread(&Nrow2, 4, 1, bStream);
  if (Nrow2 != Nsamples) {
    Rcpp::stop("Number of samples with genotype probabilities does not match the number of sample in BGEN header block.");
  }
  
  uint16_t LS; 
  ret = fread(&LS, 2, 1, bStream);
  char* snpID = new char[LS + 1];
  ret = fread(snpID, 1, LS, bStream); 
  snpID[LS] = '\0';
  String snpID_2 = string(snpID);
  
  uint16_t LR; 
  ret = fread(&LR, 2, 1, bStream);
  char* rsID = new char[LR + 1];
  ret = fread(rsID, 1, LR, bStream); 
  rsID[LR] = '\0';
  String rsID_2 = string(rsID);
  
  uint16_t LC; 
  ret = fread(&LC, 2, 1, bStream);
  char* chrStr = new char[LC + 1];
  ret = fread(chrStr, 1, LC, bStream); 
  chrStr[LC] = '\0';
  String chrStr_2 = string(chrStr);
  
  uint32_t physpos; 
  ret = fread(&physpos, 4, 1, bStream);
  
  uint32_t LA; 
  ret = fread(&LA, 4, 1, bStream);
  char* allele1 = new char[LA + 1];
  ret = fread(allele1, 1, LA, bStream); 
  allele1[LA] = '\0';
  String allele1_2 = string(allele1);
  
  uint32_t LB; 
  ret = fread(&LB, 4, 1, bStream);
  char* allele2 = new char[LB + 1];
  ret = fread(allele2, 1, LB, bStream); 
  allele2[LB] = '\0';
  String allele2_2 = string(allele2);
  
  uint16_t* probs_start;
  if (Compression == 1) {
    uint cLen; 
    ret = fread(&cLen, 4, 1, bStream);
    zBuf11.resize(cLen);
    ret = fread(&zBuf11[0], 1, cLen, bStream);
    
    if (libdeflate_zlib_decompress(decompressor, &zBuf11[0], cLen, &shortBuf11[0], destLen1, NULL) != LIBDEFLATE_SUCCESS) {
      Rcpp::stop("Decompressing " + string(rsID) + "genotype block failed with libdeflate.");
    }
    
    probs_start = &shortBuf11[0];
    
  }
  else {
    ret = fread(&zBuf11[0], 1, destLen1, bStream);
    probs_start = reinterpret_cast<uint16_t*>(&zBuf11[0]);
  }
  (void)ret;
  
  
  const double scale = 1.0 / 32768;
  double gmean = 0.0;
  double nmiss = 0.0;
  for (uint i = 0; i < Nsamples; i++) {
    double p11 = probs_start[3 * i] * scale;
    double p10 = probs_start[3 * i + 1] * scale;
    double p00 = probs_start[3 * i + 2] * scale;
    
    if (p11 == 0.0 && p10 == 0.0 && p00 == 0.0) {
      probs(i, 0) = NA_REAL;
      probs(i, 1) = NA_REAL;
      probs(i, 2) = NA_REAL;
      dosVec[i]   = NA_REAL;
      nmiss+=1.0;
      
    } else {
      double pTot = p11 + p10 + p00;
      double dosage = (2 * p00 + p10) / pTot;
      
      
      probs(i, 0) = p11;
      probs(i, 1) = p10;
      probs(i, 2) = p00;
      
      dosVec[i] = dosage;
      gmean += dosage;
      
    }
    
  }
  
  
  if (vr[0] == (nvars)){
    Rcout << "End of BGEN file has been reached. Close the file with close_bgen().\n";
  }
  
  double AF = gmean / (Nsamples-nmiss) / 2.0;
  
  delete[] snpID;
  delete[] rsID;
  delete[] chrStr;
  delete[] allele1;
  delete[] allele2;
  libdeflate_free_decompressor(decompressor);
  return(Rcpp::List::create(Named("SNPID") = snpID_2,
                            Named("RSID") = rsID_2,
                            Named("Chromosome") = chrStr_2,
                            Named("Position") = physpos,
                            Named("Allele1")  = allele1_2,
                            Named("Allele2")  = allele2_2,
                            Named("AF") = AF,
                            Named("Probabilities") = probs,
                            Named("Missing") = nmiss,
                            Named("Dosages") = dosVec));
}  

  
  
  
Rcpp::DataFrame get_vblock(SEXP bgenR_in){

    if(!bgenIsOpen) {
      Rcpp::stop("BGEN file is not open. Please reopen the BGEN file with open_bgen().");
    }
    
    Rcpp::List bgenR_list(bgenR_in);
    FILE* bStream = (FILE*)R_ExternalPtrAddr(bgenR_list["fin"]);
    
    uint offset = Rcpp::as<uint>(bgenR_list["offset"]);
    uint Layout = Rcpp::as<uint>(bgenR_list["Layout"]);
    uint nvars  = Rcpp::as<uint>(bgenR_list["M"]);
    uint Compression = Rcpp::as<uint>(bgenR_list["Compression"]);

    Rcpp::CharacterVector  vecSNPID(nvars);
    Rcpp::CharacterVector  vecRSID(nvars);
    Rcpp::CharacterVector  vecCHR(nvars);
    Rcpp::NumericVector    vecPOS(nvars);
    Rcpp::NumericVector    vecLK(nvars);
    Rcpp::CharacterVector  vecA1(nvars);
    Rcpp::CharacterVector  vecA2(nvars);
    Rcpp::NumericVector    vecBYTE(nvars);

    uint maxLA = 65536;
    char* snpID   = new char[maxLA + 1];
    char* rsID    = new char[maxLA + 1];
    char* chrStr  = new char[maxLA + 1];
    char* allele1 = new char[maxLA + 1];
    char* allele2 = new char[maxLA + 1];

    fseek(bStream, offset + 4, SEEK_SET);
    
    for (uint m = 0; m < nvars; m++) {
         int ret;
         long long unsigned int byte = ftell(bStream);
         if (Layout == 1) {
             uint Nrow;
             ret = fread(&Nrow, 4, 1, bStream);
         }

         uint16_t LS;
         ret = fread(&LS, 2, 1, bStream);
         ret = fread(snpID, 1, LS, bStream);
         snpID[LS] = '\0';
         Rcpp::String snpID_2 = string(snpID);

         uint16_t LR;
         ret = fread(&LR, 2, 1, bStream);
         ret = fread(rsID, 1, LR, bStream);
         rsID[LR] = '\0';
         Rcpp::String rsID_2 = string(rsID);

         uint16_t LC;
         ret = fread(&LC, 2, 1, bStream);
         ret = fread(chrStr, 1, LC, bStream);
         chrStr[LC] = '\0';
         Rcpp::String chrStr_2 = string(chrStr);

         uint32_t physpos;
         ret = fread(&physpos, 4, 1, bStream);

         uint16_t LKnum;
         if (Layout == 2) {
             ret = fread(&LKnum, 2, 1, bStream);

             if (LKnum != 2U) {
                 for (uint16_t a = 0; a < LKnum; a++) {
                      uint32_t LA;
                      ret = fread(&LA, 4, 1, bStream);
                      ret = fread(allele1, 1, LA, bStream);
                 }

                 if (Compression > 0) {
                     uint zLen;
                     ret = fread(&zLen, 4, 1, bStream);
                     fseek(bStream, 4 + zLen - 4, SEEK_CUR);

                 }
                 else {
                     uint zLen;
                     ret = fread(&zLen, 4, 1, bStream);
                     fseek(bStream, zLen, SEEK_CUR);
                 }

                 vecSNPID[m] = snpID_2;
                 vecRSID[m]  = rsID_2;
                 vecCHR[m]   = chrStr_2;
                 vecPOS[m]   = physpos;
                 vecLK[m]    = LKnum;
                 vecA1[m]    = NA_STRING;
                 vecA2[m]    = NA_STRING;
                 vecBYTE[m]  = byte;

                 continue;
             }
         }
         else {
             LKnum = 2U;
         }

         uint32_t LA;
         ret = fread(&LA, 4, 1, bStream);
         ret = fread(allele1, 1, LA, bStream);
         allele1[LA] = '\0';
         Rcpp::String allele1_2 = string(allele1);

         uint32_t LB;
         ret = fread(&LB, 4, 1, bStream);
         ret = fread(allele2, 1, LB, bStream);
         allele2[LB] = '\0';
         Rcpp::String allele2_2 = string(allele2);

         if (Layout == 2) {
             if (Compression > 0) {
                 uint zLen;
                 ret = fread(&zLen, 4, 1, bStream);
                 fseek(bStream, 4 + zLen - 4, SEEK_CUR);

             }
             else {
                 uint zLen;
                 ret = fread(&zLen, 4, 1, bStream);
                 fseek(bStream, zLen, SEEK_CUR);
             }
         }
         else {
             if (Compression == 1) {
                 uint zLen;
                 ret = fread(&zLen, 4, 1, bStream);
                 fseek(bStream, zLen, SEEK_CUR);
             }
             else {
                 fseek(bStream, 6 * nvars, SEEK_CUR);
             }
         }
         (void)ret;
         vecSNPID[m] = snpID_2;
         vecRSID[m]  = rsID_2;
         vecCHR[m]   = chrStr_2;
         vecPOS[m]   = physpos;
         vecLK[m]    = LKnum;
         vecA1[m]    = allele1_2;
         vecA2[m]    = allele2_2;
         vecBYTE[m] = byte;
    }
    
    fseek(bStream, offset + 4, SEEK_SET);

    delete[] snpID;
    delete[] rsID;
    delete[] chrStr;
    delete[] allele1;
    delete[] allele2;
    return(Rcpp::DataFrame::create(Named("SNPID")   = vecSNPID,
                                   Named("RSID")    = vecRSID,
                                   Named("CHR")     = vecCHR,
                                   Named("POS")     = vecPOS,
                                   Named("ALLELES") = vecLK,
                                   Named("A1")      = vecA1,
                                   Named("A2")      = vecA2));
}






vector<long long unsigned int>* get_bytes(FILE* & bStream_in, uint offset_in, uint Layout_in, uint M_in, uint Compression_in){
  
    int maxLA = 65536;
    char* snpID   = new char[maxLA + 1];
    char* rsID    = new char[maxLA + 1];
    char* chrStr  = new char[maxLA + 1];
    char* allele1 = new char[maxLA + 1];
    char* allele2 = new char[maxLA + 1];
    std::vector<long long unsigned int>* bytes = new vector<long long unsigned int>(M_in);
    vector<long long unsigned int> &vr = *bytes;
    
    fseek(bStream_in, offset_in + 4, SEEK_SET);

    for (uint m = 0; m < M_in; m++) {
      int ret;
      vr[m] = ftell(bStream_in);
      
      if (Layout_in == 1) {
        uint Nrow;
        ret = fread(&Nrow, 4, 1, bStream_in);
      }
      
      uint16_t LS;
      ret = fread(&LS, 2, 1, bStream_in);
      ret = fread(snpID, 1, LS, bStream_in);
      
      uint16_t LR;
      ret = fread(&LR, 2, 1, bStream_in);
      ret = fread(rsID, 1, LR, bStream_in);

      uint16_t LC;
      ret = fread(&LC, 2, 1, bStream_in);
      ret = fread(chrStr, 1, LC, bStream_in);
      
      uint32_t physpos;
      ret = fread(&physpos, 4, 1, bStream_in);
      
      uint16_t LKnum;
      if (Layout_in == 2) {
          ret = fread(&LKnum, 2, 1, bStream_in);
        
          for (uint16_t a = 0; a < LKnum; a++) {
               uint32_t LA;
               ret = fread(&LA, 4, 1, bStream_in);
               ret = fread(allele1, 1, LA, bStream_in);
          }
          
        
      } else {
          uint32_t LA;
          ret = fread(&LA, 4, 1, bStream_in);
          ret = fread(allele1, 1, LA, bStream_in);
          
          uint32_t LB;
          ret = fread(&LB, 4, 1, bStream_in);
          ret = fread(allele2, 1, LB, bStream_in);
      }
      
      if (Layout_in == 2) {
        if (Compression_in > 0) {
          uint zLen;
          ret = fread(&zLen, 4, 1, bStream_in);
          fseek(bStream_in, 4 + zLen - 4, SEEK_CUR);
          
        }
        else {
          uint zLen;
          ret = fread(&zLen, 4, 1, bStream_in);
          fseek(bStream_in, zLen, SEEK_CUR);
        }
      }
      else {
        if (Compression_in == 1) {
          uint zLen;
          ret = fread(&zLen, 4, 1, bStream_in);
          fseek(bStream_in, zLen, SEEK_CUR);
        }
        else {
          fseek(bStream_in, 6 * M_in, SEEK_CUR);
        }
      }
      (void)ret;
    }
    
    fseek(bStream_in, offset_in + 4, SEEK_SET);
    delete[] snpID;
    delete[] rsID;
    delete[] chrStr;
    delete[] allele1;
    delete[] allele2;
    
    
    return(bytes);
}
