#include "read_bgen.h"
#include "libdeflate-1.6/libdeflate.h"
#include "zstd-1.4.5/lib/zstd.h"
#include <Rcpp.h>
#include <stdio.h>
#include <zlib.h>




typedef unsigned int   uint;
typedef unsigned char  uchar;
typedef unsigned short ushort;


using namespace std;
using namespace Rcpp;


FILE* bStream;
BGEN bgen;
struct libdeflate_decompressor* decompressor = libdeflate_alloc_decompressor();

uLongf destLen1;
uchar* zBuf11;
uint16_t* shortBuf11;
std::vector<uchar> zBuf12;
std::vector<uchar> shortBuf12;

uint maxLA = 65536;
char* snpID   = new char[maxLA + 1];
char* rsID    = new char[maxLA + 1];
char* chrStr  = new char[maxLA + 1];
char* allele1 = new char[maxLA + 1];
char* allele0 = new char[maxLA + 1];


Rcpp::CharacterVector read_bgenSampleID();
Rcpp::List query_bgen11();
Rcpp::List query_bgen13();


Rcpp::String close_bgen(){
  
  if(bStream == NULL) {
     Rcpp::stop("BGEN file is already closed.");
  }
  
  delete[] snpID;
  delete[] rsID;
  delete[] chrStr;
  delete[] allele1;
  delete[] allele0;

  if ( bgen.Layout == 1 ) {
      free(zBuf11);
      free(shortBuf11);
  }

  fseek(bStream, 0, SEEK_END);
  fclose(bStream);
  bStream = NULL;
  
  return("BGEN file closed successfully.");
}



Rcpp::List open_bgen(SEXP bgenfile_in){
  
  string bgenfile = Rcpp::as<string>(bgenfile_in);
  
  
  bStream = fopen(bgenfile.c_str(), "rb");
  if (!bStream) { 
      Rcpp::stop("ERROR: Cannot not open BGEN file: " + bgenfile + "\n"); 
  }
  
  
  // Header Block
  fread(&bgen.offset, 4, 1, bStream);
  uint L_H; fread(&L_H, 4, 1, bStream);
  fread(&bgen.Mbgen, 4, 1, bStream);
  fread(&bgen.Nbgen, 4, 1, bStream);
  char magic[5]; fread(magic, 1, 4, bStream);
  magic[4] = '\0';
  
  fseek(bStream, L_H - 20, SEEK_CUR);
  uint flags; fread(&flags, 4, 1, bStream);
  
  
  // Header Block Flags
  bgen.Compression = flags & 3;
  bgen.Layout = (flags >> 2) & 0xf;
  bgen.SampleIdentifiers = flags >> 31;
  
  
  
  if (bgen.SampleIdentifiers == 1){
    bgen.sampleID = read_bgenSampleID();
  }
  
  if (bgen.Layout == 1) {
      destLen1 = 6 * bgen.Nbgen;
      zBuf11 = (unsigned char*)malloc(6 * bgen.Nbgen);
      shortBuf11 = (uint16_t*)malloc(6 * bgen.Nbgen);
  }
  
  fseek(bStream, bgen.offset + 4, SEEK_SET);
  bgen.Counter = 0;
  return(Rcpp::List::create(Named("offset") = bgen.offset,
                            Named("M") = bgen.Mbgen,
                            Named("N") = bgen.Nbgen,
                            Named("Compression") = bgen.Compression,
                            Named("Layout") = bgen.Layout,
                            Named("SampleIdentifier") = bgen.SampleIdentifiers,
                            Named("SampleID") = bgen.sampleID));
  
}



Rcpp::CharacterVector read_bgenSampleID(){
  
  Rcpp::CharacterVector sampleID(bgen.Nbgen);
  
  uint LS1;  fread(&LS1,  4, 1, bStream);
  uint Nrow; fread(&Nrow, 4, 1, bStream);
  
  if (Nrow != bgen.Nbgen) {
    Rcpp::Rcout << "Number of sample identifiers (" << Nrow << ") does not match the number of BGEN samples (" << bgen.Nbgen << ")\n";
    Rcpp::stop("");
  }
  
  
  char* samID = new char[maxLA + 1];
  for (uint n = 0; n < Nrow; n++) {
    ushort LSID; fread(&LSID, 2, 1, bStream);
    fread(samID, 1, LSID, bStream);
    samID[LSID] = '\0';
    sampleID[n] = samID;
  }
  
  return sampleID;
  
}




// These two functions below are being redistributed from plink2.0 
uintptr_t Bgen13GetOneVal(const unsigned char* prob_start, uint64_t prob_offset, uint32_t bit_precision, uintptr_t numer_mask) {
  const uint64_t bit_offset = prob_offset * bit_precision;
  uint64_t relevant_bits;
  memcpy(&relevant_bits, &(prob_start[bit_offset / CHAR_BIT]), sizeof(int64_t));
  return (relevant_bits >> (bit_offset % CHAR_BIT)) & numer_mask;
}


void Bgen13GetTwoVals(const unsigned char* prob_start, uint64_t prob_offset, uint32_t bit_precision, uintptr_t numer_mask, uintptr_t* first_val_ptr, uintptr_t* second_val_ptr) {
  const uint64_t bit_offset = bit_precision / 8.0;
  
  switch(bit_precision) {
  case 8:
    *first_val_ptr  = prob_start[0];
    prob_start += bit_offset;
    *second_val_ptr = prob_start[0];
    prob_start += bit_offset;
    break;
  case 16:
    *first_val_ptr  = bufAt[0]|(bufAt[1]<<8);
    prob_start += bit_offset;
    *second_val_ptr = bufAt[0]|(bufAt[1]<<8);
    prob_start += bit_offset;
    break;
  case 24:
    *first_val_ptr  = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16);
    prob_start += bit_offset;
    *second_val_ptr = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16);
    prob_start += bit_offset;
    break;
  case 32:
    *first_val_ptr  = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24);
    prob_start += bit_offset;
    *second_val_ptr = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24);
    prob_start += bit_offset;
    
  }
  
}




Rcpp::List query_bgen(){
  
  
  if(bgen.Layout == 2){
    return(query_bgen13());
  }
  else if(bgen.Layout == 1){
    return(query_bgen11());
  }

  return(NA_REAL);
} 






Rcpp::List query_bgen13(){
  
  bgen.Counter++;
  if(bgen.Counter >= (bgen.Mbgen + 1)){
     Rcpp::stop("End of BGEN file has already been reached. Please close the file with close_bgen().");
  }
  
  uint Nsamples = bgen.Nbgen;
  uint Compression = bgen.Compression;
  NumericMatrix probs(Nsamples, 2);
  NumericVector dosVec(Nsamples);

  ushort LS; fread(&LS, 2, 1, bStream);
  fread(snpID, 1, LS, bStream); snpID[LS] = '\0';
  
  
  ushort LR; fread(&LR, 2, 1, bStream);
  fread(rsID, 1, LR, bStream); rsID[LR] = '\0';
  
  ushort LC; fread(&LC, 2, 1, bStream);
  fread(chrStr, 1, LC, bStream); chrStr[LC] = '\0';
  uint physpos; fread(&physpos, 4, 1, bStream);

  ushort LKnum; fread(&LKnum,   2, 1, bStream);
  if ( LKnum != 2U ){
    Rcpp::stop("\nERROR: " + string(rsID) + " does not contain 2 alleles.\n\n");
  }
  
  uint32_t LA; fread(&LA, 4, 1, bStream);
  fread(allele1, 1, LA, bStream); allele1[LA] = '\0';
  
  uint32_t LB; fread(&LB, 4, 1, bStream);
  fread(allele0, 1, LB, bStream); allele0[LB] = '\0';
  
  
  uchar* bufAt;
  uint cLen; fread(&cLen, 4, 1, bStream);
  if (Compression == 1) {
    zBuf12.resize(cLen - 4);
    uint dLen; fread(&dLen, 4, 1, bStream);
    fread(&zBuf12[0], 1, cLen - 4, bStream);
    shortBuf12.resize(dLen);
    
    uLongf destLen = dLen;
    if (libdeflate_zlib_decompress(decompressor, &zBuf12[0], cLen - 4, &shortBuf12[0], destLen, NULL) != LIBDEFLATE_SUCCESS) {
        Rcpp::stop("ERROR: Decompressing " + string(rsID) + "genotype block failed.\n\n");
    }
    bufAt = &shortBuf12[0];
  }
  else if (Compression == 2) {
    zBuf12.resize(cLen - 4);
    uint dLen; fread(&dLen, 4, 1, bStream);
    fread(&zBuf12[0], 1, cLen - 4, bStream);
    shortBuf12.resize(dLen);
    
    uLongf destLen = dLen;
    size_t ret = ZSTD_decompress(&shortBuf12[0], destLen, &zBuf12[0], cLen - 4);
    if (ret > destLen) {
      if (ZSTD_isError(ret)) {
        cout << "ZSTD ERROR: " << ZSTD_getErrorName(ret);
      }
    }
    bufAt = &shortBuf12[0];
  }
  else {
    zBuf12.resize(cLen);
    fread(&zBuf12[0], 1, cLen, bStream);
    bufAt = &zBuf12[0];
  }
  
  
  uint32_t N; memcpy(&N, bufAt, sizeof(int32_t));
  uint16_t K; memcpy(&K, &(bufAt[4]), sizeof(int16_t));
  const uint32_t min_ploidy = bufAt[6];
  if (min_ploidy > 2) { Rcpp::stop("ERROR: Variants with ploidy > 2 is currently not supported."); }
  const uint32_t max_ploidy = bufAt[7];
  if (max_ploidy > 2) { Rcpp::stop("ERROR: Variants with ploidy > 2 is currently not supported."); }
  const unsigned char* missing_and_ploidy_info = &(bufAt[8]);
  const unsigned char* probs_start = &(bufAt[10 + N]);
  const uint32_t is_phased = probs_start[-2];

  const uint32_t B = probs_start[-1];
  const uintptr_t numer_mask = (1U << B) - 1;
  
  
  uintptr_t prob_offset = 0;
  double gmean = 0.0;
  
  if(!is_phased){
    
     for (uint32_t i = 0; i < N; i++) {
          const uint32_t missing_and_ploidy = missing_and_ploidy_info[i];

          
          if(missing_and_ploidy == 2){
             uintptr_t numer_aa;
             uintptr_t numer_ab;

             Bgen13GetTwoVals(probs_start, prob_offset, B, numer_mask, &numer_aa, &numer_ab);
             prob_offset += 2;
            
             double p11 = numer_aa / double(1.0 * (numer_mask));
             double p10 = numer_ab / double(1.0 * (numer_mask));
             double dosage = 2 * (1 - p11 - p10) + p10;
            
             probs(i, 0) = p11;
             probs(i, 1) = p10;
             dosVec[i] = dosage;
             gmean += dosage;
          
          }
          else if (missing_and_ploidy == 1) {

             const uintptr_t numer_a = Bgen13GetOneVal(probs_start, prob_offset, B, numer_mask);
             prob_offset++;

             double p11 = numer_a / double(1.0 * (numer_mask));
             double dosage = 1 - p11;

             probs(i, 0) = p11;
             probs(i, 1) = NA_REAL;
             dosVec[i]   = dosage;
             gmean += dosage;

          }
          else {
            
             probs(i, 0) = NA_REAL;
             probs(i, 1) = NA_REAL;
             dosVec[i]   = NA_REAL;
          }
     }
    
  } else if(is_phased){
       for (uint32_t i = 0; i < N; i++) {
         
            const uint32_t missing_and_ploidy = missing_and_ploidy_info[i];
         
            if(missing_and_ploidy == 2){
               uintptr_t numer_aa;
               uintptr_t numer_ab;

               Bgen13GetTwoVals(probs_start, prob_offset, B, numer_mask, &numer_aa, &numer_ab);
               prob_offset += 2;
           
               double p11 = numer_aa / double(1.0 * (numer_mask));
               double p10 = numer_ab / double(1.0 * (numer_mask));
               double dosage = 2 - (p11 + p10);
           
               probs(i, 0) = p11;
               probs(i, 1) = p10;
               dosVec[i] = dosage;
               gmean += dosage;
           
            } else if (missing_and_ploidy == 1){
                const uintptr_t numer_a = Bgen13GetOneVal(probs_start, prob_offset, B, numer_mask);
               prob_offset++;
                  
               double p11 = numer_a / double(1.0 * (numer_mask));
               double dosage = 1 - p11;
                 
               probs(i, 0) = p11;
               probs(i, 1) = NA_REAL;
               dosVec[i]   = dosage;
               gmean += dosage;
               
            } else {
               probs(i, 0) = NA_REAL;
               probs(i, 1) = NA_REAL;
               dosVec[i]   = NA_REAL;
            }
       }
      
  }
  
  if(bgen.Counter == (bgen.Mbgen)){
    Rcpp::Rcout << "End of BGEN file has been reached. Please close the file with close_bgen()." << std::endl;
  }
  
  double AF = gmean / Nsamples / 2.0;
  return(Rcpp::List::create(Named("SNPID") = string(snpID),
                            Named("RSID") = string(rsID),
                            Named("Chromosome") = string(chrStr),
                            Named("Position") = physpos,
                            Named("Allele1")  = string(allele0),
                            Named("Allele2")  = string(allele1),
                            Named("AF") = AF,
                            Named("Probabilities") = probs,
                            Named("Dosages") = dosVec));
  
}









Rcpp::List query_bgen11(){
  
    bgen.Counter++;
    if(bgen.Counter >= (bgen.Mbgen + 1)){
      Rcpp::stop("End of BGEN file has already been reached. Please close the file with close_bgen().");
    }

    
    uint Nsamples = bgen.Nbgen;
    
    uint Compression = bgen.Compression;
    NumericMatrix probs(Nsamples, 3);
    NumericVector dosVec(Nsamples);
    
    
    uint Nrow2; fread(&Nrow2, 4, 1, bStream); 
    if (Nrow2 != Nsamples) {
        Rcpp::stop("\nERROR: Number of samples with genotype probabilities for variant " + string(rsID) + " does not match the number of sample in BGEN header block.\n\n");
    }
    
    ushort LS; fread(&LS, 2, 1, bStream);
    fread(snpID, 1, LS, bStream); snpID[LS] = '\0';
    
    ushort LR; fread(&LR, 2, 1, bStream);
    fread(rsID, 1, LR, bStream); rsID[LR] = '\0';
    
    ushort LC; fread(&LC, 2, 1, bStream);
    fread(chrStr, 1, LC, bStream); chrStr[LC] = '\0';
    uint physpos; fread(&physpos, 4, 1, bStream);
    
    uint32_t LA; fread(&LA, 4, 1, bStream);
    fread(allele1, 1, LA, bStream); allele1[LA] = '\0';
    
    uint32_t LB; fread(&LB, 4, 1, bStream);
    fread(allele0, 1, LB, bStream); allele0[LB] = '\0';
    
    
  
    if (Compression == 1) {
        uint cLen; fread(&cLen, 4, 1, bStream);
        fread(zBuf11, 1, cLen, bStream);
      
        if (libdeflate_zlib_decompress(decompressor, &zBuf11[0], cLen, &shortBuf11[0], destLen1, NULL) != LIBDEFLATE_SUCCESS) {
            Rcpp::stop("ERROR: Decompressing " + string(rsID) + "genotype block failed.\n\n");
        }
    
    }
    else {
        fread(zBuf11, 1, destLen1, bStream);
        shortBuf11 = reinterpret_cast<uint16_t*>(zBuf11);
    }
  
  

    const double scale = 1.0 / 32768;
    double gmean = 0.0;
    for (uint i = 0; i < Nsamples; i++) {
         double p11 = shortBuf11[3 * i] * scale;
         double p10 = shortBuf11[3 * i + 1] * scale;
         double p00 = shortBuf11[3 * i + 2] * scale;
         double pTot = p11 + p10 + p00;
         double dosage = (2 * p00 + p10) / pTot;
        

         probs(i, 0) = p11;
         probs(i, 1) = p10;
         probs(i, 2) = p00;
         
         dosVec[i] = dosage;
         gmean += dosage;
            
    }
    
  
    if(bgen.Counter == (bgen.Mbgen)){
       Rcpp::Rcout << "End of BGEN file has been reached. Please close the file with close_bgen()." << std::endl;
    }
  
    double AF = gmean / Nsamples / 2.0;
    return(Rcpp::List::create(Named("SNPID") = string(snpID),
                              Named("RSID") = string(rsID),
                              Named("Chromosome") = string(chrStr),
                              Named("Position") = physpos,
                              Named("Allele1")  = string(allele0),
                              Named("Allele2")  = string(allele1),
                              Named("AF") = AF,
                              Named("Probabilities") = probs,
                              Named("Dosages") = dosVec));
}  


  
  
  
Rcpp::DataFrame get_variantBlock(){
  


  if(bStream == NULL) {
     Rcpp::stop("\nERROR: BGEN file is not open. Please reopen the BGEN file with open_bgen().");
  }
  
  uint Layout = bgen.Layout;
  uint Nsamples = bgen.Nbgen;
  uint Compression = bgen.Compression;
  
  Rcpp::CharacterVector  vecSNPID(Nsamples);
  Rcpp::CharacterVector  vecRSID(Nsamples);
  Rcpp::CharacterVector  vecCHR(Nsamples);
  Rcpp::NumericVector    vecPOS(Nsamples);
  Rcpp::CharacterVector  vecA1(Nsamples);
  Rcpp::CharacterVector  vecA2(Nsamples);
  
  fseek(bStream, bgen.offset + 4, SEEK_SET);
  
  
  uint maxLA = 65536;
  uint maxLB = 65536;
  char* snpID   = new char[maxLA + 1];
  char* rsID    = new char[maxLA + 1];
  char* chrStr  = new char[maxLA + 1];
  char* allele1 = new char[maxLA + 1];
  char* allele0 = new char[maxLB + 1];
  
  
  
  for (uint m = 0; m < bgen.Mbgen; m++) {
       maxLA = 65536;
       maxLB = 65536;
       
       
       if (Layout == 1) {
           uint Nrow; fread(&Nrow, 4, 1, bStream);
       }
 
       ushort LS; fread(&LS, 2, 1, bStream);
       fread(snpID, 1, LS, bStream); snpID[LS] = '\0'; 
       
       ushort LR; fread(&LR, 2, 1, bStream);
       
       fread(rsID, 1, LR, bStream); rsID[LR] = '\0'; 
       ushort LC; fread(&LC, 2, 1, bStream);
       fread(chrStr, 1, LC, bStream); chrStr[LC] = '\0';
       uint physpos; fread(&physpos, 4, 1, bStream);
       
       if (Layout == 2) {
           ushort LKnum; fread(&LKnum, 2, 1, bStream);
       }
    
       uint32_t LA; fread(&LA, 4, 1, bStream);
       fread(allele1, 1, LA, bStream); allele1[LA] = '\0';

       uint32_t LB; fread(&LB, 4, 1, bStream);
       fread(allele0, 1, LB, bStream); allele0[LB] = '\0';
    
       if (Layout == 2) {
           if (Compression > 0) {
               uint zLen;  fread(&zLen, 4, 1, bStream);
               fseek(bStream, 4 + zLen - 4, SEEK_CUR);
          
           }
           else {
               uint zLen;  fread(&zLen, 4, 1, bStream);
               fseek(bStream, zLen, SEEK_CUR);
           }
       }
       else {
           if (Compression == 1) {
               uint zLen;  fread(&zLen, 4, 1, bStream);
               fseek(bStream, zLen, SEEK_CUR);
           }
           else {
               fseek(bStream, 6 * Nsamples, SEEK_CUR);
           }
       }
    
       vecSNPID[m] = snpID;
       vecRSID[m]  = rsID;
       vecCHR[m]   = chrStr;
       vecPOS[m]   = physpos;
       vecA1[m]    = allele1;
       vecA2[m]    = allele0;
       
  }
  
  fseek(bStream, bgen.offset + 4, SEEK_SET);
  
  
  return(Rcpp::DataFrame::create(Named("SNPID") = vecSNPID,
                                 Named("RSID")  = vecRSID,
                                 Named("CHR")   = vecCHR,
                                 Named("POS")   = vecPOS,
                                 Named("A1")    = vecA1,
                                 Named("A2")    = vecA2));  
}
