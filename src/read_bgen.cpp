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
  if (!fread(&bgen.offset, 4, 1, bStream)) {
      Rcpp::stop("ERROR: Cannot read BGEN header block (offset)");
  }

  uint L_H; 
  if (!fread(&L_H, 4, 1, bStream)) {
      Rcpp::stop("ERROR: Cannot read BGEN header block (LH)");
  }

  if (!fread(&bgen.Mbgen, 4, 1, bStream)) {
      Rcpp::stop("ERROR: Cannot read BGEN header block (M)");
  }

  if (!fread(&bgen.Nbgen, 4, 1, bStream)) {
      Rcpp::stop("ERROR: Cannot read BGEN header block (N)");
  }

  char magic[5]; 
  if (!fread(magic, 1, 4, bStream)) {
      Rcpp::stop("ERROR: Cannot read BGEN header block (magic numbers)");
  }
  magic[4] = '\0';
  if (!(magic[0] == 'b' && magic[1] == 'g' && magic[2] == 'e' && magic[3] == 'n')) {
      Rcpp::stop("ERROR: BGEN file's magic number does not match 'b', 'g', 'e', 'n'.");
  }

  fseek(bStream, L_H - 20, SEEK_CUR);

  uint flags; 
  if (!fread(&flags, 4, 1, bStream)) {
      Rcpp::stop("ERROR: Cannot read BGEN header block (flags)");
  }
  
  
  // Header Block Flags
  bgen.Compression = flags & 3;
  if (bgen.Compression != 0 && bgen.Compression != 1 && bgen.Compression != 2) {
      Rcpp::stop("ERROR: BGEN compression flag should be 0, 1, or 2.");
  }

  bgen.Layout = (flags >> 2) & 0xf;
  if (bgen.Layout != 1 && bgen.Layout != 2) {
      Rcpp::stop("ERROR: BGEN layout flag should be 1 or 2.");
  }

  bgen.SampleIdentifiers = flags >> 31;
  if (bgen.SampleIdentifiers != 0 && bgen.SampleIdentifiers != 1) {
      Rcpp::stop("ERROR: BGEN layout flag should be 0 or 1.");
  }
  
  
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
  
  uint LS1;  
  if (!fread(&LS1, 4, 1, bStream)) {
      Rcpp::stop("ERROR: Cannot read BGEN sample block (LS)");
  }
  uint Nrow; 
  if (!fread(&Nrow, 4, 1, bStream)) {
      Rcpp::stop("ERROR: Cannot read BGEN sample block (N)");
  }
  
  if (Nrow != bgen.Nbgen) {
      Rcpp::stop("Number of sample identifiers (" + std::to_string(Nrow) + ") does not match the number of BGEN samples (" + std::to_string(bgen.Nbgen) + ")");
  }
  
  
  char* samID = new char[maxLA + 1];
  int ret;
  for (uint n = 0; n < Nrow; n++) {
    ushort LSID; 
    ret = fread(&LSID, 2, 1, bStream);
    ret = fread(samID, 1, LSID, bStream);
    samID[LSID] = '\0';
    
    sampleID[n] = samID;
  }
  (void)ret;
  delete[] samID;
  return sampleID;
  
}




// These two functions below are being redistributed from plink2.0 
uintptr_t Bgen13GetOneVal(const unsigned char* prob_start, uint32_t bit_precision) {
  
  switch(bit_precision) {
  case 8:
    return(prob_start[0]);
  case 16:
    return(prob_start[0]|(prob_start[1]<<8));
    break;
  case 24:
    return(prob_start[0]|(prob_start[1]<<8)|(prob_start[2]<<16));
    break;
  case 32:
    return(prob_start[0]|(prob_start[1]<<8)|(prob_start[2]<<16)|(prob_start[3]<<24));
  default:
    break;
  }
  
  return 0;
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
  int ret;

  ushort LS; 
  ret = fread(&LS, 2, 1, bStream);
  ret = fread(snpID, 1, LS, bStream); 
  snpID[LS] = '\0';
  
  ushort LR; 
  ret = fread(&LR, 2, 1, bStream);
  ret = fread(rsID, 1, LR, bStream); 
  rsID[LR] = '\0';
  
  ushort LC; 
  ret = fread(&LC, 2, 1, bStream);
  ret = fread(chrStr, 1, LC, bStream); 
  chrStr[LC] = '\0';

  uint32_t physpos; 
  ret = fread(&physpos, 4, 1, bStream);

  ushort LKnum; 
  ret = fread(&LKnum,   2, 1, bStream);
  if ( LKnum != 2U ){
    Rcpp::stop("\nERROR: " + string(rsID) + " does not contain 2 alleles.");
  }
  
  uint32_t LA; 
  ret = fread(&LA, 4, 1, bStream);
  ret = fread(allele1, 1, LA, bStream); 
  allele1[LA] = '\0';
  
  uint32_t LB; 
  ret = fread(&LB, 4, 1, bStream);
  ret = fread(allele0, 1, LB, bStream); 
  allele0[LB] = '\0';
  
  
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
        Rcpp::stop("ERROR: Decompressing " + string(rsID) + "genotype block failed with libdeflate.");
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
          Rcout << "ZSTD ERROR: " << ZSTD_getErrorName(ret);
          Rcpp::stop("\n\n");
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
  if (min_ploidy > 2) { 
      Rcpp::stop("ERROR: Variants with ploidy > 2 is currently not supported."); 
  }
  const uint32_t max_ploidy = prob_start[7];
  if (max_ploidy > 2) { 
      Rcpp::stop("ERROR: Variants with ploidy > 2 is currently not supported."); 
  }

  const unsigned char* missing_and_ploidy_info = &(prob_start[8]);
  const unsigned char* probs_start = &(prob_start[10 + N]);
  const uint32_t is_phased = probs_start[-2];
  if (is_phased != 0 && is_phased != 1) {
      Rcpp::stop("ERROR: phased value must be 0 or 1.");
  }

  const uint32_t B = probs_start[-1];
  if (B != 8 && B != 16 && B != 24 && B != 32) {
      Rcpp::stop("ERROR: Bits to store probabilities must be 8, 16, 24, or 32.");
  }
  const uintptr_t numer_mask = (1U << B) - 1;
  const uintptr_t probs_offset = B / 8;
  
  
  double gmean = 0.0;
  
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
          else if (missing_and_ploidy == 1) {

             const uintptr_t numer_a = Bgen13GetOneVal(probs_start, B);
             probs_start += probs_offset;

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

               Bgen13GetTwoVals(probs_start, B, probs_offset, &numer_aa, &numer_ab);
               probs_start += (probs_offset * 2);
           
               double p11 = numer_aa / double(1.0 * (numer_mask));
               double p10 = numer_ab / double(1.0 * (numer_mask));
               double dosage = 2 - (p11 + p10);
           
               probs(i, 0) = p11;
               probs(i, 1) = p10;
               dosVec[i] = dosage;
               gmean += dosage;
           
            } else if (missing_and_ploidy == 1){
               const uintptr_t numer_a = Bgen13GetOneVal(probs_start, B);
               probs_start += probs_offset;
                  
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
  
  if (bgen.Counter == (bgen.Mbgen)){
      Rcpp::stop("End of BGEN file has been reached. Please close the file with close_bgen().");
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
    int ret;
    
    uint Nrow2; 
    ret = fread(&Nrow2, 4, 1, bStream);
    if (Nrow2 != Nsamples) {
        Rcpp::stop("\nERROR: Number of samples with genotype probabilities for variant " + string(rsID) + " does not match the number of sample in BGEN header block.\n\n");
    }
    
    ushort LS; 
    ret = fread(&LS, 2, 1, bStream);
    ret = fread(snpID, 1, LS, bStream); 
    snpID[LS] = '\0';
    
    ushort LR; 
    ret = fread(&LR, 2, 1, bStream);
    ret = fread(rsID, 1, LR, bStream); 
    rsID[LR] = '\0';
    
    ushort LC; 
    ret = fread(&LC, 2, 1, bStream);
    ret = fread(chrStr, 1, LC, bStream); 
    chrStr[LC] = '\0';

    uint32_t physpos; 
    ret = fread(&physpos, 4, 1, bStream);
    
    uint32_t LA; 
    ret = fread(&LA, 4, 1, bStream);
    ret = fread(allele1, 1, LA, bStream); 
    allele1[LA] = '\0';
    
    uint32_t LB; ret = fread(&LB, 4, 1, bStream);
    ret = fread(allele0, 1, LB, bStream); 
    allele0[LB] = '\0';
    
    
  
    if (Compression == 1) {
        uint cLen; 
        ret = fread(&cLen, 4, 1, bStream);
        ret = fread(zBuf11, 1, cLen, bStream);
      
        if (libdeflate_zlib_decompress(decompressor, &zBuf11[0], cLen, &shortBuf11[0], destLen1, NULL) != LIBDEFLATE_SUCCESS) {
            Rcpp::stop("ERROR: Decompressing " + string(rsID) + "genotype block failed with libdeflate.\n\n");
        }
    
    }
    else if (Compression == 2) {
        uint cLen;
        ret = fread(&cLen, 4, 1, bStream);
        ret = fread(zBuf11, 1, cLen, bStream);

        size_t ret = ZSTD_decompress(&shortBuf11[0], cLen, &zBuf11[0], destLen1);
        if (ret > destLen1) {
            if (ZSTD_isError(ret)) {
                Rcout << "ZSTD ERROR: " << ZSTD_getErrorName(ret);
                Rcpp::stop("\n\n");
            }
        }
    }
    else {
        ret = fread(zBuf11, 1, destLen1, bStream);
        shortBuf11 = reinterpret_cast<uint16_t*>(zBuf11);
    }
    (void)ret;
  

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
        Rcpp::stop("End of BGEN file has been reached. Please close the file with close_bgen().");
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
       int ret;
       
       if (Layout == 1) {
           uint Nrow;
           ret = fread(&Nrow, 4, 1, bStream);
       }
 
       ushort LS; 
       ret = fread(&LS, 2, 1, bStream);
       ret = fread(snpID, 1, LS, bStream); 
       snpID[LS] = '\0';
       
       ushort LR; 
       ret = fread(&LR, 2, 1, bStream);
       ret = fread(rsID, 1, LR, bStream);
       rsID[LR] = '\0';

       ushort LC; 
       ret = fread(&LC, 2, 1, bStream);
       ret = fread(chrStr, 1, LC, bStream); 
       chrStr[LC] = '\0';

       uint32_t physpos; 
       ret = fread(&physpos, 4, 1, bStream);
       
       if (Layout == 2) {
           ushort LKnum; 
           ret = fread(&LKnum, 2, 1, bStream);
       }
    
       uint32_t LA; 
       ret = fread(&LA, 4, 1, bStream);
       ret = fread(allele1, 1, LA, bStream); 
       allele1[LA] = '\0';

       uint32_t LB; 
       ret = fread(&LB, 4, 1, bStream);
       ret = fread(allele0, 1, LB, bStream); 
       allele0[LB] = '\0';
    

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
               fseek(bStream, 6 * Nsamples, SEEK_CUR);
           }
       }
       (void)ret;
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
