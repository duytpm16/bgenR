#include "read_bgen.h"
#include "../thirdparty/libdeflate-1.6/libdeflate.h"
#include "../thirdparty/zstd-1.4.5/lib/zstd.h"
#include <Rcpp.h>
#include <stdio.h>
#include <zlib.h>




typedef unsigned int   uint;
typedef unsigned char  uchar;
typedef unsigned short ushort;


using namespace std;
using namespace Rcpp;


FILE* bgen_stream;
BGEN bgen;
struct libdeflate_decompressor* decompressor = libdeflate_alloc_decompressor();


Rcpp::CharacterVector read_bgenSampleID();
Rcpp::List query_bgen11();
Rcpp::List query_bgen13();


void close_bgen(){
  fclose(bgen_stream);
}

Rcpp::List set_bgen(SEXP bgenfile_in){
    
    string bgenfile = Rcpp::as<string>(bgenfile_in);
    
    
    bgen_stream = fopen(bgenfile.c_str(), "rb");
    if(!bgen_stream) { Rcpp::stop("ERROR: Cannot not open BGEN file: " + bgenfile + "\n"); }
    

    // Header Block
    fread(&bgen.offset, 4, 1, bgen_stream);
    uint L_H; fread(&L_H, 4, 1, bgen_stream);
    fread(&bgen.Mbgen, 4, 1, bgen_stream);
    fread(&bgen.Nbgen, 4, 1, bgen_stream);
    char magic[5]; fread(magic, 1, 4, bgen_stream);
    magic[4] = '\0';

    fseek(bgen_stream, L_H - 20, SEEK_CUR);
    uint flags; fread(&flags, 4, 1, bgen_stream);


    // Header Block Flags
    bgen.Compression = flags & 3;
    bgen.Layout = (flags >> 2) & 0xf;
    bgen.SampleIdentifiers = flags >> 31;


   
    if (bgen.SampleIdentifiers == 1){
        bgen.sampleID = read_bgenSampleID();
    }



    fseek(bgen_stream, bgen.offset + 4, SEEK_SET);

    return(Rcpp::List::create(Named("offset") = bgen.offset,
                              Named("M") = bgen.Mbgen,
                              Named("N") = bgen.Nbgen,
                              Named("Compression") = bgen.Compression,
                              Named("Layout") = bgen.Layout,
                              Named("SampleIdentifier") = bgen.SampleIdentifiers,
                              Named("SampleID") = bgen.sampleID));

}



Rcpp::CharacterVector read_bgenSampleID(){

  uint maxLA = 65536;
  Rcpp::CharacterVector sampleID(bgen.Nbgen);
  
  uint LS1;  fread(&LS1,  4, 1, bgen_stream);
  uint Nrow; fread(&Nrow, 4, 1, bgen_stream);

  if (Nrow != bgen.Nbgen) {
      Rcpp::Rcout << "Number of sample identifiers (" << Nrow << ") does not match the number of BGEN samples (" << bgen.Nbgen << ")\n";
      Rcpp::stop("");
  }


  char* samID = new char[maxLA + 1];
  for (uint n = 0; n < Nrow; n++) {
       ushort LSID; fread(&LSID, 2, 1, bgen_stream);
       fread(samID, 1, LSID, bgen_stream);
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
  const uint64_t bit_offset = prob_offset * bit_precision;
  uint64_t relevant_bits;
  // This can read slightly past the end of the buffer.
  // Note that with bit_precision=29 and variable ploidy,
  // (bit_offset % CHAR_BIT) == 7 is possible, so we may only get 57 bits when
  // we need 58; thus we don't support 29-31 bits for now.
  memcpy(&relevant_bits, &(prob_start[bit_offset / CHAR_BIT]), sizeof(int64_t));
  relevant_bits = relevant_bits >> (bit_offset % CHAR_BIT);
  *first_val_ptr = relevant_bits & numer_mask;
  *second_val_ptr = (relevant_bits >> bit_precision) & numer_mask;
}







Rcpp::List query_bgen(){
     

     if(bgen.Layout == 2){
        return(query_bgen13());
     }
     else if(bgen.Layout == 1){
        return(query_bgen11());
     }
} 


Rcpp::List query_bgen13(){
  
     uint maxLA = 65536;
     uint maxLB = 65536;
     char* snpID   = new char[maxLA + 1];
     char* rsID    = new char[maxLA + 1];
     char* chrStr  = new char[maxLA + 1];
     char* allele1 = new char[maxLA + 1];
     char* allele0 = new char[maxLA + 1];
     std::vector<uchar> zBuf12;
     std::vector<uchar> shortBuf12;
     uint Nsamples = bgen.Nbgen;
     uint Compression = bgen.Compression;
     
     
     NumericMatrix probs(Nsamples, 2);
     NumericVector dosVec(Nsamples);
     ushort LS; fread(&LS, 2, 1, bgen_stream);
     if (LS > maxLA) {
         maxLA = 2 * LS;
         delete[] snpID;
         snpID = new char[maxLA + 1];
     }
     fread(snpID, 1, LS, bgen_stream); snpID[LS] = '\0';

     
     ushort LR; fread(&LR, 2, 1, bgen_stream);
     if (LR > maxLA) {
         maxLA = 2 * LR;
         delete[] rsID;
         rsID = new char[maxLA + 1];
     }
     fread(rsID, 1, LR, bgen_stream); rsID[LR] = '\0';
     
     ushort LC; fread(&LC, 2, 1, bgen_stream);
     fread(chrStr, 1, LC, bgen_stream); chrStr[LC] = '\0';
     uint physpos; fread(&physpos, 4, 1, bgen_stream);
     ushort LKnum; fread(&LKnum,   2, 1, bgen_stream);
     if ( LKnum != 2U ){
         Rcpp::stop("\nERROR: " + string(rsID) + " does not contain 2 alleles.\n\n");
     }

     ushort LA; fread(&LA, 4, 1, bgen_stream);
     if (LA > maxLA) {
         maxLA = 2 * LA;
         delete[] allele1;
         allele1 = new char[maxLA + 1];
     }
     fread(allele1, 1, LA, bgen_stream); allele1[LA] = '\0';

     ushort LB; fread(&LB, 4, 1, bgen_stream);
     if (LB > maxLB) {
         maxLA = 2 * LB;
         delete[] allele0;
         allele0 = new char[maxLB + 1];
     }
     fread(allele0, 1, LB, bgen_stream); allele0[LB] = '\0';
     
     
     uchar* bufAt;
     uint cLen; fread(&cLen, 4, 1, bgen_stream);
     if (Compression == 1) {
         zBuf12.resize(cLen - 4);
         uint dLen; fread(&dLen, 4, 1, bgen_stream);
         fread(&zBuf12[0], 1, cLen - 4, bgen_stream);
         shortBuf12.resize(dLen);

         uLongf destLen = dLen;
         if (libdeflate_zlib_decompress(decompressor, &zBuf12[0], cLen - 4, &shortBuf12[0], destLen, NULL) != LIBDEFLATE_SUCCESS) {
           cerr << "\nERROR: Decompressing " << snpID << " block failed\n\n";
           exit(1);
         }
         bufAt = &shortBuf12[0];
      }
      else if (Compression == 2) {
          zBuf12.resize(cLen - 4);
          uint dLen; fread(&dLen, 4, 1, bgen_stream);
          fread(&zBuf12[0], 1, cLen - 4, bgen_stream);
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
          fread(&zBuf12[0], 1, cLen, bgen_stream);
          bufAt = &zBuf12[0];
      }

      
      

      uint32_t N; memcpy(&N, bufAt, sizeof(int32_t));
      if ( N != Nsamples ){
          Rcpp::stop("\nERROR: Number of samples with genotype probabilities for variant " + string(snpID) + " does not match the number of sample in BGEN header block.\n\n");
      }
      uint16_t K; memcpy(&K, &(bufAt[4]), sizeof(int16_t));
      if ( K != 2U){
          Rcpp::stop("\nERROR: " + string(snpID) + " does not contain 2 alleles.\n\n");
      }
      const uint32_t min_ploidy = bufAt[6];
      if ( min_ploidy != 2){
           Rcpp::stop("\nERROR: " + string(snpID) + "  does not have minimum ploidy value of 2.\n\n");
      }
      const uint32_t max_ploidy = bufAt[7];
      if ( max_ploidy != 2){
           Rcpp::stop("\nERROR: " + string(snpID) + "  does not maximum ploidy value of 2.\n\n");
      }
      uchar* missing_and_ploidy_iter = &(bufAt[8]);
      uchar* probs_start = &(bufAt[10 + N]);
      uint32_t  is_phased = probs_start[-2];
      if (is_phased){
         Rcpp::stop("\nERROR: " + string(snpID) + "  contains phased genotypes. Currently unsupported.\n\n");
      }
      uint32_t  bit_precision = probs_start[-1];
      uintptr_t numer_mask    = (1U << bit_precision) - 1;




     uintptr_t prob_offset = 0;
     double gmean = 0.0;
     for (uint32_t i = 0; i < N; i++) {

          const uint32_t missing_and_ploidy = *missing_and_ploidy_iter++;
          uintptr_t numer_aa;
          uintptr_t numer_ab;

          switch (missing_and_ploidy) {
          case 1:
            Bgen13GetOneVal(probs_start, prob_offset, bit_precision, numer_mask);
            prob_offset++;
            break;
          case 2:
            Bgen13GetTwoVals(probs_start, prob_offset, bit_precision, numer_mask, &numer_aa, &numer_ab);
            prob_offset += 2;
            break;
          default:
            cerr << "\nERROR: " << snpID << " contains ploidy " << missing_and_ploidy << ". Currently unsupported.\n\n";
            exit(1);
          }

          double p11 = numer_aa / double(1.0 * (numer_mask));
          double p10 = numer_ab / double(1.0 * (numer_mask));
          double dosage = 2 * (1 - p11 - p10) + p10;

          probs(i, 0) = p11;
          probs(i, 1) = p10;
          dosVec[i] = dosage;
          gmean += dosage;

     }

     double AF = gmean / Nsamples / 2.0;
     return(Rcpp::List::create(Named("SNPID") = string(snpID),
                               Named("RSID") = string(rsID),
                               Named("Chromosome") = string(chrStr),
                               Named("Position") = physpos,
                               Named("Allele1")  = string(allele0),
                               Named("Allele2")  = string(allele1),
                               Named("SampleID") = bgen.sampleID,
                               Named("AF") = AF,
                               Named("Probabilities") = probs,
                               Named("Dosages") = dosVec));

}

Rcpp::List query_bgen11(){
  return(Rcpp::List::create(Named("offset") = bgen.offset,
                                   Named("M") = bgen.Mbgen,
                                   Named("N") = bgen.Nbgen,
                                   Named("Compression") = bgen.Compression,
                                   Named("Layout") = bgen.Layout,
                                   Named("SampleIdentifier") = bgen.SampleIdentifiers,
                                   Named("SampleID") = bgen.sampleID));
}  

