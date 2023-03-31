#include "read_bgen.h"

bool bgenIsOpen = false;

Rcpp::List query_bgen(SEXP bgenR_in, SEXP seek_in)
{
  Rcpp::List bgenR_list(bgenR_in);
  
  uint Layout = Rcpp::as<uint>(bgenR_list["layout_flag"]);
  if(Layout == 2){
    return(query_bgen13(bgenR_in, seek_in));
  }
  else if(Layout == 1){
    return(query_bgen11(bgenR_in, seek_in));
  }
  
  return(NA_REAL);
} 

Rcpp::List query_bgen_zlib(SEXP bgenR_in, SEXP seek_in)
{
  Rcpp::List bgenR_list(bgenR_in);
  
  uint Layout = Rcpp::as<uint>(bgenR_list["layout_flag"]);
  if(Layout == 2){
    return(query_bgen13_zlib(bgenR_in, seek_in));
  }
  else if(Layout == 1){
    return(query_bgen11(bgenR_in, seek_in));
  }
  
  return(NA_REAL);
} 


Rcpp::String close_bgen(SEXP bgenR_in)
{
  Rcpp::List bgenR_list(bgenR_in);
  
  FILE* fstream = (FILE*)R_ExternalPtrAddr(bgenR_list["fin"]);
  if(fstream == NULL) {
     Rcpp::stop("BGEN file is already closed.");
  } else {
    fseek(fstream, 0, SEEK_END);
    fclose(fstream);
    R_ClearExternalPtr(bgenR_list["fin"]);
  }

  std::vector<llui>* findex = (std::vector<llui>*)R_ExternalPtrAddr(bgenR_list["findex"]);
  if (findex != NULL) {
    R_ClearExternalPtr(bgenR_list["findex"]);
  }
  
  R_ClearExternalPtr(bgenR_list["fcounter"]);
  
  bgenIsOpen = false;
  return("BGEN file closed successfully.");
}



Rcpp::List open_bgen(SEXP bgenfile_in, SEXP getIndices_in)
{
  if (bgenIsOpen) {
    Rcpp::stop("A BGEN file is already opened. Use close_bgen() to close the open BGEN file.");
  }
  
  std::string bgenfile = Rcpp::as<std::string>(bgenfile_in);
  FILE* fstream = fopen(bgenfile.c_str(), "rb");
  if (!fstream) { 
      Rcpp::stop("Cannot not open BGEN file: [" + bgenfile + "]."); 
  }
  bgenIsOpen = true;

  
  // Header Block
  uint offset;
  if (!fread(&offset, 4, 1, fstream)) {
      Rcpp::stop("Cannot read BGEN header block (offset).");
  }

  uint L_H; 
  if (!fread(&L_H, 4, 1, fstream)) {
      Rcpp::stop("Cannot read BGEN header block (LH).");
  }

  uint mbgen;
  if (!fread(&mbgen, 4, 1, fstream)) {
      Rcpp::stop("Cannot read BGEN header block (M).");
  }

  uint nbgen;
  if (!fread(&nbgen, 4, 1, fstream)) {
      Rcpp::stop("Cannot read BGEN header block (N).");
  }

  char magic[5]; 
  if (!fread(magic, 1, 4, fstream)) {
      Rcpp::stop("Cannot read BGEN header block (magic numbers).");
  }
  magic[4] = '\0';
  if (!(magic[0] == 'b' && magic[1] == 'g' && magic[2] == 'e' && magic[3] == 'n')) {
      Rcpp::stop("BGEN file's magic number does not match 'b', 'g', 'e', 'n'.");
  }

  fseek(fstream, L_H - 20, SEEK_CUR);

  uint flags; 
  if (!fread(&flags, 4, 1, fstream)) {
      Rcpp::stop("Cannot read BGEN header block (flags).");
  }
  

  // Header Block Flags
  uint compression_flag = flags & 3;
  if (compression_flag != 0 && compression_flag != 1 && compression_flag != 2) {
      Rcpp::stop("BGEN compression flag should be 0, 1, or 2.");
  }

  uint layout_flag = (flags >> 2) & 0xf;
  if (layout_flag != 1 && layout_flag != 2) {
      Rcpp::stop("BGEN layout flag should be 1 or 2.");
  }

  uint sampleID_flag = flags >> 31;
  if (sampleID_flag != 0 && sampleID_flag != 1) {
      Rcpp::stop("BGEN sample identifier flag should be 0 or 1.");
  }

  Rcpp::CharacterVector sampleID;
  if (sampleID_flag == 1){
      sampleID = CharacterVector(nbgen);
      uint LS1;
      if (!fread(&LS1, 4, 1, fstream)) {
          Rcpp::stop("Cannot read BGEN sample block (LS).");
      }
     
      uint nsamples;
      if (!fread(&nsamples, 4, 1, fstream)) {
          Rcpp::stop("Cannot read BGEN sample block (N).");
      }
      if (nsamples != nbgen) {
          Rcpp::stop("Number of sample identifiers (" + std::to_string(nsamples) + ") does not match the number of BGEN samples (" + std::to_string(nbgen) + ") specified in the header block.");
      }
    
    
      char* samID = new char[65537];
      int ret;
      for (uint n = 0; n < nsamples; n++) {
           uint16_t LSID;
           ret = fread(&LSID, 2, 1,    fstream);
           ret = fread(samID, 1, LSID, fstream);
           samID[LSID] = '\0';
           
           sampleID[n] = std::string(samID);
      }
      (void)ret;
      delete[] samID;

  } else {
      sampleID = CharacterVector(1);
  }
  fseek(fstream, offset + 4, SEEK_SET);
  

  // Create vector of variant position in the binary file
  bool getIndices = Rcpp::as<bool>(getIndices_in);
  std::vector<llui>* findex = get_bytes(getIndices, fstream, offset, layout_flag, mbgen, nbgen, compression_flag);
  
  // Counter
  std::vector<uint>* fcounter = new std::vector<uint>(1);
  
  
  return(Rcpp::List::create(Named("offset")           = offset,
                            Named("M")                = mbgen,
                            Named("N")                = nbgen,
                            Named("compression_flag") = compression_flag,
                            Named("layout_flag")      = layout_flag,
                            Named("sampleID_flag")    = sampleID_flag,
                            Named("sampleID")         = sampleID,
                            Named("fin")              = R_MakeExternalPtr(fstream,  R_NilValue, R_NilValue),
                            Named("fcounter")         = R_MakeExternalPtr(fcounter, R_NilValue, R_NilValue),
                            Named("findex")           = R_MakeExternalPtr(findex, R_NilValue, R_NilValue))
           );
}



void Bgen13GetTwoVals(const unsigned char* prob_start, uint32_t bit_precision, uintptr_t offset, uintptr_t* first_val_ptr, uintptr_t* second_val_ptr) 
{
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



Rcpp::List query_bgen13(SEXP bgenR_in, SEXP seek_in)
{
   Rcpp::List bgenR_list(bgenR_in);
  
   FILE* fstream = (FILE*)R_ExternalPtrAddr(bgenR_list["fin"]);
   if (fstream == NULL) {
     Rcpp::stop("BGEN file is not open. Open the BGEN file with open_bgen().");
   }
  
   
   // Initialize objects
   int ret;
   uint m = Rcpp::as<uint>(bgenR_list["M"]);
   uint n = Rcpp::as<uint>(bgenR_list["N"]);
   uint c = Rcpp::as<uint>(bgenR_list["compression_flag"]);
   
   NumericMatrix probs(n, 2);
   NumericVector dosVec(n);
   
   struct libdeflate_decompressor* decompressor = libdeflate_alloc_decompressor();
   
   std::vector<uint>* fcounter = (std::vector<uint>*)R_ExternalPtrAddr(bgenR_list["fcounter"]);
   std::vector<uint>& current_count = *fcounter;   
   
   
   // Seek
   uint seek = Rcpp::as<uint>(seek_in);
   if (seek > 0) {
     std::vector<ulli>* findex = (std::vector<ulli>*)R_ExternalPtrAddr(bgenR_list["findex"]);
     if (findex == NULL) {
       Rcpp::stop("Specify `getIndicies=TRUE` with open_bgen() to use the `seek` argument of this function.");
     }
     
     std::vector<ulli>& index  = *findex;
     fseek(fstream, index[seek-1], SEEK_SET);
     current_count[0] = seek-1;
   }
   
   
   // Check counter
   if (current_count[0] >= m){
     Rcpp::stop("End of BGEN file has been reached. Close the file or seek to previous variants.");
   }

  

   
   // Query
   llui curr = ftell(fstream);
   
   uint16_t LS;
   ret = fread(&LS, 2, 1, fstream);
   char* snpID = new char[LS + 1];
   ret = fread(snpID, 1, LS, fstream);
   snpID[LS] = '\0';

   uint16_t LR;
   ret = fread(&LR, 2, 1, fstream);
   char* rsID = new char[LR + 1];
   ret = fread(rsID, 1, LR, fstream);
   rsID[LR] = '\0';

   uint16_t LC;
   ret = fread(&LC, 2, 1, fstream);
   char* chr = new char[LC + 1];
   ret = fread(chr, 1, LC, fstream);
   chr[LC] = '\0';

   uint32_t physpos;
   ret = fread(&physpos, 4, 1, fstream);

   uint16_t LKnum;
   ret = fread(&LKnum, 2, 1, fstream);
   if (LKnum != 2U) {
       char* allele1 = new char[65537];
       for (uint16_t a = 0; a < LKnum; a++) {
            uint32_t LA;
            ret = fread(&LA, 4, 1, fstream);
            ret = fread(allele1, 1, LA, fstream);
       }

       if (c > 0) {
           uint zLen;
           ret = fread(&zLen, 4, 1, fstream);
           fseek(fstream, 4 + zLen - 4, SEEK_CUR);

       }
       else {
           uint zLen;
           ret = fread(&zLen, 4, 1, fstream);
           fseek(fstream, zLen, SEEK_CUR);
       }
       delete[] allele1;
       
       current_count[0]++;
       return(Rcpp::List::create(Named("SNPID")         = Rcpp::String(snpID),
                                 Named("RSID")          = Rcpp::String(rsID),
                                 Named("Chromosome")    = Rcpp::String(chr),
                                 Named("Position")      = physpos,
                                 Named("Alleles")       = LKnum,
                                 Named("Allele1")       = NA_REAL,
                                 Named("Allele2")       = NA_REAL,
                                 Named("AF")            = NA_REAL,
                                 Named("Probabilities") = NA_REAL,
                                 Named("Dosages")       = NA_REAL));
   }

   uint32_t LA;
   ret = fread(&LA, 4, 1, fstream);
   char* allele1 = new char[LA + 1];
   ret = fread(allele1, 1, LA, fstream);
   allele1[LA] = '\0';

   uint32_t LB;
   ret = fread(&LB, 4, 1, fstream);
   char* allele2 = new char[LB + 1];
   ret = fread(allele2, 1, LB, fstream);
   allele2[LB] = '\0';

   uint cLen;
   ret = fread(&cLen, 4, 1, fstream);
   
   uchar* prob_start;
   std::vector<uchar> zBuf12;
   std::vector<uchar> shortBuf12;
   if (c == 1) {
       zBuf12.resize(cLen - 4);
       uint dLen;
       ret = fread(&dLen, 4, 1, fstream);
       ret = fread(&zBuf12[0], 1, cLen - 4, fstream);
       shortBuf12.resize(dLen);

       uLongf destLen = dLen;
       if (libdeflate_zlib_decompress(decompressor, &zBuf12[0], cLen - 4, &shortBuf12[0], destLen, NULL) != LIBDEFLATE_SUCCESS) {
           fseek(fstream, curr, SEEK_SET);
           Rcpp::stop("Decompressing " + std::string(rsID) + " genotype block failed with libdeflate.");
       }
       prob_start = &shortBuf12[0];
   }
   else if (c == 2) {
       zBuf12.resize(cLen - 4);
       uint dLen;
       ret = fread(&dLen, 4, 1, fstream);
       ret = fread(&zBuf12[0], 1, cLen - 4, fstream);
       shortBuf12.resize(dLen);

       uLongf destLen = dLen;
       size_t ret = ZSTD_decompress(&shortBuf12[0], destLen, &zBuf12[0], cLen - 4);
       if (ret > destLen) {
           if (ZSTD_isError(ret)) {
               fseek(fstream, curr, SEEK_SET);
               Rcpp::stop("Decompressing " + std::string(rsID) + " genotype block failed with zstd.");
           }
       }
       prob_start = &shortBuf12[0];
   }
   else {
       zBuf12.resize(cLen);
       ret = fread(&zBuf12[0], 1, cLen, fstream);
       prob_start = &zBuf12[0];
   }
   (void)ret;

   uint32_t N;
   memcpy(&N, prob_start, sizeof(int32_t));
   uint16_t K;
   memcpy(&K, &(prob_start[4]), sizeof(int16_t));

   const uint32_t min_ploidy = prob_start[6];
   if (min_ploidy != 2) {
       fseek(fstream, curr, SEEK_SET);
       Rcpp::stop("Variants with ploidy != 2 is currently not supported.");
   }
   const uint32_t max_ploidy = prob_start[7];
   if (max_ploidy != 2) {
       fseek(fstream, curr, SEEK_SET);
       Rcpp::stop("Variants with ploidy != 2 is currently not supported.");
   }

   const unsigned char* missing_and_ploidy_info = &(prob_start[8]);
   const unsigned char* probs_start = &(prob_start[10 + N]);
   const uint32_t is_phased = probs_start[-2];
   if (is_phased != 0 && is_phased != 1) {
       fseek(fstream, curr, SEEK_SET);
       Rcpp::stop("phased value must be 0 or 1.");
   } 

   const uint32_t B = probs_start[-1];
   if (B != 8 && B != 16 && B != 24 && B != 32) {
       fseek(fstream, curr, SEEK_SET);
       Rcpp::stop("Bits to store probabilities must be 8, 16, 24, or 32.");
   }
   const uintptr_t numer_mask = (1U << B) - 1;
   const uintptr_t probs_offset = B / 8;


   double gmean = 0.0;
   double nmiss = 0.0;
   if (!is_phased) {
      for (uint32_t i = 0; i < N; i++) {
           const uint32_t missing_and_ploidy = missing_and_ploidy_info[i];

           if(missing_and_ploidy == 2){
              uintptr_t numer_aa;
              uintptr_t numer_ab;
              Bgen13GetTwoVals(probs_start, B, probs_offset, &numer_aa, &numer_ab);

              double p11 = numer_aa / double(1.0 * (numer_mask));
              double p10 = numer_ab / double(1.0 * (numer_mask));
              double dosage = 2 * (1 - p11 - p10) + p10;
              probs(i, 0) = p11;
              probs(i, 1) = p10;
              dosVec[i]   = dosage;
              gmean += dosage;
              
              probs_start += (probs_offset * 2);
          }
          else if (missing_and_ploidy == 130) {
              probs(i, 0) = NA_REAL;
              probs(i, 1) = NA_REAL;
              dosVec[i]   = NA_REAL;
              nmiss+=1.0;
          }
          else {
              fseek(fstream, curr, SEEK_SET);
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
                                  
                 double p11 = numer_aa / double(1.0 * (numer_mask));
                 double p10 = numer_ab / double(1.0 * (numer_mask));
                 double dosage = 2 - (p11 + p10);
                 probs(i, 0) = p11;
                 probs(i, 1) = p10;
                 dosVec[i]   = dosage;
                 gmean += dosage;
                 
                 probs_start += (probs_offset * 2);
             }
             else if (missing_and_ploidy == 130) {
                 probs(i, 0) = NA_REAL;
                 probs(i, 1) = NA_REAL;
                 dosVec[i]   = NA_REAL;
                 nmiss+=1;
             }
             else {
                 fseek(fstream, curr, SEEK_SET);
                 Rcpp::stop("Variants with ploidy != 2 is currently not supported.");
             }
        }

   }
   double AF = gmean / (n-nmiss) / 2.0;

   delete[] snpID;
   delete[] rsID;
   delete[] chr;
   delete[] allele1;
   delete[] allele2;
   libdeflate_free_decompressor(decompressor);

   current_count[0]++;
   return(Rcpp::List::create(Named("SNPID")         = Rcpp::String(snpID),
                             Named("RSID")          = Rcpp::String(rsID),
                             Named("Chromosome")    = Rcpp::String(chr),
                             Named("Position")      = physpos,
                             Named("Alleles")       = LKnum,
                             Named("Allele1")       = Rcpp::String(allele1),
                             Named("Allele2")       = Rcpp::String(allele2),
                             Named("AF")            = AF,
                             Named("Probabilities") = probs,
                             Named("Missing")       = nmiss,
                             Named("Dosages")       = dosVec));
}


Rcpp::List query_bgen13_zlib(SEXP bgenR_in, SEXP seek_in)
{
  Rcpp::List bgenR_list(bgenR_in);
  
  FILE* fstream = (FILE*)R_ExternalPtrAddr(bgenR_list["fin"]);
  if (fstream == NULL) {
    Rcpp::stop("BGEN file is not open. Open the BGEN file with open_bgen().");
  }
  
  
  // Initialize objects
  int ret;
  uint m = Rcpp::as<uint>(bgenR_list["M"]);
  uint n = Rcpp::as<uint>(bgenR_list["N"]);
  uint c = Rcpp::as<uint>(bgenR_list["compression_flag"]);
  
  NumericMatrix probs(n, 2);
  NumericVector dosVec(n);
  
  std::vector<uint>* fcounter = (std::vector<uint>*)R_ExternalPtrAddr(bgenR_list["fcounter"]);
  std::vector<uint>& current_count = *fcounter;   
  
  
  // Seek
  uint seek = Rcpp::as<uint>(seek_in);
  if (seek > 0) {
    std::vector<ulli>* findex = (std::vector<ulli>*)R_ExternalPtrAddr(bgenR_list["findex"]);
    if (findex == NULL) {
      Rcpp::stop("Specify `getIndicies=TRUE` with open_bgen() to use the `seek` argument of this function.");
    }
    
    std::vector<ulli>& index  = *findex;
    fseek(fstream, index[seek-1], SEEK_SET);
    current_count[0] = seek-1;
  }
  
  
  // Check counter
  if (current_count[0] >= m){
    Rcpp::stop("End of BGEN file has been reached. Close the file or seek to previous variants.");
  }
  
  
  
  
  // Query
  llui curr = ftell(fstream);
  
  uint16_t LS;
  ret = fread(&LS, 2, 1, fstream);
  char* snpID = new char[LS + 1];
  ret = fread(snpID, 1, LS, fstream);
  snpID[LS] = '\0';
  std::string snpID2(snpID);
  
  uint16_t LR;
  ret = fread(&LR, 2, 1, fstream);
  char* rsID = new char[LR + 1];
  ret = fread(rsID, 1, LR, fstream);
  rsID[LR] = '\0';
  std::string rsID2(rsID);
  
  uint16_t LC;
  ret = fread(&LC, 2, 1, fstream);
  char* chr = new char[LC + 1];
  ret = fread(chr, 1, LC, fstream);
  chr[LC] = '\0';
  std::string chr2(chr);
  
  uint32_t physpos;
  ret = fread(&physpos, 4, 1, fstream);
  
  uint16_t LKnum;
  ret = fread(&LKnum, 2, 1, fstream);
  if (LKnum != 2U) {
    char* allele1 = new char[65537];
    for (uint16_t a = 0; a < LKnum; a++) {
      uint32_t LA;
      ret = fread(&LA, 4, 1, fstream);
      ret = fread(allele1, 1, LA, fstream);
    }
    
    if (c > 0) {
      uint zLen;
      ret = fread(&zLen, 4, 1, fstream);
      fseek(fstream, 4 + zLen - 4, SEEK_CUR);
      
    }
    else {
      uint zLen;
      ret = fread(&zLen, 4, 1, fstream);
      fseek(fstream, zLen, SEEK_CUR);
    }
    delete[] allele1;
    
    current_count[0]++;
    return(Rcpp::List::create(Named("SNPID")         = Rcpp::String(snpID2),
                              Named("RSID")          = Rcpp::String(rsID2),
                              Named("Chromosome")    = Rcpp::String(chr2),
                              Named("Position")      = physpos,
                              Named("Alleles")       = LKnum,
                              Named("Allele1")       = NA_REAL,
                              Named("Allele2")       = NA_REAL,
                              Named("AF")            = NA_REAL,
                              Named("Probabilities") = NA_REAL,
                              Named("Dosages")       = NA_REAL));
  }
  
  uint32_t LA;
  ret = fread(&LA, 4, 1, fstream);
  char* allele1 = new char[LA + 1];
  ret = fread(allele1, 1, LA, fstream);
  allele1[LA] = '\0';
  
  uint32_t LB;
  ret = fread(&LB, 4, 1, fstream);
  char* allele2 = new char[LB + 1];
  ret = fread(allele2, 1, LB, fstream);
  allele2[LB] = '\0';
  
  uint cLen;
  ret = fread(&cLen, 4, 1, fstream);
  
  uchar* prob_start;
  std::vector<uchar> zBuf12;
  std::vector<uchar> shortBuf12;
  if (c == 1) {
    zBuf12.resize(cLen - 4);
    uint dLen;
    ret = fread(&dLen, 4, 1, fstream);
    ret = fread(&zBuf12[0], 1, cLen - 4, fstream);
    shortBuf12.resize(dLen);
    
    uLongf destLen = dLen;
    if (uncompress(&shortBuf12[0], &destLen, &zBuf12[0], cLen-4) != Z_OK || destLen != dLen) {
      fseek(fstream, curr, SEEK_SET);
      Rcpp::stop("Decompressing " + std::string(rsID) + " genotype block failed with libdeflate.");
    }
    prob_start = &shortBuf12[0];
  }
  else if (c == 2) {
    zBuf12.resize(cLen - 4);
    uint dLen;
    ret = fread(&dLen, 4, 1, fstream);
    ret = fread(&zBuf12[0], 1, cLen - 4, fstream);
    shortBuf12.resize(dLen);
    
    uLongf destLen = dLen;
    size_t ret = ZSTD_decompress(&shortBuf12[0], destLen, &zBuf12[0], cLen - 4);
    if (ret > destLen) {
      if (ZSTD_isError(ret)) {
        fseek(fstream, curr, SEEK_SET);
        Rcpp::stop("Decompressing " + std::string(rsID) + " genotype block failed with zstd.");
      }
    }
    prob_start = &shortBuf12[0];
  }
  else {
    zBuf12.resize(cLen);
    ret = fread(&zBuf12[0], 1, cLen, fstream);
    prob_start = &zBuf12[0];
  }
  (void)ret;
  
  uint32_t N;
  memcpy(&N, prob_start, sizeof(int32_t));
  uint16_t K;
  memcpy(&K, &(prob_start[4]), sizeof(int16_t));
  
  const uint32_t min_ploidy = prob_start[6];
  if (min_ploidy != 2) {
    fseek(fstream, curr, SEEK_SET);
    Rcpp::stop("Variants with ploidy != 2 is currently not supported.");
  }
  const uint32_t max_ploidy = prob_start[7];
  if (max_ploidy != 2) {
    fseek(fstream, curr, SEEK_SET);
    Rcpp::stop("Variants with ploidy != 2 is currently not supported.");
  }
  
  const unsigned char* missing_and_ploidy_info = &(prob_start[8]);
  const unsigned char* probs_start = &(prob_start[10 + N]);
  const uint32_t is_phased = probs_start[-2];
  if (is_phased != 0 && is_phased != 1) {
    fseek(fstream, curr, SEEK_SET);
    Rcpp::stop("phased value must be 0 or 1.");
  } 
  
  const uint32_t B = probs_start[-1];
  if (B != 8 && B != 16 && B != 24 && B != 32) {
    fseek(fstream, curr, SEEK_SET);
    Rcpp::stop("Bits to store probabilities must be 8, 16, 24, or 32.");
  }
  const uintptr_t numer_mask = (1U << B) - 1;
  const uintptr_t probs_offset = B / 8;
  
  
  double gmean = 0.0;
  double nmiss = 0.0;
  if (!is_phased) {
    for (uint32_t i = 0; i < N; i++) {
      const uint32_t missing_and_ploidy = missing_and_ploidy_info[i];
      
      if(missing_and_ploidy == 2){
        uintptr_t numer_aa;
        uintptr_t numer_ab;
        Bgen13GetTwoVals(probs_start, B, probs_offset, &numer_aa, &numer_ab);
        
        double p11 = numer_aa / double(1.0 * (numer_mask));
        double p10 = numer_ab / double(1.0 * (numer_mask));
        double dosage = 2 * (1 - p11 - p10) + p10;
        probs(i, 0) = p11;
        probs(i, 1) = p10;
        dosVec[i]   = dosage;
        gmean += dosage;
        
        probs_start += (probs_offset * 2);
      }
      else if (missing_and_ploidy == 130) {
        probs(i, 0) = NA_REAL;
        probs(i, 1) = NA_REAL;
        dosVec[i]   = NA_REAL;
        nmiss+=1.0;
      }
      else {
        fseek(fstream, curr, SEEK_SET);
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
        
        double p11 = numer_aa / double(1.0 * (numer_mask));
        double p10 = numer_ab / double(1.0 * (numer_mask));
        double dosage = 2 - (p11 + p10);
        probs(i, 0) = p11;
        probs(i, 1) = p10;
        dosVec[i]   = dosage;
        gmean += dosage;
        
        probs_start += (probs_offset * 2);
      }
      else if (missing_and_ploidy == 130) {
        probs(i, 0) = NA_REAL;
        probs(i, 1) = NA_REAL;
        dosVec[i]   = NA_REAL;
        nmiss+=1;
      }
      else {
        fseek(fstream, curr, SEEK_SET);
        Rcpp::stop("Variants with ploidy != 2 is currently not supported.");
      }
    }
    
  }
  double AF = gmean / (n-nmiss) / 2.0;
  
  delete[] snpID;
  delete[] rsID;
  delete[] chr;
  delete[] allele1;
  delete[] allele2;
  
  current_count[0]++;
  return(Rcpp::List::create(Named("SNPID")         = Rcpp::String(snpID),
                            Named("RSID")          = Rcpp::String(rsID),
                            Named("Chromosome")    = Rcpp::String(chr),
                            Named("Position")      = physpos,
                            Named("Alleles")       = LKnum,
                            Named("Allele1")       = Rcpp::String(allele1),
                            Named("Allele2")       = Rcpp::String(allele2),
                            Named("AF")            = AF,
                            Named("Probabilities") = probs,
                            Named("Missing")       = nmiss,
                            Named("Dosages")       = dosVec));
}


Rcpp::List query_bgen11(SEXP bgenR_in, SEXP seek_in)
{
  Rcpp::List bgenR_list(bgenR_in);
  
  FILE* fstream = (FILE*)R_ExternalPtrAddr(bgenR_list["fin"]);
  if (fstream == NULL) {
    Rcpp::stop("BGEN file is not open. Open the BGEN file with open_bgen().");
  }

  
  // Initialize objects
  int ret;
  uint m = Rcpp::as<uint>(bgenR_list["M"]);
  uint n = Rcpp::as<uint>(bgenR_list["N"]);
  uint c = Rcpp::as<uint>(bgenR_list["compression_flag"]);
  uLongf destLen1 = 6 * n;
  
  NumericMatrix probs(n, 3);
  NumericVector dosVec(n);
  
  struct libdeflate_decompressor* decompressor = libdeflate_alloc_decompressor();

  
  std::vector<uint>* fcounter = (std::vector<uint>*)R_ExternalPtrAddr(bgenR_list["fcounter"]);
  std::vector<uint>& current_count = *fcounter;   
  
  
  // Seek
  uint seek = Rcpp::as<uint>(seek_in);

  if (seek > 0) {
    std::vector<ulli>* findex = (std::vector<ulli>*)R_ExternalPtrAddr(bgenR_list["findex"]);
    if (findex == NULL) {
      Rcpp::stop("Specify `getIndicies=TRUE` with open_bgen() to use the `seek` argument of this function.");
    }

    std::vector<ulli>& index  = *findex;
    fseek(fstream, index[seek-1], SEEK_SET);
    current_count[0] = seek-1;
  }

  
  // Check counter
  if (current_count[0] >= m){
    Rcpp::stop("End of BGEN file has been reached. Close the file or seek to previous variants.");
  }

  
  
  // Query
  llui curr = ftell(fstream);

  std::vector<uchar> zBuf11;
  std::vector<uint16_t> shortBuf11;
  if (c == 0) {
    zBuf11.resize(destLen1);
  } else {
    shortBuf11.resize(destLen1);
  }


  uint nrow;
  ret = fread(&nrow, 4, 1, fstream);
  if (nrow != n) {
    fseek(fstream, curr, SEEK_SET);
    Rcpp::stop("Number of samples with genotype probabilities ()" + std::to_string(nrow) + ") does not match the number of sample in BGEN header block (" + std::to_string(n) +").");
  }

  uint16_t LS;
  ret = fread(&LS, 2, 1, fstream);
  char* snpID = new char[LS + 1];
  ret = fread(snpID, 1, LS, fstream);
  snpID[LS] = '\0';

  uint16_t LR;
  ret = fread(&LR, 2, 1, fstream);
  char* rsID = new char[LR + 1];
  ret = fread(rsID, 1, LR, fstream);
  rsID[LR] = '\0';

  uint16_t LC;
  ret = fread(&LC, 2, 1, fstream);
  char* chr = new char[LC + 1];
  ret = fread(chr, 1, LC, fstream);
  chr[LC] = '\0';

  uint32_t physpos;
  ret = fread(&physpos, 4, 1, fstream);

  uint32_t LA;
  ret = fread(&LA, 4, 1, fstream);
  char* allele1 = new char[LA + 1];
  ret = fread(allele1, 1, LA, fstream);
  allele1[LA] = '\0';

  uint32_t LB;
  ret = fread(&LB, 4, 1, fstream);
  char* allele2 = new char[LB + 1];
  ret = fread(allele2, 1, LB, fstream);
  allele2[LB] = '\0';


  uint16_t* probs_start;
  if (c == 1) {
    uint cLen;
    ret = fread(&cLen, 4, 1, fstream);
    zBuf11.resize(cLen);
    ret = fread(&zBuf11[0], 1, cLen, fstream);

    if (libdeflate_zlib_decompress(decompressor, &zBuf11[0], cLen, &shortBuf11[0], destLen1, NULL) != LIBDEFLATE_SUCCESS) {
      fseek(fstream, curr, SEEK_SET);
      Rcpp::stop("Decompressing " + std::string(rsID) + " genotype block failed with libdeflate.");
    }

    probs_start = &shortBuf11[0];

  }
  else {
    ret = fread(&zBuf11[0], 1, destLen1, fstream);
    probs_start = reinterpret_cast<uint16_t*>(&zBuf11[0]);
  }
  (void)ret;


  const double scale = 1.0 / 32768;
  double gmean = 0.0;
  double nmiss = 0.0;
  for (uint i = 0; i < n; i++) {
    double p11 = probs_start[3 * i] * scale;
    double p10 = probs_start[3 * i + 1] * scale;
    double p00 = probs_start[3 * i + 2] * scale;

    if (p11 == 0.0 && p10 == 0.0 && p00 == 0.0) {
      probs(i, 0) = NA_REAL;
      probs(i, 1) = NA_REAL;
      probs(i, 2) = NA_REAL;
      dosVec[i]   = NA_REAL;
      nmiss+=1.0;

    } 
    else {
      double pTot = p11 + p10 + p00;
      double dosage = (2 * p00 + p10) / pTot;

      probs(i, 0) = p11;
      probs(i, 1) = p10;
      probs(i, 2) = p00;

      dosVec[i] = dosage;
      gmean += dosage;

    }
  }
  double AF = gmean / (n-nmiss) / 2.0;

  
  delete[] snpID;
  delete[] rsID;
  delete[] chr;
  delete[] allele1;
  delete[] allele2;
  libdeflate_free_decompressor(decompressor);

  current_count[0]++;
  return(Rcpp::List::create(Named("SNPID")         = Rcpp::String(snpID),
                            Named("RSID")          = Rcpp::String(rsID),
                            Named("Chromosome")    = Rcpp::String(chr),
                            Named("Position")      = physpos,
                            Named("Alleles")       = 2,
                            Named("Allele1")       = Rcpp::String(allele1),
                            Named("Allele2")       = Rcpp::String(allele2),
                            Named("AF")            = AF,
                            Named("Probabilities") = probs,
                            Named("Missing")       = nmiss,
                            Named("Dosages")       = dosVec));
}  


Rcpp::DataFrame variant_block(SEXP bgenR_in)
{
  if(!bgenIsOpen) {
    Rcpp::stop("BGEN file is not open. Open the BGEN file with open_bgen().");
  }

  Rcpp::List bgenR_list(bgenR_in);
  FILE* fstream = (FILE*)R_ExternalPtrAddr(bgenR_list["fin"]);
  if (fstream == NULL) {
    Rcpp::stop("BGEN file is not open. Open the BGEN file with open_bgen().");
  }

  
  // Intialize objects
  char* snpID   = new char[65537];
  char* rsID    = new char[65537];
  char* chr     = new char[65537];
  char* allele1 = new char[65537];
  char* allele2 = new char[65537];
  
  uint offset = Rcpp::as<uint>(bgenR_list["offset"]);
  uint m = Rcpp::as<uint>(bgenR_list["M"]);
  uint n = Rcpp::as<uint>(bgenR_list["N"]);
  uint l = Rcpp::as<uint>(bgenR_list["layout_flag"]);
  uint c = Rcpp::as<uint>(bgenR_list["compression_flag"]);

  Rcpp::CharacterVector  vecSNPID(m);
  Rcpp::CharacterVector  vecRSID(m);
  Rcpp::CharacterVector  vecCHR(m);
  Rcpp::NumericVector    vecPOS(m);
  Rcpp::NumericVector    vecLK(m);
  Rcpp::CharacterVector  vecA1(m);
  Rcpp::CharacterVector  vecA2(m);

  
  int ret;
  fseek(fstream, offset + 4, SEEK_SET);
  for (uint i = 0; i < m; i++) {
       if (l == 1) {
           uint Nrow;
           ret = fread(&Nrow, 4, 1, fstream);
       }

       uint16_t LS;
       ret = fread(&LS, 2, 1, fstream);
       ret = fread(snpID, 1, LS, fstream);
       snpID[LS] = '\0';

       uint16_t LR;
       ret = fread(&LR, 2, 1, fstream);
       ret = fread(rsID, 1, LR, fstream);
       rsID[LR] = '\0';

       uint16_t LC;
       ret = fread(&LC, 2, 1, fstream);
       ret = fread(chr, 1, LC, fstream);
       chr[LC] = '\0';

       uint32_t physpos;
       ret = fread(&physpos, 4, 1, fstream);

       uint16_t LKnum;
       if (l == 2) {
           ret = fread(&LKnum, 2, 1, fstream);

           if (LKnum != 2U) {
               for (uint16_t a = 0; a < LKnum; a++) {
                    uint32_t LA;
                    ret = fread(&LA, 4, 1, fstream);
                    ret = fread(allele1, 1, LA, fstream);
               }

               if (c > 0) {
                   uint zLen;
                   ret = fread(&zLen, 4, 1, fstream);
                   fseek(fstream, 4 + zLen - 4, SEEK_CUR);

               }
               else {
                   uint zLen;
                   ret = fread(&zLen, 4, 1, fstream);
                   fseek(fstream, zLen, SEEK_CUR);
               }

               vecSNPID[i] = Rcpp::String(snpID);
               vecRSID[i]  = Rcpp::String(rsID);
               vecCHR[i]   = Rcpp::String(chr);
               vecPOS[i]   = physpos;
               vecLK[i]    = LKnum;
               vecA1[i]    = NA_STRING;
               vecA2[i]    = NA_STRING;
               continue;
           }
       }
       else {
           LKnum = 2U;
       }

       uint32_t LA;
       ret = fread(&LA, 4, 1, fstream);
       ret = fread(allele1, 1, LA, fstream);
       allele1[LA] = '\0';

       uint32_t LB;
       ret = fread(&LB, 4, 1, fstream);
       ret = fread(allele2, 1, LB, fstream);
       allele2[LB] = '\0';

       if (l == 2) {
           if (c > 0) {
               uint zLen;
               ret = fread(&zLen, 4, 1, fstream);
               fseek(fstream, 4 + zLen - 4, SEEK_CUR);

           }
           else {
               uint zLen;
               ret = fread(&zLen, 4, 1, fstream);
               fseek(fstream, zLen, SEEK_CUR);
           }
       }
       else {
           if (c == 1) {
               uint zLen;
               ret = fread(&zLen, 4, 1, fstream);
               fseek(fstream, zLen, SEEK_CUR);
           }
           else {
               fseek(fstream, 6 * n, SEEK_CUR);
           }
       }
       
       vecSNPID[i] = Rcpp::String(snpID),
       vecRSID[i]  = Rcpp::String(rsID),
       vecCHR[i]   = Rcpp::String(chr),
       vecPOS[i]   = physpos;
       vecLK[i]    = LKnum;
       vecA1[i]    = Rcpp::String(allele1);
       vecA2[i]    = Rcpp::String(allele2);
  }
  fseek(fstream, offset + 4, SEEK_SET);

  
  
  (void)ret;
  delete[] snpID;
  delete[] rsID;
  delete[] chr;
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



std::vector<llui>* get_bytes(bool getIndices_in, FILE* & fstream_in, uint offset_in, uint Layout_in, uint M_in, uint N_in, uint Compression_in)
{
  if (!getIndices_in) {
    return nullptr;
  }

  char* snpID   = new char[65537];
  char* rsID    = new char[65537];
  char* chr     = new char[65537];
  char* allele1 = new char[65537];
  char* allele2 = new char[65537];
  std::vector<llui>* bytes = new std::vector<llui>(M_in);
  std::vector<llui> &b = *bytes;

  fseek(fstream_in, offset_in + 4, SEEK_SET);

  int ret;
  for (uint m = 0; m < M_in; m++) {
    b[m] = ftell(fstream_in);

    if (Layout_in == 1) {
      uint Nrow;
      ret = fread(&Nrow, 4, 1, fstream_in);
    }

    uint16_t LS;
    ret = fread(&LS, 2, 1, fstream_in);
    ret = fread(snpID, 1, LS, fstream_in);

    uint16_t LR;
    ret = fread(&LR, 2, 1, fstream_in);
    ret = fread(rsID, 1, LR, fstream_in);

    uint16_t LC;
    ret = fread(&LC, 2, 1, fstream_in);
    ret = fread(chr, 1, LC, fstream_in);

    uint32_t physpos;
    ret = fread(&physpos, 4, 1, fstream_in);

    uint16_t LKnum;
    if (Layout_in == 2) {
        ret = fread(&LKnum, 2, 1, fstream_in);

        for (uint16_t a = 0; a < LKnum; a++) {
             uint32_t LA;
             ret = fread(&LA, 4, 1, fstream_in);
             ret = fread(allele1, 1, LA, fstream_in);
        }


    } else {
        uint32_t LA;
        ret = fread(&LA, 4, 1, fstream_in);
        ret = fread(allele1, 1, LA, fstream_in);

        uint32_t LB;
        ret = fread(&LB, 4, 1, fstream_in);
        ret = fread(allele2, 1, LB, fstream_in);
    }

    if (Layout_in == 2) {
      if (Compression_in > 0) {
        uint zLen;
        ret = fread(&zLen, 4, 1, fstream_in);
        fseek(fstream_in, 4 + zLen - 4, SEEK_CUR);

      }
      else {
        uint zLen;
        ret = fread(&zLen, 4, 1, fstream_in);
        fseek(fstream_in, zLen, SEEK_CUR);
      }
    }
    else {
      if (Compression_in == 1) {
        uint zLen;
        ret = fread(&zLen, 4, 1, fstream_in);
        fseek(fstream_in, zLen, SEEK_CUR);
      }
      else {
        fseek(fstream_in, 6 * N_in, SEEK_CUR);
      }
    }
  }

  fseek(fstream_in, offset_in + 4, SEEK_SET);
  
  (void)ret;
  delete[] snpID;
  delete[] rsID;
  delete[] chr;
  delete[] allele1;
  delete[] allele2;

  return(bytes);
}
