# bgenR

Package to read in BGEN file formats into R.  
My version that doesn't require index files. 

Currently designed to read the BGEN file iteratively and for bi-allelic unphased/phased genotypes.  
  
devtools::install_github("https://github.com/duytpm16/bgenR")  

```r
library(bgenR)

bgenInfo = open_bgen("my_bgenfile.bgen")

# Number of variants
M = bgenInfo$M


for(i in 1:M){
    # do something...  

    snpInfo <- query_bgen()  
    
    # do something...
}

close_bgen()

```
