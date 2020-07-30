# bgenR

Package to read in BGEN file formats.  
My version that doesn't require index files. 

Still in development
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