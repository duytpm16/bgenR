# bgenR

An R package to read in all versions of the [BGEN file](https://www.well.ox.ac.uk/~gav/bgen_format/spec/latest.html) format.

<br />

## Requirements  

- R (>= 3.5.0)  
- Rcpp (>= 1.0.4)  

<br />

## Installation  

```r
# Install devtools 
if (!require("devtools")) {
  install.packages("devtools")
}

# Install from Github repository
devtools::install_github("https://github.com/duytpm16/bgenR") 
```

<br /> 

## Examples

Example 1: Iterate over all variants in the BGEN file.  
```r
library(bgenR)

# Open the BGEN file
bgenfile <- system.file("extdata", "bgen12_zlib.bgen", package = "bgenR")
bgen <- open_bgen(bgenfile)

# Number of variants
M <- bgen$M

# Iterate over all variants in the BGEN file
for(i in 1:M){
  info <- query_bgen(bgen)  
  # do something...
}

# Close the BGEN file
close_bgen(bgen)
```

<br />

Example 2: Retrieve information for a variant at a specific index.  
```r
library(bgenR)

# Open the BGEN file
bgenfile <- system.file("extdata", "bgen12_zlib.bgen", package = "bgenR")
bgen <- open_bgen(bgenfile, getIndices = TRUE)

# Get information for variant at specific index
info <- query_bgen(bgen, seek = 100)
info <- query_bgen(bgen) # Information for the 101th variant
info <- query_bgen(bgen, seek = 10) # Information for the 10th variant

# Close the BGEN file
close_bgen(bgen)
```