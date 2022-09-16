\name{query_bgen}
\alias{query_bgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Query BGEN File
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns a list containing the variant and genotype block information for one variant. The query begins with the first variant in the BGEN file and the function can be called again to retrive the next variant until M variants have been read. Additionally, it is possible to seek forward and backward to a specific variant using the \code{seek} argument.
}
\usage{
query_bgen(bgenList, seek = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bgenList}{
%%     ~~Describe \code{bgenList} here~~
  The bgenR class list obtained from \code{open_bgen()}.
}
  \item{seek}{
%%     ~~Describe \code{seek} here~~
  An integer value specifying which variant in the BGEN file to retrieve information from the variant and genotype blocks. The integer value must be between [1, M], where M is the number of variants in the BGEN file. This argument can only be used if \code{getIndices=TRUE} in the \code{open_bgen()} function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
If the \code{seek} argument is used, this function will return the variant and genotype block information for the variant at the \code{seek} position. If \code{query_bgen} is called again without specifying a \code{seek} value, the function will return the variant and genotype block information for the variant after the variant at the \code{seek} position. For example, if \code{seek=100}, \code{query_bgen} will return information for the 100th variant. If \code{query_bgen} is called again without the \code{seek} argument, information from the variant and genotype blocks of the 101th variant will be returned. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list object containing the variant and genotype block information for one variant.
\item{SNPID}{The variant identifier.}
\item{RSID}{The reference variant identifier.}
\item{Chromosome}{The chromosome.}
\item{Position}{The physical position.}
\item{Alleles}{The number of alleles.}
\item{Allele1}{The first allele.}
\item{Allele2}{The second allele.}
\item{AF}{The allele frequency.}
\item{Probabilites}{For BGENv1.1, this is a (No. of Samples x 3) matrix containing the genotype probabilities. For BGENv1.2, this is a (No. of Samples x 2) matrix containing the genotype probabilties. The last probability can be inferred as one minues the sum of ther other probabilities in this BGEN file version.}
\item{Missing}{The number of samples with missing genotype probabilities.}
\item{Dosages}{A vector containing the dosage value for each sample.}
}
\references{
%% ~put references to the literature/web site here ~
See https://www.well.ox.ac.uk/~gav/bgen_format/spec/latest.html.
}
\author{
%%  ~~who you are~~
Duy T. Pham (duytpham18@gmail.com)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  bgenfile <- system.file("extdata", "bgen12_zlib.bgen", package = "bgenR")
  bgen <- open_bgen(bgenfile, getIndices = TRUE)
  q    <- query_bgen(bgen)
  q100 <- query_bgen(bgen, seek = 100)
  close_bgen(bgen)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{bgen}
