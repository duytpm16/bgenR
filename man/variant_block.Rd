\name{variant_block}
\alias{variant_block}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Variant Block of the BGEN File
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns a \code{data.frame} object containing the variant block information for each variant in the BGEN file.
}
\usage{
variant_block(bgenList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bgenList}{
%%     ~~Describe \code{bgenList} here~~
  The bgenR class list obtained from \code{open_bgen()}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A \code{data.frame} object containing the variant block information for all variants in the BGEN file.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{SNPID}{The variant identifier.}
\item{RSID}{The reference variant identifier.}
\item{Chromosome}{The chromosome.}
\item{Position}{The physical position.}
\item{Alleles}{The number of alleles.}
\item{Allele1}{The first allele.}
\item{Allele2}{The second allele.}
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
  bgen <- open_bgen(bgenfile)
  vb   <- variant_block(bgen) 
  close_bgen(bgen)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{bgen}
