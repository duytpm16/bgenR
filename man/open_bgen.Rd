\name{open_bgen}
\alias{open_bgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Open BGEN File
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Opens the BGEN file and returns a bgenR class list containing information from the header and sample identifier blocks if they exists.
}
\usage{
open_bgen(bgenFile, getIndices = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bgenFile}{
%%     ~~Describe \code{bgenFile} here~~
  Path to the BGEN file.
}
  \item{getIndices}{
%%     ~~Describe \code{getIndices} here~~
  Logical value indicating whether the starting byte position of each variant in the BGEN file should be recorded to be used in other functions of the bgenR package. This will increase the runtime of this function for large BGEN files.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A bgenR class list. All elements of the bgenR class list are used in other functions of the bgenR package. Do not modify or delete elements in this bgenR class list.
\item{offset}{The number of bytes to skip to get to the first variant.}
\item{M}{The number of variants.}
\item{N}{The number of samples.}
\item{compression_flag}{The integer value of the compression flag.}
\item{layout_flag}{The integer value of the layout flag.}
\item{sampleID_flag}{The integer value of the sample identifier flag.}
\item{sampleID}{A vector containing the sample identifier names if a sample identifier block exists in the BGEN file.}
\item{fin}{An external pointer to the FILE pointer of the BGEN file.}
\item{fcounter}{An external pointer to keep track of the current variant.}
\item{findex}{An external pointer to a vector of unsigned integers containing the starting byte position of each variant in the BGEN file.}
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
  q    <- query_bgen(bgen)
  close_bgen(bgen)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{bgen}
