\name{close_bgen}
\alias{close_bgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Close BGEN File.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Safely close the BGEN file.
}
\usage{
close_bgen(bgenList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bgenList}{The bgenR class list obtained from \code{open_bgen()}.}
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
None.
}
\references{
%% ~put references to the literature/web site here ~
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
