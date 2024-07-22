\name{WHstat}
\alias{WHstat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Phasing statistics
}
\description{
Statistics on the number of variants phased, length and number of phasing blocks.
}
\usage{
WHstat(vcf = NULL, whpath = NULL, out = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{vcf}{The phased vcf file preferably in an indexed '.vcg.gz' format}

\item{out}{The output phased .vcf.gz file name or path}

\item{whpath}{Path to WhatsHap executable. Defaults to \code{ '~/.local/bin/whatshap'} }
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
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (vcf = NULL, whpath = NULL, out = NULL)
{
    str <- paste(whpath, "stats", vcf, ">", out)
    system(str)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
