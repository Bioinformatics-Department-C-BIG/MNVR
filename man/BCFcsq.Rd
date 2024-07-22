\name{BCFcsq}
\alias{BCFcsq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict MNV consequences using bcftools csq
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BCFcsq(vcf = NULL, ref = NULL, out = NULL, gff = NULL, args = NULL, bcftpath = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcf}{
Input variant calls (.vcf.gz) to annotate consequences
}
  \item{ref}{
Reference genome file in '.fasta' format. Must be accompanied by .fai index (create with samtools faidx)
}
  \item{out}{
The output phased .vcf.gz file name or path
}
  \item{gff}{
Gene anotation file in '.gff3' format
}
  \item{args}{
%%     ~~Describe \code{args} here~~
}
  \item{bcftpath}{
Path to bcftools executable if other than the defaul ('/usr/bin/bcftools')
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
function (vcf = NULL, ref = NULL, out = NULL, gff = NULL, args = NULL,
    bcftpath = NULL)
{
    str <- paste(bcftpath, "csq", paste(args, collapse = " "),
        "-O z", "-f", ref, "-o", out, "-g", gff, vcf)
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
