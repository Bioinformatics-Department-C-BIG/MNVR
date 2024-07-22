\name{WHphase}
\alias{WHphase}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
WHphase
}
\description{
Wrapper for WhatsHap  phase function. It will perform phasing of VCF files using the corresponding alignments (.bam or .cram) and reference genome (GRCh37, hg19 or GRCh38 etc. )
}
\usage{
WHphase(vcf = NULL, bam = NULL, ref = NULL, out = NULL, args = NULL, whpath = NULL)
}
\arguments{
%- maybe also 'usage' for other objects documented here.
\item{vcf}{The unphased vcf file preferably in an indexed '.vcg.gz' format}

\item{bam}{The alignment file used to generate the variant calls in an indexed '.bam' or 'cram' format}

\item{ref}{Reference genome file in '.fasta' format. Must be accompanied by .fai index (create with samtools faidx)}

\item{out}{The output phased .vcf.gz file name or path}

\item{args}{The output phased .vcf.gz file name or path}

\item{whpath}{Path to WhatsHap executable. Defaults to \code{'~/.local/bin/whatshap'} }
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
function (vcf = NULL, bam = NULL, ref = NULL, out = NULL, args = NULL,
    whpath = NULL)
{
    str <- paste(WHpath, "phase", paste(args, collapse = " "),
        "-r", ref, "-o", out, vcf, bam)
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
