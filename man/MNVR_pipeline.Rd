\name{MNVR_pipeline}
\alias{MNVR_pipeline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wrapper function for running the full MNVR pipeline
}
\description{
Wrapper function for running the full MNVR pipeline. If installed on linux system the sample processing will proceed in parallel
}
\usage{
MNVR_pipeline(vcfdir = NULL, bamdir = NULL, whpath = "~/.local/bin/whatshap", bcftpath = "/usr/bin/bcftools", outpath = NULL, HGNC = "/home/konstantinosp/tools/hgnc_complete_set.txt", gffpath = NULL, reference = NULL, hgnc = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcfdir}{
  Vector of vcf.gz files to process
}
  \item{bamdir}{
  Vector of bam files to process. Sample names must match the vcfdir samples
}
  \item{whpath}{
  Path to WhatsHap executable. Defaults to \code{'~/.local/bin/whatshap'}
}
  \item{bcftpath}{
  Path to bcftools executable if other than the defaul ('/usr/bin/bcftools')
}
  \item{outpath}{
    Path of the output directory. This will be create if it does not exist. User must have write permissions for the specific path
}
  \item{HGNC}{
    For selecting consequences for only canonical trasncripts (MANE select). HGNC dataset in (.txt format) for identifying the canonical transcript for each gene symbol. Recommended download from ('https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt')

}
  \item{gffpath}{
Gene anotation file in '.gff3' format. Gene models and coordinates must match the selected reference genome

}
  \item{reference}{
Reference genome file in '.fasta' format. Must be accompanied by .fai index (create with samtools faidx)
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
function (vcfdir = NULL, bamdir = NULL, whpath = "~/.local/bin/whatshap",
    bcftpath = "/usr/bin/bcftools", outpath = NULL, HGNC = "/home/konstantinosp/tools/hgnc_complete_set.txt",
    gffpath = NULL, reference = NULL, hgnc = NULL)
{
    WHpath <- whpath
    BCFtpath <- bcftpath
    if (!dir.exists(outpath)) {
        dir.create(outpath)
    }
    gffpath <- "/home/mariost/DBs/ensembl/Homo_sapiens.GRCh37.87.gff3.gz"
    reference <- "/home/mariost/DBs/references/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta"
    hgnc <- read.delim(HGNC, header = T)
    hgnc$ENselect <- str_match(hgnc$mane_select, "(ENST.+)\\|.+")[,
        2]
    hgnc$ENselect <- gsub("\\..+", "", hgnc$ENselect)
    indir <- gsub("\\/\\/", "/", vcfdir)
    bamdir <- gsub("\\/\\/", "/", bamdir)
    sni <- str_match(indir, ".+/(.+)\\..+vcf.gz$")[, 2]
    bmi <- str_match(bamdir, ".+/(.+)\\..+bam$")[, 2]
    samples <- data.frame(Sample = sni, VCF = indir, BAM = bamdir[match(bmi,
        sni)])
    mclapply(c(1:nrow(samples)), function(s) {
        VCF <- samples$VCF[s]
        BAM <- samples$BAM[s]
        OUTPUT <- gsub(".vcf.gz$", ".WH.vcf.gz", VCF)
        vcfn <- gsub(".+/", "", VCF)
        OUTPUT <- paste0(outpath, gsub(".vcf.gz$", ".WH.vcf.gz",
            vcfn))
        OUTtxt <- paste0(outpath, gsub(".vcf.gz$", ".WHcGT.txt",
            vcfn))
        OUTstat <- paste0(outpath, gsub(".vcf.gz$", ".WHstats.txt",
            vcfn))
        args <- c("--distrust-genotypes --include-homozygous",
            paste("--changed-genotype-list", OUTtxt))
        WHphase(vcf = VCF, bam = BAM, ref = reference, out = OUTPUT,
            args = args, whpath = WHpath)
        system(paste0("tabix -f -p vcf ", OUTPUT))
        WHstat(vcf = OUTPUT, whpath = WHpath, out = OUTstat)
        OUTPUTpass <- gsub(".WH.vcf.gz$", ".WHPASS.vcf.gz", OUTPUT)
        OUTPUTcsqi <- gsub(".WH.vcf.gz$", ".WHcsq_UN.vcf.gz",
            OUTPUT)
        system(paste0("bcftools view -O z -o ", OUTPUTpass, " -f PASS ",
            OUTPUT))
        args <- c("--phase a --local-csq ")
        BCFcsq(vcf = OUTPUTpass, ref = reference, out = OUTPUTcsqi,
            args = args, bcftpath = BCFtpath, gff = gffpath)
        system(paste0("tabix -f -p vcf ", OUTPUTcsqi))
        OUTPUTcsqi <- gsub(".WH.vcf.gz$", ".WHcsq.vcf.gz", OUTPUT)
        args <- c("--phase a")
        BCFcsq(vcf = OUTPUTpass, ref = reference, out = OUTPUTcsqi,
            args = args, bcftpath = BCFtpath, gff = gffpath)
        system(paste0("tabix -f -p vcf ", OUTPUTcsqi))
        CSQfilter(OUTPUTcsq = OUTPUTcsqi)
    }, mc.cores = ifelse(nrow(samples) <= ceiling(0.9 * cores),
        nrow(samples), ceiling(0.9 * cores)))
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
