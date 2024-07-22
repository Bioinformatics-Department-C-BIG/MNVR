\name{WHphase}
\alias{phase}
\title{WHphase}
\usage{
WHphase(vcf=NULL,bam=NULL,ref=NULL,out=NULL,args=NULL,whpath='~/.local/bin/whatshap')
}
\description{
Wrapper for WhatsHap  phase function. It will perform phasing of VCF files using the corresponding alignments (.bam or .cram) and reference genome (GRCh37, hg19 or GRCh38 etc. )
}
\arguments{
\item{vcf}{The unphased vcf file preferably in an indexed '.vcg.gz' format}

\item{bam}{The alignment file used to generate the variant calls in an indexed '.bam' or 'cram' format}

\item{ref}{Reference genome file in '.fasta' format. Must be accompanied by .fai index (create with samtools faidx)}

\item{out}{The output phased .vcf.gz file name or path}

\item{args}{The output phased .vcf.gz file name or path}

\item{whpath}{Path to WhatsHap executable. Defaults to \code{'~/.local/bin/whatshap'} }

}
\examples{
hello()
}
