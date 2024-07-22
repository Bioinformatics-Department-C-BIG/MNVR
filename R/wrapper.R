indir<-dir('/home/mariost/depts/BG/malekou/out240704/variant_calling/deepvariant/',pattern = '.+t.vcf.gz$',full.names = T,recursive = T)
bamdir<-dir('/home/mariost/depts/BG/malekou/out240704/preprocessing/markduplicates/',pattern = '.+bam$',full.names = T,recursive = T)
bamdir<-bamdir[!grepl('edit|Joint',bamdir)]
indir<-indir[!grepl('edit|Joint',indir)]
library(tidyverse)
library(parallel)
cores<-detectCores()
MNVR_pipeline<-function(vcfdir=NULL,bamdir=NULL,whpath='~/.local/bin/whatshap',bcftpath='/usr/bin/bcftools',outpath=NULL,HGNC='/home/konstantinosp/tools/hgnc_complete_set.txt',gffpath=NULL,reference=NULL,hgnc=NULL){

WHpath<-whpath
BCFtpath<-bcftpath
# outpath<-'./outPACKTEST/'

if(!dir.exists(outpath)){
  dir.create(outpath)
  }


gffpath<-'/home/mariost/DBs/ensembl/Homo_sapiens.GRCh37.87.gff3.gz'
reference<-'/home/mariost/DBs/references/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HGNC reference transcripts
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hgnc<-read.delim(HGNC,header = T)
hgnc$ENselect<-str_match(hgnc$mane_select,'(ENST.+)\\|.+')[,2]
hgnc$ENselect <- gsub('\\..+','',hgnc$ENselect)




indir<-gsub('\\/\\/','/',vcfdir)
bamdir<-gsub('\\/\\/','/',bamdir)

sni<-str_match(indir,'.+/(.+)\\..+vcf.gz$')[,2]
bmi<-str_match(bamdir,'.+/(.+)\\..+bam$')[,2]
samples<-data.frame(Sample=sni,VCF=indir,BAM=bamdir[match(bmi,sni)])




mclapply(c(1:nrow(samples)),function(s){

VCF<-samples$VCF[s]
BAM<-samples$BAM[s]

OUTPUT<-gsub('.vcf.gz$','.WH.vcf.gz',VCF)
# vcfn<-VCF
 vcfn<-gsub('.+/','',VCF)

OUTPUT<-paste0(outpath,gsub('.vcf.gz$','.WH.vcf.gz',vcfn))
OUTtxt<-paste0(outpath,gsub('.vcf.gz$','.WHcGT.txt',vcfn))
OUTstat<-paste0(outpath,gsub('.vcf.gz$','.WHstats.txt',vcfn))

args<-c('--distrust-genotypes --include-homozygous',paste('--changed-genotype-list',OUTtxt))
#
WHphase(vcf=VCF,bam = BAM,ref=reference, out=OUTPUT,args=args,whpath = WHpath)
#
system(paste0('tabix -f -p vcf ',OUTPUT))
#
WHstat(vcf = OUTPUT,whpath = WHpath,out=OUTstat)

OUTPUTpass<-gsub('.WH.vcf.gz$','.WHPASS.vcf.gz',OUTPUT)

OUTPUTcsqi<-gsub('.WH.vcf.gz$','.WHcsq_UN.vcf.gz',OUTPUT)

system(paste0("bcftools view -O z -o ",OUTPUTpass," -f PASS ",OUTPUT))

args<-c('--phase a --local-csq ')

BCFcsq(vcf = OUTPUTpass,ref = reference,out = OUTPUTcsqi,args = args,bcftpath = BCFtpath,gff = gffpath)
system(paste0('tabix -f -p vcf ',OUTPUTcsqi))

OUTPUTcsqi<-gsub('.WH.vcf.gz$','.WHcsq.vcf.gz',OUTPUT)
args<-c('--phase a')

BCFcsq(vcf = OUTPUTpass,ref = reference,out = OUTPUTcsqi,args = args,bcftpath = BCFtpath,gff = gffpath)

system(paste0('tabix -f -p vcf ',OUTPUTcsqi))

CSQfilter(OUTPUTcsq = OUTPUTcsqi,HGNC='/home/konstantinosp/tools/hgnc_complete_set.txt',OUTpath = outpath)

 },mc.cores = ifelse(nrow(samples)<=ceiling(0.9*cores),nrow(samples),ceiling(0.9*cores)))
#},mc.cores = 1)

}



indir<-gsub('\\/\\/','/',indir)
bamdir<-gsub('\\/\\/','/',bamdir)

sni<-str_match(indir,'.+/(.+)\\..+vcf.gz$')[,2]
bmi<-str_match(bamdir,'.+/(.+)\\..+bam$')[,2]
samples<-data.frame(Sample=sni,VCF=indir,BAM=bamdir[match(bmi,sni)])

mclapply(c(1:nrow(samples)),function(s){

  VCF<-samples$VCF[s]
  BAM<-samples$BAM[s]

  OUTPUT<-gsub('.vcf.gz$','.WH.vcf.gz',VCF)
  # vcfn<-VCF
  vcfn<-gsub('.+/','',VCF)

OUTPUT<-paste0(outpath,gsub('.vcf.gz$','.WH.vcf.gz',vcfn))
OUTPUTcsqi<-gsub('.WH.vcf.gz$','.WHcsq.vcf.gz',OUTPUT)
CSQfilter(OUTPUTcsq = OUTPUTcsqi)
},mc.cores = ifelse(nrow(samples)<=ceiling(0.9*cores),nrow(samples),ceiling(0.9*cores)))

MNVR_pipeline(outpath = '/home/mariost/depts/BI/MNV/outPACK240704/',
              vcfdir =indir,
              bamdir =bamdir

              )
