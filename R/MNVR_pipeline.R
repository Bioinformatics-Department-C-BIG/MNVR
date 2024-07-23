MNVR_pipeline<-function(vcfdir=NULL,bamdir=NULL,whpath='~/.local/bin/whatshap',bcftpath='/usr/bin/bcftools',
                        outpath=NULL,hgnc=NULL,gffpath=NULL,reference=NULL,keep_temp=T){




  if(!dir.exists(outpath)){
    dir.create(outpath)
  }

  indir<-gsub('\\/\\/','/',vcfdir)
  bamdir<-gsub('\\/\\/','/',bamdir)

  sni<-str_match(indir,'.+/(.+)\\..+vcf.gz$')[,2]
  bmi<-str_match(bamdir,'.+/(.+)\\..+bam$')[,2]

  samples<-data.frame(Sample=sni,VCF=indir,BAM=bamdir[match(bmi,sni)])
if(Sys.info()[[1]]=='Linux'){
  print('Detected Linux System. Proceeding with multicore processing.')
ncores<-ifelse(nrow(samples)<=ceiling(0.9*detectCores()),nrow(samples),ceiling(0.9*detectCores()))

}else{ncores=1

}
  print(paste0('N cores: ',ncores))
  print(samples)



ALLRES<<-mclapply(c(1:nrow(samples)),function(s,refp=reference,whp=whpath,bcftp=bcftpath,gffp=gffpath,outp=outpath){

    VCF<-samples$VCF[s]
    BAM<-samples$BAM[s]

    OUTPUT<-gsub('.vcf.gz$','.WH.vcf.gz',VCF)
    # vcfn<-VCF
    vcfn<-gsub('.+/','',VCF)

    OUTPUT<-paste0(outp,gsub('.vcf.gz$','.WH.vcf.gz',vcfn))
    OUTtxt<-paste0(outp,gsub('.vcf.gz$','.WHcGT.txt',vcfn))
    OUTstat<-paste0(outp,gsub('.vcf.gz$','.WHstats.txt',vcfn))

    args<-c('--distrust-genotypes --include-homozygous',paste('--changed-genotype-list',OUTtxt))
    #
    WHphase(vcf=VCF,bam = BAM,ref=refp, out=OUTPUT,args=args,whpath = whp)
    #
    system(paste0('tabix -f -p vcf ',OUTPUT))
    #
    WHstat(vcf = OUTPUT,whpath = whp,out=OUTstat)

    OUTPUTpass<-gsub('.WH.vcf.gz$','.WHPASS.vcf.gz',OUTPUT)

    OUTPUTcsqi<-gsub('.WH.vcf.gz$','.WHcsq_UN.vcf.gz',OUTPUT)

    system(paste0("bcftools view -O z -o ",OUTPUTpass," -f PASS ",OUTPUT))

    args<-c('--phase a --local-csq ')

    BCFcsq(vcf = OUTPUTpass,ref = refp,out = OUTPUTcsqi,args = args,bcftpath = bcftp,gff = gffp)
    system(paste0('tabix -f -p vcf ',OUTPUTcsqi))

    OUTPUTcsqi<-gsub('.WH.vcf.gz$','.WHcsq.vcf.gz',OUTPUT)
    args<-c('--phase a')

    BCFcsq(vcf = OUTPUTpass,ref = refp,out = OUTPUTcsqi,args = args,bcftpath = bcftp,gff = gffp)

    system(paste0('tabix -f -p vcf ',OUTPUTcsqi))

    res<-CSQfilter(OUTPUTcsq = OUTPUTcsqi,HGNC=hgnc,OUTpath = outp,ref = refp,vcf_format = T)


   return(res)
  },mc.cores = ncores)

Allres<-as.data.frame(rbindlist(ALLRES))


if(!keep_temp){

  flist<-dir(outpath,full.names = T)
  flist<-flist[!grepl('.tsv$|MNV.vcf.gz$|.pdf$|WHstats.txt$',flist)]
  file.remove(flist)
}

return(Allres)


}
