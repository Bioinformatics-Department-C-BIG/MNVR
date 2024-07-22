WHphase<-function(vcf=NULL,bam=NULL,ref=NULL,out=NULL,args=NULL,whpath='~/.local/bin/whatshap'){

str<-paste(
  whpath,
  'phase',
  paste(args,collapse = ' '),
  '-r',
  ref,
  '-o',
  out,
  vcf,
  bam
)

system(str)

}
