
WHphase<-function(vcf=NULL,bam=NULL,ref=NULL,out=NULL,args=NULL,whpath=NULL){

str<-paste(
  WHpath,
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
