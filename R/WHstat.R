
WHstat<-function(vcf=NULL,whpath='~/.local/bin/whatshap',out=NULL){

  str<-paste(
    whpath,
    'stats',
    vcf,
    '>',
    out
  )
  system(str)

}
