
WHstat<-function(vcf=NULL,whpath='~/.local/bin/whatshap',out=NULL){

  str<-paste(
    whpath,
    'stats',
    vcf,
    '>',
    out
  )
  system(str)

  whs<-read.delim(out,sep = "")
  cin<-which(whs$statistics=='Chromosome'|whs$statistics=='ALL')
  chr<-whs$for.[cin]
  chr[length(chr)]<-'All'
  cc<-0
  sdf<-data.frame()
  for(c in cin){
    cc<-cc+1
    wh<-whs[c:(c+20),]
    # wh<-trimws((wh))
    # wh<-wh[!is.na(wh)]

    td<-data.frame(
      Sample=gsub('^X','',colnames(wh)[5]),
      Chromosome=wh$for.[1],
      Variants=as.numeric(wh$sample[2]),
      Heterozygous=as.numeric(wh$statistics[3]),
      Phased=as.numeric(wh$statistics[4]),
      Het_Phased_Perc=as.numeric(gsub('%','',wh$for.[5])),
      Blocks=as.numeric(wh$statistics[8]),
      Average_Block_size=as.numeric(wh$sample[12]),
      Largest_Block=as.numeric(wh$for.[13]),
      Smallest_Block=as.numeric(wh$for.[14]),
      Average_Block_Length=as.numeric(wh$sample[18]),
      Longest_Block_Length=as.numeric(wh$for.[19]),
      Shortest_Block_Length=as.numeric(wh$for.[20])
      )
    sdf<-rbind.data.frame(sdf,td)
    }
  return(sdf)
}
