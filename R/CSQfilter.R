CSQfilter<-function(OUTPUTcsq=NULL,HGNC=NULL,outpath=NULL,vcf_format=T,ref=NULL){



  bcf <- read.table(OUTPUTcsq)
  colnames(bcf)[-ncol(bcf)]<-c('CHR','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')

  bcf <- add_column(bcf, transcript_num = str_count(bcf$INFO, pattern = ',') + 1, .after = "INFO")

  # distance_threshold <- 2
  sn<-str_match(OUTPUTcsq,'.+/(.+)\\..+vcf.gz$')[,2]

  bcf_PASS <- bcf[bcf$FILTER =='PASS',]
  bcf_PASS$INFO<-gsub('BCSQ=','',bcf_PASS$INFO)
  mnvind<-grepl('.>.\\+.', bcf_PASS$INFO)
  bcf_mnvs_PASS <- bcf_PASS[mnvind,]

print('here')



  bcf_mnvs_tidy<-data.frame()
  for(i in 1:nrow(bcf_mnvs_PASS)){
    temp<-bcf_mnvs_PASS[i,]
    if(temp$transcript_num>1){
      temp_e<-temp
      temp_e[c(1:temp$transcript_num),]<-temp[1,]
      temp_e$INFO<-as.vector(str_split(temp$INFO,',',simplify = T))
      temp<-temp_e
    }

    csq<-data.frame(str_split(temp$INFO,'\\|',simplify = T))
    names(csq)<-c("consequence","gene","transcript","biotype","strand","aa_change","dna_change")
    temp<-cbind.data.frame(temp,csq)
    bcf_mnvs_tidy<-rbind.data.frame(bcf_mnvs_tidy,temp)

  }

  # indm<-which(diff(as.numeric(as.factor(bcf[,1])))==0 & diff(as.numeric(bcf$POS)) <= distance_threshold)
  # ind<-c(indm,(indm+1))
  ind<-grep('\\+',bcf_mnvs_tidy$dna_change)
  ind<-unique(ind[order(ind)])

  # Filtering

  bcf_mnvs_tidy<-bcf_mnvs_tidy[ind,]
  bcf_mnvs_tidy <-bcf_mnvs_tidy[grepl('missense$', bcf_mnvs_tidy$consequence) | grepl('^synonymous$', bcf_mnvs_tidy$consequence) | grepl('^stop_gained$', bcf_mnvs_tidy$consequence),]
  bcf_mnvs_tidy$consequence[ bcf_mnvs_tidy$consequence=='*missense']<-'missense'

  if(!is.null(HGNC)){
    hgnc<-read.delim(HGNC,header = T)
    hgnc$ENselect<-str_match(hgnc$mane_select,'(ENST.+)\\|.+')[,2]
    hgnc$ENselect <- gsub('\\..+','',hgnc$ENselect)
    bcf_mnvs_tidy<-bcf_mnvs_tidy[which(bcf_mnvs_tidy$transcript%in% hgnc$ENselect),]
  }

  bcf_mnvs_tidy <- na.omit(bcf_mnvs_tidy)


  temp<-cbind.data.frame(bcf_mnvs_tidy[,c(1,2,4,5,11)],as.data.frame(str_match(bcf_mnvs_tidy$dna_change,'(\\d+)(\\D+)\\>(\\D+)\\+(\\d+)(\\D+)\\>(\\D+)')[,-1]))
  temp$match<-temp$V1!=temp$POS
  temp$POS2<-as.numeric(ifelse(temp$match,temp$V1,temp$V4))
  temp$REF2<-ifelse(temp$match,temp$V2,temp$V5)
  temp$ALT2<-ifelse(temp$match,temp$V3,temp$V6)
  temp$greater<-temp$POS<temp$POS2
  gr<-GRanges(seqnames = temp$CHR,ranges=IRanges(start = ifelse(temp$greater,temp$POS,temp$POS2),end = ifelse(temp$greater,temp$POS2,temp$POS)))
  # idx<-scanFaIndex(reference)
  dna<-scanFa(ref,gr)
  temp$REFmnv<-as.vector(as.data.frame(dna))[[1]]
  temp$DIFF<-abs(temp$POS2-temp$POS)
  seqin<-GRanges(seqnames = temp$CHR,ranges=IRanges(start = ifelse(temp$greater,temp$POS,temp$POS2)+1,end = ifelse(temp$greater,temp$POS,temp$POS2)+temp$DIFF-1))
  dnain<-scanFa(ref,seqin)

  temp$DIFFseq<-as.vector(as.data.frame(dnain))[[1]]
  temp$ALTmnv<-paste0(ifelse(temp$greater,temp$ALT,temp$ALT2),temp$DIFFseq,ifelse(temp$greater,temp$ALT2,temp$ALT))


  pair<-match(paste0(bcf_mnvs_tidy$CHR,'_',temp$POS2),paste0(bcf$CHR,'_',bcf$POS))


  bcf_mnvs_tidy <- add_column(bcf_mnvs_tidy, V11 = bcf$V10[pair], .after = "V10")
if(vcf_format){
  mnvs_vcf<-bcf_mnvs_tidy[,c(1:8,10,11)]
  mnvs_vcf$POS<-ifelse(temp$greater,temp$POS,temp$POS2)
  mnvs_vcf$REF<-temp$REFmnv
  mnvs_vcf$ALT<-temp$ALTmnv


  vcf<-read.vcfR(paste0(outpath,sn,'.WH.vcf.gz'),nrows = 1)
  cn<-colnames(vcf@gt)
  vcf@gt<-cbind(mnvs_vcf$FORMAT,mnvs_vcf$V10)
  colnames(vcf@gt)<-cn

  cn<-colnames(vcf@fix)
  vcf@fix<-trimws(as.matrix(mnvs_vcf[,c(1:8)]))
  colnames(vcf@fix)<-cn
  vcf@fix[,8]<-paste0('END=',as.numeric(vcf@fix[,2])+nchar(vcf@fix[,4])-1)

  write.vcf(vcf,file = paste0(outpath,sn,'.MNV.vcf.gz'),)
}
  bcf_mnvs_tidy<-cbind.data.frame(bcf_mnvs_tidy,temp[,c(13,14,15,17,20)])
  colnames(bcf_mnvs_tidy)<-gsub('V10','GTF',colnames(bcf_mnvs_tidy))
  colnames(bcf_mnvs_tidy)<-gsub('V11','GTF2',colnames(bcf_mnvs_tidy))

  write.table(bcf_mnvs_tidy,file = paste0(outpath,sn,'_MNV_results.tsv'),quote = F,col.names = T,row.names = F)
  # bcf_mnvs_tidy<-read.table(paste0(outpath,sn,'_MNV_results.tsv'),header = T)

  sumtab<-as.data.frame.matrix(table(bcf_mnvs_tidy$CHR,bcf_mnvs_tidy$consequence))
  sumtab<-rbind.data.frame(sumtab,(colSums(sumtab)))
  rownames(sumtab)[nrow(sumtab)]<-'Total'

  write.table(sumtab,file = paste0(outpath,sn,'_MNV_Summary.tsv'),quote = F,col.names = T,row.names = T)

  ccod<-list('missense'='orange','synonymous'='royalblue','stop_gained'='orangered')

  pdf(file = paste0(outpath,sn,'_KP_plot.pdf'),width = 12,height = 8)

  plotKaryotype()%>%
    karyoploteR::kpPoints(chr = paste0('chr',bcf_mnvs_tidy$CHR),x=bcf_mnvs_tidy$POS,y=0.2, col=as.vector(unlist(ccod[bcf_mnvs_tidy$consequence])))%>%
    karyoploteR::kpPoints(chr = paste0('chr',bcf_mnvs_tidy$CHR),x=bcf_mnvs_tidy$POS,y=0.2, col=as.vector(unlist(ccod[bcf_mnvs_tidy$consequence])))

  dev.off()
  bcf_mnvs_tidy$CHR<-ordered(bcf_mnvs_tidy$CHR,levels=c(1:22,'X','Y','MT',paste0('chr',c(1:22,'X','Y','M'))))
  p1<-ggplot(bcf_mnvs_tidy,aes(x=consequence,fill=consequence))+geom_bar(show.legend = F)+theme_bw()
  p2<-ggplot(bcf_mnvs_tidy,aes(x=CHR,fill=consequence))+geom_bar()+theme_bw()+facet_grid(.~consequence)

  pdf(file = paste0(outpath,sn,'_Consequences_plot.pdf'),width = 12,height = 8)
  gridExtra::grid.arrange(grobs=list(p1,p2),nrow=2)
  dev.off()
  bcf_mnvs_tidy$sample<-sn
  return(bcf_mnvs_tidy)
}
