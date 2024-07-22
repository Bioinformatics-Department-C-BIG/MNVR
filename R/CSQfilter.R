CSQfilter<-function(OUTPUTcsq=NULL,HGNC='/home/konstantinosp/tools/hgnc_complete_set.txt',OUTpath=NULL){


  outpath<-OUTpath
  hgnc<-read.delim(HGNC,header = T)
  hgnc$ENselect<-str_match(hgnc$mane_select,'(ENST.+)\\|.+')[,2]
  hgnc$ENselect <- gsub('\\..+','',hgnc$ENselect)


  bcf <- read.table(OUTPUTcsq)
  colnames(bcf)[-ncol(bcf)]<-c('CHR','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')

  bcf <- add_column(bcf, transcript_num = str_count(bcf$INFO, pattern = ',') + 1, .after = "INFO")

  distance_threshold <- 2
  sn<-str_match(OUTPUTcsq,'.+/(.+)\\..+vcf.gz$')[,2]

  bcf_mnvs_PASS <- bcf[bcf$FILTER =='PASS',]
  bcf_mnvs_PASS$INFO<-gsub('BCSQ=','',bcf_mnvs_PASS$INFO)
  bcf_mnvs_PASS <- bcf_mnvs_PASS[grepl('.>.\\+.', bcf_mnvs_PASS$INFO),]

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

  bcf_mnvs_tidy<-bcf_mnvs_tidy[which(bcf_mnvs_tidy$transcript%in% hgnc$ENselect),]

  bcf_mnvs_tidy <- na.omit(bcf_mnvs_tidy)

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
  gridExtra::grid.arrange(grobs=list(p1,p2),ncol=2)
  dev.off()
}
