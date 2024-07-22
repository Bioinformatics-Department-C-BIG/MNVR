# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


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
