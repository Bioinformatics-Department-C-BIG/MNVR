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


BCFcsq<-function(vcf=NULL,ref=NULL,out=NULL,gff=NULL,args=NULL,bcftpath='/usr/bin/bcftools'){

  str<-paste(
    bcftpath,
    'csq',
    paste(args,collapse = ' '),
    '-O z',
    '-f',
    ref,
    '-o',
    out,
    '-g',
    gff,
    vcf
  )
  print(str)
  system(str)

}



