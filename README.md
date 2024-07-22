# Multi-Nucleotide Analysis in R (MNVR)

Wrapper functions and pipeline for phasing and MNV identification, consequence annotation and analaysis from VCF/BAM files

## Getting Started

You can install this package in R using the following devtools command:

```         
require('devtools')
install_github("Bioinformatics-Department-C-BIG/MNVR")
```

## External dependencies

This package requires WhatsHap and bcftools already installed on a linux system.

### Linux

WhatsHap can be installed either using pip

```         
pip3 install --user whatshap
```

or conda

```         
conda create -n whatshap-env whatshap
```

bcftools can be installed using:

```         

sudo apt install bcftools
```

or

build from source

<https://samtools.github.io/bcftools/>


## Authors

Marios Tomazou, PhD CING Bioinformatics Department (mariost@cing.ac.cy)
Konstantinos Panagiotou, BSc CING Bioinformatics Department (konstantinosp@cing.ac.cy)
Athina Theodosiou, PhD CING Cytogenetics Department (athinat@cing.ac.cy)
Styliana Menelaou, MSc CING Bioinformatics Department (stylianam@cing.ac.cy)
