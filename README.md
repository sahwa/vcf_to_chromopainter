# vcf_to_chromopainter

This is a basic tool which can be used to convert phased genotypes in vcf format into the format used by chromopainter. Unfortunately it requires a few libraries, which can be easily installed using ``install.packages(c("stringr", "data.table", "Rcpp", "optparse"))``. I've only tested it on Linux, so apologies if it does not work on mac or windows. 
