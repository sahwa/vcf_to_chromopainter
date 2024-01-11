# vcf_to_chromopainter

This is a basic tool which can be used to convert phased genotypes in vcf format into the format used by chromopainter. It requires a few libraries, which can be easily installed using ``install.packages(c("stringr", "data.table", "Rcpp", "optparse"))``. 

Tested on macOS Catalina, Ubuntu 20+ and CentOS 7.9.2009 (Core).

This tool hasn't been tested thoroughly, so please email/message me if you use it and have any issues.

# usage

```
Usage: vcf_to_chromopainter_main.R [options]

Options:
        -g CHARACTER, --genotypes=CHARACTER
                genotype file name. no default - required

        -l CHARACTER, --genotypelikelihoods=CHARACTER
                genotype likelihoods file name

        -u LOGICAL, --uncertaintyMode=LOGICAL
                Run in uncertainty mode

        -o CHARACTER, --output=CHARACTER
                Output file path. no default - required

        -h, --help
                Show this help message and exit
```

Example to make a ChromoPainter Uncertainty input using phased genotypes and genotype likelihoods:

```
Rscript vcf_to_chromopainter_main.R \
  -g a1.solved.vcf \ ## this file should contain phased genotypes (e.g. '0|1')
  -l a1.imputed.vcf \ ## this file should contain genotype likelihoods (e.g. gl='0.01,0.98,0.01')
  -u TRUE \ 
  -o a1
  ```
This format is for an unreleased version of ChromoPainter. If you want to make the input for ChromoPainter V2 (which is very likely), then run the above script like this:

```
Rscript vcf_to_chromopainter_main.R \
  -g input.vcf \ ## this file should contain phased genotypes (e.g. '0|1')
  -u FALSE \ 
  -o output
```

And then run:

```
sed -i '1,3!s/ //g' output.chromompainter.inp # or replace this with whatever the filename is
```


