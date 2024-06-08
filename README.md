# vcf_to_chromopainter

This is a basic tool written which can be used to convert phased genotypes in vcf format into the format used by chromopainter. It requires you to have R installed and `Rscript` in your path, as well a few libraries, which can be easily installed using ``install.packages(c("stringr", "data.table", "Rcpp", "optparse"))``. 

The C++ code compiles on on macOS Catalina, Ubuntu 20+ and CentOS 7.9.2009 (Core).

This tool hasn't been tested thoroughly, so please email/message me if you use it and have any issues.

# usage

```
Usage: src/vcf_to_chromopainter_main.R [options]


Options:
	-g CHARACTER, --genotypes=CHARACTER
		genotype file name. no default - required

	-l CHARACTER, --genotypelikelihoods=CHARACTER
		genotype likelihoods file name

	-u LOGICAL, --uncertaintyMode=LOGICAL
		Run in uncertainty mode

	-o CHARACTER, --output=CHARACTER
		Output file stem. no default - required

	-m CHARACTER, --genmap=CHARACTER
		plink formatted genetic map

	-h, --help
		Show this help message and exit
```


If you want to make the input for ChromoPainter V2 (which is very likely), then run the above script like this:

```
Rscript vcf_to_chromopainter_main.R \
  -g input.vcf \ ## this file should contain phased genotypes (e.g. '0|1')
  -u FALSE \ 
  -o output \
  -m genetic_map.txt
```

You will need to download a genetic map for the relevant build and chromosome you are analysing. These can be downloaded from here `https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/`.


