#!/usr/bin/env Rscript
#### NEED PACKAGES 'data.table', 'Rcpp', 'stringr', optparse' FOR THIS PROGRAM TO WORK ####
#### IF RUNNING IN UNCERTAINTY MODE, THERE NEEDS TO BE 2 DIFFERENT VCF FILES - ONE OF THEM CONTAINING 
#### If you want to print a recomrates file, you must have the "CM" field in the INFO tag of the genotypes file - SHAPEIT4 does this automaticall the tag must also be exactly named  CM=x and the field seperated by ";" ##


option_list = list(
	optparse::make_option(c("-g", "--genotypes"), 
		type="character", 
		default=NULL, 
		help="genotype file name. no default - required", 
		metavar="character"), 
	optparse::make_option(c("-l", "--genotypelikelihoods"), 
		type="character", 
		help="genotype likelihoods file name", 
		metavar="character"), 
	optparse::make_option(c("-u", "--uncertaintyMode"), 
		type="logical", 
		default="F", 
		help="Run in uncertainty mode", 
		metavar="logical"), 
	optparse::make_option(c("-o", "--output"), 
		type="character", 
		default=NULL, 
		help="Output file path. no default - required", 
		metavar="character")
    );


opt_parser = optparse::OptionParser(option_list=option_list);
opt = optparse::parse_args(opt_parser);

# opt$genotypes = "allAncients.allModerns.chr22.PHASED2.highAccuracySNPs.subset.subsetPos.vcf.gz"
##opt$output = "test.out"

if (is.null(opt$output)){
  optparse::print_help(opt_parser)
  stop("Must supply output file name\n", call.=FALSE)
}

if (is.null(opt$genotypes)){
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (input file)\n", call.=FALSE)
}

if (opt$uncertaintyMode == TRUE) {
	cat("Running in Uncertainty mode!\n")
}

if (opt$uncertaintyMode == TRUE && is.null(opt$genotypelikelihoods) == TRUE) {
	optparse::print_help(opt_parser)
	stop("Must supply genotype likelihoods if running in uncertainty mode\n", call.=FALSE)
}

chromopainteroutput = paste0(opt$output, ".chromopainter.inp")
recomapoutput = paste0(opt$output, ".recomrates.txt")
idfileoutput = paste0(opt$output, ".idfile.txt")

Rcpp::sourceCpp("/cluster/project8/hellenthal/SamMorris/program_files/functions.cpp")
cat("Successfully souced and compiled source code!\n")

if (opt$uncertaintyMode == TRUE) {
	command = paste0("zgrep -o -a -m 1 -h -n '#CHROM' ",opt$genotypelikelihoods, " | cut -d':' -f1")
	lineskip = as.numeric(system(command, intern=T)) - 1 
	start.time = Sys.time()
	cat("Reading in genotype likelihood Data\n")
	likelihoods = data.table::fread(opt$genotypelikelihoods, skip=lineskip)
	likelihoods = as.matrix(likelihoods)
	end.time = Sys.time()
	cat(paste("Read in genotype likelihoods in", round(end.time - start.time, 2), "seconds!\n"))

	if (colnames(likelihoods)[1] != "#CHROM") {
		stop("First field isn't chromosome - skipped the wrong number of lines? Exiting....")
	}

	field = grep("GP", likelihoods[1,])
	format = as.character(likelihoods[1,field])
	dsField = which(unlist(stringr::str_split(format, ":")) == "DS")
	gtField = which(unlist(stringr::str_split(format, ":")) == "GT")
	
	if (length(dsField) != 1) {
		stop("Can't find DS field in VCF! Exiting...\n", call.=FALSE)
	} 

	if (length(gtField) != 1) {
		stop("Can't find GT field in VCF! Exiting...\n", call.=FALSE)
	} 

	subsetStart = which(colnames(likelihoods) == "FORMAT")
	likelihoods = likelihoods[,-(1:subsetStart)]
}

command = paste0("zgrep -o -a -m 1 -h -n '#CHROM' ",opt$genotypes, " | cut -d':' -f1")
lineskip = as.numeric(system(command, intern=T)) - 1 

cat("Reading in genotype Data\n")
genotypes = data.table::fread(opt$genotypes, skip=lineskip)


if (colnames(genotypes)[1] != "#CHROM") {
	stop("First field isn't chromosome - skipped the wrong number of lines? Exiting....")
}

posCol = which(colnames(genotypes) == "POS")
positions = genotypes[,..posCol]
positions = positions$POS

#### make the recombination map ###

if (stringr::str_detect(genotypes$INFO[1], "CM")) {
	cat("Found CM field - generating recombination rate file\n")
	cM = as.numeric(stringr::str_split(stringr::str_split(string=genotypes$INFO, pattern=";", simplify=T)[,stringr::str_which(unlist(stringr::str_split(string=genotypes$INFO, pattern=";")[1]), "CM")], pattern="=", simplify=T)[,2])
	recomap = as.matrix(data.frame(pos = positions, cM = cM))
	rates = ReturnGenMap2(recomap)
	recomapout = as.matrix(data.table::data.table(start.pos = positions, recom.rate.perbp=rates))
	format(recomapout, scientific=F)
	data.table::fwrite(data.table::as.data.table(recomapout), recomapoutput, col.names=TRUE, row.names=FALSE, sep=" ", quote=FALSE, scipen=999)
} else {
	cat("Can't find CM field in INFO tag - not generating recombination rate file\n")
}

positions[1] = paste("P", positions[1])
positions = as.character(positions)

subsetStart = which(colnames(genotypes) == "FORMAT")
genotypes = genotypes[,-(1:subsetStart)]
genotypes = as.matrix(genotypes)


if (stringr::str_detect(as.character(genotypes[1,1]), "/") == TRUE)  {
	stop("Looks like input isn't phased... Exiting...")
}

if (opt$uncertaintyMode == TRUE) {
	cat("Combining genotypes and likelihoods\n")
	ChromoPainterOutput = ReturnChromopainterUncerainty(genotypes, likelihoods, dsField-1)
} else {
	ChromoPainterOutput = ReturnChromopainter(genotypes) 	
}

if(opt$uncertaintyMode == TRUE) {
	ChromoPainterOutput = format(round(ChromoPainterOutput, 3), nsmall = 3)
}

ChromoPainterOutput = data.table::as.data.table(ChromoPainterOutput)
colnames(ChromoPainterOutput) = positions
data.table::fwrite(ChromoPainterOutput, chromopainteroutput, sep=" ", col.names=T, row.names=F, quote=F)

command = paste0("sed -i \"1i ", ncol(ChromoPainterOutput), "\" ", chromopainteroutput)
system(command)
command = paste0("sed -i \"1i ", nrow(ChromoPainterOutput), "\" ", chromopainteroutput)
system(command)

##### makes an idfile ######

idfile = data.table::data.table(V1 = colnames(genotypes), V2 = colnames(genotypes), V3 = 1)
data.table::fwrite(idfile, idfileoutput, col.names=F, row.names=F, sep=" ", quote=F)

