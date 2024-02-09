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
		metavar="character"),
    optparse::make_option(c("-m", "--genmap"),
    	type="character",
    	default=NULL,
    	help="plink formatted genetic map",
   		metavar="character")
    );

parse_args = function() {

	opt_parser = optparse::OptionParser(option_list=option_list);
	opt = optparse::parse_args(opt_parser);

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
}


source_cpp_f = function() {
	Rcpp::sourceCpp("vcf_to_chromopainter_functions.cpp")
	cat("Successfully souced and compiled source code!\n")
}

make_recombination_map = function(geneticmap, positions) {
	positions = sort(positions)
	data.table(start.pos = pos, recom.rate.perbp = approx(x = geneticmap$V4, y = geneticmap$V3, xout = positions)$y)
}

readVCF = function(likelihoods_file) {
	fread(cmd=paste("zgrep -v '^##'", likelihoods_file))
}

getDSField = function(formatField){
	str_which(str_split_1(formatField, ":"), "DS")
}

getGTField = function(formatField){
	str_which(str_split_1(formatField, ":"), "GT")
}

getGPField = function(formatField){
	str_which(str_split_1(formatField, ":"), "GP")
}

processGenotypeLikelihoods = function(filename) {
	## read in genotype likelihoods 
	likelihoods = readVCF(filename)
	if (colnames(likelihoods)[1] != "#CHROM") {
		stop("First field isn't chromosome - skipped the wrong number of lines? Exiting....")
	}
	field = grep("GP", likelihoods[1,])
	format = as.character(likelihoods[1,..field])
	dsField = getDSField(format)
	gtField = getGTField(format)
	gpField = getGPField(format)

	GLpositions = likelihoods$POS
	subsetStart = which(colnames(likelihoods) == "FORMAT")
	dropcols = colnames(likelihoods)[1:subsetStart]
	likelihoods[, (dropcols) := NULL]

}

processGenotypes = function(filename) {
	## read in genotype likelihoods 
	likelihoods = readVCF(filename)
	if (colnames(likelihoods)[1] != "#CHROM") {
		stop("First field isn't chromosome - skipped the wrong number of lines? Exiting....")
	}
	subsetStart = which(colnames(likelihoods) == "FORMAT")
	dropcols = colnames(likelihoods)[1:subsetStart]
	likelihoods[, (dropcols) := NULL]

}

processGenotypeLikelihoods = function(filename) {
	## read in genotype likelihoods 
	likelihoods = readGLs(filename)
	if (colnames(likelihoods)[1] != "#CHROM") {
		stop("First field isn't chromosome - skipped the wrong number of lines? Exiting....")
	}
	field = grep("GP", likelihoods[1,])
	format = as.character(likelihoods[1,..field])
	dsField = getDSField(format)
	gtField = getGTField(format)
	gpField = getGPField(format)

	GLpositions = likelihoods$POS
	subsetStart = which(colnames(likelihoods) == "FORMAT")
	dropcols = colnames(likelihoods)[1:subsetStart]
	likelihoods[, (dropcols) := NULL]

}






########################################################################################

if (opt$uncertaintyMode) {
	GLs = processGenotypeLikelihoods(opt$genotypelikelihoods)
	GTs = processGenotypes(opt$genotypes)
	CPout = ReturnChromopainterUncerainty(GTs, GLs, )
	recomap = make_recombination_map(geneticmap, positions)
}

if (!opt$uncertaintyMode) {
	GTs = processGenotypes(opt$genotypes)
	CPout = ReturnChromopainterUncerainty(GTs, GLs, )
	recomap = make_recombination_map(geneticmap, positions)
}
