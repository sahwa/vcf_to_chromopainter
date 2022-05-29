#include <Rcpp.h>
#include <string>
#include <sstream>
#include <cmath>
#include <../include/vcf_to_chromopainter_functions.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix ReturnChromopainter(Rcpp::StringMatrix vcfGenotypes) {

	int nInds = vcfGenotypes.ncol();
	int nSnps = vcfGenotypes.nrow();
	Rcpp::NumericMatrix chromopainterOutput(nInds*2, nSnps);
	printf("Processing VCF with %d individuals at %d SNPS \n", nInds, nSnps);
	for (int i=0; i<nSnps; i++) {
		for (int j=0; j<nInds; j++) {
			char a = vcfGenotypes(i,_)[j][0];
			int aa = a - '0';
			char b = vcfGenotypes(i,_)[j][2];
			int bb = b - '0';
			if (j==0) {
				chromopainterOutput(0, i) = aa;
				chromopainterOutput(1, i) = bb;
			} else { 
				chromopainterOutput((j*2), i) = aa;
				chromopainterOutput(((j*2)+1), i) = bb;
			}
		}
	}	

	return chromopainterOutput;
}

// [[Rcpp::export]]
double ReturnUncertainty(String vcfield, int DSfield, int aa, int bb) {

	std::string stdfield = vcfield; 
	std::vector <std::string> tokensMain;
	std::stringstream checkMain(stdfield);
	std::string intermediate;  

	while (getline(checkMain, intermediate, ':')) { 
		tokensMain.push_back(intermediate); 
	} 
	
	std::string DS = tokensMain[DSfield];
	double DSf;
	std::istringstream(DS) >> DSf;

	double aaD = static_cast<double>(aa);
	double bbD = static_cast<double>(bb);

	double genoSum = aaD + bbD;
	return abs(genoSum - DSf);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix ReturnChromopainterUncerainty(Rcpp::StringMatrix genotypes, Rcpp::StringMatrix likelihoods, int GPfield) {

	int nInds = genotypes.ncol();
	int nSnps = genotypes.nrow();

	Rcpp::NumericMatrix chromopainterOutput(nInds*2, nSnps);
	printf("Processing VCF with %d individuals at %d SNPS \n", nInds, nSnps);
	printf("Running in uncertainty mode \n");

	for (int i=0; i<nSnps; i++) {
		// printf("Processing SNP %d\n", i);
		for (int j=0; j<nInds; j++) {
			char a = genotypes(i,_)[j][0];
			int aa = a - '0';
			char b = genotypes(i,_)[j][2];
			int bb = b - '0';
			int genosum = aa + bb;
					   
			std::string GPstring = split_string_n(likelihoods(i,j), ':', GPfield);
			std::vector<double> GPvector = split_string_to_vector(GPstring, ';');
			double max_GP = *max_element(GPvector.begin(), GPvector.end());

			double uncertainty = 1 - max_GP;

			double genoUncertaintyA;
			if (aa == 0) {
				genoUncertaintyA = 0 + uncertainty;
			} else {
				genoUncertaintyA = 1 - uncertainty;
			} 
			
			double genoUncertaintyB;
			if (bb == 0) {
				genoUncertaintyB = 0 + uncertainty;
			} else {
				genoUncertaintyB = 1 - uncertainty;
			}        
			
			if (j==0) {
				chromopainterOutput(0, i) = genoUncertaintyA;
				chromopainterOutput(1, i) = genoUncertaintyB;
			} else { 
				chromopainterOutput((j*2), i) = genoUncertaintyA;
				chromopainterOutput(((j*2)+1), i) = genoUncertaintyB;
			}
		}
	}

	return chromopainterOutput;
}

// [[Rcpp::export]]
Rcpp::NumericVector ReturnGenMap2(Rcpp::NumericMatrix recomap) {
	Rcpp::NumericVector GenMap(recomap.nrow());
	for (int i=0; i<recomap.nrow()-1; i++) {
		long double physicaldist = recomap(i+1,0) - recomap(i,0);
		long double gendist = recomap(i+1,1) - recomap(i,1);
		long double cmdist = (gendist / physicaldist) / 100;
		GenMap(i) = cmdist;
	}
	return GenMap;
}

// [[Rcpp::export]]
double ReturnUncertainty2(Rcpp::String vcfield, int DSfield, int aa, int bb, int gpField) {
	std::string stdfield = vcfield;
	std::stringstream checkMain(stdfield);
	double dose = 0;
	int i = 0;
	for (std::string intermediate; std::getline(checkMain, intermediate, ':'); i++) {
		if (i == DSfield) {
			dose = std::stod(intermediate);
			break;
		}
	}
	return std::fabs(aa + bb - dose);
}

// [[Rcpp::export]]
std::string split_string_n(Rcpp::String original1, char separator, int gpField) {	
	std::string original = original1;
	std::vector<std::string> results;
	std::string::const_iterator start = original.begin();
	std::string::const_iterator end = original.end();
	std::string::const_iterator next = std::find(start, end, separator);
	while (next != end) {
		results.push_back(std::string(start, next));
		start = next + 1;
		next = std::find(start, end, separator);
	}
	results.push_back(std::string(start, next));
	
	return results[gpField];
}

// [[Rcpp::export]]
std::vector<double> split_string_to_vector(Rcpp::String original1, char separator) {
	
	std::string original = original1;
	std::vector<std::string> results;
	std::string::const_iterator start = original.begin();
	std::string::const_iterator end = original.end();
	std::string::const_iterator next = std::find( start, end, separator );
	while (next != end) {
		results.push_back(std::string(start, next));
		start = next + 1;
		next = std::find(start, end, separator);
	}
	results.push_back(std::string(start, next));

	std::vector<double> results_double(results.size());
	std::transform(results.begin(), results.end(), results_double.begin(), [](const std::string& val) {
		return std::stod(val);
	});
	
	return results_double;
}


// // [[Rcpp::export]]
// double getMaxGP(const std::vector<std::string>& stringVector) {
// 	std::vector<double> doubleVector(stringVector.size());
// 	std::transform(stringVector.begin(), stringVector.end(), doubleVector.begin(), [](const std::string& val) {
//         return stod(val);
//     });
//     double max = *max_element(doubleVector.begin(), doubleVector.end());
// 	return max;
// }

// double strparse(string_view s, int field, int a, int b) {
//     double candidates[3];
//     auto format = std::chars_format::fixed;
//     auto start_point = s.data();
//     auto end_point = s.data() + s.size();
//     switch (field) {
//     case 0:
//         start_point = std::from_chars(start_point, end_point, candidates[0], format);
//     case 1:
//         start_point = std::from_chars(start_point, end_point, candidates[1], format);
//     case 2:
//         start_point = std::from_chars(start_point, end_point, candidates[2], format);
//     }

//     double ad = a;
//     double bd = b;
//     return std::abs(ad + bd - candidates[field]);
// }
