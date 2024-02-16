#include <../include/vcf_to_chromopainter_functions.h>
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <utility>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix ReturnChromopainter(Rcpp::StringMatrix vcfGenotypes) {

  int nInds = vcfGenotypes.ncol();
  int nSnps = vcfGenotypes.nrow();
  Rcpp::NumericMatrix chromopainterOutput(nInds * 2, nSnps);
  printf("Processing VCF with %d individuals at %d SNPS \n", nInds, nSnps);
  for (int i = 0; i < nSnps; i++) {
    for (int j = 0; j < nInds; j++) {
      char a = vcfGenotypes(i, _)[j][0];
      int aa = a - '0';
      char b = vcfGenotypes(i, _)[j][2];
      int bb = b - '0';
      if (j == 0) {
        chromopainterOutput(0, i) = aa;
        chromopainterOutput(1, i) = bb;
      } else {
        chromopainterOutput((j * 2), i) = aa;
        chromopainterOutput(((j * 2) + 1), i) = bb;
      }
    }
  }

  return chromopainterOutput;
}

// dont need this any more but keep just in case
// [[Rcpp::export]]
double ReturnUncertainty(String vcfield, int DSfield, int aa, int bb) {

  std::string stdfield = vcfield;
  std::vector<std::string> tokensMain;
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
Rcpp::NumericVector GetHaploidDosage(std::vector<double> GPvector,
                                     std::string MaxPhased) {
  // haploid dosages according to Ringbauer
  double p00, p01, p10, p11;
  double a1, a2;
	
  p00 = GPvector[0];
  p11 = GPvector[2];
  auto maxIt = std::max_element(GPvector.begin(), GPvector.end());
  double maxValue = *maxIt;
  auto maxIndex = std::distance(GPvector.begin(), maxIt);
  if (maxIndex == 1) { // i.e. if max likelihood unphased geno is het
    if (MaxPhased == "0|1") {
      p01 = GPvector[1];
      p10 = 0;
    } else {
      p01 = 0;
      p10 = GPvector[1];
    } // Closing bracket for the else block
  } else {
    p01 = GPvector[1] / 2;
    p10 = GPvector[1] / 2;
  }
  a1 = p11 + p10;
  a2 = p11 + p01;
	// sometimes there is weird case where the sum of GPvector is greater than 1. 
	// simple and dumb fix is just to reduce either a1 or a2 to 1 if they exceed 1
	if (a1 > 1) a1 = 1;
  if (a2 > 1) a2 = 1;
	// similarly if they are less than zero, just make them zero
	// not sure why this would ever happen but just in case
	if (a1 < 0) a1 = 0;	
	if (a2 < 0) a2 = 0;
  return Rcpp::NumericVector::create(a1, a2);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix
ReturnChromopainterUncerainty(Rcpp::StringMatrix genotypes,
                              Rcpp::StringMatrix likelihoods, int GPfield, char seperator) {
  int nInds = genotypes.ncol();
  int nSnps = genotypes.nrow();

  Rcpp::NumericMatrix chromopainterOutput(nInds * 2, nSnps);
  printf("Processing VCF with %d individuals at %d SNPS \n", nInds, nSnps);
  printf("Running in uncertainty mode \n");

  for (int i = 0; i < nSnps; i++) {
    // printf("Processing SNP %d\n", i);
    for (int j = 0; j < nInds; j++) {
      std::string PhasedGP = Rcpp::as<std::string>(genotypes(i, j));
      std::string GPstring = split_string_n(likelihoods(i, j), ':', GPfield);
      std::vector<double> GPvector = split_string_to_vector(GPstring, seperator);
      double max_GP = *max_element(GPvector.begin(), GPvector.end());
      auto dosage = GetHaploidDosage(GPvector, PhasedGP);
      chromopainterOutput((j * 2), i) = dosage[0];
      chromopainterOutput(((j * 2) + 1), i) = dosage[1];
    }
  }
  return chromopainterOutput;
}

// [[Rcpp::export]]
Rcpp::NumericVector ReturnGenMap2(Rcpp::NumericMatrix recomap) {
  Rcpp::NumericVector GenMap(recomap.nrow());
  for (int i = 0; i < recomap.nrow() - 1; i++) {
    long double physicaldist = recomap(i + 1, 0) - recomap(i, 0);
    long double gendist = recomap(i + 1, 1) - recomap(i, 1);
    long double cmdist = (gendist / physicaldist) / 100;
    GenMap(i) = cmdist;
  }
  return GenMap;
}

// also dont need this any more
// [[Rcpp::export]]
double ReturnUncertainty2(Rcpp::String vcfield, int DSfield, int aa, int bb,
                          int gpField) {
  std::string stdfield = vcfield;
  std::stringstream checkMain(stdfield);
  double dose = 0;
  int i = 0;
  for (std::string intermediate; std::getline(checkMain, intermediate, ':');
       i++) {
    if (i == DSfield) {
      dose = std::stod(intermediate);
      break;
    }
  }
  return std::fabs(aa + bb - dose);
}

// [[Rcpp::export]]
std::string split_string_n(Rcpp::String original1, char separator,
                           int gpField) {
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
std::vector<double> split_string_to_vector(Rcpp::String original1,
                                           char separator) {
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

  std::vector<double> results_double(results.size());
  std::transform(results.begin(), results.end(), results_double.begin(),
                 [](const std::string &val) { return std::stod(val); });

  return results_double;
}

// dont need either of these either
// // [[Rcpp::export]]
// double getMaxGP(const std::vector<std::string>& stringVector) {
// 	std::vector<double> doubleVector(stringVector.size());
// 	std::transform(stringVector.begin(), stringVector.end(),
// doubleVector.begin(), [](const std::string& val) {
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
//         start_point = std::from_chars(start_point, end_point, candidates[0],
//         format);
//     case 1:
//         start_point = std::from_chars(start_point, end_point, candidates[1],
//         format);
//     case 2:
//         start_point = std::from_chars(start_point, end_point, candidates[2],
//         format);
//     }

//     double ad = a;
//     double bd = b;
//     return std::abs(ad + bd - candidates[field]);
// }
