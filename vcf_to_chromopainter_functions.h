#include <Rcpp.h>

double ReturnUncertainty2(Rcpp::String vcfield, int DSfield, int aa, int bb) ;
Rcpp::NumericVector ReturnGenMap2(Rcpp::NumericMatrix recomap);
std::string split_string_n(Rcpp::String original1, char separator, int gpField);
std::vector<double> split_string_to_vector(Rcpp::String original1, char separator);

