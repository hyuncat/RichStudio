#include "KappaMap.h"

// KappaMap constructor: Initializes an indexable, flattened list of kappa scores
//
// @param n: int, the size of the KappaMap (number of terms)
KappaMap::KappaMap(int n) : n_terms(n) {
  map.resize(n * n); // Vector size is n_terms^2
}

// getKappa: Get kappa score for terms with indices (i, j) from the KappaMap, where
//           i and j correspond to a particular term's index in the original term_vec
//
// @param i: int, row index
// @param j: int, column index
// @return: The kappa score at indices (i, j)
double KappaMap::getKappa(int i, int j) const {
  int kappaIndex = (i * n_terms) + j;
  return map[kappaIndex];
}

// setKappa: Set the kappa score for each term pair with indices (i, j) in the KappaMap
//
// @param kappa: double, the kappa score to be set
// @param i: int, row index
// @param j: int, column index
void KappaMap::setKappa(double kappa, int i, int j) {
  int kappaIndex = (i * n_terms) + j;
  map[kappaIndex] = kappa;

  kappaIndex = (j * n_terms) + i;
  map[kappaIndex] = kappa;
}


//' Convert distance vector to 2-dimensional DataFrame
//'
//' @param myTerms CharacterVector, R vector with names of all terms in dataset
//' @return kappa_df DataFrame, the kappa similarity matrix as an R dataframe
Rcpp::DataFrame KappaMap::KappaSimilarityMatrix(Rcpp::CharacterVector myTerms) const {

  Rcpp::NumericMatrix kappa_matrix(n_terms, n_terms);

  // Fill matrix with kappa scores
  int k = 0;
  for (int i=0; i<n_terms; ++i) {
    for (int j=0; j<n_terms; ++j) {
      if (i == j) {
        kappa_matrix(i, j) = -99.0; // dummy value when term is the same
      } else {
        kappa_matrix(i, j) = map[k];
        kappa_matrix(j, i) = map[k];
      }
      k++;
    }
  }

  // Convert to R dataframe
  Rcpp::DataFrame kappa_df = Rcpp::internal::convert_using_rfunction(kappa_matrix, "as.data.frame");

  // Set row/colnames to all term names in input dataset
  kappa_df.attr("names") = myTerms;
  kappa_df.attr("row.names") = myTerms;

  return kappa_df;
}
