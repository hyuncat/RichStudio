#ifndef KAPPAMAP_H
#define KAPPAMAP_H

#include <vector>
#include <Rcpp.h>


class KappaMap {
public:
  KappaMap(int n);
  double getKappa(int i, int j) const;
  void setKappa(double kappa, int i, int j);
  Rcpp::DataFrame KappaSimilarityMatrix(Rcpp::CharacterVector myTerms) const;
  Rcpp::DataFrame KappaVector() const { return Rcpp::DataFrame::create(Rcpp::_["TermIndex"] = map); }

private:
  std::vector<double> map;
  int n_terms;
};

#endif
