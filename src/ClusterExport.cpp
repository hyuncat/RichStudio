#include <Rcpp.h>
#include "ClusterManager.h"
#include "DistanceMetric.h"
#include <string>

//' @export
// [[Rcpp::export]]
Rcpp::List RichCluster(std::string distanceMetric, double distanceCutoff,
                       Rcpp::CharacterVector termNameColumn,
                       Rcpp::CharacterVector geneIDColumn,
                       Rcpp::NumericVector PvalueColumn) {

  ClusterManager CM(termNameColumn,geneIDColumn, PvalueColumn);
  DistanceMetric DM(distanceMetric, distanceCutoff);

  CM.calculateDistanceScores(DM);
  return Rcpp::List::create(
    Rcpp::_["DistanceMatrix"] = CM.exportR_DistanceMatrix(),
    Rcpp::_["ClusterMap"] = CM.exportR_ClusterMap()
  );
}
