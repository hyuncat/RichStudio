#include <Rcpp.h>
#include "ClusterManager.h"
#include "DistanceMetric.h"
#include "MergeStrategy.h"
#include <string>

//' @export
// [[Rcpp::export]]
Rcpp::List RichCluster(std::string distanceMetric, double distanceCutoff,
                       std::string mergeStrategy, double membershipCutoff,
                       Rcpp::CharacterVector termNameColumn,
                       Rcpp::CharacterVector geneIDColumn,
                       Rcpp::NumericVector PvalueColumn) {

  ClusterManager CM(termNameColumn,geneIDColumn, PvalueColumn);
  DistanceMetric DM(distanceMetric, distanceCutoff);

  CM.calculateDistanceScores(DM);
  
  /* MergeStrategy requires the following parameters:
      std::string mergeStrategy, (ex: "DAVID")
      double mergeCutoff, (0-1)
      std::string membershipStrategy, (ex: "DAVID")
      double membershipCutoff, (0-1)
  */
  MergeStrategy MS("DAVID", 0.5, "DAVID", membershipCutoff);
  CM.filterSeeds(MS);
  Rcpp::DataFrame FilteredSeedMap = CM.exportR_SeedMap();
  CM.mergeSeeds(MS);
  return Rcpp::List::create(
    Rcpp::_["DistanceMatrix"] = CM.exportR_DistanceMatrix(),
    Rcpp::_["SeedMap"] = CM.exportR_SeedMap(),
    Rcpp::_["FilteredSeeds"] = FilteredSeedMap,
    Rcpp::_["MergedSeeds"] = CM.exportR_ClusterList()
  );
}
