#ifndef CLUSTERMANAGER_H
#define CLUSTERMANAGER_H

#include <Rcpp.h>
#include <string>
#include "DistanceMatrix.h"
#include "ClusterMap.h"
#include "DistanceMetric.h"


// ClusterManager: For efficiently clustering a single dataset in multiple ways
class ClusterManager {
public:
  static constexpr double SAME_TERM_DISTANCE = -99;

  ClusterManager(Rcpp::CharacterVector termNameColumn,
                 Rcpp::CharacterVector geneIDColumn,
                 Rcpp::NumericVector PvalueColumn);
  void calculateDistanceScores(DistanceMetric distanceMetric);
  // void filterInitialClusters(double minClusterSize);
  // void mergeClusters(double mergeThreshold);
  Rcpp::NumericMatrix exportR_DistanceMatrix() {
    return distanceMatrix.export_RMatrix();
  };
  Rcpp::DataFrame exportR_ClusterMap() {
    return clusterMap.export_RDataFrame();
  };

private:
  std::vector<std::string> _termNames;
  std::vector<std::string> _geneIDstrings; // Each string contains all geneIDs for a single term
  std::vector<double> _Pvalues;
  int _nterms;
  DistanceMatrix distanceMatrix;
  ClusterMap clusterMap;
};

#endif
