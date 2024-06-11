#ifndef CLUSTERLIST_H
#define CLUSTERLIST_H

#include "SeedMap.h"
#include "MergeStrategy.h"
#include "DistanceMatrix.h"
#include <list>
#include <unordered_set>

class ClusterList {
public:
  // Constructor initializes and filters at same time
  // Maybe in future, can pass in R dataframe (maybe method of seedMap tho)
  ClusterList(SeedMap& seedMap, DistanceMatrix& DM)
    : _seedMap(seedMap), 
      _distanceMatrix(DM),
      _distanceFunction([&DM](int a, int b) { return DM.getDistance(a, b); }
  ) {};
  
  void filterSeeds(MergeStrategy MS);
  // void mergeSeeds(MergeStrategy MS, const DistanceMatrix& DM);
  
  struct MergePartner {
    int clusterNumber;
    double mergeScore;
  };
  
  MergePartner getBestMergePartner(const MergePartner& partner);
  Rcpp::DataFrame export_RDataFrame(); // R conversion util
  
private:
  std::list<std::unordered_set<int>> clusterList;
  SeedMap& _seedMap;
  DistanceMatrix& _distanceMatrix;
  std::function<double(int, int)> _distanceFunction;
};

#endif