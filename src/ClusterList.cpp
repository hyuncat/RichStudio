#include "ClusterList.h"
#include "MergeStrategy.h"
#include "DistanceMatrix.h"
#include "StringUtils.h"
#include <functional>




// void ClusterMap::mergeSeeds(MembershipStrategy MS, const DistanceMatrix& DM) 
// {
//   auto distanceFunction = [&DM](int a, int b) { return DM.getDistance(a, b); };
//   for (const auto& cluster1 : clusterList) {
//     for (const auto& cluster2 : clusterList) 
//     {
//       std::unordered_set<int> clusterGroup = cluster1;
//       clusterGroup.insert(cluster2.begin(), cluster2.end());
//       
//       double membership = MS.calculateMembership(clusterGroup, distanceFunction);
//       if (membership >= MS.getMembershipCutoff()) {
//         // merge the two clusters
//       }
//     }
//     
//   }
// }

void ClusterList::filterSeeds(MergeStrategy MS) {
  const auto& seedMapInstance = _seedMap.getSeedMap();
  clusterList.clear();
  
  for (const auto& termPair : seedMapInstance) {
    std::unordered_set<int> clusterGroup;
    clusterGroup.insert(termPair.first);
    clusterGroup.insert(termPair.second.begin(), termPair.second.end());
    
    double membership = MS.calculateMembership(clusterGroup, _distanceFunction);
    if (membership >= MS.getMembershipCutoff()) {
      clusterList.push_back(clusterGroup);
    }
  }
}


// Export ClusterList as a R dataframe with termNames and index values
Rcpp::DataFrame ClusterList::export_RDataFrame() {
  std::vector<std::string> clusterIndicesColumn;
  std::vector<std::string> clusterColumn;
  
  std::vector<std::string>& termNames = *(_seedMap._termNames); // dereference pointer
  
  // Iterate through ClusterList and convert to vectors
  for (const auto& clusterGroup : clusterList) {
    // Turn the ints into a string
    std::string clusterIndicesString = StringUtils::unorderedSetToString(clusterGroup, ", ");
    clusterIndicesColumn.push_back(clusterIndicesString);
    
    std::vector<std::string> clusterGroupTerms;
    for (const auto& term_index : clusterGroup) {
      std::string term = termNames[term_index];
      clusterGroupTerms.push_back(term);
    }
    // Convert vector/unordered_set to one comma-delimited string
    std::string clusterGroupTerms_string = StringUtils::vectorToString(clusterGroupTerms, ", ");
    clusterColumn.push_back(clusterGroupTerms_string);
  }
  
  // Create and return a DataFrame using Rcpp
  return Rcpp::DataFrame::create(Rcpp::Named("Cluster") = clusterColumn,
                                 Rcpp::Named("ClusterIndices") = clusterIndicesColumn);
}
