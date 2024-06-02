#include "ClusterManager.h"
#include <Rcpp.h>
#include "StringUtils.h"
#include "DistanceMetric.h"
#include <stdexcept>

ClusterManager::ClusterManager(Rcpp::CharacterVector termNameColumn,
                               Rcpp::CharacterVector geneIDColumn,
                               Rcpp::NumericVector PvalueColumn)
  // Store as C++ vectors
  : _termNames(
    Rcpp::as<std::vector<std::string>>(termNameColumn)
  ),
  _geneIDstrings(
    Rcpp::as<std::vector<std::string>>(geneIDColumn)
  ),
  _Pvalues(
    Rcpp::as<std::vector<double>>(PvalueColumn)
  ),
  // Initialize distanceMatrix and clusterMap with vector info
  _nterms(_termNames.size()),
  distanceMatrix(_nterms, _termNames),
  clusterMap(&_termNames, &_geneIDstrings, &_Pvalues)
{
  // Ensure all vectors are of the same size
  if (_termNames.size() != _geneIDstrings.size() || _termNames.size() != _Pvalues.size()) {
    throw std::invalid_argument("All input columns must have the same size.");
  }
}

// For calculating distance scores using geneID similarity
void ClusterManager::calculateDistanceScores(DistanceMetric distanceMetric) {

  int totalGeneCount = StringUtils::countUniqueElements(_geneIDstrings);

  // Fill distanceMatrix and initialize clusterMap with calculated distances
  for (int i=0; i<_nterms; ++i)
  {
    std::unordered_set<std::string> term1_genes = StringUtils::splitStringToUnorderedSet(_geneIDstrings[i], ",");

    // Calculate distance between each possible pair
    for (int j=0; j<_nterms; ++j)
    {
      if (i == j) {
        distanceMatrix.setDistance(ClusterManager::SAME_TERM_DISTANCE,
                                   i, j);
      }
      else {
        std::unordered_set<std::string> term2_genes = StringUtils::splitStringToUnorderedSet(_geneIDstrings[j], ",");
        double distanceScore = distanceMetric.calculateDistance(term1_genes, term2_genes, totalGeneCount);
        distanceMatrix.setDistance(distanceScore, i, j);

        if (distanceScore >= distanceMetric.getDistanceCutoff()) {
          clusterMap.addTermPair(i, j);

          // Print distance score information
          std::cout << distanceMetric.getDistanceMethod() << " for (" << _termNames[i] << ", ";
          std::cout << _termNames[j] << "): " << distanceScore << std::endl;

          std::cout << _termNames[i] << " genes: " << _geneIDstrings[i] << ", ";
          std::cout <<  _termNames[j] << " genes: " << _geneIDstrings[j] << std::endl;
        }
      }
    }
  }
}
