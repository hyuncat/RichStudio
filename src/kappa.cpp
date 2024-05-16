#include <Rcpp.h>

#include <iostream>
#include <algorithm>
#include <utility>
#include <functional>
#include <cmath>

#include <unordered_map>
#include <vector>
#include <set>

#include "utils.h"
#include "KappaMap.h"


//' Function to calculate kappa score between two terms
//' @param t1_genes, vector containing term1's geneIDs
//' @param t2_genes, vector containing term2's geneIDs
//' @param totalGeneCount, double representing count of all unique genes in dataset of interest
//' @return kappaScore, double
//' @name getKappa
double getKappa(const std::unordered_set<std::string>& t1_genes, const std::unordered_set<std::string>& t2_genes, double totalGeneCount) {

  // Calculate the intersection of t1_genes and t2_genes
  std::unordered_set<std::string> intersection;
  std::set_intersection(t1_genes.begin(), t1_genes.end(), t2_genes.begin(), t2_genes.end(), std::inserter(intersection, intersection.begin()));

  double common = static_cast<double>(intersection.size()); // Number of common genes

  if (common == 0) {
    return 0.0; // return 0 if no overlapping genes
  }

  double t1_only = t1_genes.size() - common; // Genes unique to t1_genes
  double t2_only = t2_genes.size() - common; // Genes unique to t2_genes

  double unique = totalGeneCount - common - t1_only - t2_only; // Count of all genes not found in either term

  double relative_observed_agree = (common + unique) / totalGeneCount;
  double chance_yes = ((common + t1_only) / totalGeneCount) * ((common + t2_only) / totalGeneCount);
  double chance_no = ((unique + t1_only) / totalGeneCount) * ((unique + t2_only) / totalGeneCount);
  double chance_agree = chance_yes + chance_no;

  if (chance_agree == 1) {
    return 0.0; // Prevent divide by zero
  } else {
    return (relative_observed_agree - chance_agree) / (1 - chance_agree); // Return kappa!
  }
}


// Returns true if term2 is present in sigTermIndices[term1]'s unordered_set
bool existsInSet(const std::unordered_map<int, std::unordered_set<int>>& sigTermIndices, int term1_idx, int term2_idx) {
  auto it = sigTermIndices.find(term1_idx);
  if (it != sigTermIndices.end()) {
    const auto& set = it->second;
    return set.find(term2_idx) != set.end(); // O(1) complexity!
  }
  return false;
}


std::unordered_map<int, std::unordered_set<int>> createInitialSeeds(const std::unordered_map<int, std::unordered_set<int>>& sigTermIndices,
                                                                    const KappaMap map, double membershipThreshold = 0.50) {
  std::unordered_map<int, std::unordered_set<int>> initialSeeds;
  std::cout << "Creating initial seeds..." << std::endl;

  int count = 0;
  for (const auto& pair : sigTermIndices) {
    int key = pair.first;
    const auto& values = pair.second;

    double passedPairs = values.size();
    double totalPairs = values.size();
    double membershipOverThresholdPercent = 0;

    for (auto it1 = values.begin(); it1 != values.end(); ++it1) {
      if (values.size() == 1) {
        std::cout << "(" << count << " | case 1): only 1 term in unordered_set" << std::endl;
        break;
      }
      int term1_idx = *it1;

      for (auto it2 = std::next(it1); it2 != values.end(); ++it2) {
        int term2_idx = *it2;
        if (existsInSet(sigTermIndices, term1_idx, term2_idx)) {
          std::cout << "(" << count << " | case 2): term pair (" << term1_idx << ", " << term2_idx << ") is above 0.50" << std::endl;
          passedPairs++;
        } else {
          std::cout << "(" << count << " | case 3): term pair (" << term1_idx << ", " << term2_idx << ") is not above 0.50" << std::endl;
        }
        totalPairs++;
      }
    }

    membershipOverThresholdPercent = passedPairs / totalPairs;
    std::cout << count << " | membershipOverThresholdPercent = " << membershipOverThresholdPercent << std::endl;

    if (membershipOverThresholdPercent >= membershipThreshold) {
      initialSeeds.insert(pair);
    }
    count++;

  }

  return initialSeeds;
}


double calculateOverlap(int currentTermIdx, int otherTermIdx, std::unordered_map<int, std::unordered_set<int>>& remainingSeeds) {
  double common = 0, total;
  double currentSeedSize = remainingSeeds[currentTermIdx].size() + 1;
  double newSeedSize = remainingSeeds[otherTermIdx].size() + 1;

  if (existsInSet(remainingSeeds, currentTermIdx, otherTermIdx)) {
    common++; // The other term itself exists
  }

  for (int term_idx : remainingSeeds[otherTermIdx]) {
    if (existsInSet(remainingSeeds, currentTermIdx, term_idx)) {
      common++;
    }
  }

  total = currentSeedSize + newSeedSize;
  return (common * 2) / total;
}

std::unordered_map<int, std::unordered_set<int>> singleMergeBestSeeds(int currentTermIdx, std::unordered_map<int, std::unordered_set<int>>& remainingSeeds) {
  double bestOverlapping = 0;
  int bestOverlapIndex = -1;

  for (const auto& pair : remainingSeeds) {
    int key = pair.first;
    if (key == currentTermIdx) continue; // Skip the current term itself

    double overlapping = calculateOverlap(currentTermIdx, key, remainingSeeds);
    if (overlapping >= 0.50 && overlapping > bestOverlapping) {
      bestOverlapping = overlapping;
      bestOverlapIndex = key;
    }
  }

  if (bestOverlapIndex != -1) {
    // Merge and erase the best mergable seed
    remainingSeeds[currentTermIdx].insert(remainingSeeds[bestOverlapIndex].begin(), remainingSeeds[bestOverlapIndex].end());
    remainingSeeds.erase(bestOverlapIndex);
  }

  return remainingSeeds;
}

std::unordered_map<int, std::unordered_set<int>> mergeSeeds(std::unordered_map<int, std::unordered_set<int>>& initialSeeds) {
  std::unordered_map<int, std::unordered_set<int>> remainingSeeds = initialSeeds;

  auto it = remainingSeeds.begin();
  while (it != remainingSeeds.end()) {
    remainingSeeds = singleMergeBestSeeds(it->first, remainingSeeds);
    ++it;
  }

  return remainingSeeds;
}





//' Create a Kappa similarity matrix
//'
//' This function computes kappa scores between pairs of biological terms
//' based on shared gene IDs and creates data frames for analysis in R.
//'
//' @param myTerms, A character vector (R) containing biological terms
//' @param myGenes, A character vector (R) where each item is a string of comma-separated geneIDs corresponding to a single term
//' @param myPvalues, A numeric vector (R) containing p-values corresponding to each term.
//' @param kappaCutoff, Numeric value specifying the cutoff threshold for kappa scores. Default is 0.5.
//' @return A list containing data frames for significant term indices, all kappa scores, and the kappa score matrix.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List KappaSimilarityMatrix(Rcpp::CharacterVector myTerms, Rcpp::CharacterVector myGenes,
                                 Rcpp::NumericVector myPvalues, double kappaCutoff = 0.5) {

   // Convert R vectors to C++
   std::vector<std::string> all_terms = Rcpp::as<std::vector<std::string>>(myTerms);
   std::vector<std::string> all_genes = Rcpp::as<std::vector<std::string>>(myGenes);
   std::vector<double> all_pvalues = Rcpp::as<std::vector<double>>(myPvalues);

   std::unordered_map<int, std::unordered_set<int>> sigTermIndices;
   std::unordered_map<int, std::unordered_set<int>> initialSeeds, mergeSeedGroups;

   // all_kappas: Flattened distance matrix
   // Index a given pair (term_i, term_j) with all_kappas[(i*j)+j]
   int n_terms = all_terms.size();
   KappaMap kappa_map(n_terms);

   // Calculate kappa score between each pair of terms
   double kappa;
   double totalGeneCount = countUniqueGenes(all_genes);
   std::unordered_set<std::string> term1_genes, term2_genes;
   int count = 0;

   // Create all_kappas distance matrix
   for (int i=0; i<n_terms; ++i) {
     term1_genes = splitString(all_genes[i], ',');
     for (int j=i; j<n_terms; ++j) {
       // kappaIndex = (i*n_terms)+j;
       if (j == i) {
         kappa_map.setKappa(-99, i, j);
       } else {
         term2_genes = splitString(all_genes[j], ',');
         kappa = getKappa(term1_genes, term2_genes, totalGeneCount);
         kappa_map.setKappa(kappa, i, j);



         if (kappa >= kappaCutoff) {
           addTermToMap(sigTermIndices, i, j);
           std::cout << "kappa for (" << all_terms[i] << ", " << all_terms[j] << "): " << kappa << std::endl;
           std::cout << all_terms[i] << " genes: " << all_genes[i] << ", " << all_terms[j] << " genes: " << all_genes[j] << std::endl;
         }
       }
     }
   }

   initialSeeds = createInitialSeeds(sigTermIndices, kappa_map);
   mergeSeedGroups = mergeSeeds(initialSeeds);

   Rcpp::DataFrame SigTermIndicesDf = mapToDataFrame(sigTermIndices);
   Rcpp::DataFrame AllKappasVecDf = kappa_map.KappaVector();
   Rcpp::DataFrame KappaSimilarityMatrixDf = kappa_map.KappaSimilarityMatrix(myTerms);

   Rcpp::DataFrame initialSeedsDf = mapToDataFrame(initialSeeds);
   Rcpp::DataFrame mergeSeedsDf = mapToDataFrame(mergeSeedGroups);


   // DataFrame dfClusters = makeClusterDf(seedResult.clusters);

   return Rcpp::List::create(
     Rcpp::_["SigTermIndices"] = SigTermIndicesDf,
     Rcpp::_["AllKappaScores_Vector"] = AllKappasVecDf,
     Rcpp::_["KappaSimilarityMatrix"] = KappaSimilarityMatrixDf,
     Rcpp::_["InitialSeeds"] = initialSeedsDf,
     Rcpp::_["MergeSeeds"] = mergeSeedsDf
   );
 }



