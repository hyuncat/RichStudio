#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <algorithm>
#include <utility>
#include <functional>

#include <unordered_map>
#include <vector>


// A hash function used to hash a pair of any kind
struct hash_pair {
  template <class T1, class T2>
  size_t operator()(const std::pair<T1, T2>& p) const
  {
    std::hash<T1> hash1;
        std::hash<T2> hash2;

        auto hashValue1 = hash1(p.first);
        auto hashValue2 = hash2(p.second);

        if (hashValue1 != hashValue2) {
            return hashValue1 ^ hashValue2;
        }

        // If hashValue1 == hashValue2, their XOR is zero.
        return hashValue1;
  }
};

// Rcpp function to split a comma-separated string into a vector
std::vector<std::string> splitString(const std::string& input, char delimiter) {
  std::vector<std::string> result;
  std::istringstream tokenStream(input);
  std::string token;
  while (std::getline(tokenStream, token, delimiter)) {
    result.push_back(token);
  }
  return result;
}

// For counting the # of unique GeneIDs in an vector (like an R data.frame column)
double countUniqueGenes(std::vector<std::string> geneVec) {
  std::unordered_set<std::string> uniqueGenes;
  std::vector<std::string> currentGenes;

  // For each comma-delimited geneID cell
  for (std::string& gene : geneVec) {
    currentGenes = splitString(gene, ','); // Split
    // For each individual geneID
    for (std::string& g : currentGenes) {
      uniqueGenes.insert(g); // if not present in uniqueGenes, add
    }
  }
  return static_cast<double>(uniqueGenes.size()); // Return number of unique genes
}


// Function to calculate kappa score between two terms
// t1_genes: Vector containing term1's geneIDs
// t2_genes: Vector containing term2's geneIDs
double calKappa(const std::vector<std::string> t1_genes, const std::vector<std::string> t2_genes, double totalGeneCount) {

  // Return 0 if no overlapping genes
  if (std::find_first_of(t1_genes.begin(), t1_genes.end(), t2_genes.begin(), t2_genes.end()) == t1_genes.end()) {
    return 0.0;
  }

  // Else, calculate kappa of overlap
  double common, t1_only, t2_only = 0;

  for (const std::string& gene : t1_genes) {
    // If t1_gene found in t2_genes, increment count of common genes
    if (std::find(t2_genes.begin(), t2_genes.end(), gene) != t2_genes.end()) {
      common++;
    } else {  // Else, gene unique to t1
      t1_only++;
    }
  }
  for (const std::string& gene : t2_genes) {
    if (std::find(t2_genes.begin(), t2_genes.end(), gene) == t2_genes.end()) {
      t2_only++;
    }
  }

  //t2_only = t2_genes.size() - common; // Genes unique to term2
  double unique = totalGeneCount - common - t1_only - t2_only; // Count of all genes not found in either term
  //double total = common + t1_only + t2_only + totalGeneCount;

  // Get probability of relative agreement and chance agreement
  // double oab = (tmp[0][0] + tmp[1][1]) / static_cast<double>(tmp[0][0] + tmp[0][1] + tmp[1][0] + tmp[1][1]);
  // double aab = ((tmp[0][0] + tmp[1][0]) * (tmp[0][0] + tmp[0][1]) + (tmp[0][1] + tmp[1][1]) * (tmp[1][0] + tmp[1][1])) /
  //   (static_cast<double>(tmp[0][0] + tmp[0][1] + tmp[1][0] + tmp[1][1]) * (tmp[0][0] + tmp[0][1] + tmp[1][0] + tmp[1][1]));

  double relative_agree = (common + unique)/totalGeneCount;
  double chance_yes = ((common + t1_only)/totalGeneCount) + ((common + t2_only)/totalGeneCount);
  double chance_no = ((unique + t1_only)/totalGeneCount) + ((unique + t2_only)/totalGeneCount);
  double chance_agree = chance_yes + chance_no;

  if (chance_agree == 1) {
    return 0.0; // Prevent divide by zero
  } else {
    return (relative_agree - chance_agree) / (1 - chance_agree);// Kappa
  }
}


// Helper function to add an element to the end of the vector associated with a key
void addTermToMap(std::unordered_map<std::string, std::vector<std::string>>& map, const std::string& key, const std::string& element) {
  auto it = map.find(key);
  if (it != map.end()) {
    // Key exists, add the element to the existing vector
    it->second.push_back(element);
  } else {
    // Key does not exist, create a new vector with the element and insert into the map
    std::vector<std::string> newVector = {element};
    map.insert({key, newVector});
  }
}


// Temporary functions to test functionality of code
// Function to convert std::unordered_map<std::string, std::vector<std::string>> to DataFrame
DataFrame mapToDataFrame(const std::unordered_map<std::string, std::vector<std::string>>& map) {
  std::vector<std::string> keys;
  std::vector<std::vector<std::string>> values;

  for (const auto& entry : map) {
    keys.push_back(entry.first);
    values.push_back(entry.second);
  }

  return DataFrame::create(_["Key"] = keys, _["Values"] = values);
}

// Function to convert std::unordered_map<std::pair<std::string, std::string>, double> to DataFrame
DataFrame pairMapToDataFrame(const std::unordered_map<std::pair<std::string, std::string>, double, hash_pair>& map) {
  std::vector<std::string> keys1;
  std::vector<std::string> keys2;
  std::vector<double> values;

  for (const auto& entry : map) {
    keys1.push_back(entry.first.first);
    keys2.push_back(entry.first.second);
    values.push_back(entry.second);
  }

  return DataFrame::create(_["Key1"] = keys1, _["Key2"] = keys2, _["Values"] = values);
}


//' @export
// [[Rcpp::export]]
List RcppKappaCluster(CharacterVector myTerms, CharacterVector myGenes,
                      double kappaCutoff = 0.5) {
  std::unordered_map<std::pair<std::string, std::string>, double, hash_pair> kappaScores;
  std::unordered_map<std::string, std::vector<std::string>> signifKappaTerms;

  std::vector<std::string> term_vec = Rcpp::as<std::vector<std::string>>(myTerms);
  std::vector<std::string> gene_vec = Rcpp::as<std::vector<std::string>>(myGenes);

  // Calculate kappa score between each pair of terms
  double kappa;
  double totalGeneCount = countUniqueGenes(gene_vec);
  std::vector<std::string> term1_genes, term2_genes;

  for (auto term1 = term_vec.begin(); term1 != term_vec.end(); ++term1) {
    term1_genes = splitString(*term1, ',');
    for (auto term2 = term1 + 1; term2 != term_vec.end(); ++term2) {
      term2_genes = splitString(*term2, ',');
      kappa = calKappa(term1_genes, term2_genes, totalGeneCount);
      kappaScores[std::make_pair(*term1, *term2)] = kappa; // Update kappaScores
      // Add term pairs with kappa > kappaCutoff to signifKappaTerms
      if (kappa >= kappaCutoff) {
        addTermToMap(signifKappaTerms, *term1, *term2);
      }
    }
  }

  // Convert std::unordered_map to DataFrame
  DataFrame dfSignifKappaTerms = mapToDataFrame(signifKappaTerms);
  DataFrame dfKappaScores = pairMapToDataFrame(kappaScores);

  // Return a List containing DataFrames
  return List::create(_["SignifKappaTerms"] = dfSignifKappaTerms, _["KappaScores"] = dfKappaScores);
}
