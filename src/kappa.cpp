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

// Rcpp function to conjoin a vector into a comma-separated string
template<typename T>
std::string vecToString(const std::vector<T>& vec, char delimiter) {
  if (vec.empty()) {
    return "";
  }
  std::ostringstream oss;
  oss << vec[0];
  for (size_t i = 1; i < vec.size(); ++i) {
    oss << delimiter << vec[i];
  }
  return oss.str();
}

// For counting the # of unique GeneIDs in an vector (like an R data.frame column)
double countUniqueGenes(std::vector<std::string> geneVec) {
  std::unordered_set<std::string> uniqueGenes;
  std::vector<std::string> currentGenes;

  // For each comma-delimited geneID cell
  int initCount = 0;
  for (std::string& gene : geneVec) {
    currentGenes = splitString(gene, ','); // Split
    initCount += currentGenes.size();
    // For each individual geneID
    for (std::string& g : currentGenes) {
      uniqueGenes.insert(g); // if not present in uniqueGenes, add
    }
  }
  std::string oldSize = "original # geneIds: "+std::to_string(initCount);
  std::string newSize = "filtered # unique geneIds: "+std::to_string(uniqueGenes.size());
  std::cout << oldSize << std::endl;
  std::cout << newSize << std::endl;
  return static_cast<double>(uniqueGenes.size()); // Return number of unique genes
}


// Function to calculate kappa score between two terms
// t1_genes: Vector containing term1's geneIDs
// t2_genes: Vector containing term2's geneIDs
//' @export
 // [[Rcpp::export]]
double calKappa(const std::vector<std::string> t1_genes, const std::vector<std::string> t2_genes, double totalGeneCount) {

  // Return 0 if no overlapping genes
  if (std::find_first_of(t1_genes.begin(), t1_genes.end(), t2_genes.begin(), t2_genes.end()) == t1_genes.end()) {
    return 0.0;
  }

  // Else, calculate kappa of overlap
  double common = 0, t1_only = 0, t2_only = 0;

  for (const std::string& gene : t1_genes) {
    // t1_gene found in t2_genes
    if (std::find(t2_genes.begin(), t2_genes.end(), gene) != t2_genes.end()) {
      common++;
    } else {  // Gene unique to t1
      t1_only++;
    }
  }
  for (const std::string& gene : t2_genes) {
    if (std::find(t2_genes.begin(), t2_genes.end(), gene) == t2_genes.end()) {
      t2_only++;
    }
  }

  double unique = totalGeneCount - common - t1_only - t2_only; // Count of all genes not found in either term

  double relative_agree = (common + unique)/totalGeneCount;
  double chance_yes = ((common + t1_only)/totalGeneCount) + ((common + t2_only)/totalGeneCount);
  double chance_no = ((unique + t1_only)/totalGeneCount) + ((unique + t2_only)/totalGeneCount);
  double chance_agree = chance_yes * chance_no;

  if (chance_agree == 1) {
    return 0.0; // Prevent divide by zero
  } else {
    return (relative_agree - chance_agree) / (1 - chance_agree); // Return kappa!
  }
}

std::unordered_map<std::string, std::vector<std::string>> termMembershipOverThreshold(
          std::string term,
          std::vector<std::string> termVec,
          std::unordered_map<std::pair<std::string, std::string>, double, hash_pair> kappaScores) {
  double totalPairs = 0, goodPairs = 0;
  for (std::vector<std::string>::iterator t=termVec.begin(); t!=termVec.end(); ++t)
  {
    std::cout<<*t<<std::endl;
  }
}

std::unordered_map<std::string, std::vector<std::string>> createInitSeeds(
                    std::unordered_map<std::string, std::vector<std::string>> signifKappaTerms,
                    std::unordered_map<std::pair<std::string, std::string>, double, hash_pair> kappaScores,
                    int initGroupMembership) {
  std::unordered_map<std::string, std::vector<std::string>> qualifiedSeeds;
  for (auto term : signifKappaTerms) {
    if (term.second.size() > initGroupMembership) {


    }
  }
  return qualifiedSeeds;
}


// Helper function to add an element to the end of the vector associated with a key
// "SignifKappaTerms" = strMap, strElement are not null
// "SignifKappaScores" = intMap, intElement are not null
template <typename T>
void addTermToMap(std::unordered_map<std::string, std::vector<T>>& map,
                  const std::string& key, const T& element) {

  // Adding term to SignifKappaTerms
  auto it = map.find(key);
  std::vector<T> keyValue;
  if (it != map.end()) {
    // Key exists, add the element to the existing vector
    keyValue = it->second;
    keyValue.push_back(element);
    it->second = keyValue;
  } else {
    // Key does not exist, create a new vector with the element and insert into the map
    keyValue = {element};
    map[key] = keyValue;
  }

}


// Temporary functions to test functionality of code
// Function to convert std::unordered_map<std::string, std::vector<std::string>> to DataFrame
template <typename T>
DataFrame mapToDataFrame(const std::unordered_map<std::string, std::vector<T>>& map=NULL,
                         bool isTerm=true) {
  std::vector<std::string> keys;
  std::vector<std::string> values;
  std::string valueString;

  if (isTerm == true) {
    for (const auto& entry : map) {
      keys.push_back(entry.first);
      valueString = vecToString(entry.second, ',');
      values.push_back(valueString);
    }
    return DataFrame::create(_["Term"] = keys, _["Signif_TermPairs"] = values);
  } else {
    for (const auto& entry : map) {
      keys.push_back(entry.first);
      valueString = vecToString(entry.second, ',');
      values.push_back(valueString);
    }
    return DataFrame::create(_["Term"] = keys, _["Pair_Kappas"] = values);
  }
}

// Function to convert std::unordered_map<std::pair<std::string, std::string>, double> to DataFrame
// "KappaScores" dataframe
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
  std::unordered_map<std::string, std::vector<double>> signifKappaScores;

  std::vector<std::string> term_vec = Rcpp::as<std::vector<std::string>>(myTerms);
  std::vector<std::string> gene_vec = Rcpp::as<std::vector<std::string>>(myGenes);

  // Calculate kappa score between each pair of terms
  double kappa;
  double totalGeneCount = countUniqueGenes(gene_vec);
  std::vector<std::string> term1_genes, term2_genes;

  for (int i=0; i<term_vec.size(); ++i) {
    term1_genes = splitString(gene_vec[i], ','); // Split current term genes
    for (int j=i+1; j<term_vec.size(); ++j) {
      term2_genes = splitString(gene_vec[j], ','); // Split compared term genes

      // Calculate the Kappa
      kappa = calKappa(term1_genes, term2_genes, totalGeneCount);
      kappaScores[std::make_pair(term_vec[i], term_vec[j])] = kappa; // Update kappaScores

      // Add term pairs with kappa > kappaCutoff to signifKappaTerms
      if (kappa >= kappaCutoff) {
        addTermToMap(signifKappaTerms, term_vec[i], term_vec[j]);
        addTermToMap(signifKappaScores, term_vec[i], kappa);
        //std::string result = "term 1: " + *term1 + ", term 2: " + *term2 + ", kappa: " + std::to_string(kappa);
        //std::cout << result << std::endl;
      }
    }
  }

  // Convert std::unordered_map to DataFrame
  DataFrame dfSignifKappaTerms = mapToDataFrame(signifKappaTerms, true);
  DataFrame dfSignifKappaScores = mapToDataFrame(signifKappaScores, false);
  DataFrame dfKappaScores = pairMapToDataFrame(kappaScores);

  // Return a List containing DataFrames
  return List::create(_["SignifKappaTerms"] = dfSignifKappaTerms, _["SignifKappaScores"] = dfSignifKappaScores, _["KappaScores"] = dfKappaScores);
}
