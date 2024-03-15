#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <algorithm>
#include <utility>
#include <functional>

#include <unordered_map>
#include <vector>


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

// Convert distance vector to 2-dimensional DataFrame
DataFrame makeAllKappaDf(std::vector<double>& all_kappas, CharacterVector myTerms, int terms_length) {
  NumericMatrix kappa_matrix(terms_length, terms_length);

  // Fill matrix with kappa scores
  int i = 0, j = 0;
  int all_kappas_size = all_kappas.size();
  for (int k=0; k<all_kappas_size; ++k) {
    i = k/terms_length;
    j = k%terms_length;
    kappa_matrix(i, j) = all_kappas[k];
    kappa_matrix(j, i) = all_kappas[k];
  }
  //DataFrame kappa_df = DataFrame::create(kappa_matrix);
  DataFrame kappa_df = internal::convert_using_rfunction(kappa_matrix, "as.data.frame");
  kappa_df.attr("names") = myTerms;
  kappa_df.attr("row.names") = myTerms;

  return kappa_df;
}


// Helper function adds term pair indices to associated term in sigTermIndices
void addTermToMap(std::unordered_map<double, std::vector<double>>& map,
                  double key, double element) {

  // Adding term to SignifKappaTerms
  auto it = map.find(key);
  std::vector<double> keyValue;
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

  // Do the same for the other term!
  auto it2 = map.find(element);
  std::vector<double> keyValue2;
  if (it2 != map.end()) {
    // Key exists, add the element to the existing vector
    keyValue2 = it2->second;
    keyValue2.push_back(key);
    it2->second = keyValue2;
  } else {
    // Key does not exist, create a new vector with the element and insert into the map
    keyValue2 = {key};
    map[element] = keyValue2;
  }

}

// Temporary functions to test functionality of code
// Function to convert std::unordered_map<double, std::vector<double>> to DataFrame
DataFrame mapToDataFrame(const std::unordered_map<double, std::vector<double>>& map) {
  std::vector<double> keys;
  std::vector<std::string> values;
  std::string valueString;

  for (const auto& entry : map) {
    keys.push_back(entry.first);
    valueString = vecToString(entry.second, ',');
    values.push_back(valueString);
  }
  return DataFrame::create(_["TermIndex"] = keys, _["Pair_KappaIndices"] = values);
}


//' @export
// [[Rcpp::export]]
List RcppKappaCluster2(CharacterVector myTerms, CharacterVector myGenes,
                      NumericVector myPvalues, double kappaCutoff = 0.5) {

  // Convert R vectors to C++
  std::vector<std::string> all_terms = Rcpp::as<std::vector<std::string>>(myTerms);
  std::vector<std::string> all_genes = Rcpp::as<std::vector<std::string>>(myGenes);
  std::vector<double> all_pvalues = Rcpp::as<std::vector<double>>(myPvalues);

  std::unordered_map<double, std::vector<double>> sigTermIndices;

  // all_kappas: Flattened distance matrix
  // Index a given pair (term_i, term_j) with all_kappas[(i*j)+j]
  int terms_length = all_terms.size();
  std::vector<double> all_kappas(terms_length*terms_length);
  //std::vector<double> all_kappas;

  // Calculate kappa score between each pair of terms
  double kappa, kappaIndex;
  double totalGeneCount = countUniqueGenes(all_genes);
  std::vector<std::string> term1_genes, term2_genes;

  // Create all_kappas distance matrix
  for (int i=0; i<terms_length; ++i) {
    term1_genes = splitString(all_genes[i], ',');
    for (int j=i; j<terms_length; ++j) {
      kappaIndex = (i*terms_length)+j;
      if (j == i) {
        all_kappas[kappaIndex] = -99.0;
      } else {
        term2_genes = splitString(all_genes[j], ',');
        kappa = calKappa(term1_genes, term2_genes, totalGeneCount);
        all_kappas[kappaIndex] = kappa;

        if (kappa >= kappaCutoff) {
          addTermToMap(sigTermIndices, i, j);
        }
      }
    }
  }
  DataFrame dfSigTermIndices = mapToDataFrame(sigTermIndices);
  DataFrame dfAllKappasVec = DataFrame::create(_["TermIndex"] = all_kappas);
  DataFrame dfAllKappasMat = makeAllKappaDf(all_kappas, myTerms, terms_length);

  return List::create(_["dfSigTermIndices"] = dfSigTermIndices, _["all_kappas"] = dfAllKappasVec, _["AllKappaMatrix"] = dfAllKappasMat);
}


