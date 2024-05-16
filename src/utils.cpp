#include <Rcpp.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "utils.h"


// utils.h contains helper functions/classes for clustering
//  -> splitString
//  -> helper fcts for SigTermIndices: addTermToMap, mapToDataFrame


// Explicit instantiation for std::string
// template std::string setToString(const std::unordered_set<std::string>& set, char delimiter);


//' Rcpp function to split a comma-separated string into an unordered_set
//'
//' @param input: string, the input string to be split
//' @param delimiter: char, the character separating desired values in input
//' @return result: vector, of strings containing the split values
std::unordered_set<std::string> splitString(const std::string& input, char delimiter) {
  std::unordered_set<std::string> result;
  std::istringstream tokenStream(input);
  std::string token;
  while (std::getline(tokenStream, token, delimiter)) {
    result.insert(token); // also, set insertion removes duplicates!
  }
  return result;
}

// Rcpp function to conjoin a vector into a comma-separated string
template <typename T>
std::string setToString(const std::unordered_set<T>& set, char delimiter) {
  if (set.empty()) {
    return "";
  }
  std::ostringstream oss;
  auto it = set.begin();
  oss << *it; // Output the first element directly
  ++it;
  for (; it != set.end(); ++it) {
    oss << delimiter << *it;
  }
  return oss.str();
}


// For counting the # of unique GeneIDs in an vector (like an R data.frame column)
double countUniqueGenes(const std::vector<std::string>& genes) {
  std::unordered_set<std::string> uniqueGenes;

  // For each comma-delimited geneID cell
  // int initCount = 0;
  for (const std::string& gene : genes) {
    std::unordered_set<std::string> currentGenes = splitString(gene, ','); // Split
    // For each individual geneID
    for (const std::string& g : currentGenes) {
      uniqueGenes.insert(g); // if not present in uniqueGenes, add
    }
  }
  return static_cast<double>(uniqueGenes.size()); // Return number of unique genes
}


// Helper function adds term pair indices to associated term in sigTermIndices
void addTermToMap(std::unordered_map<int, std::unordered_set<int>>& map,
                  int key, int element) {

  // Adding term to SignifKappaTerms
  auto it = map.find(key);
  std::unordered_set<int> keyValue;
  if (it != map.end()) {
    // Key exists, add the element to the existing set
    keyValue = it->second;
    keyValue.insert(element);
    it->second = keyValue;
  } else {
    // Key does not exist, create a new vector with the element and insert into the map
    keyValue = {element};
    map[key] = keyValue;
  }

  // Do the same for the other term!
  auto it2 = map.find(element);
  std::unordered_set<int> keyValue2;
  if (it2 != map.end()) {
    // Key exists, add the element to the existing set
    keyValue2 = it2->second;
    keyValue2.insert(key);
    it2->second = keyValue2;
  } else {
    // Key does not exist, create a new vector with the element and insert into the map
    keyValue2 = {key};
    map[element] = keyValue2;
  }

}

// Temporary functions to test functionality of code
// Function to convert std::unordered_map<int, std::vector<int>> to DataFrame
Rcpp::DataFrame mapToDataFrame(const std::unordered_map<int, std::unordered_set<int>>& map) {
  std::vector<int> keys;
  std::vector<std::string> values;

  for (const auto& entry : map) {
    keys.push_back(entry.first);
    std::string valueString = setToString(entry.second, ',');
    values.push_back(valueString);
  }
  return Rcpp::DataFrame::create(Rcpp::_["TermIndex"] = keys, Rcpp::_["Pair_KappaIndices"] = values);
}


