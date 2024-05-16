#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>


// Rcpp function to split a comma-separated string into a vector
std::unordered_set<std::string> splitString(const std::string& input, char delimiter);

// Rcpp function to conjoin a vector into a comma-separated string
template <typename T>
std::string setToString(const std::unordered_set<T>& set, char delimiter);

// For counting the # of unique GeneIDs in an vector (like an R data.frame column)
double countUniqueGenes(const std::vector<std::string>& geneVec);

// Helper function adds term pair indices to associated term in sigTermIndices
void addTermToMap(std::unordered_map<int, std::unordered_set<int>>& map,
                  int key, int element);

// Function to convert std::unordered_map<int, std::vector<int>> to DataFrame
Rcpp::DataFrame mapToDataFrame(const std::unordered_map<int, std::unordered_set<int>>& map);


#endif
