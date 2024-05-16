// #include <Rcpp.h>
// using namespace Rcpp;
//
// #include <iostream>
// #include <algorithm>
// #include <utility>
// #include <functional>
// #include <cmath>
//
// #include <unordered_map>
// #include <vector>
// #include <set>
//
//
// struct KappaMap {
//   std::vector<double> map;
//   int n_terms;
//   KappaMap(int n) : n_terms(n){
//     map.resize(n*n); // Vector size is n_terms^2
//   }
//   // Get the encoded kappa index
//   double getKappa(int i, int j) const {
//     int kappaIndex = (i*n_terms)+j;
//     return map[kappaIndex];
//   }
//   // Encode a kappa score with a transformed index
//   void setKappa(double kappa, int i, int j) {
//     int kappaIndex = (i*n_terms)+j;
//     map[kappaIndex] = kappa;
//
//     kappaIndex = (j*n_terms)+i;
//     map[kappaIndex] = kappa;
//   }
// };
//
// struct KappaCluster {
//   std::set<int> pair_indices;
//   double avg_kappa;
// };
//
// struct InitSeedResult {
//   std::vector<KappaCluster> clusters;
//   std::vector<int> clusterIndices;
// };
//
// // Rcpp function to split a comma-separated string into a vector
// std::vector<std::string> splitString(const std::string& input, char delimiter) {
//   std::vector<std::string> result;
//   std::istringstream tokenStream(input);
//   std::string token;
//   while (std::getline(tokenStream, token, delimiter)) {
//     result.push_back(token);
//   }
//   return result;
// }
//
// // Rcpp function to conjoin a vector into a comma-separated string
// template<typename T>
// std::string vecToString(const std::vector<T>& vec, char delimiter) {
//   if (vec.empty()) {
//     return "";
//   }
//   std::ostringstream oss;
//   oss << vec[0];
//   for (size_t i = 1; i < vec.size(); ++i) {
//     oss << delimiter << vec[i];
//   }
//   return oss.str();
// }
//
// // Rcpp function to conjoin a set into a comma-separated string
// template<typename T>
// std::string setToString(const std::set<T>& mySet, char delimiter) {
//   if (mySet.empty()) {
//     return "";
//   }
//   std::ostringstream oss;
//   auto it = mySet.begin();
//   oss << *it;
//   ++it;
//
//   for (; it != mySet.end(); ++it) {
//     oss << delimiter << *it;
//   }
//   return oss.str();
// }
//
// // For counting the # of unique GeneIDs in an vector (like an R data.frame column)
// double countUniqueGenes(const std::vector<std::string>& geneVec) {
//   std::unordered_set<std::string> uniqueGenes;
//   std::vector<std::string> currentGenes;
//
//   // For each comma-delimited geneID cell
//   // int initCount = 0;
//   for (const std::string& gene : geneVec) {
//     currentGenes = splitString(gene, ','); // Split
//     // initCount += currentGenes.size();
//     // For each individual geneID
//     for (const std::string& g : currentGenes) {
//       uniqueGenes.insert(g); // if not present in uniqueGenes, add
//     }
//   }
//   // std::string oldSize = "original # geneIds: "+std::to_string(initCount);
//   // std::string newSize = "filtered # unique geneIds: "+std::to_string(uniqueGenes.size());
//   // std::cout << oldSize << std::endl;
//   // std::cout << newSize << std::endl;
//   return static_cast<double>(uniqueGenes.size()); // Return number of unique genes
// }
//
//
//
// // Function to calculate kappa score between two terms
// // t1_genes: Vector containing term1's geneIDs
// // t2_genes: Vector containing term2's geneIDs
// double calKappa(const std::vector<std::string>& t1_genes, const std::vector<std::string>& t2_genes, double totalGeneCount) {
//
//    // Return 0 if no overlapping genes
//    if (std::find_first_of(t1_genes.begin(), t1_genes.end(), t2_genes.begin(), t2_genes.end()) == t1_genes.end()) {
//      return 0.0;
//    }
//
//    // Else, calculate kappa of overlap
//    double common = 0, t1_only = 0, t2_only = 0;
//
//    for (const std::string& gene : t1_genes) {
//      // t1_gene found in t2_genes
//      if (std::find(t2_genes.begin(), t2_genes.end(), gene) != t2_genes.end()) {
//        common++;
//      } else {  // Gene unique to t1
//        t1_only++;
//      }
//    }
//    for (const std::string& gene : t2_genes) {
//      if (std::find(t2_genes.begin(), t2_genes.end(), gene) == t2_genes.end()) {
//        t2_only++;
//      }
//    }
//
//    double unique = totalGeneCount - common - t1_only - t2_only; // Count of all genes not found in either term
//
//    double relative_agree = (common + unique)/totalGeneCount;
//    double chance_yes = ((common + t1_only)/totalGeneCount) + ((common + t2_only)/totalGeneCount);
//    double chance_no = ((unique + t1_only)/totalGeneCount) + ((unique + t2_only)/totalGeneCount);
//    double chance_agree = chance_yes * chance_no;
//
//    if (chance_agree == 1) {
//      return 0.0; // Prevent divide by zero
//    } else {
//      return (relative_agree - chance_agree) / (1 - chance_agree); // Return kappa!
//    }
// }
//
//
// double getKappaDelta2(const KappaMap& map, KappaCluster& seed, int term_j) {
//   double delKappaSum = 0;
//   for (const auto term_i : seed.pair_indices) {
//     delKappaSum += map.getKappa(term_i, term_j);
//   }
//   return delKappaSum/seed.pair_indices.size();
// }
//
//
// void testInitSeeds(const KappaMap& map, const std::unordered_map<int, std::vector<int>>& sigTermIndices, double kappaCutoff=0.5) {
//   std::vector<int> j_terms;
//
//   for (const auto&termPair : sigTermIndices) {
//     int i = termPair.first;
//     j_terms = termPair.second;
//     for (const int& j : j_terms) {
//       std::cout << "{(" << i <<", " << j << "): " <<  map.getKappa(i, j) << "}, ";
//     }
//     std::cout << std::endl;
//   }
// }
//
//
// // Get initial seed clusters for all term pairs with kappa=1
// InitSeedResult getInitSeeds(const KappaMap& map, const std::unordered_map<int, std::vector<int>>& sigTermIndices, double kappaCutoff=0.5) {
//   InitSeedResult result;
//   std::vector<int> clusterIndices(map.n_terms, -1); // Stores cluster index if term is in cluster, else -1
//   double k_delta;
//   int n_clusters = 0;
//
//   std::vector<KappaCluster> clusters;
//   KappaCluster initSeed;
//
//   std::vector<int> j_terms;
//   for (const auto& termPair : sigTermIndices) {
//
//     int i = termPair.first;
//     j_terms = termPair.second;
//
//     for (const int& j : j_terms) {
//       initSeed.pair_indices = {}; // Reset pair indices
//       initSeed.avg_kappa = -1;
//
//       // Case 1: both terms not yet in cluster group
//       if (clusterIndices[i] == -1 && clusterIndices[j] == -1) {
//         std::cout << "case 1: both terms not yet in cluster" << std::endl;
//         k_delta = map.getKappa(i, j);
//         std::cout << "kappa (" << i << ", " << j << "): " << k_delta << std::endl; //for debug
//         if (k_delta >= kappaCutoff) {
//           // Set dummy variable with indices and avg_kappa
//           initSeed.pair_indices.insert(i);
//           initSeed.pair_indices.insert(j);
//           initSeed.avg_kappa = k_delta;
//
//           // Mark cluster index in clusterIndices
//           clusterIndices[i] = n_clusters;
//           clusterIndices[j] = n_clusters;
//           n_clusters++;
//
//           clusters.push_back(initSeed); // Append dummy to list of clusters
//         }
//       }
//       // Case 2: both terms in clusters
//       else if (clusterIndices[i] != -1 && clusterIndices[j] != -1) {
//         std::cout << "case 2" << std::endl;
//         // Both in same cluster (skip)
//         if (clusterIndices[i] == clusterIndices[j]) {
//           continue;
//         }
//         // Different clusters
//         else {
//           int i_idx = clusterIndices[i];
//           int j_idx = clusterIndices[j];
//
//           int k_delta1 = getKappaDelta2(map, clusters[i_idx], j);
//           int k_delta2 = getKappaDelta2(map, clusters[j_idx], i);
//
//           if (k_delta1 >= k_delta2) {
//             if (k_delta1 >= kappaCutoff) {
//               // merge j into clusters[i_idx]
//               clusterIndices[j] = i_idx; // Set j cluster# to same as i
//               clusters[j_idx].pair_indices.erase(j); // Erase from original j cluster
//
//               double prev_CS = clusters[i_idx].pair_indices.size(); // prev cluster size
//               double prev_k = (1/2)*prev_CS*(prev_CS-1); // nCr combinations given n=prev_CS, r=2
//               double prev_k_delta = clusters[i_idx].avg_kappa / prev_k;
//               double new_kappa_avg = (prev_k_delta + k_delta1) / (prev_k + prev_CS); // 1*prev_CS unique combinations added
//
//               clusters[i_idx].pair_indices.insert(j);
//               clusters[i_idx].avg_kappa = new_kappa_avg;
//               std::cout << "merging term j into cluster i" << k_delta << std::endl; //for debug
//             }
//           } else if (k_delta2 > k_delta1) {
//             if (k_delta2 >= kappaCutoff) {
//               // merge i into clusters[j_idx]
//               clusterIndices[i] = j_idx; // Set i cluster# to same as j
//               clusters[i_idx].pair_indices.erase(i); // Erase from original i cluster
//
//               double prev_CS = clusters[j_idx].pair_indices.size(); // prev cluster size
//               double prev_k = (1/2)*prev_CS*(prev_CS-1); // nCr combinations given n=prev_CS, r=2
//               double prev_k_delta = clusters[j_idx].avg_kappa / prev_k;
//               double new_kappa_avg = (prev_k_delta + k_delta2) / (prev_k + prev_CS); // 1*prev_CS unique combinations added
//
//               clusters[j_idx].pair_indices.insert(i);
//               clusters[j_idx].avg_kappa = new_kappa_avg;
//             }
//           }
//         //   if (k_delta >= kappaCutoff) {
//         //     // Merge into the cluster with lower index
//         //     int low_cluster_idx = std::min(i_idx, j_idx);
//         //     int high_cluster_idx = std::max(i_idx, j_idx);
//         //     for (const auto& term2 : clusters[high_cluster_idx].pair_indices) {
//         //       if (clusterIndices[term2] == low_cluster_idx) {
//         //         continue;
//         //       } else {
//         //         clusterIndices[term2] = low_cluster_idx;
//         //
//         //         double prev_CS = clusters[low_cluster_idx].pair_indices.size(); // prev cluster size
//         //         double prev_k = (1/2)*prev_CS*(prev_CS-1); // nCr combinations given n=prev_CS, r=2
//         //         double prev_k_delta = clusters[low_cluster_idx].avg_kappa / prev_k;
//         //         double new_kappa_avg = (prev_k_delta + k_delta) / (prev_k + prev_CS); // 1*prev_CS unique combinations added
//         //
//         //         clusters[low_cluster_idx].pair_indices.insert(term2);
//         //         clusters[low_cluster_idx].avg_kappa = new_kappa_avg;
//         //
//         //       }
//         //     }
//         //     clusters.erase(clusters.begin() + high_cluster_idx);
//         //   }
//         }
//       } else if (clusterIndices[i] != -1 && clusterIndices[j] == -1) {
//         // case 3: i is cluster, j is term
//         std::cout << "case 3" << std::endl;
//         int cluster_i_index = clusterIndices[i];
//         k_delta = getKappaDelta2(map, clusters[cluster_i_index], j);
//         std::cout << "cluster_i, term_j KappaDelta: " << k_delta << std::endl; //for debug
//
//         if (k_delta >= kappaCutoff) {
//           clusterIndices[j] = cluster_i_index; // Set j cluster# to same as i
//
//           double prev_CS = clusters[cluster_i_index].pair_indices.size(); // prev cluster size
//           double prev_k = (1/2)*prev_CS*(prev_CS-1); // nCr combinations given n=prev_CS, r=2
//           double prev_k_delta = clusters[cluster_i_index].avg_kappa / prev_k;
//           double new_kappa_avg = (prev_k_delta + k_delta) / (prev_k + prev_CS); // 1*prev_CS unique combinations added
//
//           clusters[cluster_i_index].pair_indices.insert(j);
//           clusters[cluster_i_index].avg_kappa = new_kappa_avg;
//         }
//       } else {
//         // case 4: i is term, j is cluster
//         std::cout << "case 4" << std::endl;
//         k_delta = getKappaDelta2(map, clusters[j], i);
//         std::cout << "cluster_i, term_j KappaDelta: " << k_delta << std::endl; //for debug
//
//         if (k_delta >= kappaCutoff) {
//           int cluster_j_index = clusterIndices[j];
//           clusterIndices[i] = cluster_j_index; // Set i cluster# to same as j
//
//           double prev_CS = clusters[cluster_j_index].pair_indices.size(); // prev cluster size
//           double prev_k = (1/2)*prev_CS*(prev_CS-1); // nCr combinations given n=prev_CS, r=2
//           double prev_k_delta = clusters[cluster_j_index].avg_kappa / prev_k;
//           double new_kappa_avg = (prev_k_delta + k_delta) / (prev_k + prev_CS); // 1*prev_CS unique combinations added
//
//           clusters[cluster_j_index].pair_indices.insert(i);
//           clusters[cluster_j_index].avg_kappa = new_kappa_avg;
//         }
//       }
//
//     }
//   }
//   result.clusters = clusters;
//   result.clusterIndices = clusterIndices;
//   return result;
// }
//
//
// // Convert distance vector to 2-dimensional DataFrame
// DataFrame makeAllKappaDf(const KappaMap& map, CharacterVector myTerms) {
//   NumericMatrix kappa_matrix(map.n_terms, map.n_terms);
//   // Fill matrix with kappa scores
//   int k = 0;
//   for (int i=0; i<map.n_terms; ++i) {
//     for (int j=0; j<map.n_terms; ++j) {
//       if (i == j) {
//         kappa_matrix(i, j) = -99.0;
//       } else {
//         kappa_matrix(i, j) = map.map[k];
//         kappa_matrix(j, i) = map.map[k];
//       }
//       k++;
//     }
//   }
//   DataFrame kappa_df = internal::convert_using_rfunction(kappa_matrix, "as.data.frame");
//   kappa_df.attr("names") = myTerms;
//   kappa_df.attr("row.names") = myTerms;
//   return kappa_df;
// }
//
//
// // Helper function adds term pair indices to associated term in sigTermIndices
// void addTermToMap(std::unordered_map<int, std::vector<int>>& map,
//                   int key, int element) {
//
//   // Adding term to SignifKappaTerms
//   auto it = map.find(key);
//   std::vector<int> keyValue;
//   if (it != map.end()) {
//     // Key exists, add the element to the existing vector
//     keyValue = it->second;
//     keyValue.push_back(element);
//     it->second = keyValue;
//   } else {
//     // Key does not exist, create a new vector with the element and insert into the map
//     keyValue = {element};
//     map[key] = keyValue;
//   }
//
//   // Do the same for the other term!
//   auto it2 = map.find(element);
//   std::vector<int> keyValue2;
//   if (it2 != map.end()) {
//     // Key exists, add the element to the existing vector
//     keyValue2 = it2->second;
//     keyValue2.push_back(key);
//     it2->second = keyValue2;
//   } else {
//     // Key does not exist, create a new vector with the element and insert into the map
//     keyValue2 = {key};
//     map[element] = keyValue2;
//   }
//
// }
//
// // Temporary functions to test functionality of code
// // Function to convert std::unordered_map<int, std::vector<int>> to DataFrame
// DataFrame mapToDataFrame(const std::unordered_map<int, std::vector<int>>& map) {
//   std::vector<int> keys;
//   std::vector<std::string> values;
//   std::string valueString;
//
//   for (const auto& entry : map) {
//     keys.push_back(entry.first);
//     valueString = vecToString(entry.second, ',');
//     values.push_back(valueString);
//   }
//   return DataFrame::create(_["TermIndex"] = keys, _["Pair_KappaIndices"] = values);
// }
//
// // Convert vector of KappaClusters into a dataframe
// DataFrame makeClusterDf(const std::vector<KappaCluster>& clusters) {
//
//   std::vector<double> avgKappaVec;
//   std::vector<std::string> clusterIndexVec;
//   std::string index;
//
//   for (const auto& cluster : clusters) {
//     avgKappaVec.push_back(cluster.avg_kappa);
//     index = setToString(cluster.pair_indices, ',');
//     clusterIndexVec.push_back(index);
//   }
//   return DataFrame::create(_["AvgKappa"] = avgKappaVec, _["ClusterPairIndices"] = clusterIndexVec);
// }
//
//
// //' @export
// // [[Rcpp::export]]
// List RcppKappaCluster2(CharacterVector myTerms, CharacterVector myGenes,
//                       NumericVector myPvalues, double kappaCutoff = 0.5) {
//
//   // Convert R vectors to C++
//   std::vector<std::string> all_terms = Rcpp::as<std::vector<std::string>>(myTerms);
//   std::vector<std::string> all_genes = Rcpp::as<std::vector<std::string>>(myGenes);
//   std::vector<double> all_pvalues = Rcpp::as<std::vector<double>>(myPvalues);
//
//   std::unordered_map<int, std::vector<int>> sigTermIndices;
//
//   // all_kappas: Flattened distance matrix
//   // Index a given pair (term_i, term_j) with all_kappas[(i*j)+j]
//   int n_terms = all_terms.size();
//   KappaMap kappa_map(n_terms);
//   kappa_map.n_terms = n_terms;
//
//   // Calculate kappa score between each pair of terms
//   double kappa;
//   double totalGeneCount = countUniqueGenes(all_genes);
//   std::vector<std::string> term1_genes, term2_genes;
//
//   // Create all_kappas distance matrix
//   for (int i=0; i<n_terms; ++i) {
//     term1_genes = splitString(all_genes[i], ',');
//     for (int j=i; j<n_terms; ++j) {
//       // kappaIndex = (i*n_terms)+j;
//       if (j == i) {
//         kappa_map.setKappa(-99, i, j);
//       } else {
//         term2_genes = splitString(all_genes[j], ',');
//         kappa = calKappa(term1_genes, term2_genes, totalGeneCount);
//         kappa_map.setKappa(kappa, i, j);
//
//         if (kappa >= kappaCutoff) {
//           addTermToMap(sigTermIndices, i, j);
//         }
//       }
//     }
//   }
//
//
//   std::cout << "this works 0" << std::endl;
//   // testInitSeeds(kappa_map, sigTermIndices, kappaCutoff);
//   // InitSeedResult seedResult = getInitSeeds(kappa_map, sigTermIndices, kappaCutoff);
//
//   DataFrame dfSigTermIndices = mapToDataFrame(sigTermIndices);
//   DataFrame dfAllKappasVec = DataFrame::create(_["TermIndex"] = kappa_map.map);
//   DataFrame dfAllKappasMat = makeAllKappaDf(kappa_map, myTerms);
//   // DataFrame dfClusters = makeClusterDf(seedResult.clusters);
//   //
//   // return List::create(_["dfSigTermIndices"] = dfSigTermIndices, _["all_kappas"] = dfAllKappasVec,
//   //                     _["AllKappaMatrix"] = dfAllKappasMat, _["Clusters"] = dfClusters);
//   return List::create(
//     _["dfSigTermIndices"] = dfSigTermIndices,
//     _["all_kappas"] = dfAllKappasVec,
//     _["AllKappaMatrix"] = dfAllKappasMat
//   );
// }
//
//
