#ifndef FUZZYCLUSTERING_H
#define FUZZYCLUSTERING_H

#include "SeedMap.h"

class FuzzySeeding {
public:
  FuzzySeeding(SeedMap& seedMap) : _seedMap(seedMap) {}
  
  
private:
  SeedMap& _seedMap;
};

#endif