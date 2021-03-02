/* C3D: Computational protein Design by Duplication and Divergence (2021)
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#ifndef DEF_UTILS
#define DEF_UTILS

#include <vector>
#include <string>
#include <cstdlib>
#include <assert.h>
#include <iostream>
#include <ctime>
#include <random>

typedef unsigned long uInt;

struct Samples{
  uInt* data;
  uInt Nsamples;
  uInt N;
  uInt q;
};

struct PottsModel{
  uInt N;
  uInt q;
  float* P;
};

struct Parameters{
  std::string outPrefix;
  std::string prmInputFile;
  std::string samplesInputFile;
  std::string mutantsListFile;
  std::string mutatePositionsFile;
  int randomSeed;
  uInt Nsweeps;
  uInt sweepsPerSample;
  uInt blmIterations;
  uInt nMutants;
  uInt nPointMutations;
  uInt nPointMutationsDom2;
  uInt nDomainSplit;
  float learningRate;
  float regularizationLambda;
  float T;
  float probThreshold;
};

std::vector<std::string> splitString(const std::string& s, const char& delimiter);

inline uInt h_at(uInt i, uInt A,uInt q){
  return q*i +A;
}

inline uInt J_at(uInt i,uInt j, uInt A, uInt B,uInt N, uInt q){
  assert (i!=j);
  if(i<j)
    return q*N + (i*N - (i+2)*(i+1)/2 +j)*q*q + q*A + B;
  else
    return q*N + (j*N - (j+2)*(j+1)/2 +i)*q*q + q*B + A;
}
#endif
