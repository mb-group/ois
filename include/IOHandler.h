/* OIS: Orthogonal Interacting Sequences (2022)  
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#ifndef DEF_IOHANDLER
#define DEF_IOHANDLER

#include <cstring>
#include <fstream>
#include <algorithm>
#include "Utils.h"

class IOHandler{

 public:
  IOHandler();
  ~IOHandler();
  void setOutPrefix(std::string out);
  Parameters parseInputArgs(int argc, char** argv);
  Parameters loadArgsTrain(int argc, char** argv);
  Parameters loadArgsMCGenerate(int argc, char** argv);
  Parameters loadArgsOrtho(int argc, char** argv);
  Parameters loadArgsSelect(int argc, char** argv); 
  
  PottsModel loadModel(std::string prmFile);
  Samples loadSamples(std::string sampleFile);
  
  void writeConfigurations(uInt* configurations, uInt Nsamples, uInt N);
  void writeEnergies(float* Es, uInt Nsamples,uInt nFields=1);
  void writePrm(float* P,uInt N, uInt q);
  void writeFrequencies(float* F,uInt N, uInt q);
  
 private:
  std::string outPrefix;
};
#endif
