/* OrthoSeq: Computational protein Design by Duplication and Divergence (2021)
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#include "Selector.h"

void SEL_Select(std::string mutantsFile, float probThreshold, std::string outPrefix){

  // Load Orthogonal energies dInter_A*B and dInter_AB* from input file
  std::ifstream in(mutantsFile.c_str());
  std::string tmp;
  std::vector<std::string> splitLine;

  std::vector<double> dEsAsB;
  std::vector<double> dEsABs;
  while(!in.eof()){
    getline(in,tmp);
    splitLine=splitString(tmp,'\t');
    if(tmp[0]!='#' && splitLine.size()>5){
      dEsAsB.push_back(atof((splitLine[splitLine.size()-2]).c_str()));
      dEsABs.push_back(atof((splitLine[splitLine.size()-1]).c_str()));
    }
  }
  in.close();
  unsigned int Nmutants=dEsAsB.size();
    
  // Find threshold values by iterative line search
  float xm=0;
  float ym=0;
  for(unsigned int i=0;i<Nmutants;i++){
    xm+=dEsAsB[i];
    ym+=dEsABs[i];
  }
  xm/=Nmutants;
  ym/=Nmutants;
  for(unsigned int i=0;i<Nmutants;i++){
    dEsAsB[i]-=xm;
    dEsABs[i]-=ym;
  }
  double m=*std::max_element(std::begin(dEsAsB),std::end(dEsAsB));
 
  unsigned int nSteps=10000;
  float t;
  float frac=0;
  float Threshold;
  for(unsigned int s=0;s<nSteps;s++){
    t=s*(m/(nSteps-1));
    frac=0;
    for(unsigned int i=0;i<Nmutants;i++){
      frac+=(dEsAsB[i]>=t)*(dEsABs[i]>=t);
    }
    frac/=Nmutants;
    if(frac<=probThreshold){
      Threshold=t;
      break;
    }
  }
  float t1=xm+Threshold;
  float t2=ym+Threshold;

  // Extract mutants above the thresholds
  in.open(mutantsFile.c_str());
  std::ofstream out(outPrefix+"_selected.dat");
  while(!in.eof()){
    getline(in,tmp);
    splitLine=splitString(tmp,'\t');
    if(tmp[0]=='#')
      out<<tmp<<std::endl;
    else if(splitLine.size()>5 and atof((splitLine[splitLine.size()-2]).c_str())>t1 && atof((splitLine[splitLine.size()-1]).c_str())>t2)
      out<<tmp<<std::endl;
  }
  in.close();
  out.close();
}
