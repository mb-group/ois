/* OIS: Orthogonal Interacting Sequences (2022) 
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#include <iostream>
#include <time.h>  
#include "IOHandler.h"
#include "Selector.h"

using namespace std;

int mainSelect(int argc, char** argv){

  // Parse parameters and setup run
  IOHandler io;
  Parameters params = io.loadArgsSelect(argc,argv);
  if(params.samplesInputFile==""){
    cout<<"No mutants file provided"<<endl;
    return -1;
  }
  else if (params.probThreshold<0 || params.probThreshold>1){
    cout<<"Invalid probability threshold (or not provided)."<<endl;
    return -1;
  }
    
  SEL_Select(params.samplesInputFile,params.probThreshold,params.outPrefix);
  
  return 0;
}
