/* C3D: Computational protein Design by Duplication and Divergence (2021)
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#include "Modes.h"
#include <cstring>
#include <iostream>
#include <mpi.h>

int main(int argc, char** argv){
  
  // Setup MPI 
  MPI_Init(&argc,&argv);
  
  // Select the running mode.
  int retVal=1;
  if (argc==1 || std::string(argv[1])=="-h" || std::string(argv[1])=="--help"){
    std::cout<<"Usage: c3d mode [options]"<<std::endl;
    std::cout<<"Available modes are: train, generate, select, ortho"<<std::endl;
    std::cout<<"No execution mode provided. Exiting."<<std::endl;
    return -1;
  }
  else if(!strcmp(argv[1],"train")){
    retVal=mainTrain(argc,argv);
  }
  else if(!strcmp(argv[1],"generate")){
    retVal=mainMCGenerate(argc,argv);
  }
  else if(!strcmp(argv[1],"select")){
    retVal=mainSelect(argc,argv);
  }
  else if(!strcmp(argv[1],"ortho")){
    retVal=mainOrtho(argc,argv);
  }
  else{
    std::cout<<"Unrecognized execution mode"<<std::endl;
  }

  MPI_Finalize();
  return retVal;
}
