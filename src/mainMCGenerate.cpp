/* C3D: Computational protein Design by Duplication and Divergence (2021)
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#include <iostream>
#include <mpi.h>
#include <time.h>  
#include "IOHandler.h"
#include "MCGenerator.h"

using namespace std;

int mainMCGenerate(int argc, char** argv){

  // Setup MPI 
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

  // Parse parameters and setup run
  IOHandler io;
  Parameters params = io.loadArgsMCGenerate(argc,argv);
  if(params.prmInputFile==""){
    cout<<"No prm file provided."<<endl;
    return -1;
  }
  else if (params.samplesInputFile==""){
    cout<<"No native sample provided."<<endl;
    return -1;
  }
  else if (params.nPointMutations==0 && params.nPointMutationsDom2==0){
    cout<<"Number of requested mutations on a domain is 0. Exiting."<<endl;
      return -1;
  }

  io.setOutPrefix(params.outPrefix+std::to_string(mpi_rank));
  if(params.randomSeed==0)
    srand(mpi_rank+time(NULL));
  else
    srand(params.randomSeed+mpi_rank);

  // Load Potts model and native sample
  PottsModel model = io.loadModel(params.prmInputFile);
  if(params.nDomainSplit==0)
    params.nDomainSplit=model.N-1;

  
  Samples samples = io.loadSamples(params.samplesInputFile);
  if(samples.Nsamples >1){
    cout<<"Warning: Multiple samples provided, first sample in file used as native."<<endl;
  }
  
  // Main routine: Generate the mutants starting from the native sample
  MC_generate(samples,params,model); 
  
  return 0;
}
