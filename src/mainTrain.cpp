/* OrthoSeq: Computational protein Design by Duplication and Divergence (2021)
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#include <iostream>
#include <mpi.h>
#include <time.h>  
#include "IOHandler.h"
#include "Training.h"


using namespace std;

int mainTrain(int argc, char** argv){

  // Setup MPI 
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  
  // Parse parameters and setup run
  IOHandler io;
  Parameters params = io.loadArgsTrain(argc,argv);
  if(params.samplesInputFile==""){
    cout<<"No sample file provided"<<endl;
    return -1;
  }
  if(params.randomSeed==0)
    srand(mpi_rank+time(NULL));
  else
    srand(params.randomSeed+mpi_rank);
  
  // Load Samples to fit
  Samples samples = io.loadSamples(params.samplesInputFile);
  
  // Setup initial parameters: Initiate Potts parameters to 0 if not supplied (i.e. no restart)
  uInt N = samples.N;
  uInt q = samples.q;
  uInt nParams = q*N + q*q*N*(N-1)/2;
  float* P;
  if(params.prmInputFile!=""){
    PottsModel model = io.loadModel(params.prmInputFile);
    P = model.P;
  }
  else{
    P = new float[nParams];
    std::fill_n(P,nParams,0.);
  }

  // Perform the Boltzman learning procedure
  TRN_train(P,samples,params);

  if(mpi_rank==0)
    io.writePrm(P,N,q);

  // Clean Up
  delete[] P;

  return 0;
}
