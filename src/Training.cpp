/* C3D: Computational protein Design by Duplication and Divergence (2021)
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#include "Training.h"
#include "IOHandler.h"

void TRN_train(float* P,Samples samples,Parameters params){

  // Setup MPI variables
  int mpi_rank,mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);

  // Setup model variables
  uInt N = samples.N;
  uInt q = samples.q;
  uInt nParams = q*N + q*q*N*(N-1)/2;
  uInt Nsamples = params.Nsweeps/params.sweepsPerSample;  
  float* sampledEs = new float[Nsamples];
  uInt* sampledData = new uInt[Nsamples*N];
  uInt* allSampledData;
  float* allSampledEs;
  uInt* startSequence = new uInt[N];
  uInt startId=0;
  std::ofstream log;
  if(mpi_rank==0){
    allSampledData = new uInt[Nsamples*N*mpi_size];
    allSampledEs = new float[Nsamples*mpi_size];
    log.open(params.outPrefix+".log");
    std::cout<<"# Training by Boltzman Learning"<<std::endl;
    std::cout<<"# N="<<N<<" , q="<<q<<" , B="<<samples.Nsamples<<std::endl;
    std::cout<<"# BLM steps="<<params.blmIterations<<" , Nsweeps="<<
      params.Nsweeps<<" , Nupdate="<<params.sweepsPerSample<<std::endl;
    std::cout<<"# Number of replicas="<<mpi_size<<std::endl;
    std::cout<<"# Random seed:"<<params.randomSeed<<std::endl;
    std::cout<<"# ============================================"<<std::endl;
    std::cout<<"#it"<<"   "<<"likelihood"<<"   "<<"meanError"<<"   "<<"time"<<std::endl;
    log<<"# Training by Boltzman Learning"<<std::endl;
    log<<"# N="<<N<<" , q="<<q<<" , B="<<samples.Nsamples<<std::endl;
    log<<"# BLM steps="<<params.blmIterations<<" , Nsweeps="<<
      params.Nsweeps<<" , Nupdate="<<params.sweepsPerSample<<std::endl;
    log<<"# Number of replicas="<<mpi_size<<" , "<<std::endl;
    log<<"# Random seed:"<<params.randomSeed<<std::endl;
    log<<"# ============================================"<<std::endl;
    log<<"#it"<<"   "<<"likelihood"<<"   "<<"meanError"<<"   "<<"time"<<std::endl;
  }

  // Setup optimization variables
  float* F_data = new float[nParams];
  TRN_computeFrequencies(F_data,samples.data,samples.Nsamples,samples.N,samples.q);
  float* F_model = new float[nParams];

  float likelihood,rmse,gdi,dt;
  clock_t t1;

  // Optimize parameters by simple gradient descent
  for(uInt it=0;it<params.blmIterations;it++){
    t1 = std::clock();

    // Compute the Log-likelihood gradient by MC sampling
    // Set the starting sequence
    startId = rand()%samples.Nsamples;
    for(uInt i=0;i<N;i++)
      startSequence[i]=samples.data[N*startId+i];
    
    // Perform the MC sampling
    MC_sample(startSequence,P,N,q,params.Nsweeps,params.sweepsPerSample,sampledEs,sampledData);
    
    // Gather all the samples on process 0
    MPI_Gather(sampledData,Nsamples*N,MPI_UNSIGNED_LONG,allSampledData,
	       Nsamples*N,MPI_UNSIGNED_LONG,0,MPI_COMM_WORLD);
    MPI_Gather(sampledEs,Nsamples,MPI_FLOAT,allSampledEs,Nsamples,
	       MPI_FLOAT,0,MPI_COMM_WORLD);
    // Compute the gradient update step on process 0
    if(mpi_rank==0){
      // Compute the model frequencies
      TRN_computeFrequencies(F_model,allSampledData,Nsamples*mpi_size,N,q);
    }
    
    if(mpi_rank==0){
	// Update the parameters
	likelihood=0.;
	rmse=0.;
	for(uInt i=0;i<nParams;i++){
	  gdi = params.learningRate*(F_data[i] - F_model[i] + params.regularizationLambda*P[i]);
	  P[i]-= gdi;
	  likelihood+= P[i]*F_data[i]+0.5*params.regularizationLambda*P[i]*P[i];
	  rmse+=gdi*gdi/(params.learningRate*params.learningRate);
	}
	
	// Print log
	likelihood/=N;
	rmse = sqrt(rmse/nParams);
	dt = float(std::clock()-t1)/CLOCKS_PER_SEC;
	std::cout<<it<<"   "<<likelihood<<"   "<<rmse<<"   "<<dt<<std::endl;
	log<<it<<"   "<<likelihood<<"   "<<rmse<<"   "<<dt<<std::endl;
    }
 
    // Broadcast the updated parameters P to all processes
    MPI_Bcast(P,nParams,MPI_FLOAT,0,MPI_COMM_WORLD);
  }
  
  //Clean up

  delete[] F_data;
  delete[] F_model;
  delete[] sampledEs;
  delete[] sampledData;
  delete[] startSequence;
  if(mpi_rank==0){
    delete[] allSampledData;
    delete[] allSampledEs;
  }
}

void TRN_computeFrequencies(float* F, uInt* data, uInt Nsamples,uInt N, uInt q){
  uInt xi=0;
  std::fill_n(F,q*N+ N*(N-1)*q*q/2,0.);
  for(uInt b=0;b<Nsamples;b++)
    for(uInt i=0;i<N;i++){
      xi = data[b*N+i];
      F[h_at(i,xi,q)]++;
      for(uInt j=i+1;j<N;j++)
	F[J_at(i,j,xi,data[b*N+j],N,q)]++;
    }
  for(uInt i=0;i<q*N+q*q*N*(N-1)/2;i++)
    F[i]/=Nsamples;
}
