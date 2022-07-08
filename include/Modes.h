/* C3D: Computational protein Design by Duplication and Divergence (2021)
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#ifndef DEF_MODES
#define DEF_MODES

// Training mode: Learn the Potts model parameters from data in alignment.
int mainTrain(int argc, char** argv);

// Generative mode: Generate mutants, evaluating energies with respect to endogenous partners.
int mainMCGenerate(int argc, char** argv);

// Selection mode: Select generated mutants for orthogonality properties.
int mainSelect(int argc, char** argv);

// Orthogonality mode: Directly sample from low E(A*,B*) and high E(A*,B) and E(A,B*)
int mainOrtho(int argc, char** argv);
#endif
