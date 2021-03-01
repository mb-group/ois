# Computational protein Design by Duplication and Divergence (C3D)

  1. [Installation & Requirements](#installation&requirements)	
  2. [Basic usages](#basic-usages)
  3. [File formats](#file-formats)
  4. [Utility functions](#utility-functions)
  5. [Tutorial](#tutorial)
  
  
# Installation & Requirements
To compile and run c3d, a working MPI library must be installed. A C++ compiler supporting at least the C++17 standard is required.
C3D has been tested using [OpenMPI](https://www.open-mpi.org/) on MacOS and Ubuntu.

To compile the software run

```shell
	$ make
```
from the root C3D folder. This creates a single executable c3d, which is used to perform the different generation and analysis steps (see [Full example](#full-example)). 

To use the utility python scripts in the utilities/ folder, python dependencies can be isntalled by
```shell
	$ pip install -r requirements.txt
```

# Basic usages
C3D is built around three exection modes: Model training, sample generation and sample selection. The -h options displays the available options for any mode.

```shell
	$ c3d -h
	
	Usage: c3d mode [options]
	Available modes are: train, generate, select
	No execution mode provided. Exiting.
```
## Training mode
Train a generalized potts model (the global statistical model) to sequences in a multiple-sequence alignment. To convert the sequencse in the apropriate format, use the fastaToMatrix.py utility (see [Utility Functions](#utility-functions)). The training is performed by approximate minimization of a L2-regularized log-likelihood function. The gradient is estimated by Markov-Chain Monte-Carlo sampling at each iteration (Boltzman-learning approach). To run **X** parallel replicas to estimate the gradient, run

```shell
	$ mpirun -N X c3d train -f samplesFile [options]
```
Available options are 
```shell
	$ c3d train -h

	  Usage: c3d train -f samplesFile [options]
                 -f       : Input file for learning (space delimited raw format)
                 -N       : Total number of sweeps to perform [default 10000]
                 -n       : Number of sweeps between recording two samples [default 10]
                 -p       : Potts model starting parameters file in prm format
                 -o       : Output prefix for saving files [default "output"]
                 -b       : Number of iterations of the BLM optimization [default 1000]
                 -l       : Regularization parameter lambda [default 0.01]
                 -r       : Learning rate eta [default 0.01]
                 --seed   : Random number generator seed [default 0 == time].
```

## Generation mode
Generate sample from the global statistical model by constraint Markov-Chain Monte-Carlo sampling. The sampling can be restriced to using only a subset of all positions (-pi option). The -m and -m2 options define how many mutations are to be made on the first (resp. second) protein. The  -ns option indicates where the first protein end in the concatenated alignment. The constraint sampling ensures that each generated sample has m (resp. m2) mutations on the first (resp. second) protein.

Available options are 

```shell
	$ c3d generate -h
	
	  Usage: c3d generate -f nativeFile -p prmFile [options]
                 -f       : Native sample file (space delimited one-line sample file)
                 -p       : Potts model parameters file in prm format
                 -M       : Number of mutants to compute. [default 100]
                 -m       : Number of point mutations per mutant [default 0]
                 -m2      : Number of point mutations per mutant on the second domain [default 0]
                 -T       : Virtual scaling temperature in the MCMC sampling [default 1].
                 -o       : Output prefix for saving files [default "output"]
                 -pi      : List of positions (0-based indexes) onto which to restrict the mutations [default None]
                 -ns      : Last index (inclusive, 0-based) of first domain,
                            used for domain split E computation.  [default None]
                 -n       : Number of sweeps between recording two mutants [default 10]
                 --seed   : Random number generator seed [default 0 == time].
```
## Selection mode
Select candidate mutants with potential for being orthogonal. The selection criteria is the probability threshold, i.e. the fraction of mutants lying in the upper-right quadrant of the non-cognate energies (see paper). The absolute threshold are found by iterative line search.

Available options are 
```shell
	$ c3d select -h

	  Usage: c3d select -f mutantsFile -t probThreshold [options]]
                 -f       : Mutants file, as output by c3d generate (comprising mutations and energies).
                 -t       : Probability threshold (in [0,1]) to select orthogonal mutants.
                 -o       : Output prefix for saving files [default "output"]
```

# File formats
# Utility functions

# Tutorial
