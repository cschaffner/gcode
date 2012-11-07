MATLAB code for genetic code optimizations. 
These tools have been used to obtain the numerical results and plots in a forthcoming article
by Harry Buhrman, Peter van der Gulik, Dave Spijer, Gunnar Klau, Leen Stougie, Christian Schaffner

created: Oct 21, 2011
last modified: Nov 7, 2012
by Christian Schaffner, c.schaffner@uva.nl

Installation:
* Make sure the lib/ subdirectory is part of the path


This directory contains the following files:
VariousGraphs.m : main program module to set parameters and run computations
makegraph.m : draws histograms from computed values

PaperGraphs.m : used to generate the data for histograms in the paper
makePapergraph.m : used to draw the histograms that appear in main body of the paper
makeAppendixGraph.m: used to draw the histograms that appear in the appendix of the paper

output/ : output directory where .fig and .pdf files are saved to
PaperOutput/ : output directory with the graphs from the paper

Eppstein/ : contains MATLAB code for the inverse parametric optimization
  Eppstein_objective.m  :  gives the objective function           
  Eppstein_q20_constraint.m : gives the constraints
  Eppstein_steer.m : The main procedure
  Eppstein_steer_all1weights_popt.mat : stores intermediate results
  Eppstain_comparison.m : contains several sets of 20 numerical values which make the SGC optimal
                          and plots them together with polar requirements.

min_assignments/ : directory with matlab code for checking which are the minimal assignments to 
                   fix such that the sgc is optimal when using the subset approach
  permutecode_allsubsets.m : checks all permutations allowed by the subset approach and keeps the ones
                             that are better than the sgc
  minimal_fixed_assignment.m : based on the output of permutecode_allsubsets.m , finds out which are the 
                               minimal subsets one can fix such that the sgc is optimal
  th_polar_FH_weights_all_good_subset_permutations.mat : results of permutecode_allsubsets.m for FH weights
  th_polar_all1weights_all_good_subset_permutations.mat : results of permutecode_allsubsets.m for all 1 weights

aaindex/ : directory with Japanese amino acid property database and procedures to check them
  aaindex1 : database of amino acid properties from http://www.genome.jp/aaindex/
  aaindex2 : matrices from http://www.genome.jp/aaindex/
  aaindex/BatchCheck.m : checking all amino-acid values from aaindex1 
  aaindex/BatchCheckAnalysis.m : analysis of results produced by BatchCheck.m
  and more

lib/ : this directory contains lots of helper functions, it has to be part of the MATLAB path
  lib/CreateQAP.m : exports QAP problem to solver-readable format in subdirectory solver
  lib/displaygcode.m :  displays gen code given a permutation and a set of 20 amino acid values
  lib/dsxy2figxy.m : converts data coordinates into figure coordinates
  lib/geneticcode.m : clears workspace and reads in necessary data in order to play around and perform optimizations
  lib/invertp.m : inverts a permutation
  lib/istransit.m : helper function for geneticcode.m 
  lib/mypdist.m : helper function for geneticcode.m
  lib/mypdistweights.m : helper function for geneticcode.m
  lib/nonzeropermute.m : permutes non-zero entries of a vector
  lib/permutecode_random.m : generates random permutations and computes scores
  lib/permutecode_subsets.m : checks all permutations among fixed subsets (specified in the code)
  lib/randfixperm.m : samples a random permutation with some fixed positions
  lib/randperm.m : samples a random permutation
  lib/timing.m : timing and status functions for loops

solver/ : contains non-MATLAB solvers
  gqapd.tar.gz : original source package of GRASP heuristic
  gqapd : compiled version, might have to be recompiled depending on the platform
  gqapd_dir/ : contains FORTRAN files and Makefile for GRASP heuristic
    driver_mine.f : slightly modified original file. Output repeated at the end, so that only the last two lines 
                    have to be read from MATLAB
  mysimqap : compiled simulated annealing heuristic, might have to be recompiled depending on the platform
  qapsim_mine.f : slightly modified FORTRAN code for simulated anealing from QAPLIB
  qapbb.f : branch & bound FORTRAN solver from QAPLIB http://www.seas.upenn.edu/qaplib/codes.html
  qapbb : compiled version, might have to be recompiled depending on the platform
  had12.dat : small test data set (for testing and to illustrate the expected format)

molstruct/ : contains code for computing molecular distances between amino-acids
  /amino : contains .mol files for all 20 amino acids
  /results : contains .png pictures of all pairwise comparisons
  *.out : contains the isomorphisms to the maximal common subgraph of all pairwise comparisons
  *.mol : contains the maximal common subgraphs of all pairwise comparisons
  struc.py : Python program to actually compute the pairwise distances
  	     in order to find the maximal common subgraphs, it runs the Java SMSD from 
	     http://www.ebi.ac.uk/thornton-srv/software/SMSD/
	     which needs to be downloaded and installed separately
