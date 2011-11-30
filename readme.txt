MATLAB code for genetic code optimizations. 
These tools have been used to obtain the numerical results and plots in a forthcoming article
by Harry Buhrman, Peter van der Gulik, Dave Spijer, Gunnar Klau, Leen Stougie, Christian Schaffner

created: Oct 21, 2011
by Christian Schaffner, c.schaffner@uva.nl

Installation:
* Make sure the lib/ subdirectory is part of the path

This directory contains the following files:
VariousGraphs.m : main program module to set parameters and run computations
makegraph.m : draws histograms from computed values

output/ : output directory where .fig and .pdf files are saved to

aaindex/ : directory with Japanese amino acid property database and procedures to check them
  aaindex1 : database of amino acid properties from http://www.genome.jp/aaindex/
  aaindex2 : matrices from http://www.genome.jp/aaindex/
  aaindex/BatchCheck.m : checking all amino-acid values from aaindex1 
  aaindex/BatchCheckAnalysis.m : analysis of results produced by BatchCheck.m
  and more

lib/ : this directory contains lots of helper functions, it has to be part of the MATLAB path
  lib/geneticcode.m : clears workspace and reads in necessary data in order to play around and perform optimizations
  lib/permutecode_random.m : generates random permutations and computes scores
  lib/permutecode_subsets.m : checks all permutations among fixed subsets (specified in the code)
  lib/CreateQAP.m : exports QAP problem to solver-readable format in subdirectory solver

  lib/randfixperm.m : samples a random permutation with some fixed positions
  lib/istransit.m : helper function for geneticcode.m 
  lib/mypdist.m : helper function for geneticcode.m
  lib/mypdistweights.m : helper function for geneticcode.m
  lib/randperm.m : samples a random permutation
  lib/nonzeropermute.m : permutes non-zero entries of a vector
  lib/timing.m : timing and status functions for loops
  lib/displaygcode.m :  displays gen code given a permutation and a set of 20 amino acid values
  lib/invertp.m : inverts a permutation

solver/ : contains non-MATLAB solvers
  solver/qapbb.f : branch & bound FORTRAN solver from QAPLIB http://www.seas.upenn.edu/qaplib/codes.html
  solver/qapbb : compiled version
  solver/had12.dat : small test data set (for testing and to illustrate the expected format)
