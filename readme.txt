MATLAB code for genetic code optimizations. 
These tools have been used to obtain the numerical results and plots in a forthcoming article
by Harry Buhrman, Peter van der Gulik, Dave Spijer, Gunnar Klau, Leen Stougie, Christian Schaffner

created: Oct 21, 2011
by Christian Schaffner, c.schaffner@uva.nl


This directory contains the following files:
* geneticcode.m : clears workspace and reads in necessary data in order to play around and perform optimizations
* VariousGraphs.m : main program module to set parameters and run computations
* makegraph.m : draws histograms from computed values

* output/ : output directory where .fig and .pdf files are saved to

* aaindex1 : database of amino acid properties from http://www.genome.jp/aaindex/
* permutecode_random.m : generates random permutations and computes scores
* randfixperm.m : samples a random permutation with some fixed positions
* istransit.m : helper function for geneticcode.m 
* mypdist.m : helper function for geneticcode.m
* mypdistweights.m : helper function for geneticcode.m
* randperm.m : samples a random permutation
