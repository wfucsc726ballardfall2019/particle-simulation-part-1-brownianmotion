This problem is borrowed from UC Berkeley's Parallel Computing class (CS 267).

Here is the most recent website with full description (and helpful resources):
https://sites.google.com/lbl.gov/cs267-spr2019/hw-2
This assignment focuses on Part 1 (sequential and OpenMP).

The goal of the assignment is to design and implement efficient sequential code for particle simulation and then design and implement efficient parallel code using OpenMP.  You can complete this by rewriting sequential.cpp and openmp.cpp codes.  The rest of the files can be used as is (but you'll probably want to edit them as you develop).  The Makefile should work to compile on the DEAC cluster, and the batch scripts auto-deac-serial and auto-deac-openmp should work as benchmarking scripts (submit using "sbatch auto-deac-serial" on a head node, for example).  The autograder will help you analyze your results.  All executables can be run with a help option (./serial -h) to see usage instructions.

I've also included executables to visualize your simulation in the visualize directory.  This needs to be done on your laptop (not on the cluster).  There are two executables, one for Mac and one for Windows.  I make no guarantees, but if you specify the output file from your simulation as a command line argument, you should be able to see what happens to your particles as they step through time (correctness bugs usually show up in the visualization).  For Mac users, you might have to download a framework called SDL, see the readme in the visualize directory.

Testing git