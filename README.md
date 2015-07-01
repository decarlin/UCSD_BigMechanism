# UCSD_BigMechanism

This is the repo for development of the ucsd big mechanism contribution.

The basic heat diffusion for a query is executed in automateHeatKernel.py script.  See the diffused.txt rule in examples/Makefile for execution syntax.  By default, the output of this script is filtered.sif (which contains the top N=30 by default most relevent terms in it) and diffused.txt, which contains all first order neighbors in the network and the diffused query score for each.