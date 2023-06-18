# WangLandauMC
This repository contains a barebones Monte Carlo (MC) code written in C++ that performance Wang-Landau (WL) sampling to calculate the density of states. Included in this repository are versions to handle single component and binary Lennard-Jones systems (for both 1D and 2D density of states. This code performs parallel sampling to speed convegence, although I note this is not a parallel replica-exchange implementation (see https://journals.aps.org/pre/abstract/10.1103/PhysRevE.90.023302 for replica-exchange WL algorithm).  This code was initially developed during research into creating a hybrid WLMC/statistical temperature MD (STMD) approach, allowing STMD to be used to calculated the 2D density of states  (and thus allowing examination of the density landscape; this was of particular interest given our prior work looking at phase transitions in lipid systems, see https://pubs.aip.org/aip/jcp/article-abstract/139/5/054505/351963/Examining-the-phase-transition-behavior-of?redirectedFrom=fulltext ).  This code was originally written June 2013; I found it in my code archive and felt it was worthwhile uploading so I don't lose it.   

## For more info on the algorithm, see:

Efficient, Multiple-Range Random Walk Algorithm to Calculate the Density of States

Fugao Wang and D. P. Landau
Phys. Rev. Lett. 86, 2050 – Published 5 March 2001

https://doi.org/10.1103%2FPhysRevLett.86.2050
