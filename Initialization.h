//
//  Initialization.h
//  SingleChainAverage
//
//  Created by Linling Miao on 6/8/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#ifndef Initialization_h
#define Initialization_h

void Initialization(int argc, char * argv[]);

/* Function used to read intial traj from a file */
void readChains();

/* Read the binary iteration file */
void readIteration();

/* Read pre-averaged matrix file */
void readMatrixAvg();
void readZimmMatrix();

void generateOutput();

long initRan();

/* Initialize the parallel environment */
void initParaEvir();

void initChain();

/* Initial ring structure, add it after the function initChain() */
void initRingStructure();

float ran1(long *idum);

void printInitial();

#endif /* Initialization_h */
