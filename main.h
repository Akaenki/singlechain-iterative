//
//  main.h
//  SingleChainAverage
//
//  Created by Linling Miao on 6/8/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#ifndef main_h
#define main_h

#include "Parameters.h"

int main(int argc, char * argv[]);


Vector3D_t getNID(int i,int j);
void applyPBC();

void initForce();

/* Bonding forces */
void forceBond();
void forceSpringRing();

/* Exclusive volume interations */
void forceLJ();

void updateChain();

/* Calculation CoM */
Vector3D_t CenterOfMass(int i);

void printTrajectory(int t);

float gasdev(long *idum);
long long timer();

#endif /* main_h */
