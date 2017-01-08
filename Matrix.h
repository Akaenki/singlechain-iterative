 //
//  Matrix.h
//  SingleChainAverage
//
//  Created by Linling Miao on 6/8/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#ifndef Matrix_h
#define Matrix_h

void updateMatrix();

void regularRPY();
//void EwaldSumRPY();

/* print the matrix in binary or .txt file */
/* In Binary file the head contains
void printMatrix(int isBinary);

/* print the binary iteration file 
 * The binary matrix contains a 32 bit (4 bytes) header of number of beads
 * The natrix body contains 3Nc*3Nc single precision floating points (4 bytes each)
 * The whole matrix is read row-wise */
void printIteration();

//void regularRPYcheck();
//int choleskycheck(float *mm);
//void printErrorTraj();
//void printErrorMatrix(float *mm);

/* Will print the dagonal components of the mobility matrix */
void DiagComponents();

#endif /* Matrix_h */
