//
//  Parameters.h
//  SingleChainAverage
//
//  Created by Linling Miao on 6/8/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#ifndef Parameters_h
#define Parameters_h

#if defined(HYBRID_PARALLEL)
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <stdbool.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cblas.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define NDIV (1+IMM1/NTAB)
#define EP 0.9624
#define THETA0 31.7

#ifndef M_PI
#define M_PI 3.1415926535
#endif


/////////////////////////////////////////////////////////////////////////////
// Types
/////////////////////////////////////////////////////////////////////////////

typedef struct {
    double x,y,z;
} Vector3D_t;

typedef struct {
    double rx,ry,rz;
    double fx,fy,fz;
} Chain_t;

/* Mobility Matrix */
typedef struct {
    float *DiffMatrix; //Matrix calculated in current iteration
    float *DiffMatrixAvg; //Averaged matrix read from the binary file
} Matrix_t;

/* Geyer-Winter(TEA) Parameters */
typedef struct {
    double beta_ij;
    float* C_i;
} gwParm_t;

/////////////////////////////////////////////////////////////////////////////
// Gloable
/////////////////////////////////////////////////////////////////////////////

long *idum;

uint16_t iteration;

int HImode; //0 for FD, 1 for HI, 2 for HI as first iteration
int AVGmode;
double AVGfactor;

uint32_t N;
int numThread;

int flowType;
double flowRate;

char* flowratec;

Vector3D_t *R_CoM;
Vector3D_t CoMstored;
double MSD_CoM;
Vector3D_t ext;

Chain_t* Chain;

float *Ddecomp;//Cholesky Decomposition
float ****m2;

gwParm_t* gwParm;

//Matrix
Matrix_t* matrix;

/////////////////////////////////////////////////////////////////////////////
// Others
/////////////////////////////////////////////////////////////////////////////

int trajStep;
int step;

double Rg,Ree;

char* outputName;
char* trajName;
char* directory;

bool ReadInitFromFile;
bool CheckOverlap;
bool GeyerWinter;
bool EwaldSum;
bool usePBC;
bool printNonBinary;
bool CoilStart;

#endif /* Parameters_h */
