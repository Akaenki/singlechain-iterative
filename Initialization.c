//
//  Initialization.c
//  SingleChainAverage
//
//  Created by Linling Miao on 6/8/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#include "Initialization.h"
#include "Parameters.h"
#include "GeyerWinterCalculation.h"
#include "Properties.h"
#include "util.h"

void Initialization(int argc, char * argv[]){
    ParseInput(argc, argv);
    
    matrix = (Matrix_t*)malloc(sizeof(Matrix_t));
    matrix->DiffMatrix = (float*)calloc(9*N*N, sizeof(float));
    gwParm = (gwParm_t*)malloc(sizeof(gwParm_t));
    gwParm->C_i = (float*)calloc(3*N, sizeof(float));
    
    trajStep = 1000;
    
    ReadInitFromFile = false; //True for all special start (not necessary read from a file)
    
    CoilStart = true; //only if ReadInitFromFile is true, false for fully stretched start
    
    GeyerWinter = false;
    
    directory = "./output";
    
    if(HImode){
        if(HImode==2){
            readZimmMatrix();
            iteration = 1;
        } else if(HImode==1){
            readIteration();
            readMatrixAvg();
            iteration++;
        }
        
        if(GeyerWinter) GeyerWinterCalculation();
        else Cholesky();
        printf("Now Running iteration:%d ...\n",iteration);
    } else{
        iteration = 1;
        printf("Now Running iteration:%d ...\n",iteration);
    }
    
    outputName = malloc(100*sizeof(char));
    sprintf(outputName,"%s/sc_N%d_hsr%s.txt",directory,N,flowratec);
    trajName = malloc(100*sizeof(char));
    sprintf(trajName,"%s/sc_N%d_hsr%s_%d.xyz",directory,N,flowratec,iteration);
////Don't Change These//////////////////////////////////////////////////////////////////////////////////////
    CheckOverlap = true;
    
    printNonBinary = false;
    
    EwaldSum = false;
    
    if(EwaldSum) usePBC = true;
    else usePBC = false;
    
    MSD_CoM = 0.0; step = 0;
    
    ext.x = 0.0; ext.y = 0.0, ext.z = 0.0;
////End of Changable Variables//////////////////////////////////////////////////////////////////////////////
    idum = malloc(sizeof(long));
    *idum = -1;
    ran1(idum);
    *idum = -1*initRan();
    
    initParaEvir();
    
    if(ReadInitFromFile) readChains();
    else initChain();
    
    //generateOutput();

    //printInitial();
    
    AVGmode = 0;
    /* Average Methods: *
     * use AVGmode = 0 if don't want to average the Matrix calculated in current iteration with any of the previous matrices
     * use AVGmode = 1 if want to averaged the current matrix with the previous matrices
            * The averaging algorism:
                    ** mm'[n] = AVGfactor * mm[n] + (1 - AVGfactor) * mm'[n-1] **
            * Specify AVGfactor below or use initial value 0.5 *
     */
    
    AVGfactor = 0.5;
}

void readChains(){
    Chain = calloc(N,sizeof(Chain_t));
    
    if(CoilStart){
        FILE *trajec;
        char* str = malloc(20*sizeof(char));
        sprintf(str,"Coil_N%d.xyz",N);
        trajec = fopen(str,"r");
        for(int i = 0; i<N; ++i){
            fscanf(trajec, "A %lf %lf %lf\n", &Chain[i].rx, &Chain[i].ry, &Chain[i].rz);
        }
        fclose(trajec);
        free(str);
    }
    else{
        Chain[0].rx = 0.0;
        Chain[0].ry = 0.0;
        Chain[0].rz = 0.0;
        
        for(int i = 1; i<N; ++i){
            Chain[i].rx = Chain[i-1].rx + 2.05;
            Chain[i].ry = 0.0;
            Chain[i].rz = 0.0;
        }
    }
}

void generateOutput(){
    FILE *outputfile;
    outputfile = fopen(outputName, "w");
    fprintf(outputfile, "epsilon = %lf\n", EPSILON);
    fprintf(outputfile, "kappa = %lf\n", KAPPA);
    fprintf(outputfile, "N = %d\n", N);
    fprintf(outputfile, "dt = %lf\n", DT);
    fprintf(outputfile, "tmax = %d\n", TMAX);
    fprintf(outputfile, "\n\n");
    fprintf(outputfile, "SEED %ld\n", *idum);
    fclose(outputfile);
}

void initParaEvir(){
    bool isOpenMP = false;
#ifdef _OPENMP
    isOpenMP = true;
    omp_set_num_threads(numThread);
#endif
    if(isOpenMP){
        openblas_set_num_threads(1);
    } else{
 
    }
    openblas_set_num_threads(numThread);
}

void initChain(){
    Chain = calloc(N,sizeof(Chain_t));
    
    Chain[0].rx = 0.0;
    Chain[0].ry = 0.0;
    Chain[0].rz = 0.0;
    
    int i;
    for(i = 1; i<N; ++i){
        int test = 0;
        while(test==0){
            test = 1;
            double theta = ran1(idum)*2.0*3.14158;
            double phi = acos(2.0*ran1(idum)-1.0);
            Chain[i].rx = Chain[i-1].rx + 2.05*cos(theta)*sin(phi);
            Chain[i].ry = Chain[i-1].ry + 2.05*sin(theta)*sin(phi);
            Chain[i].rz = Chain[i-1].rz + 2.05*cos(phi);
            
            if(usePBC){
                Chain[i].rx -= round(Chain[i].rx/L)*L;
                Chain[i].ry -= round(Chain[i].ry/L)*L;
                Chain[i].rz -= round(Chain[i].rz/L)*L;
            }
            
            if(CheckOverlap){
                int j;
                for(j = 0; j<i; ++j){
                    double dx = Chain[j].rx - Chain[i].rx;
                    double dy = Chain[j].ry - Chain[i].ry;
                    double dz = Chain[j].rz - Chain[i].rz;
                    
                    if(usePBC){
                        dx -= round(dx/L)*L;
                        dy -= round(dy/L)*L;
                        dz -= round(dz/L)*L;
                    }
                    
                    if(dx*dx+dy*dy+dz*dz<4.0){
                        test = 0;
                    }
                }
            }
        }
    }
}

void initRingStructure(){
    int test = 0;
    double instKappa = 0.5;
    
    while(test==0){
        for(int i = 0; i<N; ++i){
            int j;
            if(i==0) j = N-1;
            else j = i-1;
            
            double dx = Chain[i].rx - Chain[j].rx;
            double dy = Chain[i].ry - Chain[j].ry;
            double dz = Chain[i].rz - Chain[j].rz;
            
            double r = sqrt(dx*dx+dy*dy+dz*dz);
            double Fs = -instKappa*(r-2.0);
            
            Chain[i].fx += Fs*dx/r;
            Chain[i].fy += Fs*dy/r;
            Chain[i].fz += Fs*dz/r;
            Chain[j].fx -= Fs*dx/r;
            Chain[j].fy -= Fs*dy/r;
            Chain[j].fz -= Fs*dz/r;
            
        }
        
        for(int i = 0; i<N; ++i){
            Chain[i].rx += DT*Chain[i].fx;
            Chain[i].ry += DT*Chain[i].fy;
            Chain[i].rz += DT*Chain[i].fz;
        }
        
        double dx = Chain[N-1].rx - Chain[0].rx;
        double dy = Chain[N-1].ry - Chain[0].ry;
        double dz = Chain[N-1].rz - Chain[0].rz;
        
        if(dx*dx+dy*dy+dz*dz <= 5.0){
            test = 1;
        }
    }
}

void readIteration(){
    char* str = malloc(100*sizeof(char));
    sprintf(str,"%s/iteration_N%d_hsr%s.bin",directory,N,flowratec);
    FILE *iterat = fopen(str,"rb");
    fread(&iteration,sizeof(uint16_t),1,iterat);
    fclose(iterat);
    free(str);
}

void readMatrixAvg(){
    FILE *Matrix;
    char* str = malloc(50*sizeof(char));
    //sprintf(str,"%s/Matrix_N%d_hsr%s_%d.bin",directory,N,flowratec,iteration);
    sprintf(str,"%s/Matrix_N%d_hsr%s.bin",directory,N,flowratec);
    Matrix = fopen(str,"rb");
    
    uint32_t Nc;
    
    /* Header */
    fread(&Nc,sizeof(uint32_t),1,Matrix);
    
    uint32_t mtotal = 9*Nc*Nc;
    
    matrix->DiffMatrixAvg = (float*)malloc(mtotal*sizeof(float));
    
    /* Matrix */
    fread(matrix->DiffMatrixAvg,sizeof(float),mtotal,Matrix);
    
    fclose(Matrix);
    
    free(str);
}

void readZimmMatrix(){
    FILE *Matrix;
    char* str = malloc(50*sizeof(char));
    if(CoilStart){
        sprintf(str,"Zimm_N%d.bin",N);//MARK:name
    } else{
        sprintf(str,"Stretch_N%d.bin",N);//MARK:name
    }
    
    Matrix = fopen(str,"rb");
    
    uint32_t Nc;
    
    /* Header */
    fread(&Nc,sizeof(uint32_t),1,Matrix);
    
    uint32_t mtotal = 9*Nc*Nc;
    
    matrix->DiffMatrixAvg = (float*)malloc(mtotal*sizeof(float));
    
    /* Matrix */
    fread(matrix->DiffMatrixAvg,sizeof(float),mtotal,Matrix);
    
    fclose(Matrix);
    
    free(str);
}

void printInitial(){
    FILE *Trajectory;
    Trajectory = fopen(trajName,"a");
    fprintf(Trajectory, "%d\n-1\n", N);
    for(int i = 0; i<N; ++i){
        fprintf(Trajectory, "A %lf %lf %lf\n", Chain[i].rx, Chain[i].ry, Chain[i].rz);
    }
    fclose(Trajectory);
}


long initRan(){
    unsigned long a = clock();
    unsigned long b = time(NULL);
    unsigned long c = getpid();
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c%1000000000+1000000000;
}

float ran1(long *idum){
    long j;
    long k;
    static long idum2 = 123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    
    if(*idum <= 0){
        if(-(*idum)<1) *idum=1;
        else *idum = -(*idum);
        idum2 = (*idum);
        for(j=NTAB+7;j>=0;--j)
        {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if(*idum<0) *idum+=IM1;
            if(j<NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if(*idum<0) *idum += IM1;
    k=idum2/IQ2;
    if(*idum<0) idum2+= IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if(iy<1) iy += IMM1;
    if((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}



