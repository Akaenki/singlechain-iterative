//
//  GeyerWinterCalculation.c
//  SingleChainAverage
//
//  Created by Linling Miao on 6/8/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#include "GeyerWinterCalculation.h"
#include "Parameters.h"

void GeyerWinterCalculation(){
    double eps = 0.0;
    int nt = 3*N;/*dimension*/
    
    int i;
    for(i = 0; i<nt*nt; ++i){
        eps += matrix->DiffMatrixAvg[i];
    }
    
    eps -= nt; /*exclude diagonal(self) components*/
    
    eps /= nt*nt - nt;
    
    gwParm->beta_ij = (1-sqrt(1-((nt-1)*eps*eps-(nt-2)*eps)))/((nt-1)*eps*eps-(nt-2)*eps);
    //printf("beta=%f\n", gwParm->beta_ij);
    
    gwParm->C_i = calloc(nt,sizeof(float));
    
    for(i = 0; i<nt; ++i){
        int j;
        for(j = 0; j<nt; ++j){
            gwParm->C_i[i] += matrix->DiffMatrixAvg[i*nt+j]*matrix->DiffMatrixAvg[i*nt+j];
        }
    }
    
    for(i = 0; i<nt; ++i){
        gwParm->C_i[i] -= 1.0;
        gwParm->C_i[i] = gwParm->C_i[i]*gwParm->beta_ij*gwParm->beta_ij+1.0;
        gwParm->C_i[i] = 1.0/sqrt(gwParm->C_i[i]);
        //printf("C%d=%f\n",i,gwParm->C_i[i]);
    }
    
    //??
    /*D_ij = calloc(nt*nt,sizeof(float));
    for(i = 0; i<nt*nt; ++i){
        D_ij[i] = matrix->DiffMatrixAvg[i];
    }
    
    int ii;
    for(ii = 0; ii<nt; ++ii){
        D_ij[ii*nt+ii] /= gwParm->beta_ij;
    }*/
}

void Cholesky(){
    int n = 3*N;
    int error = 0;
    
    Ddecomp = calloc(n*n,sizeof(float));
    
    for (int i = 0; i<n; i++){
        for (int j = 0; j<(i+1); j++){
            double s = 0;
            for (int k = 0; k < j; k++){
                s += Ddecomp[i*n+k]*Ddecomp[j*n+k];
            }
            
            if(i == j){
                if((matrix->DiffMatrixAvg[i*n+i]-s)<=0){
                    error = 1;
                }
                Ddecomp[i*n+j] = sqrt(matrix->DiffMatrixAvg[i*n+i]-s);
            } else{
                Ddecomp[i*n+j] = (matrix->DiffMatrixAvg[i*n+j]-s)/Ddecomp[j*n+j];
            }
        }
    }
    
    if(error){
        printf("ERROR: Matrix Not Positive Definete\nFollowing Run Will Be Free Draining\n");
        for(int i = 0; i<n; ++i){
            for(int j = 0; j<n; ++j){
                Ddecomp[i*n+j] = 0.0;
                if(i == j) Ddecomp[i*n+j] = 1.0;
            }
        }
    }
}