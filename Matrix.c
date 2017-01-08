//
//  Matrix.c
//  SingleChainAverage
//
//  Created by Linling Miao on 6/8/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#include "Matrix.h"
#include "Parameters.h"

void updateMatrix(){
    regularRPY();
    
    step++;

    //printMatrix(1);
    
    //DiagComponents();
    
    //printIteration();
}

void regularRPY(){
    int nn = 3*N;
    
//#pragma omp parallel for private(ii) schedule(static)
    for(int ii = 0; ii<N; ++ii){
        for(int jj = 0; jj<N; ++jj){
            if(ii == jj){
                matrix->DiffMatrix[ii*3*nn+jj*3] += 1.0;
                matrix->DiffMatrix[(ii*3+1)*nn+jj*3+1] += 1.0;
                matrix->DiffMatrix[(ii*3+2)*nn+jj*3+2] += 1.0;
            } else{
                double dx = Chain[ii].rx - Chain[jj].rx;
                double dy = Chain[ii].ry - Chain[jj].ry;
                double dz = Chain[ii].rz - Chain[jj].rz;
                double rr = dx*dx+dy*dy+dz*dz;
                double r = sqrt(rr);
                
                double mm1,mm2;
                if(r>=2.0){
                    mm1 = 3.0/(4.0*r)*(1.0+2.0/(3.0*rr));
                    mm2 = 3.0/(4.0*r)*(1.0-2.0/rr)/rr;
                }
                else{
                    mm1 = 1.0-9.0*r/32.0;
                    mm2 = 3.0/(32.0*r);
                }
                
                matrix->DiffMatrix[ii*3*nn+jj*3] += mm1+mm2*dx*dx;
                matrix->DiffMatrix[ii*3*nn+jj*3+1] += mm2*dx*dy;
                matrix->DiffMatrix[ii*3*nn+jj*3+2] += mm2*dx*dz;
                matrix->DiffMatrix[(ii*3+1)*nn+jj*3] += mm2*dy*dx;
                matrix->DiffMatrix[(ii*3+1)*nn+jj*3+1] += mm1+mm2*dy*dy;
                matrix->DiffMatrix[(ii*3+1)*nn+jj*3+2] += mm2*dy*dz;
                matrix->DiffMatrix[(ii*3+2)*nn+jj*3] += mm2*dz*dx;
                matrix->DiffMatrix[(ii*3+2)*nn+jj*3+1] += mm2*dz*dy;
                matrix->DiffMatrix[(ii*3+2)*nn+jj*3+2] += mm1+mm2*dz*dz;
                
            }
        }
    }
}

/*void regularRPYcheck(){
    int nn = 3*N;
    float *mm = calloc(nn*nn,sizeof(float));
    
    int ii;
//#pragma omp parallel for private(ii) schedule(dynamic)
    for(ii = 0; ii<N; ++ii){
        int jj;
        for(jj = 0; jj<N; ++jj){
            if(ii == jj){
                mm[ii*3*nn+jj*3] += 1.0;
                mm[(ii*3+1)*nn+jj*3+1] += 1.0;
                mm[(ii*3+2)*nn+jj*3+2] += 1.0;
            } else{
                double dx = Chain[ii].rx - Chain[jj].rx;
                double dy = Chain[ii].ry - Chain[jj].ry;
                double dz = Chain[ii].rz - Chain[jj].rz;
                double rr = dx*dx+dy*dy+dz*dz;
                double r = sqrt(rr);
                
                double mm1,mm2;
                if(r>=2.0){
                    mm1 = 3.0/(4.0*r)*(1.0+2.0/(3.0*rr));
                    mm2 = 3.0/(4.0*r)*(1.0-2.0/rr)/rr;
                }
                else{
                    mm1 = 1.0-9.0*r/32.0;
                    mm2 = 3.0/(32.0*r);
                }
                
                mm[ii*3*nn+jj*3] += mm1+mm2*dx*dx;
                mm[ii*3*nn+jj*3+1] += mm2*dx*dy;
                mm[ii*3*nn+jj*3+2] += mm2*dx*dz;
                mm[(ii*3+1)*nn+jj*3] += mm2*dy*dx;
                mm[(ii*3+1)*nn+jj*3+1] += mm1+mm2*dy*dy;
                mm[(ii*3+1)*nn+jj*3+2] += mm2*dy*dz;
                mm[(ii*3+2)*nn+jj*3] += mm2*dz*dx;
                mm[(ii*3+2)*nn+jj*3+1] += mm2*dz*dy;
                mm[(ii*3+2)*nn+jj*3+2] += mm1+mm2*dz*dz;
            }
        }
    }
    
    int error = choleskycheck(mm);
    
    if(error) {
        printErrorTraj();
        printErrorMatrix(mm);
    }
    
    int i;
    for(i=0; i<nn*nn; ++i){
        matrix->DiffMatrix[i] += mm[i];
    }
    
    free(mm);
}*/

/*void EwaldSumRPY(){
    int nn = 3*N;
    double vol = L*L*L;
    double alpha = 6.0/L;
    
    int ii;
//#pragma omp parallel for private(ii) schedule(dynamic)
    for(ii = 0; ii<N; ++ii){
        int jj;
        for(jj = 0; jj<N; ++jj){
            double dx = Chain[ii].rx - Chain[jj].rx;
            double dy = Chain[ii].ry - Chain[jj].ry;
            double dz = Chain[ii].rz - Chain[jj].rz;
            
            //Reciprocal Part
            int kx;
            for(kx = -kmax; kx<kmax+1; ++kx){
                double rkx = 2.0*M_PI*kx/L;
                int ky;
                for(ky = -kmax; ky<kmax+1; ++ky){
                    double rky = 2.0*M_PI*ky/L;
                    int kz;
                    for(kz = -kmax; kz<kmax+1; ++kz){
                        double rkz = 2.0*M_PI*kz/L;
                        double kk = kx*kx+ky*ky+kz*kz;
                        if(kk != 0){
                            double mm4 = cos(rkx*dx+rky*dy+rkz*dz);
                            matrix->DiffMatrix[ii*3*nn+jj*3] += mm4*m2[kx+kmax][ky+kmax][kz+kmax][0]/vol;
                            matrix->DiffMatrix[ii*3*nn+jj*3+1] += mm4*m2[kx+kmax][ky+kmax][kz+kmax][1]/vol;
                            matrix->DiffMatrix[ii*3*nn+jj*3+2] += mm4*m2[kx+kmax][ky+kmax][kz+kmax][2]/vol;
                            matrix->DiffMatrix[(ii*3+1)*nn+jj*3] += mm4*m2[kx+kmax][ky+kmax][kz+kmax][3]/vol;
                            matrix->DiffMatrix[(ii*3+1)*nn+jj*3+1] += mm4*m2[kx+kmax][ky+kmax][kz+kmax][4]/vol;
                            matrix->DiffMatrix[(ii*3+1)*nn+jj*3+2] += mm4*m2[kx+kmax][ky+kmax][kz+kmax][5]/vol;
                            matrix->DiffMatrix[(ii*3+2)*nn+jj*3] += mm4*m2[kx+kmax][ky+kmax][kz+kmax][6]/vol;
                            matrix->DiffMatrix[(ii*3+2)*nn+jj*3+1] += mm4*m2[kx+kmax][ky+kmax][kz+kmax][7]/vol;
                            matrix->DiffMatrix[(ii*3+2)*nn+jj*3+2] += mm4*m2[kx+kmax][ky+kmax][kz+kmax][8]/vol;
                        }
                    }
                }
            }
            
            //Real Part (Cutoff 0.5L)
            dx -= round(dx/L)*L; dy -= round(dy/L)*L; dz -= round(dz/L)*L;
            double rr = dx*dx+dy*dy+dz*dz;
            double r = sqrt(rr);
            
            if(ii == jj){
                double diag = 1.0-6.0/sqrt(M_PI)*alpha+40.0/3.0/sqrt(M_PI)*alpha*alpha*alpha;
                matrix->DiffMatrix[ii*3*nn+jj*3] += diag;
                matrix->DiffMatrix[(ii*3+1)*nn+jj*3+1] += diag;
                matrix->DiffMatrix[(ii*3+2)*nn+jj*3+2] += diag;
            }
            else{
                double c1 = 0.75/r+0.5*pow(r,-3.0);
                double c2 = 4.0*pow(alpha,7.0)*pow(r,4.0)+3.0*pow(alpha,3.0)*rr-20.0*pow(alpha,5.0)*rr-4.5*alpha+14.0*pow(alpha,3.0)+alpha/rr;
                double c3 = 0.75/r-1.5*pow(r,-3.0);
                double c4 = -4.0*pow(alpha,7.0)*pow(r,4.0)-3.0*pow(alpha,3.0)*rr+16.0*pow(alpha,5.0)*rr+1.5*alpha-2.0*pow(alpha,3.0)-3.0*alpha/rr;
                
                double mm1,mm2;
                if(r>=2.0){
                    mm1 = c1*erfc(alpha*r)+c2*exp(-alpha*alpha*rr)/sqrt(M_PI);
                    mm2 = c3*erfc(alpha*r)+c4*exp(-alpha*alpha*rr)/sqrt(M_PI);
                }
                else{
                    mm1 = c1*erfc(alpha*r)+c2*exp(-alpha*alpha*rr)/sqrt(M_PI)+(1.0-9.0*r/32.0-3.0/4.0/r-0.5*pow(r,-3.0));
                    mm2 = c3*erfc(alpha*r)+c4*exp(-alpha*alpha*rr)/sqrt(M_PI)+(3.0*r/32.0-3.0/4.0/r+1.5*pow(r,-3.0));
                }
                matrix->DiffMatrix[ii*3*nn+jj*3] += mm1+mm2*dx*dx/rr;
                matrix->DiffMatrix[ii*3*nn+jj*3+1] += mm2*dx*dy/rr;
                matrix->DiffMatrix[ii*3*nn+jj*3+2] += mm2*dx*dz/rr;
                matrix->DiffMatrix[(ii*3+1)*nn+jj*3] += mm2*dy*dx/rr;
                matrix->DiffMatrix[(ii*3+1)*nn+jj*3+1] += mm1+mm2*dy*dy/rr;
                matrix->DiffMatrix[(ii*3+1)*nn+jj*3+2] += mm2*dy*dz/rr;
                matrix->DiffMatrix[(ii*3+2)*nn+jj*3] += mm2*dz*dx/rr;
                matrix->DiffMatrix[(ii*3+2)*nn+jj*3+1] += mm2*dz*dy/rr;
                matrix->DiffMatrix[(ii*3+2)*nn+jj*3+2] += mm1+mm2*dz*dz/rr;
            }
        }
    }
}*/

void DiagComponents(){
    FILE *comp;
    char* str = malloc(50*sizeof(char));
    sprintf(str,"%s/DMatrix_N%d_hsr%s_%d.txt",directory,N,flowratec,iteration);
    
    int nn = 3*N;
    
    int i;
    float **rematrix = calloc(N,sizeof(rematrix[0]));
    for(i = 0; i<N; ++i){
        rematrix[i] = calloc(3,sizeof(rematrix[0][0]));
    }
    
    for(i = 0; i<N; ++i){
        int j;
        for(j = 0; j<3; ++j){
            //printf("index:%d comp:%f\n",j*nn+i*3+i,matrix[j*nn+i*3+j]);
            rematrix[i][j] = matrix->DiffMatrix[j*nn+i*3+j]/step;
        }
    }
    
    comp = fopen(str,"w");
    
    for(i = 0; i<N; ++i){
        fprintf(comp,"%f %f %f\n",rematrix[i][0],rematrix[i][1],rematrix[i][2]);
    }
    
    fclose(comp);
    
    for(i = 0; i<N; ++i){
        free(rematrix[i]);
    }
    free(rematrix);
    
    free(str);
}

void printMatrix(int isBinary){
    int nn = 3*N;
    int mtotal = 9*N*N;
    
    
    if(isBinary){
        FILE *DiffusionB;
        char* str = malloc(50*sizeof(char));
        //sprintf(str,"%s/Matrix_N%d_hsr%s_%d.bin",directory,N,hsrchar,iteration);
        sprintf(str,"%s/Matrix_N%d_hsr%s.bin",directory,N,flowratec);//MARK: name
        DiffusionB = fopen(str, "wb");
        //Header
        fwrite(&N,sizeof(uint32_t),1,DiffusionB);
        
        //Matrix
        if(AVGmode==0 || HImode==0){
            int i;
            for(i = 0; i<mtotal; ++i){
                float value = matrix->DiffMatrix[i]/step;
                fwrite(&value,sizeof(float),1,DiffusionB);
            }
        } else if(AVGmode==1){
            int i;
            for(i = 0; i<mtotal; ++i){
                float value = AVGfactor*(matrix->DiffMatrix[i]/step) + (1 - AVGfactor)*matrix->DiffMatrixAvg[i];
                fwrite(&value,sizeof(float),1,DiffusionB);
            }
        }
        
        fclose(DiffusionB);
        
        free(str);
        
    } else{
        FILE *Diffusion;
        char* str = malloc(50*sizeof(char));
        sprintf(str,"%s/DMatrix_N%d_hsr%s_%d.txt",directory,N,flowratec,iteration);
        
        Diffusion = fopen(str,"w");
        
        //fprintf(Diffusion,"step = %d\n",step);
        
        int i;
        for(i = 0; i<nn; ++i){
            int j;
            for(j = 0; j<nn; ++j){
                float value = matrix->DiffMatrix[i*nn+j]/step;
                fprintf(Diffusion,"%f ",value);
            }
            fprintf(Diffusion,"\n");
        }
        
        fclose(Diffusion);
        
        free(str);
        
    }
}

void printIteration(){
    char* str = malloc(100*sizeof(char));
    sprintf(str,"%s/iteration_N%d_hsr%s.bin",directory,N,flowratec);
    FILE *iterat = fopen(str,"wb");
    fwrite(&iteration,sizeof(uint16_t),1,iterat);
    fclose(iterat);
    free(str);
}


/*int choleskycheck(float *mm){
    int n = 3*N;
    int error = 0;
    
    float *lm = calloc(n*n,sizeof(float));
    
    int i;
    for (i = 0; i<n; i++){
        int j;
        for (j = 0; j<(i+1); j++){
            double s = 0;
            int k;
            for (k = 0; k < j; k++){
                s += lm[i*n+k]*lm[j*n+k];
            }
            
            if(i == j){
                if((mm[i*n+i]-s)<=0){
                    error = 1;
                    i = n;
                    j = i+1;
                    break;
                }
                lm[i*n+j] = sqrt(mm[i*n+i]-s);
            } else{
                lm[i*n+j] = (mm[i*n+j]-s)/lm[j*n+j];
            }
        }
    }
    
    free(lm);
    
    return error;
}*/

/*void printErrorTraj(){
    FILE *Trajectory;
    Trajectory = fopen("./output/ErrorTrajectory.txt","a");
    //fprintf(Trajectory, "%d\n-1\n", N);
    int i;
    for(i = 0; i<N; ++i){
        fprintf(Trajectory, "%lf %lf %lf\n", Chain[i].rx, Chain[i].ry, Chain[i].rz);
    }
    fclose(Trajectory);
}*/

/*void printErrorMatrix(float *mm){
    int nn = 3*N;
    
    FILE *Diffusion = fopen("./output/ErrorMatrix.txt","w");
    
    int i;
    for(i = 0; i<nn; ++i){
        int j;
        for(j = 0; j<nn; ++j){
            float value = mm[i*nn+j];
            fprintf(Diffusion,"%f ",value);
        }
        fprintf(Diffusion,"\n");
    }
    
    fclose(Diffusion);

}*/




