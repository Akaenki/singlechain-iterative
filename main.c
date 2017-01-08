//
//  main.c
//  SingleChainExtensionalAvg
//
//  Created by Linling Miao on 6/14/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#include "main.h"
#include "Initialization.h"
#include "Matrix.h"
#include "Properties.h"

int main(int argc, char * argv[]) {
    Initialization(argc, argv);
    
    int t;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  BEGAIN OF TIME LOOP
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(t = 0; t<TMAX; ++t){
        storeCoM();
        
        initForce();
        
        forceBond();
        forceLJ();
        
        if(t%100==0) {
            updateMatrix();
        }
        
        updateChain();
        
        if(usePBC) applyPBC();
        
        printTrajectory(t);
        
        //R_gyration(t);
        //R_endtoend(t);
        //Diffusivity(t);
        Extension(t);
        //exit(1);
    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  END OF TIME LOOP
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    printIteration();
    printMatrix(1);
    char* command = malloc(200*sizeof(char));
    if(iteration < MAXITER){
        sprintf(command,"bash run.sh %s %u 1",flowratec,N);
        system(command);
        //exit(1);
    }
    free(command);
    return 0;
}

Vector3D_t getNID(int i,int j){
    Vector3D_t NID;
    
    double dx = Chain[i].rx - Chain[j].rx;
    double dy = Chain[i].ry - Chain[j].ry;
    double dz = Chain[i].rz - Chain[j].rz;
    
    dx -= L*round(dx/L);
    dy -= L*round(dy/L);
    dz -= L*round(dz/L);
    
    NID.x = dx; NID.y = dy; NID.z = dz;
    
    return NID;
}

void applyPBC(){
    for(int i = 0; i<N*NP; ++i){
        Chain[i].rx -= L*round(Chain[i].rx/L);
        Chain[i].ry -= L*round(Chain[i].ry/L);
        Chain[i].rz -= L*round(Chain[i].rz/L);
    }
}

void initForce(){
    for(int i = 0; i<N; ++i){
        Chain[i].fx = 0.0;
        Chain[i].fy = 0.0;
        Chain[i].fz = 0.0;
    }
}

void forceBond(){
    for(int i = 1; i<N; ++i){
        Vector3D_t NID;
        if(usePBC){
            NID = getNID(i,i-1);
        } else{
            NID.x = Chain[i].rx - Chain[i-1].rx;
            NID.y = Chain[i].ry - Chain[i-1].ry;
            NID.z = Chain[i].rz - Chain[i-1].rz;
        }
        double rr = NID.x*NID.x+NID.y*NID.y+NID.z*NID.z;
        double r = sqrt(rr);
        double Fs = -KAPPA*(r-2.0);
        
        Chain[i].fx += Fs*NID.x/r;
        Chain[i].fy += Fs*NID.y/r;
        Chain[i].fz += Fs*NID.z/r;
        Chain[i-1].fx -= Fs*NID.x/r;
        Chain[i-1].fy -= Fs*NID.y/r;
        Chain[i-1].fz -= Fs*NID.z/r;
    }
}

void forceSpringRing(){
    //#pragma omp parallel for private(i) schedule(dynamic)
    for(int i = 0; i<N; ++i){
        int j;
        if(i==0) j = N-1;
        else j = i-1;
        
        Vector3D_t NID;
        
        NID.x = Chain[i].rx - Chain[j].rx;
        NID.y = Chain[i].ry - Chain[j].ry;
        NID.z = Chain[i].rz - Chain[j].rz;
        
        double r = sqrt(NID.x*NID.x+NID.y*NID.y+NID.z*NID.z);
        double Fs = -KAPPA*(r-2.0);
        
        Chain[i].fx += Fs*NID.x/r;
        Chain[i].fy += Fs*NID.y/r;
        Chain[i].fz += Fs*NID.z/r;
        Chain[j].fx -= Fs*NID.x/r;
        Chain[j].fy -= Fs*NID.y/r;
        Chain[j].fz -= Fs*NID.z/r;
    }
}

void forceLJ(){
    for(int i = 0; i<N; ++i){
        for(int j = i+1; j<N; ++j){
            Vector3D_t NID;
            if(usePBC){
                NID = getNID(i,j);
            } else{
                NID.x = Chain[i].rx - Chain[j].rx;
                NID.y = Chain[i].ry - Chain[j].ry;
                NID.z = Chain[i].rz - Chain[j].rz;
            }
            
            double rr = NID.x*NID.x+NID.y*NID.y+NID.z*NID.z;
            if(rr<25.0){
                double ratio = 4.00/rr;
                double r6 = ratio*ratio*ratio;
                if(r6>3) r6 = 3;
                double coeff = (12.0*EPSILON/rr)*(r6*r6-r6);
                Chain[i].fx += coeff*NID.x;
                Chain[i].fy += coeff*NID.y;
                Chain[i].fz += coeff*NID.z;
                Chain[j].fx -= coeff*NID.x;
                Chain[j].fy -= coeff*NID.y;
                Chain[j].fz -= coeff*NID.z;
                //printf("%d:%f %f %f %lf\n",i,Chain[i].fx,Chain[i].fy,Chain[i].fz,rr);
            }
        }
    }
}

void updateChain(){
    int i;
    int nn = 3*N;
    double p = sqrt(2*DT);
    
    float *RR = calloc(nn,sizeof(float));
    for(i = 0; i<nn; ++i){
        RR[i] = gasdev(idum);
    }
    
    if(HImode){
        float *force = calloc(nn,sizeof(float));
        for(i = 0; i<N; ++i){
            force[i*3] = Chain[i].fx;
            force[i*3+1] = Chain[i].fy;
            force[i*3+2] = Chain[i].fz;
        }
        
        float *D_force = calloc(nn,sizeof(float));
        float *D_noise = calloc(nn,sizeof(float));
        
        if(GeyerWinter){
            float *D_ij = calloc(nn*nn,sizeof(float));
            for(i = 0; i<nn*nn; ++i){
                D_ij[i] = matrix->DiffMatrixAvg[i];
            }
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nn,1,nn,1,matrix->DiffMatrixAvg,nn,force,1,1,D_force,1);
            int ii;
            for(ii = 0; ii<nn; ++ii){
                D_ij[ii*nn+ii] /= gwParm->beta_ij;
            }
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nn,1,nn,gwParm->beta_ij,D_ij,nn,RR,1,1,D_noise,1);
            
            if(flowType==0){ //Static
//#pragma omp parallel for private(j) schedule(static)
                for(int j = 0; j<N; ++j){
                    Chain[j].rx += DT*D_force[j*3] + p*gwParm->C_i[j*3]*D_noise[j*3];
                    Chain[j].ry += DT*D_force[j*3+1] + p*gwParm->C_i[j*3+1]*D_noise[j*3+1];
                    Chain[j].rz += DT*D_force[j*3+2] + p*gwParm->C_i[j*3+2]*D_noise[j*3+2];
                }
            } else if(flowType==1){ //Elongational
//#pragma omp parallel for private(j) schedule(static)
                for(int j = 0; j<N; ++j){
                    Chain[j].rx += DT*D_force[j*3] + DT*flowRate*(Chain[j].rx - CoMstored.x) + p*gwParm->C_i[j*3]*D_noise[j*3];
                    Chain[j].ry += DT*D_force[j*3+1] - DT*flowRate*(Chain[j].ry - CoMstored.y) + p*gwParm->C_i[j*3+1]*D_noise[j*3+1];
                    Chain[j].rz += DT*D_force[j*3+2] + p*gwParm->C_i[j*3+2]*D_noise[j*3+2];
                }
            } else if(flowType==2){ //shear
//#pragma omp parallel for private(j) schedule(static)
                for(int j = 0; j<N; ++j){
                    Chain[j].rx += DT*D_force[j*3] + DT*flowRate*(Chain[j].ry - CoMstored.y) + p*gwParm->C_i[j*3]*D_noise[j*3];
                    Chain[j].ry += DT*D_force[j*3+1] + p*gwParm->C_i[j*3+1]*D_noise[j*3+1];
                    Chain[j].rz += DT*D_force[j*3+2] + p*gwParm->C_i[j*3+2]*D_noise[j*3+2];
                }
            }

            
            free(D_ij);
            
        } else{
            
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nn,1,nn,1,matrix->DiffMatrixAvg,nn,force,1,1,D_force,1);
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nn,1,nn,1,Ddecomp,nn,RR,1,1,D_noise,1);
            
            if(flowType==0){
//#pragma omp parallel for private(j) schedule(dynamic)
                for(int j = 0; j<N; ++j){ //Static
                    Chain[j].rx += DT*D_force[j*3] + p*D_noise[j*3];
                    Chain[j].ry += DT*D_force[j*3+1] + p*D_noise[j*3+1];
                    Chain[j].rz += DT*D_force[j*3+2] + p*D_noise[j*3+2];
                }
                
            } else if(flowType==1){
//#pragma omp parallel for private(j) schedule(dynamic)
                for(int j = 0; j<N; ++j){ //Elongational
                    Chain[j].rx += DT*D_force[j*3] + DT*flowRate*(Chain[j].rx - CoMstored.x) + p*D_noise[j*3];
                    Chain[j].ry += DT*D_force[j*3+1] - DT*flowRate*(Chain[j].ry - CoMstored.y) + p*D_noise[j*3+1];
                    Chain[j].rz += DT*D_force[j*3+2] + p*D_noise[j*3+2];
                }
                
            } else if(flowType==2){
//#pragma omp parallel for private(j) schedule(dynamic)
                for(int j = 0; j<N; ++j){ //shear
                    Chain[j].rx += DT*D_force[j*3] + DT*flowRate*(Chain[j].rx - CoMstored.x) + p*D_noise[j*3];
                    Chain[j].ry += DT*D_force[j*3+1] + p*D_noise[j*3+1];
                    Chain[j].rz += DT*D_force[j*3+2] + p*D_noise[j*3+2];
                }
            }
        }
        
        free(D_force);
        free(D_noise);
        free(force);
        
    } else{
        if(flowType==0){
            //#pragma omp parallel for private(j) schedule(dynamic)
            for(int j = 0; j<N; ++j){
                Chain[j].rx += DT*Chain[j].fx + p*RR[j*3];
                Chain[j].ry += DT*Chain[j].fy + p*RR[j*3+1];
                Chain[j].rz += DT*Chain[j].fz + p*RR[j*3+2];
            }
        } else if(flowType==1){
            //#pragma omp parallel for private(j) schedule(dynamic)
            for(int j = 0; j<N; ++j){
                Chain[j].rx += DT*Chain[j].fx + DT*flowRate*(Chain[j].rx - CoMstored.x) + p*RR[j*3];
                Chain[j].ry += DT*Chain[j].fy - DT*flowRate*(Chain[j].ry - CoMstored.y) + p*RR[j*3+1];
                Chain[j].rz += DT*Chain[j].fz + p*RR[j*3+2];
            }
        } else if(flowType==2){
            //#pragma omp parallel for private(j) schedule(dynamic)
            for(int j = 0; j<N; ++j){
                Chain[j].rx += DT*Chain[j].fx + DT*flowRate*(Chain[j].ry - CoMstored.y) + p*RR[j*3];
                Chain[j].ry += DT*Chain[j].fy + p*RR[j*3+1];
                Chain[j].rz += DT*Chain[j].fz + p*RR[j*3+2];
            }
        }

    }
    
    free(RR);
}

Vector3D_t CenterOfMass(int i){
    Vector3D_t CoM;
    if(usePBC){
        CoM.x = Chain[i*N].rx;
        CoM.y = Chain[i*N].ry;
        CoM.z = Chain[i*N].rz;
        for(int j = 1; j<N; ++j){
            Vector3D_t NID = getNID(i*N+j,i*N+j-1);
            CoM.x += (N-j)*NID.x/N;
            CoM.y += (N-j)*NID.y/N;
            CoM.z += (N-j)*NID.z/N;
        }
        CoM.x -= round(CoM.x/L)*L;
        CoM.y -= round(CoM.y/L)*L;
        CoM.z -= round(CoM.z/L)*L;
    } else{
        CoM.x = Chain[i*N].rx/N;
        CoM.y = Chain[i*N].ry/N;
        CoM.z = Chain[i*N].rz/N;
        for(int j = 1; j<N; ++j){
            CoM.x += Chain[i*N+j].rx/N;
            CoM.y += Chain[i*N+j].ry/N;
            CoM.z += Chain[i*N+j].rz/N;
        }
    }
    return CoM;
}

void printTrajectory(int t){
    FILE *Trajectory;
    if(t%trajStep==0){
        Trajectory = fopen(trajName, "a");
        fprintf(Trajectory, "%d\n%d\n", N, t);
        for(int i = 0; i<N; ++i){
            fprintf(Trajectory, "A %lf %lf %lf\n", Chain[i].rx, Chain[i].ry, Chain[i].rz);
        }
        fclose(Trajectory);
    }
}

float gasdev(long *idum){
    float ran1(long *idum);
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;
    
    if (*idum < 0) iset = 0;
    if (iset == 0){
        do{
            v1 = 2.0*ran1(idum)-1.0;
            v2 = 2.0*ran1(idum)-1.0;
            rsq = v1*v1+v2*v2;
        }
        while(rsq >= 1.0 || rsq == 0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;
        return v2*fac;
    }
    else{
        iset = 0;
        return gset;
    }
}

/*long long timer(){
    struct timespec ts;
#ifdef __MACH__
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), SYSTEM_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts.tv_sec = mts.tv_sec;
    ts.tv_nsec = mts.tv_nsec;
#else
    clock_gettime(CLOCK_MONOTONIC, &ts);
#endif
    return 1000000000*ts.tv_sec + ts.tv_nsec;
}*/
