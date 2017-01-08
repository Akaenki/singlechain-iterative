//
//  Properties.c
//  SingleChainAverage
//
//  Created by Linling Miao on 6/8/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

/* This file contains functions calculating different properties */

#include "Properties.h"
#include "Parameters.h"

void R_gyration(int t){
    if(t==0) {
        Rg = 0.0;
        Ree = 0.0;
    }
    
    Vector3D_t CoM;
    if(usePBC){
        CoM.x = Chain[0].rx;
        CoM.y = Chain[0].ry;
        CoM.z = Chain[0].rz;
        for(int i = 1; i<N; ++i){
            double dx = Chain[i].rx - Chain[i-1].rx;
            double dy = Chain[i].ry - Chain[i-1].ry;
            double dz = Chain[i].rz - Chain[i-1].rz;
            
            dx -= round(dx/L)*L;
            dy -= round(dy/L)*L;
            dz -= round(dz/L)*L;
            
            CoM.x += (N-i)*dx/N;
            CoM.y += (N-i)*dy/N;
            CoM.z += (N-i)*dz/N;
        }
        
        CoM.x -= round(CoM.x/L)*L;
        CoM.y -= round(CoM.y/L)*L;
        CoM.z -= round(CoM.z/L)*L;
        
    } else{
        CoM.x = Chain[0].rx/N;
        CoM.y = Chain[0].ry/N;
        CoM.z = Chain[0].rz/N;
        for(int j = 1; j<N; ++j){
            CoM.x += Chain[j].rx/N;
            CoM.y += Chain[j].ry/N;
            CoM.z += Chain[j].rz/N;
        }
    }
    
    for(int i = 0; i<N; ++i){
        double dx = Chain[i].rx - CoM.x;
        double dy = Chain[i].ry - CoM.y;
        double dz = Chain[i].rz - CoM.z;
        
        Rg += dx*dx+dy*dy+dz*dz;
    }
    
    
    if(usePBC){
        double dx = 0;
        double dy = 0;
        double dz = 0;

        for(int i = 1; i<N; ++i){
            Vector3D_t NID;
            
            NID.x = Chain[i].rx - Chain[i-1].rx;
            NID.y = Chain[i].ry - Chain[i-1].ry;
            NID.z = Chain[i].rz - Chain[i-1].rz;
            
            NID.x -= round(NID.x/L)*L;
            NID.y -= round(NID.y/L)*L;
            NID.y -= round(NID.y/L)*L;
            
            dx += NID.x;
            dy += NID.y;
            dz += NID.z;
        }
        
        Ree += dx*dx+dy*dy+dz*dz;
        
    } else{
        double dx = Chain[0].rx - Chain[N-1].rx;
        double dy = Chain[0].ry - Chain[N-1].ry;
        double dz = Chain[0].rz - Chain[N-1].rz;
        
        Ree += dx*dx+dy*dy+dz*dz;
    }
    
    /*if(t>0 && t%trajStep==0){
        
        FILE *statscale = fopen(RgName,"a");
        fprintf(statscale,"%lf %lf\n",Rg/trajStep/N,Ree/trajStep);
        fclose(statscale);
    }*/
    
    if(t%trajStep==0){
        Ree = 0.0;
        Rg = 0.0;
    }
}

void Diffusivity(int t){
    Vector3D_t CoM;
    if(usePBC){
        CoM.x = Chain[0].rx;
        CoM.y = Chain[0].ry;
        CoM.z = Chain[0].rz;
        int i;
        for(i = 1; i<N; ++i){
            double dx = Chain[i].rx - Chain[i-1].rx;
            double dy = Chain[i].ry - Chain[i-1].ry;
            double dz = Chain[i].rz - Chain[i-1].rz;
            
            dx -= round(dx/L)*L;
            dy -= round(dy/L)*L;
            dz -= round(dz/L)*L;
            
            CoM.x += (N-i)*dx/N;
            CoM.y += (N-i)*dy/N;
            CoM.z += (N-i)*dz/N;
        }
        
        CoM.x -= round(CoM.x/L)*L;
        CoM.y -= round(CoM.y/L)*L;
        CoM.z -= round(CoM.z/L)*L;
        
    } else{
        CoM.x = Chain[0].rx/N;
        CoM.y = Chain[0].ry/N;
        CoM.z = Chain[0].rz/N;
        int j;
        for(j = 1; j<N; ++j){
            CoM.x += Chain[j].rx/N;
            CoM.y += Chain[j].ry/N;
            CoM.z += Chain[j].rz/N;
        }
    }
    
    Vector3D_t dCoM;
    dCoM.x = CoM.x - CoMstored.x;
    dCoM.y = CoM.y - CoMstored.y;
    dCoM.z = CoM.z - CoMstored.z;
    
    if(usePBC){
        dCoM.x -= round(dCoM.x/L)*L;
        dCoM.y -= round(dCoM.y/L)*L;
        dCoM.z -= round(dCoM.z/L)*L;
    }
    
    MSD_CoM += dCoM.x*dCoM.x+dCoM.y*dCoM.y+dCoM.z*dCoM.z;
    
    /*if(t>0 && t%trajStep==0){
        
        FILE *dg = fopen(DiffName,"a");
        fprintf(dg,"%lf\n",MSD_CoM/6.0/(t*dt));
        fclose(dg);
    }*/
    
    storeCoM();
}


void storeCoM(){
    if(usePBC){
        CoMstored.x = Chain[0].rx;
        CoMstored.y = Chain[0].ry;
        CoMstored.z = Chain[0].rz;

        for(int i = 1; i<N; ++i){
            double dx = Chain[i].rx - Chain[i-1].rx;
            double dy = Chain[i].ry - Chain[i-1].ry;
            double dz = Chain[i].rz - Chain[i-1].rz;
            
            dx -= round(dx/L)*L;
            dy -= round(dy/L)*L;
            dz -= round(dz/L)*L;
            
            CoMstored.x += (N-i)*dx/N;
            CoMstored.y += (N-i)*dy/N;
            CoMstored.z += (N-i)*dz/N;
        }
        
        CoMstored.x -= round(CoMstored.x/L)*L;
        CoMstored.y -= round(CoMstored.y/L)*L;
        CoMstored.z -= round(CoMstored.z/L)*L;
        
    } else{
        CoMstored.x = Chain[0].rx/N;
        CoMstored.y = Chain[0].ry/N;
        CoMstored.z = Chain[0].rz/N;

        for(int j = 1; j<N; ++j){
            CoMstored.x += Chain[j].rx/N;
            CoMstored.y += Chain[j].ry/N;
            CoMstored.z += Chain[j].rz/N;
        }
    }
}


void R_endtoend(int t){
    double dx = 0;
    double dy = 0;
    double dz = 0;
    if(usePBC) {
        for(int i = 1; i<N; ++i){
            double rx = Chain[i].rx - Chain[i-1].rx;
            double ry = Chain[i].ry - Chain[i-1].rx;
            double rz = Chain[i].rz - Chain[i-1].rx;
            
            rx -= round(rx/L)*L;
            ry -= round(ry/L)*L;
            rz -= round(rz/L)*L;
            
            dx += rx;
            dy += ry;
            dz += rz;
        }
    } else{
        dx = Chain[N-1].rx - Chain[0].rx;
        dy = Chain[N-1].ry - Chain[0].ry;
        dz = Chain[N-1].rz - Chain[0].rz;
    }
    
    /*if(t>0 && t%1000==0){
        FILE *ree = fopen(ReeName,"a");
        fprintf(ree,"%lf %lf %lf\n", dx,dy,dz);
        fclose(ree);
    }*/
}


void Extension(int t){
    
    double dxmax = 0.0;
    double dxmin = 0.0;
    double ddx = 0.0;
    
    int i;
    for(i = 1; i<N; ++i){
        double dx = Chain[i].rx - Chain[i-1].rx;
        
        if(usePBC) dx -= round(dx/L)*L;
        
        ddx += dx;
        if(ddx<dxmin) dxmin = ddx;
        else if(ddx>dxmax) dxmax = ddx;
    }
    ext.x += dxmax - dxmin;
    
    /*if(t>0 && t%1000==0){
        FILE *extension = fopen(extName,"a");
        fprintf(extension,"%lf\n",ext/1000);
        fclose(extension);
    }*/
    
    if(t%1000==0) ext.x = 0.0;
    
}

void ExtensionRing(int t){
    for(int i = 0; i<NP; ++i){
        double dxmax = 0.0; double dymax = 0.0; double dzmax = 0.0;
        double dxmin = 0.0; double dymin = 0.0; double dzmin = 0.0;
        double ddx = 0.0; double ddy = 0.0; double ddz = 0.0;
        
        for(int j = 1; j<N; ++j){
            Vector3D_t NID;
            NID.x = Chain[i*N+j].rx - Chain[i*N+j-1].rx;
            NID.y = Chain[i*N+j].ry - Chain[i*N+j-1].ry;
            NID.z = Chain[i*N+j].rz - Chain[i*N+j-1].rz;
            
            ddx += NID.x;
            ddy += NID.y;
            ddz += NID.z;
            if(ddx<dxmin) dxmin = ddx;
            else if(ddx>dxmax) dxmax = ddx;
            if(ddy<dymin) dymin = ddy;
            else if(ddy>dymax) dymax = ddy;
            if(ddz<dzmin) dzmin = ddz;
            else if(ddz>dzmax) dzmax = ddz;
        }
        ext.x += dxmax - dxmin;
        ext.y += dymax - dymin;
        ext.z += dzmax - dzmin;
    }
    
    /*if(t>0 && t%1000==0){
     FILE *Extensionout;
     Extensionout = fopen(extName, "a");
     fprintf(Extensionout, "%lf %lf %lf\n", ext.dx/1000, ext.dy/1000, ext.dz/1000);
     fclose(Extensionout);
     }*/
    
    if(t%1000==0){
        ext.x = 0.0; ext.y = 0.0; ext.z = 0.0;
    }
}


