//
//  Properties.h
//  SingleChainAverage
//
//  Created by Linling Miao on 6/8/16.
//  Copyright © 2016 Linling Miao. All rights reserved.
//

#ifndef Properties_h
#define Properties_h

void R_gyration(int t); //Calculate Rg and Ree
void R_endtoend(int t); //Calculate end-to-end vectors
void Diffusivity(int t); //self-diffusivity
void Extension(int t); //Extension in x direction
void ExtensionRing(int t); //Extesion in all directions for ring polymers
void storeCoM(); //Store the CoM of chains

#endif /* Properties_h */
