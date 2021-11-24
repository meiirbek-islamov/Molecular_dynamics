//
//  main.cpp
//  MD code
//
//  Created by Meiirbek Islamov  on 14.02.2020.
//  Copyright © 2020 Meiirbek Islamov . All rights reserved.
//


#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <cmath>
using namespace std;
//*****************

// declaration of input and output stream

ifstream data_in("10.txt"); // input data
ofstream data_out1("position.txt"); // output data: position
ofstream data_out2("velocity.txt"); // output data: velocity
ofstream data_out3("force.txt");    // outpud data: force
ofstream data_out4("total_U.txt");  // outpud data: Total Potential Energy
ofstream data_out5("total_K.txt");  // outpud data: Total kinetic Energy
ofstream data_out6("total_P.txt");  // outpud data: Total momentum
ofstream data_out7("position.xyz");  // outpud data: Total momentum

// declaration of function

double potential(double a, double b, double c, double d, double e, double f); // a pair potential energy

int main (){
    
    // variable declarations
    
    int c=10; // number of particles
    int i,j,k,n,m,p,t; //loop counters
    int N = 1000;   // Number of timestep
    double position[c][3]; //initial positions
    double position1[c][3]; // updated position at each timestep
    double v[c][3]; // initial velocity
    double v1[c][3]; // updated velocity
    double rij[3]; //pair distance vector
    double F[c][3]; // force
    double F1[c][3]; // updated force
    double U; //potential energy
    double U0;
    double K; //kinetic energy
    double K0;
    double ts = 0.002; //timestep
    double v12[c][3]; // velocity at half timestep
    double r2,r4,r8;
    double forcebyr; //force divided by r
    double Px,Py,Pz;
    
    
    
// read in data from file
    
    for(n=0;n<c;n++)
    {
           for (m=0; m<3; m++)
           {
           data_in>>position[n][m]; //reading in the initial positions
               v[n][m] = 0; //Initialisation of particle velocity to zero
               F[n][m] = 0; //Initialisation of the force to zero;
               
           }
    }
    
// calculating the total potential energy from initial positions
    
    U=0;
    for (i=0;i<c;i++)
    {
        for (j=i+1;j<c;j++)
        {
        U0=potential(position[i][0], position[i][1], position[i][2],
                     position[j][0], position[j][1], position[j][2]);
        U+=U0;
        }
    }
    
   
// calculating the force from initial positions
    
    for (i=0;i<c;i++)
    {
        for (j=i+1;j<c;j++)
        {
            for (k=0;k<3;k++)
            {
                rij[k] =position[i][k]-position[j][k];
            }
                r2 = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
                r4 = r2*r2;
                r8 = r4*r4;
            forcebyr = (48/(r8*r4*r2))-(24/(r8));
            
            for (k=0;k<3;k++)
            {
                F[i][k] += forcebyr*rij[k];
                F[j][k] -= forcebyr*rij[k];
            }
            }
    }
    

// Let's start the loop over number of timestep N
    
    for(p=0;p<N;p++) // for loop for timestep
    {
        
// velocity Verlet algorithm
        
        for (i=0;i<c;i++)
        {

            for (k=0;k<3;k++){
                v12[i][k] = v[i][k]+F[i][k]*ts/2; // velocity at half timestep
                position1[i][k] = position[i][k]+ts*v12[i][k]; // updated position
                data_out1<< position1[i][k]<<"       ";
                
                
                F1[i][k] = 0; // Initialize the updated force to zero
                }
                data_out1<<endl;
            }
                data_out1<<endl;
//outputting the xyz positions readible using vmd
        t=0;
        for (i=0;i<c;i++){
            t+=1;
        }
            data_out7<<t<<endl<<endl;
            for (i=0;i<c;i++){
                data_out7<<1<<"     ";
            for (k=0;k<3;k++){
            data_out7<<position1[i][k]<<"   ";
            }
            data_out7<<endl;
        }
        
        
        
//Calculating the updated force
        
        for (i=0;i<c;i++){
                for (j=i+1;j<c;j++){
                    for (k=0;k<3;k++){
                        rij[k] =position1[i][k]-position1[j][k];
                    }
                        r2 = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
                        r4 = r2*r2;
                        r8 = r4*r4;
                        forcebyr = (48/(r8*r4*r2))-(24/(r8));
                    for (k=0;k<3;k++){
                        F1[i][k] += forcebyr*rij[k];
                        F1[j][k] -= forcebyr*rij[k];
                    }
                    }
        }
// Output the updated force
        
        for (i=0;i<c;i++)
        {
        for (k=0;k<3;k++)
        {
            data_out3<<F1[i][k]<<"       ";
            
        }
            data_out3<<endl;
        }
        
//calculating the potential energy
        
        U=0;
        for (i=0;i<c;i++)
        {
            for (j=i+1;j<c;j++)
            {
            U0=potential(position1[i][0], position1[i][1], position1[i][2],
                         position1[j][0], position1[j][1], position1[j][2]); // use function
            U+=U0;
            }
        }
            data_out4<<U<<endl;
        
        
// calculate the updated velocity at timestep and output
        
        for (i=0;i<c;i++){
        for (k=0;k<3;k++){
                v1[i][k] = v12[i][k]+F1[i][k]*ts/2;
                data_out2<< v1[i][k]<<"            ";
        }
                data_out2<<endl; // output the velocity
        }
        
// calculating the kinetic energy
        
        K=0;
        for (i=0;i<c;i++){
            K0 =(v1[i][0]*v1[i][0]+v1[i][1]*v1[i][1]+v1[i][2]*v1[i][2])/2;
            K+=K0;
        }
            data_out5<<K<<endl; // output the total kinetic energy at each timestep

// calculating the momentum
        
        Px=0; Py=0; Pz=0;
        for (i=0;i<c;i++){
            Px+=v1[i][0];
            Py+=v1[i][1];
            Pz+=v1[i][2];
            
        }
// output the momentum
        
            data_out6<<Px<<"               "<<Py<<"               "<<Pz<<"               "<<endl;
        
        
// Update the initial force, position and velocity for the next run
        
        for (i=0;i<c;i++){
        for (k=0;k<3;k++){
            F[i][k]=F1[i][k];
            position[i][k]=position1[i][k];
            v[i][k]=v1[i][k];
        }
        }
        
    }
   
    }


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double potential(double a, double b, double c, double d, double e, double f)
{
    int k;
    double rij[3];
    double p[2][3];
    double r2,r6,r12;
    double U;
    
    p[0][0]=a; p[0][1]=b; p[0][2]=c; p[1][0]=d; p[1][1]=e; p[1][2]=f;
    
    for (k=0;k<3;k++)
    {
        rij[k] =p[0][k]-p[1][k];
    }
        r2 = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
        r6 =1/(r2*r2*r2);
        r12 = r6*r6;
        U = 4*(r12-r6);
    
     return U;
}

