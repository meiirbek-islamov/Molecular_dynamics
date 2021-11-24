//
//  main.cpp
//  HW#4_1c
//
//  Created by Meiirbek Islamov  on 16.03.2020.
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

ifstream data_in("liquid256.txt"); // input data
ofstream data_out4("total_U.txt");  // outpud data: Total Potential Energy
ofstream data_out5("total_K.txt");  // outpud data: Total kinetic Energy
ofstream data_out6("total_P.txt");  // outpud data: Total momentum
ofstream data_out7("position.xyz");  // outpud data: position.xyz
ofstream data_out8("temperature.txt");  // outpud data: Temperature
ofstream data_out9("pressure.txt");  // outpud data: Pressure

// declaration of function

double potential(double a, double b, double c, double d, double e, double f); // a pair potential energy
double randVel();
//  Initialize velocities based on the Initial Temperature


int main (){
    
    srand(time(NULL)); //
    
    // variable declarations
    
    int c=256; // number of particles
    int i,j,k,n,m,p,t; //loop counters
    int N = 100000;   // Number of timestep
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
    double r,r2,r4,r8;
    double forcebyr; //force divided by r
    double Px,Py,Pz;
    double L = 7.4; // size of the box
    double Rc = 2.5; //cut-off radius
    double dudrRc; //Force at Rc
    double T; //Temperature
    double kb = 1;
    double Ts = 121; // Temperature scale for Argon in K
    double Tstar; // dimensionless temperature
    double V = L*L*L; // volume of the box
    double r6,r12;
    double forcer; //pressure virial term, force times by position
    double forcer0;
    double P; // dimensionless pressure;
    double Pstar; // dimensionless pressure
    double Ps = 425e5; // Pressure scale for Argon in P
    double Tinit=1; // Initial dimensionless Temperature (different in each run)
    double mu = 0; // initial conidtion for thermostat constant
    double Tset = 0.826446;
       
    
    
    
    
// read in data from file
    
    for(n=0;n<c;n++)
    {
           for (m=0; m<3; m++)
           {
           data_in>>position[n][m]; //reading in the initial positions
               
               F[n][m] = 0; //Initialisation of the force to zero;
               
           }
    }
    
// Initializing the particle velocities
    
    for (i=0; i<c; i++) {
        for (k=0; k<3; k++) {
            //  Pull a number from a Random Velocity generation function
            v[i][k] = randVel();
        }
    }
//compute a center of mass
    double vc[3] = {0, 0, 0};

    for (i=0; i<c; i++) {
      for (k=0; k<3; k++) {

        vc[k] += v[i][k];

      }
    }
    
    for (k=0; k<3; k++){
        vc[k] /= c; //velocity of center of mass
        
    }
    
//  Subtract out the center-of-mass velocity from the velocity of each particle
    
    for (i=0; i<c; i++) {
      for (k=0; k<3; k++) {

        v[i][k] -= vc[k];

      }
    }
// scaling the average velocity of the system by a factor according to our intial temperature
    
    double vs; //sum of square of velocities
    double eta; // scaling factor
    vs=0.;
    for (i=0; i<c; i++) {
      for (k=0; k<3; k++) {

        vs += v[i][k]*v[i][k];

      }
    }

    eta = sqrt(3*(c-1)*Tinit/vs);

    for (i=0; i<c; i++) {
      for (k=0; k<3; k++) {

        v[i][k] *= eta;

      }
    }
    
    
// Check if an atom still in the simulation cell
    
    
    for (i=0;i<c;i++){
        for (k=0;k<3;k++){
            if (position[i][k] < 0){
                position[i][k]=position[i][k]+L;
            }
            if (position[i][k] > L){
                position[i][k]=position[i][k]-L;
            }
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
                
                if (rij[k]>L/2){
                    
                    rij[k] = rij[k] - L;
                }
                
                if (rij[k]<-L/2){
                    
                    rij[k] = rij[k] + L;
                }
            }
                r = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
            
              if (r<=Rc){
                    r2 = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
                    r4 = r2*r2;
                    r8 = r4*r4;
                    dudrRc = (-48/pow(Rc,13))+(24/pow(Rc,7));
                    forcebyr = (48/(r8*r4*r2))-(24/(r8))+dudrRc/r;
                    
                    for (k=0;k<3;k++){
                    F[i][k] += forcebyr*rij[k];
                    F[j][k] -= forcebyr*rij[k];
                }
            }
                else {
                    for (k=0;k<3;k++){
                    F[i][k] += 0;
                    F[j][k] -= 0;
                }
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
                v12[i][k] = v[i][k]+(F[i][k]-mu*v[i][k])*ts/2; // velocity at half timestep
                position1[i][k] = position[i][k]+ts*v12[i][k]; // updated position
             // data_out1<< position1[i][k]<<"       ";
                
                
                F1[i][k] = 0; // Initialize the updated force to zero
                }
            //  data_out1<<endl;
            }
            //  data_out1<<endl;
        
// Check if an atom still in the simulation cell
        
        
        for (i=0;i<c;i++){
            for (k=0;k<3;k++){
                if (position1[i][k] < 0){
                    position1[i][k]=position1[i][k]+L;
                }
                if (position1[i][k] > L){
                    position1[i][k]=position1[i][k]-L;
                }
            }
            
        }
        
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
        
        forcer=0;
        for (i=0;i<c;i++){
                for (j=i+1;j<c;j++){
                    for (k=0;k<3;k++){
                        rij[k] =position1[i][k]-position1[j][k];
                        
                        if (rij[k]>L/2){
                            
                            rij[k] = rij[k] - L;
                        }
                        
                        if (rij[k]<-L/2){
                            
                            rij[k] = rij[k] + L;
                        }
                    }
                    
                    r = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
                    
                    if (r<=Rc){
                        r2 = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
                        r4 = r2*r2;
                        r8 = r4*r4;
                        r6 =1/(r2*r2*r2);
                        r12 = r6*r6;
                        dudrRc = (-48/pow(Rc,13))+(24/pow(Rc,7));
                        forcebyr = (48/(r8*r4*r2))-(24/(r8))+dudrRc/r;
                        forcer0 = 48*r12-24*r6+dudrRc*r; // pressure virial term, force times by distance
                        
                        for (k=0;k<3;k++){
                        F1[i][k] += forcebyr*rij[k];
                        F1[j][k] -= forcebyr*rij[k];
                    }
                }
                    else {
                        for (k=0;k<3;k++){
                        F1[i][k] += 0;
                        F1[j][k] -= 0;
                        forcer0 = 0;
                    }
                }
                    forcer+=forcer0;
            }
        }
        
        
// calculating the thermostat constant
        
        // Calculate the initial kinetic energy for thermostat constant equation of motion
        double K1,K01;
        K1=0;
        for (i=0;i<c;i++){
            K01 =(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2])/2;
            K1+=K01;
        }

        // Calculating the Temperature
        double T1;
        T1 = 2*K1/(3*(c-1)*kb); // dimensionless temperature for use in the equation of motion of thermostat constant
// calculating the updated thermostat constant
        
        double mu1; // updated thermostat constant
        double tau=0.05; // thermostat time constant
        
        mu1 = mu + (ts/(tau*tau))*((T1/Tset)-1);
        
        
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
                v1[i][k] = (v12[i][k]+F1[i][k]*ts/2)/(1+mu1*ts/2);
     //           data_out2<< v1[i][k]<<"            ";
        }
    //         data_out2<<endl; // output the velocity
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
        
         //   data_out6<<Px<<"               "<<Py<<"               "<<Pz<<"               "<<endl;
        
// calculating the Temperature
        
        Tstar = 2*K/(3*(c-1)*kb); // dimensionless temperature
        T = Tstar*Ts; // temperature in K
        data_out8<<T<<endl; // output the temperature at each timestep

// calculating the pressure
        
        Pstar = (c*kb*Tstar/V) + forcer/(3*V); //dimensionless pressure
        P = Pstar*Ps; // pressure in Pa
        
        data_out9<<P<<endl; // output the pressure at each timestep
        cout<<p+1<<"        "<<T<<endl;
        
// Update the initial force, position and velocity for the next run
        
        for (i=0;i<c;i++){
        for (k=0;k<3;k++){
            F[i][k]=F1[i][k];
            position[i][k]=position1[i][k];
            v[i][k]=v1[i][k];
        }
        }
        mu = mu1;
        
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
    double dudrRc;
    double r;
    double Rc = 2.5;
    double URc;
    double L = 7.4;
    
    p[0][0]=a; p[0][1]=b; p[0][2]=c; p[1][0]=d; p[1][1]=e; p[1][2]=f;
    
    for (k=0;k<3;k++)
    {
        rij[k] =p[0][k]-p[1][k];
        
        if (rij[k]>L/2){
            
            rij[k] = rij[k] - L;
        }
        
        if (rij[k]<-L/2){
            
            rij[k] = rij[k] + L;
        }
    }
        r = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
        
    if (r<=Rc){
        r2 = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
        r6 =1/(r2*r2*r2);
        r12 = r6*r6;
        URc = 4*((1/pow(Rc,12))-(1/pow(Rc,6)));
        dudrRc = (-48/pow(Rc,13))+(24/pow(Rc,7));
        U = 4*(r12-r6)-URc-(r-Rc)*dudrRc;
    }
    else {
        U = 0;
    }
    
     return U;
}

double randVel() {
    
    double step;
    int a=0;
    int b=10;
    double v;
    step = (b-a)/double(RAND_MAX);
    v=a+rand()*step;
    return v;
    
}

