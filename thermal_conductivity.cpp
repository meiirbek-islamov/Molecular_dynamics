//
//  main.cpp
//  conductivity
//
//  Created by Meiirbek Islamov  on 27.04.2020.
//  Copyright © 2020 Meiirbek Islamov . All rights reserved.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <math.h>
#include <stdlib.h>
// #include <stdio.h>
using namespace std;

ifstream data_in1("Heat_flux.txt");
ifstream data_in2("Heat_flux[1].txt");
ifstream data_in3("Heat_flux[2].txt");
//ofstream data_out1("hcacf.txt");
ofstream data_out2("k_data.txt");
//ofstream data_out3("k_data1.txt");


int main()
{
    double dt = 0.002*5;
    double L = 6.83;
    int N = 199999;
    int M = 2000; //correlation length
    double T_d = 70/121.0;
    double V = L*L*L;
    double heatf1[N][3], heatf2[N][3], heatf3[N][3];
    double heatx1, heaty1, heatz1, heatx2, heaty2, heatz2, heatx3, heaty3, heatz3;
    double hcacfx1[M], hcacfy1[M], hcacfz1[M], hcacfx2[M], hcacfy2[M], hcacfz2[M], hcacfx3[M], hcacfy3[M], hcacfz3[M]; // Heat current autocorrelation function
    double t[M]; //correlation time
    double hcacf[M]; // average heat current autocorrelation function
    double k_b = 1; //Boltzmann's constant
    double kt = 0;

  
    
    for(int i=0; i<N; ++i)
    {
        for(int j=0; j<3;++j)
        {
            data_in1 >> heatf1[i][j];
            
        }
    }

    for(int i=0; i<N; ++i)
    {
        for(int j=0; j<3;++j)
        {
          
            data_in2 >> heatf2[i][j];
            
            
        }
    }
    
    for(int i=0; i<N; ++i)
    {
        for(int j=0; j<3;++j)
        {
          
            data_in3 >> heatf3[i][j];
            
        }
    }
    
 
    for (int m=0; m<M; m++){
        
        heatx1 = 0;
        heaty1 = 0;
        heatz1 = 0;
        
        for (int n=0; n<(N-m); n++){
            
            heatx1 = heatx1 + heatf1[n][0] * heatf1[m+n][0];
            heaty1 = heaty1 + heatf1[n][1] * heatf1[m+n][1];
            heatz1 = heatz1 + heatf1[n][2] * heatf1[m+n][2];
            
        }

        
        hcacfx1[m] = heatx1/(N-m);
        hcacfy1[m] = heaty1/(N-m);
        hcacfz1[m] = heatz1/(N-m);
        

    }
    
    for (int m=0; m<M; m++){
            
            heatx2 = 0;
            heaty2 = 0;
            heatz2 = 0;
            
      
            for (int n=0; n<(N-m); n++){
                
                heatx2 = heatx2 + heatf2[n][0] * heatf2[m+n][0];
                heaty2 = heaty2 + heatf2[n][1] * heatf2[m+n][1];
                heatz2 = heatz2 + heatf2[n][2] * heatf2[m+n][2];
                
            }
            
            hcacfx2[m] = heatx2/(N-m);
            hcacfy2[m] = heaty2/(N-m);
            hcacfz2[m] = heatz2/(N-m);

        }
    
    for (int m=0; m<M; m++){
               
               heatx3 = 0;
               heaty3 = 0;
               heatz3 = 0;
               
               for (int n=0; n<(N-m); n++){
                   
                   heatx3 = heatx3 + heatf3[n][0] * heatf3[m+n][0];
                   heaty3 = heaty3 + heatf3[n][1] * heatf3[m+n][1];
                   heatz3 = heatz3 + heatf3[n][2] * heatf3[m+n][2];
                   
               }
               
               hcacfx3[m] = heatx3/(N-m);
               hcacfy3[m] = heaty3/(N-m);
               hcacfz3[m] = heatz3/(N-m);
               
           }
       
    
    
    // Thermal conductivity calculation
    
    
    for (int i=0; i<M; ++i){ // correlation time list
        t[i] = i*dt*2.14;
        
    }
    //First approach to integrate hcacf (Trapezoidal rule)
    kt = 0;
    for (int i=0; i<M; ++i){
        
        hcacf[i] = (hcacfx1[i]+hcacfy1[i]+hcacfz1[i]+hcacfx2[i]+hcacfy2[i]+hcacfz2[i]+hcacfx3[i]+hcacfy3[i]+hcacfz3[i])/9;
        
        if (i == 0 || i == (M-1)) {
            
            kt =  kt+ (hcacf[i] * dt * 0.0190)/(2 * V * k_b * T_d * T_d);
            data_out2<<t[i]<<"    "<<hcacf[i]<<"    "<<kt<<endl;
        }
        
        else{
            kt =  kt + (hcacf[i] * dt * 0.0190)/(V * k_b * T_d * T_d);
            data_out2<<t[i]<<"    "<<hcacf[i]<<"    "<<kt<<endl;
        }
    }
    
//        Second approach ti integrate hcacf (Trapezoidal rule)
//
//        double kt1=0;
//        for (int i=1; i<M; ++i){
//
//                kt1 =  kt1 + ((hcacf[i-1]+hcacf[i])/2 * dt * 0.0190)/(V * k_b * T_d * T_d);
//                data_out3<<t[i]<<"    "<<hcacf[i]<<"    "<<kt1<<endl;
//
//        }
//

    // LJ thermal conductivity scale - 0.0189, 20K-1.2W/mk
            
}
