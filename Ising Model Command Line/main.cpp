//
//  main.cpp
//  Ising Model Command Line
//
//  Created by Borun Chowdhury on 4/26/15.
//  Copyright (c) 2015 Borun Chowdhury and Mouchumi Bhuyan. All rights reserved.
//

#include <iostream>


#include "IsingModel.h"

using namespace std;

float cv_Wolff(int N, float T,IsingModel &im);
float cv_Wolff_2(int N, float T,IsingModel &im);
float cv_HeatBath(int N, float T,IsingModel &im);
float cv_Metropolis(int N, float T,IsingModel &im);



int main() {

    int N;
    float T;

    
    /*for(T=1.0;T<3.0;T+=.1) {
        IsingModel im(N,T);
        //cout << "N,T, cv(W), cv(HB), cv(M): "<< N <<", " << T << ", "<< cv_Wolff(N,T,im) << ", "  << cv_HeatBath(N,T,im) << ", " << cv_Metropolis(N, T,im)  << endl;
        cout << "N,T, cv(W): "<< N <<", " << T << ", "<< cv_Wolff(N,T,im) << endl;
    }*/
    ;
    cout << "Sampling energies only when N*N spins have been flipped" << endl;
    
    for(int i=0;i<=5;i++) {
        N=2<<i;
        T=2.27;
        IsingModel im(N,T);
        //cout << "N,T, cv(W), cv(HB), cv(M): "<< N <<", " << T << ", "<< cv_Wolff(N,T,im) << ", "  << cv_HeatBath(N,T,im) << ", " << cv_Metropolis(N, T,im)  << endl;
        cout << "N,T, cv(W): "<<N << ", "<< T << ", " <<cv_Wolff(N, T,im) << endl;
    }
    
    
    cout << endl << "Sampling energies on every flip of a cluster" << endl;

    
    for(int i=0;i<=5;i++) {
        N=2<<i;
        T=2.27;
        IsingModel im(N,T);
        //cout << "N,T, cv(W), cv(HB), cv(M): "<< N <<", " << T << ", "<< cv_Wolff(N,T,im) << ", "  << cv_HeatBath(N,T,im) << ", " << cv_Metropolis(N, T,im)  << endl;
        cout << "N,T, cv(W): "<<N << ", "<< T << ", " <<cv_Wolff_2(N, T,im) << endl;
    }
    
    
    cout << endl << "Results for Heat Bath and Metropolis " << endl;
    
    
    for(int i=0;i<=5;i++) {
        N=2<<i;
        T=2.27;
        IsingModel im(N,T);
        cout << "N,T, cv(HB), cv(M): "<< N <<", " << T << ", "<< cv_HeatBath(N,T,im) << ", " << cv_Metropolis(N, T,im)  << endl;
        //cout << "N,T, cv(W): "<<N << ", "<< T << ", " <<cv_Wolff_2(N, T,im) << endl;
    }
    
    
    
}



float cv_HeatBath(int N, float T,  IsingModel &im) {

    float E,e=0.0;
    float ESq=0.0;
    
 
    int nsteps=10000;
    
    // Drop the first 10 runs for equilibriation
    for(int i=0;i<1000;i++) {
        im.SweepHeatBath(1);
    }
    
    for(int i=0;i<nsteps;i++) {
        im.SweepHeatBath(1);
        e=im.Energy();
        E+=e;
        ESq+=e*e;
    }
    
    return (ESq/nsteps-(E/nsteps)*(E/nsteps))/(N*N)/(T*T);
}


float cv_Metropolis(int N, float T, IsingModel &im) {

    float E,e=0.0;
    float ESq=0.0;
    
    
    int nsteps=10000;
    
    // Drop the first 10 runs for equilibriation
    for(int i=0;i<1000;i++) {
        im.SweepHeatBath(1);
    }
    
    for(int i=0;i<nsteps;i++) {
        im.SweepMetropolis(1);
        e=im.Energy();
        E+=e;
        ESq+=e*e;
    }
    return (ESq/nsteps-(E/nsteps)*(E/nsteps))/(N*N)/(T*T);
}

float cv_Wolff(int N, float T,IsingModel &im) {

    float E,e=0.0;
    float ESq=0.0;
    
   
    
    int p_sweep=0;

    int nsteps=10000;
    
    
    // Drop the first 10 runs for equilibriation
    for(int i=0;i<1000;i++){
        p_sweep=im.SweepWolff(1,p_sweep);
    }
    
    for(int i=0;i<nsteps;i++) {
        p_sweep=im.SweepWolff(1,p_sweep);
        e=im.Energy();
        E+=e;
        ESq+=e*e;
    }
    
    return (ESq/nsteps-(E/nsteps)*(E/nsteps))/(N*N)/(T*T);
}


float cv_Wolff_2(int N, float T,IsingModel &im) {
    
    float E,e=0.0;
    float ESq=0.0;
    
    
    
    int p_sweep=0;
    
    int nsteps=10000;
    
    
    // Drop the first 10 runs for equilibriation
    for(int i=0;i<1000;i++){
        p_sweep=im.WolffMove();
    }
    
    for(int i=0;i<nsteps;i++) {
        p_sweep=im.WolffMove();
        e=im.Energy();
        E+=e;
        ESq+=e*e;
    }
    
    return (ESq/nsteps-(E/nsteps)*(E/nsteps))/(N*N)/(T*T);
}
