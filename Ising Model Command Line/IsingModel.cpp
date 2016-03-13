//
//  IsingModel.cpp
//  Ising Model Command Line
//
//  Created by Borun Chowdhury on 4/26/15.
//  Copyright (c) 2015 Borun Chowdhury and Mouchumi Bhuyan. All rights reserved.
//

#include <stdio.h>

#include <stack>
#include <iostream>
#include <fstream>
#include "IsingModel.h"
#include <time.h>

using namespace std;

IsingModel::IsingModel(int n, float t, int d) {
    N=n;
    T=t;
    D=d;
    H=0;
    
    // Initialize random N*N lattice
    srand (time(0));
    
    Lattice=new bool *[N];
    for (int i=0;i<N;i++){
        Lattice[i]=new bool[N];
        for(int j=0;j<N;j++){
            Lattice[i][j]=rand() % 2;
        }
    }

 
    
    
    // Initialize probabilities for heat bath algorithm
    
    initializeHeatBath();
    
    // Initialize probabilities for Metropolis algorithm
    
     initializeMetropolis();
    
    // Initialize probability for Wolff algorithm
    
    initializeWolf();
}



void IsingModel::initializeHeatBath(){
    heat_bath_prob_up = new float[5];
    for(int nUp=0;nUp<5;nUp++){
        int sum_neighbors=2*(nUp-2);
        
        float eUp=-sum_neighbors-H;
        float eDown=sum_neighbors-H;
        if(T != 0){
            float boltzUp = exp(-eUp/T);
            float boltzDown = exp(-eDown/T);
            heat_bath_prob_up[nUp]=  boltzUp/(boltzUp+boltzDown);
        } /*else {
            if(eUp>0) heat_bath_prob_up[sum_neighbors]=0;
            else if (eUp<0)  heat_bath_prob_up[sum_neighbors]=1;
            else  heat_bath_prob_up[sum_neighbors]=.5;
            
        }*/
    }
}




void IsingModel::initializeMetropolis(){
    int sum_neighbors;
    metropolis_prob_up[0]=new float[5];
    metropolis_prob_up[1]=new float[5];
    float deltaEUp,deltaEDown;
    for(int nUp=0;nUp<5;nUp++){
        sum_neighbors=2*(nUp-2);
        // if spin is up then
        deltaEUp=2*sum_neighbors;
        // if spin is down then
        deltaEDown=-2*sum_neighbors;
        
        if(T != 0){
            metropolis_prob_up[0][nUp]=exp(-deltaEDown/T);
            metropolis_prob_up[1][nUp]=exp(-deltaEUp/T);
            /*if(eDown>eUp) {
                // down spin is unstable so flip it
                metropolis_prob_up[0][nUp]=1;
                metropolis_prob_up[1][nUp]=1-exp(-(eDown-eUp)/T);
            } else {
                // up spin is unstable so SHOULD flip it but Sethna has RHS=0 in the second line!
                metropolis_prob_up[0][nUp]=1-exp(-(eUp-eDown)/T);
                metropolis_prob_up[1][nUp]=1;
            }*/
        } /*else {
            if(eDown>eUp) {
                metropolis_prob_up[0][nUp]=1;
                metropolis_prob_up[1][nUp]=0;
            } else {
                metropolis_prob_up[0][nUp]=0;
                metropolis_prob_up[1][nUp]=1;
            }
        }*/
    }
}






int *IsingModel::Neighbors(int *neigbor_rows_cols, int i, int j ){
    int ip = (i+1) % N;
    int im = (i-1) % N;
    int jp = (j+1) % N;
    int jm = (j-1) % N;
    if (im<0) {
        im+=N;
    }
    if (jm<0){
        jm+=N;
    }
    *neigbor_rows_cols=ip;
    *(neigbor_rows_cols+1)=im;
    *(neigbor_rows_cols+2)=jp;
    *(neigbor_rows_cols+3)=jp;
    return neigbor_rows_cols;
}

void IsingModel::RandomLocations(int *rows, int *cols){
    for(int i=0;i<N*N;i++){
            rows[i]=rand() % N;
            cols[i]=rand() % N;
    }
}

void IsingModel::RandomRealUniforms(float *randomArr){
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, 1);
    for(int i=0;i<N*N;i++){
            randomArr[i]=dis(gen);
    }
}




void IsingModel::SweepHeatBath(int nTimes){
    
    int iArr[N*N];
    int jArr[N*N];
    float randomArr[N*N];
        
    for(int i=0;i<nTimes;i++){
        RandomLocations(iArr,jArr);

        RandomRealUniforms(randomArr);
        for(int k=0;k<N*N;k++){
            if(randomArr[k]< heat_bath_prob_up[NeighborsUp(iArr[k],jArr[k])]){
                setSpin(iArr[k],jArr[k],1);
            } else {
                setSpin(iArr[k],jArr[k],-1);
            }
        }
    }
}

void IsingModel::SweepMetropolis(int nTimes){

    
    int iArr[N*N];
    int jArr[N*N];
    float randomArr[N*N];
    
    for(int i=0;i<nTimes;i++){
        RandomLocations(iArr,jArr);
        RandomRealUniforms(randomArr);
        
        for(int k=0;k<N*N;k++){
            if(randomArr[k]< metropolis_prob_up[getSpin(iArr[k],jArr[k])][NeighborsUp(iArr[k],jArr[k])]){
                flipSpin(iArr[k],jArr[k]);
            }
        }
    }
}

void IsingModel::SaveLattice(){
    // Yet to implement
}


int IsingModel::WolffMove() {
    
    // Set up generator for uniform distribution
    
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, 1);
    
    
    // Number of spins flipped in this call
    int spinsFlipped=0;
    
    // Start with a random spin
    int i=rand() % N;
    int j=rand() % N;
    Spin s;
    s.row=i;
    s.col=j;
    
    int oldSpin=getSpin(i,j);
    

    
    // Create a stack to store members of the "cluster"
    stack<Spin> mystack;
    
    // Push the original random spin in
    mystack.push(s);
    
    // While the stack is not empty
    while(!mystack.empty()){
        
        // take out the top element
        s=mystack.top();
        mystack.pop();
        
        i=s.row;
        j=s.col;
        
        // Only do something if the element has not already been visited and therefore has not been flipped
        if(oldSpin==getSpin(i,j)){
            
            // Flip the spin of this element
            flipSpin(i,j);
            spinsFlipped++;
            
            int neighbors_rows_cols[4];
            Neighbors(neighbors_rows_cols, i,j);
            int ip,im,jp,jm;
            ip=neighbors_rows_cols[0];
            im=neighbors_rows_cols[1];
            jp=neighbors_rows_cols[2];
            jm=neighbors_rows_cols[3];
            
            

            if ((getSpin(ip,j) == oldSpin) && (dis(gen) < WolffP)) {
                s.row=ip;
                s.col=j;
                mystack.push(s);
            }
            if ((getSpin(im,j) == oldSpin) && (dis(gen) < WolffP) ) {
                s.row=im;
                s.col=j;
                mystack.push(s);
            }
            if ((getSpin(i,jp) == oldSpin) && (dis(gen) < WolffP)) {
                s.row=i;
                s.col=jp;
                mystack.push(s);
            }
            if ((getSpin(i,jm) == oldSpin)  && (dis(gen) < WolffP)) {
                s.row=i;
                s.col=jm;
                mystack.push(s);
            }
        }
        
    }
    return spinsFlipped;
    
}


// Setup probability to add a spin to a cluster
void IsingModel::initializeWolf() {
    WolffP=1.0-exp(-2.0/T);
}

void IsingModel::flipSpin(int i,int j){
    Lattice[i][j]= ! Lattice[i][j];
}






// can make this inline

inline int IsingModel::getSpin(int row, int col) {
    return (Lattice[row][col]==0)? -1:1;
}

void IsingModel::setSpin(int row, int col, int val) {
    Lattice[row][col]=(val==-1)? 0: 1;
};

int IsingModel::sumNeighbors(int i, int j) {
    int neighbors_rows_cols[4];
    Neighbors(neighbors_rows_cols, i,j);
    int ip,im,jp,jm;
    ip=neighbors_rows_cols[0];
    im=neighbors_rows_cols[1];
    jp=neighbors_rows_cols[2];
    jm=neighbors_rows_cols[3];
    return getSpin(ip,j)+getSpin(im,j)+getSpin(i,jp)+getSpin(i,jm);
};


int IsingModel::NeighborsUp(int i, int j) {
    return sumNeighbors(i,j)/2+2;
}

// Does one sweep of Wolff Algorithm i.e. N*N spins are flipped. Any excess of that is returned.

int IsingModel::SweepWolff(int nTimes, int partialSweep) {
    for (int time=0; time<nTimes; time++) {
        while (partialSweep < N * N) {
            partialSweep += WolffMove();
        }
        partialSweep = partialSweep - N*N;
    }
    return partialSweep;
}

float IsingModel::getMagnetization() {
    float sum=0;
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            sum+= getSpin(i,j);
        }
        
    }
    return sum/(N*N);
}

void IsingModel::display() {
    cout << endl << endl;
    char p[]="+";
    char n[]="-";
    
    for(int i=0;i<N;i++) {
        for (int j=0;j<N;j++) {
            (getSpin(i,j)==1) ? cout << p << " ":cout << n << " " ;
        }
        cout << endl;
    }
    
    cout << endl << endl;
}


float IsingModel::Energy() {
    float E = 0.0;
    for(int i=0;i<N;i++) {
        for (int j=0;j<N;j++) {
            E -= getSpin(i,j) * sumNeighbors(i,j);
        }
    }
    return 0.5*E;
}


IsingModel::~IsingModel() {
    delete heat_bath_prob_up;
    delete metropolis_prob_up[0];
    delete metropolis_prob_up[1];
    delete Lattice;
}



