//
//  IsingModel.h
//  Ising Model Command Line
//
//  Created by Borun Chowdhury on 4/26/15.
//  Copyright (c) 2015 Borun Chowdhury and Mouchumi Bhuyan. All rights reserved.
//

#ifndef Ising_Model_Command_Line_IsingModel_h
#define Ising_Model_Command_Line_IsingModel_h


#endif

#include <iostream>
#include <random>
using namespace std;

class Spin{
public:
    int row, col;
};


class IsingModel{
    
    
private:
    bool **Lattice;
    int D;
    int N;
    float T;
    float H;
    float *heat_bath_prob_up;
    float *metropolis_prob_up[2];
    float WolffP;
    
    random_device rd;
        
    
    void initializeWolf();
    void initializeHeatBath();
    void initializeMetropolis();
    
    void flipSpin(int i,int j);
    
    
  
    
    //int NeighborsUp(int i, int j);
    
    int *Neighbors(int *neigbor_rows_cols, int i, int j );
    
    void RandomRealUniforms(float *randomArr);
    
    void RandomLocations(int *rows, int *cols);
    
    
    int NeighborsUp(int i, int j);
    
    int sumNeighbors(int row, int col);
    
    inline int getSpin(int row, int col);
    
    void setSpin(int row, int col, int val);
    
public:
    IsingModel(int n, float t, int d=2);
    
    
    
    
    
    int SweepWolff(int nTimes, int partialSweep);
    
    int WolffMove();
    
    void SweepHeatBath(int nTimes);
    void SweepMetropolis(int nTimes=1);
    
    float getMagnetization();
    
    void display();
    
    float Energy();
    void SaveLattice();

    ~IsingModel();

};
