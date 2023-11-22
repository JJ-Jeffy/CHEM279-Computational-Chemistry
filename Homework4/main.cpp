// Problem set 3 
// University of California, Berkeley 

// Created by : Jeffy Jeffy 
// date: 10/16/23

// Description: This file contains the main program for this project 

#include "molecule_test.h"
#include "utils.h"
#include <iostream>

using namespace std;

int main(){
    string filename = "C2H4.txt";
    Molecule mol(filename);
    mol.print_geom_info();
    cout << "Molecule was created" << endl;

    //calculate the gamma matrix
    cout << "Calculating the gamma matrix" << endl;
    cout << mol.calcGammaMatrix() << endl;

    cout << "Calculating the overlap integral" << endl;
    cout << mol.overlapMatrix() << endl;

   // calculate the hcore matrix 
    cout << "Calculating the hcore matrix" << endl;
    cout << mol.calcHCoreMatrix() << endl;

    double convergence = 1e-6;
    int maxIterations = 100; 
    mol.runSCF(convergence, maxIterations);

    return 0;
}