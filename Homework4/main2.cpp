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
    cout << "Let's evaluate N2 first" << endl;
    string filename = "N2.txt"; 
    Molecule mol(filename);
    mol.print_geom_info();
    cout << "Molecule " << filename << " was created" << endl;

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

    cout << "Let's evaluate N now" << endl; 
    filename = "N.txt";
    Molecule mol2(filename);
    mol2.print_geom_info();
    cout << "Molecule " << filename << " was created" << endl;

    //calculate the gamma matrix
    cout << "Calculating the gamma matrix" << endl;
    cout << mol2.calcGammaMatrix() << endl;

    cout << "Calculating the overlap integral" << endl;
    cout << mol2.overlapMatrix() << endl;

    // calculate the hcore matrix
    cout << "Calculating the hcore matrix" << endl;
    cout << mol2.calcHCoreMatrix() << endl;

    mol2.runSCF(convergence, maxIterations);

    cout << " The total energy of N2 and N are as follows: " << endl;
    cout <<  "N2: " << mol.totalEnergy << endl;
    cout << "N: " << mol2.totalEnergy << endl;
    cout << "The difference in energy is: " << mol.totalEnergy - (2*mol2.totalEnergy) << endl;

    return 0; 
}