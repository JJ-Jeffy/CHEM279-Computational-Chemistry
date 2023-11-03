// Problem set 3 
// University of California, Berkeley 

// Created by : Jeffy Jeffy 
// date: 10/16/23

// Description: This file contains the main program for this project 

#include "molecule.h"
#include "utils.h"
#include <iostream>

using namespace std;

int main(){
    string filename = "C2H4.txt";
    Molecule mol(filename);
    cout << "Molecule was created" << endl;
    mol.print_geom_info();
    
    mat overlap = mol.overlapMatrix();

    cout << "The various matrices for "<< filename << "is " << mol.calcEnergy() << endl;

    return 0;
}