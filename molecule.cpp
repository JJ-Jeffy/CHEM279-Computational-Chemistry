// Problem set 1 
// University of California, Berkeley 

// Creator: Jeffy Jeffy 
// Date: 9/11/23

// Description: This file contains the class definition for a molecule.

// molecule.h contains the class headers for the Molecule class.

#include "molecule.h"
#include <iostream>
#include <fstream>

using namespace std; 

Molecule::Molecule(const char* filename) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        cout << "Error opening file" << endl;
        exit(1);
    }
    infile >> numAtoms;
    zvals.resize(numAtoms);
    coordinates.resize(numAtoms, vector<double>(3));

    for (int i = 0; i < numAtoms; i++) {
        infile >> zvals[i] >> coordinates[i][0] >> coordinates[i][1] >> coordinates[i][2];
        cout << zvals[i] << endl;

        // check if zvals is anything other than gold 
        if (zvals[i] != 79) {
            cout << "Error: zvals is not gold" << endl;
            exit(1);
        }
    }
       
    infile.close();
}