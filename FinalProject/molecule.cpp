// Problem set 1 
// University of California, Berkeley 

// Creator: Jeffy Jeffy 
// Date: 12/07/23 

// Description: This file contains the class definition for a molecule.

// molecule.h contains the class headers for the Molecule class.

#include "molecule.h"
#include <fstream>

// Constructor for the molecule class 
Molecule::Molecule(string filename){
    ifstream infile(filename); 
    if (!infile){
        cout << "File does not exist" << endl; 
        exit(1); 
    }

    // map for the elemts 
    map<string, int> element = {{"H",1}, {"C",6}, {"N",7}, {"O",8}, {"F",9}, {"Au", 79}};
    map<string, double> atom_mass = {{"H",1.0079}, {"C",12.0107}, {"N",14.0067}, {"O",15.9994}, {"F",18.9984}, {"Au", 196.9665}};

    infile >> numAtoms; 
    zvals = zeros<vec>(numAtoms);
    coords = zeros<mat>(numAtoms, 3);
    masses = zeros<vec>(numAtoms);

    for(int i = 0; i < numAtoms; i++){
        infile >> atomName >> coords(i,0) >> coords(i,1) >> coords(i,2);
        atoms.push_back(atomName);
        zvals(i) = element[atomName];
        masses(i) = atom_mass[atomName];
    }
    infile.close();

    if (filename == "Au3.txt"){
        // fill all the values in sigm matrix with 2.951
        sigma = zeros<mat>(numAtoms, numAtoms);
        sigma.fill(2.951);

        // fill all the values in epsilon matrix with 0.159
        epsilon = zeros<mat>(numAtoms, numAtoms);
        epsilon.fill(5.29); 
    }

    if (filename == "Au4.txt"){
        // fill all the values in sigm matrix with 2.951
        sigma = zeros<mat>(numAtoms, numAtoms);
        sigma.fill(2.951);

        // fill all the values in epsilon matrix with 0.159
        epsilon = zeros<mat>(numAtoms, numAtoms);
        epsilon.fill(5.29); 
    }

    if (filename == "C2H4.txt"){
        // fill all the values in sigm matrix with 2.951
        sigma = zeros<mat>(numAtoms, numAtoms);
        sigma = {{3.52, 3.52, 2.72, 2.72, 2.72, 2.72},
                 {3.52, 3.52, 2.72, 2.72, 2.72, 2.72},
                 {2.72, 2.72, 2.72, 2.72, 2.72, 2.72},
                 {2.72, 2.72, 2.72, 2.72, 2.72, 2.72},
                 {2.72, 2.72, 2.72, 2.72, 2.72, 2.72},
                 {2.72, 2.72, 2.72, 2.72, 2.72, 2.72}};

        // fill all the values in epsilon matrix with 0.159
        epsilon = zeros<mat>(numAtoms, numAtoms);
        epsilon = {{0.066, 0.066, 0.05, 0.05, 0.05, 0.05},
                   {0.066, 0.066, 0.05, 0.05, 0.05, 0.05},
                   {0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
                   {0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
                   {0.05, 0.05, 0.05, 0.05, 0.05, 0.05},
                   {0.05, 0.05, 0.05, 0.05, 0.05, 0.05}};
    }

    if (filename == "C2H2.txt"){
        // fill all the values in sigm matrix with 2.951
        sigma = zeros<mat>(numAtoms, numAtoms);
        sigma = {{2.55, 2.72, 2.72, 2.55},
                 {2.72, 3.52, 3.52, 2.72},
                 {2.72, 3.52, 3.52, 2.72},
                 {2.55, 2.72, 2.72, 2.55}};

        // fill all the values in epsilon matrix with 0.159
        epsilon = zeros<mat>(numAtoms, numAtoms);
        epsilon = {{0.015, 0.05, 0.05, 0.015},
                   {0.05, 0.066, 0.066, 0.05},
                   {0.05, 0.066, 0.066, 0.05},
                   {0.015, 0.05, 0.05, 0.015}};
    }

    if (filename == "H2.txt"){
        // fill all the values in sigm matrix with 2.951
        sigma = zeros<mat>(numAtoms, numAtoms);
        sigma.fill(2.55);

        // fill all the values in epsilon matrix with 0.159
        epsilon = zeros<mat>(numAtoms, numAtoms);
        epsilon.fill(0.015); 
    }
}
 

