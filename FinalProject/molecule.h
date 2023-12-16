// Problem set 1 
// University of California, Berkeley 

// Creator: Jeffy Jeffy 
// Date: 9/11/23

// Description: This file contains the class definition for a molecule.

#pragma once
#include <iostream> 
#include <armadillo>
#include <vector>

using namespace std; 
using namespace arma; 

class Molecule{
    public:
        int numAtoms;
        string atomName;
        vector<string> atoms; 
        vec zvals; 
        mat coords; 
        vec masses;
        mat sigma; 
        mat epsilon;


    // Constructor
    Molecule(string filename);
};