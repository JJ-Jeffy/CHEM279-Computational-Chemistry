// Problem set 1 
// University of California, Berkeley 

// Creator: Jeffy Jeffy 
// Date: 9/11/23

// Description: This file contains the class definition for a molecule.

#pragma once 
#include <iostream> 
#include <cmath> 
#include <string>
#include <vector>
#include <fstream>

using namespace std;

class Molecule {
public:
    int numAtoms;
    vector<vector<double>> coordinates; 
    vector<int> zvals; 
    
    //constructor 
    Molecule(const char* filename);
};