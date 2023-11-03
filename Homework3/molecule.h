// Problem set 3 
// University of California, Berkeley 

// Created by : Jeffy Jeffy 
// date: 10/16/23 

// Description: This header file contains the header for the class molecule along with structs that define shells and orbitals.

#pragma once
#include <armadillo>
#include <iostream>
#include <vector>
#include <string>
#include "utils.h"

using namespace std;
using namespace arma;

struct basisFunction{
    string name; 
    vec center; 
    vec quantumNumbers;
    vec exponents;
    vec coefficients;
    vec normalization;
};

class Molecule{
    public:
        int numAtoms;
        vec z_val;
        mat coords;
        int numBasisSets; 
        int numElectrons;
        vector<basisFunction> basisFunctionList;
        int countBasisSets();
        int countElectrons();

        // Function: Calculate the Normalization Constants for a given basis function 
        // Input: A basis function 
        // Output: The normalization constant for that basis function
        void calculateNormalization(basisFunction &basisFunction);

        // Function: Calculate the Normalization Constants for all of the basis functions
        // Input: A vector of basis functions
        // Output: The normalization constants for all of the basis functions
        void calcNormalizationAll(vector<basisFunction>& basisFunctionList);

        // Constructor for the molecule class
        Molecule(string filename);

        // print the geom info
        void print_geom_info();

        // quantum numbers for the basis functions
        vec s_lmn = {0, 0, 0};
        vec px_lmn = {1, 0, 0};
        vec py_lmn = {0, 1, 0};
        vec pz_lmn = {0, 0, 1};

        // exponents for the basis functions
        vec s_exponents = {3.42525091, 0.62391373, 0.16885540};
        vec p_exponents = {2.9412494, 0.6834831, 0.2222899};

        // coefficients for the basis functions
        vec s_coefficients = {0.15432897, 0.53532814, 0.44463454};
        vec twos_coefficients = {-0.09996723, 0.39951283, 0.70011547};
        vec p_coefficients = {0.15591627, 0.60768372, 0.39195739};

        // build the basis functions
        vector<basisFunction> buildBasisFunctions();

        //Function: Calculate the contracted overlap integral for two primitive Gaussian shells
        //Input: Two basis functions 
        //Output: The contracted overlap integral for two primitive Gaussian shells
        double contractedOverlapIntegral(basisFunction basisFunction1, basisFunction basisFunction2);

        // Function: Calculate the entire overlap matrix 
        // Input: a molecule 
        // Output: The overlap matrix
        mat overlapMatrix();

        // create a map that stores the hamiltonian value for the respective atom type 
        map<string, double> hamiltonianVal = {{"1s", -13.6},{"2s", -21.4},{"2px", -11.4},{"2py", -11.4},{"2pz", -11.4}};

        // Function: Create the hamiltonian matrix based on the basis functions and the overlap matrix
        double calcEnergy();

};