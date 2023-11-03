// Problem set 3
// University of California, Berkeley

// Created by: Jeffy Jeffy
// Date: 10/16/23

// Description: This file contains the code for the molecule.h header file 

#include "molecule.h"
#include <fstream>

// Contructor for the molecule class 
Molecule::Molecule(string filename){
    ifstream infile(filename);
    if(!infile){
        cout << "Error opening file" << endl;
        exit(1);
    }
    infile >> numAtoms;
    z_val = zeros<vec>(numAtoms);
    coords = zeros<mat>(numAtoms, 3);

    for(int i = 0; i < numAtoms; i++){
        infile >> z_val(i) >> coords(i, 0) >> coords(i, 1) >> coords(i, 2);
    }

    infile.close();

    // count the number of basis sets and electrons
    numBasisSets = countBasisSets();
    numElectrons = countElectrons();

    // build the basis functions
    basisFunctionList = buildBasisFunctions();

    // normalize the basis functions 
    calcNormalizationAll(basisFunctionList);
}

// Count the number of basis sets
int Molecule::countBasisSets(){
    int count = 0;
    for(int i = 0; i < numAtoms; i++){
        if(z_val(i) == 1){
            count += 1;
        }
        else if(z_val(i) == 6){
            count += 4;
        }
    }
    return count;
}

// Count the number of electrons 
int Molecule::countElectrons(){
    int count = 0;
    for(int i = 0; i < numAtoms; i++){
        if(z_val(i) == 1){
            count += 1;
        }
        else if(z_val(i) == 6){
            count += 4;
        }
    }

    if(count % 2 == 0){
        return count / 2;
    }
    else{
        cout << "Error: Odd number of electrons" << endl;
        exit(1);
    }
}

// Build basis functions 
vector<basisFunction> Molecule::buildBasisFunctions(){
    vector<basisFunction> basisFunctionList;

    // Loop over all the atoms and build the basis functions for hydrogen and carbon 
    for(int i = 0; i < numAtoms; i++){
        if(z_val(i) == 1){
            struct basisFunction s = {"1s", coords.row(i).t(), s_lmn, s_exponents, s_coefficients, zeros<vec>(3)};
            basisFunctionList.push_back(s);
        } else if(z_val(i) == 6){
            struct basisFunction twos = {"2s", coords.row(i).t(), s_lmn, s_exponents, twos_coefficients, zeros<vec>(3)};
            struct basisFunction px = {"2px", coords.row(i).t(), px_lmn, p_exponents, p_coefficients, zeros<vec>(3)};
            struct basisFunction py = {"2py", coords.row(i).t(), py_lmn, p_exponents, p_coefficients, zeros<vec>(3)};
            struct basisFunction pz = {"2pz", coords.row(i).t(), pz_lmn, p_exponents, p_coefficients, zeros<vec>(3)};
            basisFunctionList.push_back(twos);
            basisFunctionList.push_back(px);
            basisFunctionList.push_back(py);
            basisFunctionList.push_back(pz);
        }
    }
    return basisFunctionList;
}

// Function that prints the geometry of the molecule
void Molecule::print_geom_info(){
    cout << "The molecule has " << numAtoms << " atoms" << endl;
    cout << "The coordinates of the atoms are: " << endl;
    cout << coords << endl;
    cout << "The atomic numbers of the atoms are: " << endl;
    cout << z_val << endl;
    cout << "The number of basis sets is: " << numBasisSets << endl;
    cout << "The number of electrons is: " << numElectrons << endl;
}

// Function: Calculate the Normalization Constants for a given basis function 
// Input: A basis function 
// Output: The normalization constant for that basis function
void Molecule::calculateNormalization(basisFunction &basisFunction){
    // Calculate the normalization constant for the basis function
    double norm1 = overlapIntegral(basisFunction.center, basisFunction.center, basisFunction.exponents(0), basisFunction.exponents(0), basisFunction.quantumNumbers, basisFunction.quantumNumbers);
    double norm2 = overlapIntegral(basisFunction.center, basisFunction.center, basisFunction.exponents(1), basisFunction.exponents(1), basisFunction.quantumNumbers, basisFunction.quantumNumbers);
    double norm3 = overlapIntegral(basisFunction.center, basisFunction.center, basisFunction.exponents(2), basisFunction.exponents(2), basisFunction.quantumNumbers, basisFunction.quantumNumbers);

    norm1 = pow(norm1, -0.5);
    norm2 = pow(norm2, -0.5);
    norm3 = pow(norm3, -0.5);

    basisFunction.normalization = {norm1, norm2, norm3};
}

// Function: Calculate the Normalization Constants for all of the basis functions
// Input: A vector of basis functions
// Output: The normalization constants for all of the basis functions
void Molecule::calcNormalizationAll(vector<basisFunction>& basisFunctionList){
    for(int i = 0; i < basisFunctionList.size(); i++){
        calculateNormalization(basisFunctionList[i]);
    }
}

//Function: Calculate the contracted overlap integral for two primitive Gaussian shells
//Input: Two basis functions 
//Output: The contracted overlap integral for two primitive Gaussian shells
double Molecule::contractedOverlapIntegral(basisFunction basisFunction1, basisFunction basisFunction2){
    double overlap = 0; 

    // loop over all the exponent combinations 
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            overlap += basisFunction1.coefficients(i)*basisFunction2.coefficients(j)*basisFunction1.normalization(i)*basisFunction2.normalization(j)*
                overlapIntegral(basisFunction1.center, basisFunction2.center, basisFunction1.exponents(i), basisFunction2.exponents(j), basisFunction1.quantumNumbers, basisFunction2.quantumNumbers);
        }
    }

    return overlap;
}

// Function: Calculate the entire overlap matrix 
// Input: a molecule 
// Output: The overlap matrix
mat Molecule::overlapMatrix(){
    mat overlap = arma::zeros<mat>(numBasisSets, numBasisSets);
    for(int i = 0; i < numBasisSets; i++){
        for(int j = 0; j < numBasisSets; j++){
            overlap(i, j) = contractedOverlapIntegral(basisFunctionList[i], basisFunctionList[j]);
        }
    }
    return overlap;
}

// Function: Calculate the energy 
// Input: a molecule
double Molecule::calcEnergy(){

    // calculate the overlap matrix 
    mat overlapMat = overlapMatrix();
    cout << "OV Matrix: " << endl;
    cout << overlapMat << endl;

    // First calculate the hamiltonian matrix 
    mat hamiltonian = arma::zeros<mat>(numBasisSets, numBasisSets);
    // loop over all the basis functions
    for(int i=0; i<numBasisSets; i++){
        for(int j=0; j<numBasisSets; j++){
            // get the hamiltonian matrix value based on the atom type 
            string name1 = basisFunctionList[i].name;
            string name2 = basisFunctionList[j].name;

            double h_ij = hamiltonianVal[name1];
            double h_ji = hamiltonianVal[name2];


            // assign diagonal elements
            if(i == j){
                hamiltonian(i, j) = h_ij;
            }
            else{
                hamiltonian(i, j) = 1.75/2 * (h_ij + h_ji) * overlapMat(i,j);
            }
        }
    }

    cout << "Hamiltonian Matrix: " << endl;
    cout << hamiltonian << endl;

    // calculate the X transformation matrix 
    mat X = arma::zeros<mat>(numBasisSets, numBasisSets);

    // diagnolize the overlap matrix
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, overlapMat);

    // get inverse of square root of eigenvalues
    vec eigval_inv = 1/sqrt(eigval);

    // diagonalize the eigval_inv matrix
    mat eigval_inv_mat = diagmat(eigval_inv);

    // calculate the X matrix
    X = eigvec * eigval_inv_mat * eigvec.t();

    cout << "X Matrix: " << endl;
    cout << X << endl;

    // Calculate the Prime Matrix 
    mat prime = arma::zeros<mat>(numBasisSets, numBasisSets);

    prime = X.t() * hamiltonian * X;

    cout << "Prime Matrix: " << endl;
    cout << prime << endl;

    // Calculate the energy 
    // diagnolize the hamiltonian prime matrix 
    vec eigval_prime;
    mat eigvec_prime;
    eig_sym(eigval_prime, eigvec_prime, prime);

    // form the MO coefficients matrix
    mat C = X * eigvec_prime;

    cout << "MO Coefficients Matrix: " << endl;
    cout << C << endl;

    cout << "Energy Vec is: " << endl;
    cout << eigval_prime << endl;

    // calculate the energy
    double energy = 0;
    for(int i = 0; i < numElectrons; i++){
        energy += eigval_prime(i);
    }

    cout << "and the energy is: " << endl;
    return energy*2;
}

