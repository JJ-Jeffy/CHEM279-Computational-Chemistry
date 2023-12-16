// Problem set 1 
// University of California, Berkeley 

// Creator: Jeffy Jeffy 
// Date: 9/11/23

// Description: This file contains the main class for the program.

#include <iostream>
#include <array>
#include <fstream>
#include "molecule.h"
#include "utils.h"

using namespace std;

int main() {
    // Create a molecule object
    string filename; 
    cout << "Enter the name of the file: " << endl;
    cin >> filename;

    Molecule mol(filename);
    // // Print the molecule
    // print_geom(mol);

    // // Calculate the energy of the molecule
    // double energy = calculate_energy(mol);
    // cout << "The energy of the molecule is: " << energy << endl;

    // // Calculate the forces on the molecule
    // mat forces = calculate_total_force(mol);
    // cout << "The forces on the molecule are: " << endl;
    // cout << forces << endl;

    // // Calculate the forward difference
    // double h = 0.0001;
    // mat forward_difference = calculate_forward_difference(mol, h);
    // cout << "The forward difference is: " << endl;
    // cout << forward_difference << endl;

    // // Calculate the central difference
    // mat central_difference = calculate_central_difference(mol, h);
    // cout << "The central difference is: " << endl;
    // cout << central_difference << endl;

    // // print out the geometry of the molecule
    // cout << "The starting geometry of the molecule is: " << endl;
    // cout << mol.coords << endl;

    // // Optimize the geometry of the molecule
    // double tol = 0.0001;
    // int maxiter = 100;
    // optimize_geometry(mol, tol, maxiter, h);

    // // Print out the optimized geometry of the molecule
    // cout << "The optimized geometry of the molecule is: " << endl;
    // cout << mol.coords << endl;

    // // // Calculate the hessian matrix
    // mat hessian = calculate_hessian(mol, h);
    // cout << "The hessian matrix is: " << endl;
    // cout << hessian << endl;
        arma::mat hessian = {
        {4.91131634e-01, 0.00000000e+00, 0.00000000e+00, -4.91875163e-01, 0.00000000e+00, 0.00000000e+00, 7.35149743e-04, 0.00000000e+00, 0.00000000e+00, 8.37914092e-06, 0.00000000e+00, 0.00000000e+00},
        {0.00000000e+00, -1.45091454e-04, 0.00000000e+00, 0.00000000e+00, 2.51938938e-04, 0.00000000e+00, 0.00000000e+00, -1.05649914e-04, 0.00000000e+00, 0.00000000e+00, -1.19756821e-06, 0.00000000e+00},
        {0.00000000e+00, 0.00000000e+00, -1.45091454e-04, 0.00000000e+00, 0.00000000e+00, 2.51938938e-04, 0.00000000e+00, 0.00000000e+00, -1.05649914e-04, 0.00000000e+00, 0.00000000e+00, -1.19756821e-06},
        {-4.91875163e-01, 0.00000000e+00, 0.00000000e+00, 8.82958470e-01, 0.00000000e+00, 0.00000000e+00, -3.91834606e-01, 0.00000000e+00, 0.00000000e+00, 7.51299118e-04, 0.00000000e+00, 0.00000000e+00},
        {0.00000000e+00, 2.51938938e-04, 0.00000000e+00, 0.00000000e+00, -5.40836323e-04, 0.00000000e+00, 0.00000000e+00, 3.96882650e-04, 0.00000000e+00, 0.00000000e+00, -1.07985265e-04, 0.00000000e+00},
        {0.00000000e+00, 0.00000000e+00, 2.51938938e-04, 0.00000000e+00, 0.00000000e+00, -5.40836323e-04, 0.00000000e+00, 0.00000000e+00, 3.96882650e-04, 0.00000000e+00, 0.00000000e+00, -1.07985265e-04},
        {7.35149743e-04, 0.00000000e+00, 0.00000000e+00, -3.91834606e-01, 0.00000000e+00, 0.00000000e+00, 8.82974619e-01, 0.00000000e+00, 0.00000000e+00, -4.91875163e-01, 0.00000000e+00, 0.00000000e+00},
        {0.00000000e+00, -1.05649914e-04, 0.00000000e+00, 0.00000000e+00, 3.96882650e-04, 0.00000000e+00, 0.00000000e+00, -5.43171673e-04, 0.00000000e+00, 0.00000000e+00, 2.51938935e-04, 0.00000000e+00},
        {0.00000000e+00, 0.00000000e+00, -1.05649914e-04, 0.00000000e+00, 0.00000000e+00, 3.96882650e-04, 0.00000000e+00, 0.00000000e+00, -5.43171673e-04, 0.00000000e+00, 0.00000000e+00, 2.51938935e-04},
        {8.37914092e-06, 0.00000000e+00, 0.00000000e+00, 7.51299118e-04, 0.00000000e+00, 0.00000000e+00, -4.91875163e-01, 0.00000000e+00, 0.00000000e+00, 4.91115485e-01, 0.00000000e+00, 0.00000000e+00},
        {0.00000000e+00, -1.19756821e-06, 0.00000000e+00, 0.00000000e+00, -1.07985265e-04, 0.00000000e+00, 0.00000000e+00, 2.51938935e-04, 0.00000000e+00, 0.00000000e+00, -1.42756104e-04, 0.00000000e+00},
        {0.00000000e+00, 0.00000000e+00, -1.19756821e-06, 0.00000000e+00, 0.00000000e+00, -1.07985265e-04, 0.00000000e+00, 0.00000000e+00, 2.51938935e-04, 0.00000000e+00, 0.00000000e+00, -1.42756104e-04}
    };


    //    Get the dimensions of the numpy matrix
    // const size_t rows = sizeof(numpyMatrix) / sizeof(numpyMatrix[0]);
    // const size_t cols = sizeof(numpyMatrix[0]) / sizeof(numpyMatrix[0][0]);

    // Create an Armadillo matrix and copy data from the Numpy matrix
    // arma::mat hessian(reinterpret_cast<double*>(numpyMatrix), rows, cols, false, true);


    // Generate the mass weighted hessian matrix 
    mat mass_weighted_hessian = generate_mass_weighted_hessian(hessian, mol);
    cout << "The mass weighted hessian matrix is: " << endl;
    cout << mass_weighted_hessian << endl;

    // Calculate the eigenvalues and eigenvectors of the mass weighted hessian matrix
    vec eigenvalues;
    mat eigenvectors;
    eig_sym(eigenvalues, eigenvectors, hessian );
    cout << "The eigenvalues of the mass weighted hessian matrix are: " << endl;
    cout << eigenvalues << endl;
    cout << "The eigenvectors of the mass weighted hessian matrix are: " << endl;
    cout << eigenvectors << endl;

    return 0;
}