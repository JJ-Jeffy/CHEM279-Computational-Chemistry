// Problem set 1 
// University of California, Berkeley 

// Creator: Jeffy Jeffy 
// Date: 9/11/23

// Description: This file contains the main class for the program.

#include <iostream>
#include <array>
#include <fstream>
#include "utils.h"
#include "molecule.h"

using namespace std;

int main() {
    // Create a molecule object
    Molecule mol("geom_1.txt");

    // print the geometry of the molecule
    print_geom(mol);

    // calculate the energy of the molecule
    double energy = calculate_energy(mol);

    std::cout << "The energy is " << energy << " hartrees." << std::endl;

    // calculate the forces on each atom
    vector<vector<double>> forces = calculate_total_force(mol);

    cout << "The forces on each atom are: " << endl;
    print_forces(forces);

    // calculate the forward difference approximation of the force on each atom
    vector<vector<double>> forward_difference = calculate_forward_difference(mol, 0.1);

    cout << "The forward difference approximation of the forces on each atom are: " << endl;
    print_forces(forward_difference);

    // calculate the central difference approximation of the force on each atom
    vector<vector<double>> central_difference = calculate_central_difference(mol, 0.1);

    cout << "The central difference approximation of the forces on each atom are: " << endl;
    print_forces(central_difference);

    // loop through the given h and calculate the error between the analytical and numerical forces 
    // store the log of the error and the log of the step size in a file
    double step_size[5] = { 0.1, 0.01, 0.001, 0.0001, 0.00001 };

    ofstream myfile;
    myfile.open("error.txt");
    myfile << "log(h) " << "log(error_forward) " << "log(error_central) " << endl;
    for (int i = 0; i < 5; i++) {
        vector<vector<double>> forward_difference = calculate_forward_difference(mol, step_size[i]);
        vector<vector<double>> central_difference = calculate_central_difference(mol, step_size[i]);
        double norm_forward = calculate_norm(forward_difference);
        double norm_central = calculate_norm(central_difference);
        double norm_force = calculate_norm(forces);
        double error_forward = abs(norm_forward - norm_force); 
        double error_central = abs(norm_central - norm_force);
        myfile << log(step_size[i]) << " " << log(error_forward) << " " << log(error_central) << endl;
    }
    myfile.close();

    // loop through the molecule and optimize it using line search and steepest descent
    optimize_mol(mol, 0.01);

    return 0;
}