// Problem set 1 
// University of California, Berkeley 

// Creator: Jeffy Jeffy 
// Date: 9/11/23

// Description: This file contains the class definition for utility functions needed for Molecule class.

#include "utils.h"
#include <iostream>

using namespace std;

// Function to print the geometry of a molecule
void print_geom(Molecule mol) {
    cout << "There are " << mol.numAtoms << " atoms." << endl;
    cout << "Coordinates of this molecule are:" << endl;
    for (int i = 0; i < mol.numAtoms; i++) {
        cout << mol.zvals[i] << " " << mol.coordinates[i][0] << " " << mol.coordinates[i][1] << " " << mol.coordinates[i][2] << endl;
    }
}

// Function to print the forces on each atom
void print_forces(vector<vector<double>> forces) {
    for (int i = 0; i < forces.size(); i++) {
        cout << forces[i][0] << " " << forces[i][1] << " " << forces[i][2] << endl;
    }
}

// Function to calculate the distance between two particles 
double calculate_distance(vector<double> p1, vector<double> p2) {
    double distance = 0.0;
    for (int i = 0; i < 3; i++) {
        distance += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    }
    return sqrt(distance);
}

// Function to calculate the lennard jones interaction between two particles 
double calculate_LJ(double rij){
    double sigma = 2.951; 
    double epsilon = 5.29; 
    double r12 = pow(sigma/rij, 12);
    double r6 = pow(sigma/rij, 6);
    return epsilon * (r12 - (2*r6));
}

// Function to calculate the total energy of a molecule
double calculate_energy(Molecule mol) {
    double energy = 0.0;
    for (int i = 0; i < mol.numAtoms; i++) {
        for (int j = i + 1; j < mol.numAtoms; j++) {
            double rij = calculate_distance(mol.coordinates[i], mol.coordinates[j]);
            energy += calculate_LJ(rij);
        }
    }
    return energy;
}

// Function to calculate the forces on each x,y,z in the atom 
double calculate_force(double rij, double x1, double x2){
    double sigma = 2.951;
    double epsilon = 5.29;
    double r12_term = -12*pow(sigma,12)*(pow(1/rij, 13));
    double r6_term = 12*pow(sigma,6)*(pow(1/rij, 7));
    double dist = (x1 - x2)/rij;
    return -1*(epsilon * (r12_term + r6_term) * dist);
}

// Function to calculate the forces on each atom in a molecule
vector<vector<double>> calculate_total_force(Molecule mol){
    vector<vector<double>> forces;
    for (int i = 0; i < mol.numAtoms; i++){
        vector<double> force;
        for (int j = 0; j < 3; j++){
            force.push_back(0.0);
        }
        forces.push_back(force);
    }
    
    for (int i = 0; i <mol.numAtoms; i++){
        for (int j = i + 1; j < mol.numAtoms; j++){
            double rij = calculate_distance(mol.coordinates[i], mol.coordinates[j]);
            double x1 = mol.coordinates[i][0];
            double x2 = mol.coordinates[j][0];
            double y1 = mol.coordinates[i][1];
            double y2 = mol.coordinates[j][1];
            double z1 = mol.coordinates[i][2];
            double z2 = mol.coordinates[j][2];
            forces[i][0] += calculate_force(rij, x1, x2);
            forces[j][0] += calculate_force(rij, x2, x1);
            forces[i][1] += calculate_force(rij, y1, y2);
            forces[j][1] += calculate_force(rij, y2, y1);
            forces[i][2] += calculate_force(rij, z1, z2);
            forces[j][2] += calculate_force(rij, z2, z1);
        }
    }

    return forces;
}

// Calculate the norm of the force 
double calculate_norm(vector<vector<double>> force){
    double norm = 0.0;
    for (int i = 0; i < force.size(); i++){
        for (int j = 0; j < force[i].size(); j++){
            norm += force[i][j] * force[i][j];
        }
    }
    return sqrt(norm);
}

// Calculate the forward difference approximation of force 
vector<vector<double>> calculate_forward_difference(Molecule mol, double h){
    // calculate the distance between two points
    double energy = calculate_energy(mol);
    vector<vector<double>> forces;
    for (int i = 0; i < mol.numAtoms; i++){
        vector<double> force;
        for (int j = 0; j < 3; j++){
            force.push_back(0.0);
        }
        forces.push_back(force);
    }

    // calculate the forces on each atom
    for (int i = 0; i < mol.numAtoms; i++){
        for (int j = 0; j < 3; j++){
            mol.coordinates[i][j] += h;
            double new_energy = calculate_energy(mol);
            forces[i][j] = -(new_energy - energy)/h;
            mol.coordinates[i][j] -= h;
        }
    }

    return forces;
}


// Calculate the central difference approximation of force
vector<vector<double>> calculate_central_difference(Molecule mol, double h){
    vector<vector<double>> forces;
    for (int i = 0; i < mol.numAtoms; i++){
        vector<double> force;
        for (int j = 0; j < 3; j++){
            force.push_back(0.0);
        }
        forces.push_back(force);
    }

    // calculate the forces on each atom
    for (int i = 0; i < mol.numAtoms; i++){
        for (int j = 0; j < 3; j++){
            mol.coordinates[i][j] += h;
            double add_energy = calculate_energy(mol);
            mol.coordinates[i][j] -= 2*h;
            double sub_energy = calculate_energy(mol);
            forces[i][j] = -(add_energy - sub_energy)/(2*h);
            mol.coordinates[i][j] += h;
        }
    }

    return forces;
}

void optimize_mol(Molecule mol, double h){
    vector<vector<double>> forces = calculate_total_force(mol);
    double norm_force = calculate_norm(forces);
    double energy = calculate_energy(mol);
    int iterations = 0; 
    while(norm_force > 0.00001 && iterations < 100000){
        for (int i = 0; i < mol.numAtoms; i++){
            for (int j = 0; j < 3; j++){
                mol.coordinates[i][j] += h * forces[i][j];
            }
        }

        double new_energy = calculate_energy(mol);
        if (new_energy < energy){
            energy = new_energy;
            forces = calculate_total_force(mol);
            norm_force = calculate_norm(forces);
            h = h*1.2; 
            iterations += 1;
        }
        else{
            h = h/2;
            iterations += 1;
        }
    }

    cout << "The optimized geometry is: " << endl;
    print_geom(mol);

    cout << "The optimized energy is: " << endl;
    cout << energy << endl;

    cout << "The optimized forces are: " << endl;
    print_forces(forces);

    cout << "The number of iterations is: " << endl;
    cout << iterations << endl;
}

