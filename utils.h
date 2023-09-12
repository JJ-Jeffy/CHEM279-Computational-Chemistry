// Problem set 1 
// University of California, Berkeley 

// Creator: Jeffy Jeffy 
// Date: 9/11/23

// Description: This file contains the class definition for utility functions needed for Molecule class.

#pragma once
#include "molecule.h"
#include <iostream>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

// Function to print the geometry of a molecule
void print_geom(Molecule mol); 

// Function to print the forces on each atom
void print_forces(vector<vector<double>> forces);

// Function to calculate the distance between two particles 
double calculate_distance(vector<double> p1, vector<double> p2);

// Function to calculate the lennard jones interaction between two particles
double calculate_LJ(double rij);

// Function to calculate the total energy of a molecule
double calculate_energy(Molecule mol);

// Function to calculate the forces on each x,y,z in the atom
double calculate_force(double rij, double x1, double x2);

// Function to calculate the forces on each atom
vector<vector<double>> calculate_total_force(Molecule mol);

// Calculate the norm of the force 
double calculate_norm(vector<vector<double>> force);

// Calculate the forward difference approximation of force 
vector<vector<double>> calculate_forward_difference(Molecule mol, double h);

// Calculate the central difference approximation of force
vector<vector<double>> calculate_central_difference(Molecule mol, double h);

// Calculate the difference between the analytical and numerical forces
vector<vector<double>> calculate_difference(vector<vector<double>> analytical, vector<vector<double>> numerical);

// loop through the given convergence criteria and update the coordinates of the molecule
void optimize_mol(Molecule mol, double h);
