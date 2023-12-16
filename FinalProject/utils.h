// Problem set 1 

// University of California, Berkeley

// Creator: Jeffy Jeffy
// Date: 12/11/23

#pragma once 
#include "molecule.h"
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <armadillo>

using namespace std; 
using namespace arma;

// Function to print the geometry of the molecule
void print_geom(Molecule mol); 

// Function to calculate the distance between two vectors
double calculate_distance(const vec& p1, const vec& p2);

// Function to calculate the lennard jones potential 
double calculate_LJ(double r_ij, double sigma, double epsilon);

// Function to calculate the energy given a molecule object
double calculate_energy(Molecule mol);

// Function to calculate the energy given a matrix of coordinates
double calculate_energy(mat coords, mat sigma, mat epsilon);

// Function to calculate the forces on each atom
double calculate_forces(double rij, double x1, double x2); 

// Function to calculate the forces on the whole molecule 
mat calculate_total_force(Molecule mol);

// Function to calculate the forces given a matrix of coordinates
mat calculate_total_force(mat coords, mat sigma, mat epsilon);

// Function to calculate the forward difference 
mat calculate_forward_difference(Molecule mol, double h); 

// Function to calculate the forward difference given a matrix of coordinates
mat calculate_forward_difference(mat coords, double h);

// Function to calculate the central difference
mat calculate_central_difference(Molecule mol, double h);

// Function to calculate the central difference given a matrix of coordinates
mat calculate_central_difference(mat coords, double h);

// Function to optimize the geometry of the molecule
void optimize_geometry(Molecule& mol, double tol, int maxiter, double h);

// Function to calculate the gradient of the energy
mat calculate_gradient_energy(Molecule& mol); 

// Function to calculate the hessian element 
double calculate_hessian_element(Molecule& mol, double h, int i, int j);

// Function to calculate the hessian element of multiple coordinates
double calculate_mixed_hessian_element(Molecule& mol, double h, int i, int j, int k, int l);

// Function to calculate the hessian matrix 
mat calculate_hessian(Molecule mol, double h);

// Generate mass weighted hessian matrix 
mat generate_mass_weighted_hessian(mat hessian, Molecule& mol); 

// Function to calculate the vibrational frequencies 
vec calculate_vibrational_frequencies(mat mass_weighted_hessian);
