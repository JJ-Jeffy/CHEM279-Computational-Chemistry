//Problem set 2
//University of California, Berkeley 

// Created by: Jeffy Jeffy 
// Date: 9/25/23 


// Description: This file contains the function prototypes for the utility functions used in this project 

#pragma once
#include <iostream>
#include <armadillo>
#include "gaussian.h"

using namespace arma;
using namespace std;

// Function: Evaluate gaussian overlap at a given position 
// Input: Two gaussian objects, and a position
// Output: The value of the overlap integral at that position
double overlap(Gaussian g1, Gaussian g2, double position);

// Function: Evaluate the gaussian using trapzoidal rule
// Input: Two gaussian objects, and a starting and ending position and the step size
// Output: The value of the overlap integral using trapzoidal rule
double trapz_overlap(Gaussian g1, Gaussian g2, double start, double end, int n);

// Function: Evaluate the double factorial of a number
// Input: A number
// Output: The double factorial of that number
int double_factorial(int n);

// Function: Evaluate the center of the product of the Gaussian
// Input: Two gaussian objects
// Output: The center of the product of the Gaussian
vec product_center(Gaussian g1, Gaussian g2);

// Function: Evaluate the exponential prefactor of the product of the Gaussian
// Input: The center of the two gaussian objects, and the alpha of the two gaussian objects
// Output: The exponential prefactor of the product of the Gaussian
double product_exponential_prefactor(double xa, double xb, double alpha1, double alpha2);

// Function: Evaluate the factorial of a given number 
// Input: A number
// Output: A number 
int factorial(int n); 

// Function: Evaluate the binomial prefactor of the product of the Gaussian
// Input: Two inputs numbers (e.g. l, m, and n)
// Output: The binomial prefactor of the product of the Gaussian
int product_binomial_prefactor(int i, int x); 

// Function: Evaluate the overlap integral at a specific dimension 
// Input: alpha1, alpha2, vec at a specific position
// Output: The overlap at a particular position 
double overlap_at_dimension(double alpha1, double alpha2, double xa, double xb, int l1, int l2, double xp); 

// Function: Evaluate the overlap integral for two given gaussians
// Input: Two gaussian objects
// Output: The overlap integral for two given gaussians
void overlap_integral(Gaussian g1, Gaussian g2);