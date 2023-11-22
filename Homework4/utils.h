// Problem set 3 
// University of California, Berkeley 
// Created by: Jeffy Jeffy 

// Date: 10/16/23

// Description: this header file contains the utility functions 
//needed for this project along with some functions needed for the molecule class

#pragma once 
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

// Function: Evaluate the double factorial of a number
// Input: A number
// Output: The double factorial of that number
int double_factorial(int n);

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
double overlap_at_1D(double alpha1, double alpha2, double xa, double xb, int l1, int l2);

//Function: Evaluate the overlap integral for two primitive Gaussian shells
//Input: Two Gaussian objects
//Output: The overlap integral for two primitive Gaussian shells
double overlapIntegral(vec center1, vec center2, double alpha1, double alpha2, vec lmn1, vec lmn2);

// Function: Helper function calculate the 2e integral at a specific dimension
// Input: alpha1, alpha2, vec at a specific position
// Output: The 2e integral at a particular position
double I2e_pG(vec &Ra, vec Rb, double sigmaA, double sigmaB);