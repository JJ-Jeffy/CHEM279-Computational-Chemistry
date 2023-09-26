// Problem set 2 
// University of California, Berkeley

// Created by: Jeffy Jeffy
// Date: 9/25/23

// Description: This file contains the function prototypes for the utility functions used in this project

#include <iostream>
#include <armadillo>
#include <cmath>
#include "gaussian.h"
#include "utils.h"

using namespace std;
using namespace arma;

// Function: Evaluate gaussian overlap at a given position
// Input: Two gaussian objects, and a position
// Output: The value of the overlap integral at that position
double overlap(Gaussian g1, Gaussian g2, double position){
    double alpha1 = g1.alpha;
    double alpha2 = g2.alpha;
    int l1 = g1.angular_momentum;
    int l2 = g2.angular_momentum;
    double A = g1.center;
    double B = g2.center;
    double S = 0;
    S += pow((position - A),l1)*pow((position - B), l2)*exp(pow((position - A),2)*(-alpha1))*exp(pow((position - B),2)*(-alpha2));
    return S;
}

// Function: Evaluate gaussian overlap using trapzoidal rule
// Input: Two gaussian objects, and a starting and ending position, and the step size 
// Output: The value of the overlap integral using trapzoidal rule
double trapz_overlap(Gaussian g1, Gaussian g2, double start, double end, int n){
    double h = (end - start)/(n);
    double S = 0;
    double f_a = overlap(g1, g2, start);
    double f_b = overlap(g1, g2, end);
    S = (f_a + f_b)/2;
    for(int i = 1; i < n; i++){
        double f_c = overlap(g1, g2, start + i*h);
        S += f_c;
    }

    S = S*h;
    return S;
}

// Function: Evaluate the double factorial of a number
// Input: A number
// Output: The double factorial of that number
int double_factorial(int n){
    if(n <= 1){
        return 1;
    }
    else{
        return n*double_factorial(n-2);
    }
}

// Function: Evaluate the center of the product of the Gaussian
// Input: Two gaussian objects
// Output: The center of the product of the Gaussian
vec product_center(Gaussian g1, Gaussian g2){
    vec center1 = g1.center_vectors;
    vec center2 = g2.center_vectors;
    vec center = (g1.alpha*center1 + g2.alpha*center2)/(g1.alpha + g2.alpha);
    return center;
}

// Function: Evaluate the exponential prefactor of the product of the Gaussian
// Input: The center of the two gaussian objects, and the alpha of the two gaussian objects
// Output: The exponential prefactor of the product of the Gaussian
double product_exponential_prefactor(double xa, double xb, double alpha1, double alpha2){
    double product = exp((-alpha1*alpha2/(alpha1 + alpha2))*pow((xa - xb),2));
    return product*(sqrt(M_PI/(alpha1 + alpha2)));
}

// Function: Evaluate the factorial of a given number
// Input: A number
// Output: A number
int factorial(int n){
    if(n == 0 || n == 1){
        return 1;
    }
    else{
        return n*factorial(n-1);
    }
}

// Function: Evaluate the binomial prefactor of the product of the Gaussian
// Input: Two inputs numbers (e.g. l, m, and n)
// Output: The binomial prefactor of the product of the Gaussian
int product_binomial_prefactor(int i, int x){
    return (factorial(i))/(factorial(x)*factorial(i - x));
}

// Function: Evaluate the overlap integral at a specific dimension 
// Input: alpha1, alpha2, vec at a specific position
// Output: The overlap at a particular position 
double overlap_at_dimension(double alpha1, double alpha2, double xa, double xb, int l1, int l2, double xp){
    double exponential_prefactor = product_exponential_prefactor(xa, xb, alpha1, alpha2);
    double overlap = 0;
    for(int i = 0; i <= l1; i++){
        for (int j = 0; j <= l2; j++){
            if ((i + j) % 2 == 0){
                double binominal = product_binomial_prefactor(l1, i)*product_binomial_prefactor(l2, j);
                double numerator = double_factorial(i+j-1)*pow(xp-xa,l1-i)*pow(xp-xb,l2-j);
                double denominator = pow(2*(alpha1+alpha2),(i+j)/2);
                overlap += binominal*(numerator/denominator);
            }
        }
    }
    return exponential_prefactor*overlap;
}

// Function: Evaluate the overlap integral for two given gaussians
// Input: Two gaussian objects
// Output: The overlap integral for two given gaussians
void overlap_integral(Gaussian g1, Gaussian g2){
    // product vector of the two centers
    vec prd_center = product_center(g1, g2);

    // get the shell info of the two gaussians
    mat shell_info1 = g1.shell_info;
    mat shell_info2 = g2.shell_info;

    // loop through the matrix and evaluate the overlap integral at each position
    if (g1.angular_momentum == 0 && g2.angular_momentum == 0){
        for(int i=0; i<shell_info1.n_rows; i++){
            double overlap = 1;
            for(int j=0; j<shell_info2.n_cols;j++){
                overlap *= overlap_at_dimension(g1.alpha, g2.alpha, g1.center_vectors(j), g2.center_vectors(j), shell_info1(i,j), shell_info2(i,j), prd_center(j));
            }
            cout << overlap; 
        }
        cout << endl;
    }else if(g1.angular_momentum == 1 && g2.angular_momentum == 0){
        for(int i=0; i<shell_info2.n_rows; i++){
            for(int j=0; j<shell_info1.n_cols;j++){
                double overlap = 1;
                for(int k=0; k<shell_info2.n_cols; k++){
                    overlap *= overlap_at_dimension(g1.alpha, g2.alpha, g1.center_vectors(k), g2.center_vectors(k), shell_info1(j,k), shell_info2(i,k), prd_center(k));
                }
                cout << overlap << "    "; 
            }
            cout << endl;
        }
    }else {
        for(int i=0; i<shell_info1.n_rows; i++){
            for(int j=0; j<shell_info2.n_cols;j++){
                double overlap = 1; 
                for(int k=0; k<shell_info1.n_cols; k++){
                    overlap  *= overlap_at_dimension(g1.alpha, g2.alpha, g1.center_vectors(k), g2.center_vectors(k), shell_info1(i,k), shell_info2(j,k), prd_center(k));
                }
                cout << overlap << "    ";
            }
            cout << endl;
        }
    }
}






