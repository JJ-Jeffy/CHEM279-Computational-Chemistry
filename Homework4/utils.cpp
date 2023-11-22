// Problem set 3
// University of California, Berkeley

// Created by: Jeffy Jeffy
// Date: 10/16/23

// Description: This file contains the code for the utils.h header file

#include "utils.h"

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
double overlap_at_1D(double alpha1, double alpha2, double xa, double xb, int l1, int l2){
    double product_exp = product_exponential_prefactor(xa, xb, alpha1, alpha2);
    double xp = (alpha1*xa + alpha2*xb)/(alpha1 + alpha2);
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
    return product_exp*overlap;
}

//Function: Evaluate the overlap integral for two primitive Gaussian shells
//Input: Two Gaussian objects
//Output: The overlap integral for two primitive Gaussian shells
double overlapIntegral(vec center1, vec center2, double alpha1, double alpha2, vec lmn1, vec lmn2){
    double overlap = 1;
    for(int i = 0; i < 3; i++){
        overlap *= overlap_at_1D(alpha1, alpha2, center1(i), center2(i), lmn1(i), lmn2(i));
    }
    return overlap;
}

// Function: Helper function calculate the 2e integral at a specific dimension
// Input: alpha1, alpha2, vec at a specific position
// Output: The 2e integral at a particular position
double I2e_pG(vec &Ra, vec Rb, double sigmaA, double sigmaB){
    double U = pow(M_PI * sigmaA, 1.5) * pow(M_PI * sigmaB, 1.5);
    double V2 = 1.0/(sigmaA + sigmaB);

    double Rd = norm(Ra - Rb, 2);

    if (Rd == 0){
        return U * sqrt(2*V2)*sqrt(2/M_PI);
    }

    double srT = sqrt(V2)*Rd;

    double result = U/Rd * erf(srT);

    return result;
}