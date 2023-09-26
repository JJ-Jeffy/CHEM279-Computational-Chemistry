// Problem set 2 
// University of California, Berkeley

// Created by: Jeffy Jeffy
// Date: 9/25/23

#include <iostream>
#include <armadillo>
#include "gaussian.h"
#include "utils.h"

using namespace std; 
using namespace arma;

int main(){
    cout << "3D Gaussian Overlap Integral Calculator" << endl;
    cout << "Please enter the information for your first gaussian: " << endl;
    double alpha1, center1;
    int angular_momentum1;
    vec center_vectors1(3);
    cout << "Alpha: ";
    cin >> alpha1;
    cout << "Angular momentum: ";
    cin >> angular_momentum1;
    cout << "Center: ";
    cin >> center_vectors1(0) >> center_vectors1(1) >> center_vectors1(2);

    cout << "Please enter the information for your second gaussian: " << endl;
    double alpha2, center2;
    int angular_momentum2;
    vec center_vectors2(3);
    cout << "Alpha: ";
    cin >> alpha2;
    cout << "Angular momentum: ";
    cin >> angular_momentum2;
    cout << "Center: ";
    cin >> center_vectors2(0) >> center_vectors2(1) >> center_vectors2(2);

    Gaussian g1(alpha1, angular_momentum1, center_vectors1);
    Gaussian g2(alpha2, angular_momentum2, center_vectors2);

    cout << "The result for the overlap integral is: " << endl;
    overlap_integral(g1, g2);

    return 0;
}
