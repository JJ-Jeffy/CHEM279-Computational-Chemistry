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
    cout << "1D Gaussian Overlap Integral Calculator" << endl;
    cout << "Please enter the information for your first gaussian: " << endl;
    double alpha1, center1;
    int angular_momentum1;
    cout << "Alpha: ";
    cin >> alpha1;
    cout << "Angular momentum: ";
    cin >> angular_momentum1;
    cout << "Center: ";
    cin >> center1;

    cout << "Please enter the information for your second gaussian: " << endl;
    double alpha2, center2;
    int angular_momentum2;
    cout << "Alpha: ";
    cin >> alpha2;
    cout << "Angular momentum: ";
    cin >> angular_momentum2;
    cout << "Center: ";
    cin >> center2;

    Gaussian g1(alpha1, angular_momentum1, center1);
    Gaussian g2(alpha2, angular_momentum2, center2);

    cout << "Please enter the starting and ending position for the overlap integral: " << endl;
    double start, end;
    cout << "Start: ";
    cin >> start;
    cout << "End: ";
    cin >> end;

    cout << "Please enter the number of steps for the trapzoidal rule: " << endl;
    int n;
    cin >> n;

    double S = trapz_overlap(g1, g2, start, end, n);
    cout << "The value of the overlap integral is: " << S << endl;

    return 0;
}
