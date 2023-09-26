// Problem set 2
// University of California, Berkeley 

// Created by: Jeffy Jeffy
// Date: 9/25/23

// Description: This is a header file for the Gaussian class. It contains properties for various gaussian shells 

#pragma once
#include <iostream>
#include <armadillo>

using namespace arma;

class Gaussian{
    public:
        double alpha;
        int angular_momentum;
        double center;
        vec center_vectors;
        mat shell_info; 
        Gaussian(double alpha, int angular_momentum, double center);
        Gaussian(double alpha, int angular_momentum, vec center_vectors);
}; 