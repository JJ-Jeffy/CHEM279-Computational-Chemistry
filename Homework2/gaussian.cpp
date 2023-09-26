// Problem set 2
// University of California, Berkeley

// Created by: Jeffy Jeffy
// Date: 9/25/23

// Description: This is a cpp file for the Gaussian class. It contains properties for various gaussian shells

#include <armadillo>
#include "gaussian.h"

using namespace arma;

Gaussian::Gaussian(double alpha, int angular_momentum, double center){
    this->alpha = alpha;
    this->angular_momentum = angular_momentum;
    this->center = center;
}

Gaussian::Gaussian(double alpha, int angular_momentum, vec center){
    this->alpha = alpha;
    this->angular_momentum = angular_momentum;
    this->center_vectors = center;

    // Construct a vec to store l, m, and n values 
    if (angular_momentum == 0){
        this-> shell_info = {0, 0, 0};
    }

    if (angular_momentum == 1){
        this-> shell_info = {{1, 0, 0},
                             {0, 1, 0},
                             {0, 0, 1}};
    }
}
