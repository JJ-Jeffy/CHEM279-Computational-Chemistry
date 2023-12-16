// Problem set 1
// University of California, Berkeley

// Creator: Jeffy Jeffy
// Date: 9/11/23

#include "utils.h"

// Function to print the geometry of the molecule
void print_geom(Molecule mol){
    cout << "The number of atoms is " << mol.numAtoms << endl;
    cout << "The coordinates of this molecule are: " << endl;
    for (int i=0; i < mol.numAtoms; i++){
        cout << mol.atoms[i] << " " << mol.coords(i,0) << " " << mol.coords(i,1) << " " << mol.coords(i,2) << endl;
        cout << "The atomic number of this atom is: " << mol.zvals(i) << endl;
        cout << "The mass of this atom is: " << mol.masses(i) << endl;
    }
}

// Function to calculate the distance between two vectors
double calculate_distance(const vec& p1, const vec& p2) {
    return norm(p1 - p2, 2);
}

// Function to calculate the lennard jones potential 
double calculate_LJ(double r_ij, double sigma, double epsilon){
    double r12 = pow(sigma/r_ij, 12);
    double r6 = pow(sigma/r_ij, 6);
    return epsilon * (r12 - (2*r6));
};

// Function to calculate the energy given a molecule object
double calculate_energy(Molecule mol){
    double energy = 0.0; 
    for (int i=0; i < mol.numAtoms; i++){
        for (int j=i+1; j < mol.numAtoms; j++){
            double rij = norm(mol.coords.row(i) - mol.coords.row(j), 2);
            double sigma = mol.sigma(i,j);
            double epsilon = mol.epsilon(i,j);
            energy += calculate_LJ(rij, sigma, epsilon);
        }
    }
    return energy;
}

// Function to calculate the energy given a matrix of coordinates
double calculate_energy(mat coords, mat sigma, mat epsilon){
    double energy = 0.0; 
    int numAtoms = coords.n_rows; 
    for (int i=0; i < numAtoms; i++){
        for (int j=i+1; j < numAtoms; j++){
            double rij = norm(coords.row(i) - coords.row(j), 2);
            double sigma_val = sigma(i,j);
            double epsilon_val = epsilon(i,j);
            energy += calculate_LJ(rij, sigma_val, epsilon_val);
        }
    }
    return energy;
}

// Function to calculate the forces on each atom
double calculate_forces(double rij, double x1, double x2, double sigma, double epsilon){
    double r12_term = -12*pow(sigma,12)*(pow(1/rij, 13));
    double r6_term = 12*pow(sigma,6)*(pow(1/rij, 7));
    double dist = (x1 - x2)/rij;
    return -1*(epsilon * (r12_term + r6_term) * dist);  
}

// Function to calculate the forces on the whole molecule 
mat calculate_total_force(Molecule mol){
    mat forces = zeros<mat>(mol.numAtoms, 3);
    cout << "we are calculating the forces" << endl;
    for (int i=0; i < mol.numAtoms; i++){
        for (int j=i+1; j < mol.numAtoms; j++){
            double rij = norm(mol.coords.row(i) - mol.coords.row(j), 2);
            double x1 = mol.coords(i,0);
            double x2 = mol.coords(j,0);
            double y1 = mol.coords(i,1);
            double y2 = mol.coords(j,1);
            double z1 = mol.coords(i,2);
            double z2 = mol.coords(j,2);
            forces(i,0) += calculate_forces(rij, x1, x2, mol.sigma(i,j), mol.epsilon(i,j));
            forces(j,0) += calculate_forces(rij, x2, x1, mol.sigma(i,j), mol.epsilon(i,j));
            forces(i,1) += calculate_forces(rij, y1, y2, mol.sigma(i,j), mol.epsilon(i,j));
            forces(j,1) += calculate_forces(rij, y2, y1, mol.sigma(i,j), mol.epsilon(i,j));
            forces(i,2) += calculate_forces(rij, z1, z2, mol.sigma(i,j), mol.epsilon(i,j));
            forces(j,2) += calculate_forces(rij, z2, z1, mol.sigma(i,j), mol.epsilon(i,j));
        }
    }
    return forces;
}

// Function to calculate the forces given a matrix of coordinates
mat calculate_total_force(mat coords, mat sigma, mat epsilon){
    int numAtoms = coords.n_rows;
    mat forces = zeros<mat>(numAtoms, 3);
    for (int i=0; i < numAtoms; i++){
        for (int j=i+1; j < numAtoms; j++){
            double rij = norm(coords.row(i) - coords.row(j), 2);
            double x1 = coords(i,0);
            double x2 = coords(j,0);
            double y1 = coords(i,1);
            double y2 = coords(j,1);
            double z1 = coords(i,2);
            double z2 = coords(j,2);
            double sigma_val = sigma(i,j);
            double epsilon_val = epsilon(i,j);
            forces(i,0) += calculate_forces(rij, x1, x2, sigma_val, epsilon_val);
            forces(j,0) += calculate_forces(rij, x2, x1, sigma_val, epsilon_val);
            forces(i,1) += calculate_forces(rij, y1, y2, sigma_val, epsilon_val);
            forces(j,1) += calculate_forces(rij, y2, y1, sigma_val, epsilon_val);
            forces(i,2) += calculate_forces(rij, z1, z2, sigma_val, epsilon_val);
            forces(j,2) += calculate_forces(rij, z2, z1, sigma_val, epsilon_val);
        }
    }
    return forces;
}

// Function to calculate the forward difference 
mat calculate_forward_difference(Molecule mol, double h){
    mat forces = zeros<mat>(mol.numAtoms, 3);
    mat coords = mol.coords;
    double energy = calculate_energy(mol);
    for (int i=0; i < mol.numAtoms; i++){
        for (int j=0; j < 3; j++){
            coords(i,j) += h;
            double new_energy = calculate_energy(coords, mol.sigma, mol.epsilon);
            forces(i,j) = -(new_energy - energy)/h;
            coords(i,j) -= h;
        }
    }
    return forces;
} 

// Function to calculate the central difference
mat calculate_central_difference(Molecule mol, double h){
    mat forces = zeros<mat>(mol.numAtoms, 3);
    mat coords = mol.coords;
    for (int i = 0; i < mol.numAtoms; i++){
        for(int j = 0; j < 3; j++){
            coords(i,j) += h;
            double new_energy = calculate_energy(coords, mol.sigma, mol.epsilon); 
            coords(i,j) -= 2*h;
            double new_energy2 = calculate_energy(coords, mol.sigma, mol.epsilon);
            coords(i,j) += h;
            forces(i,j) = -(new_energy - new_energy2)/(2*h);
        }
    }
    return forces;
}

// Function to optimize the geometry of the molecule
void optimize_geometry(Molecule& mol, double tol, int maxiter, double h){
    double search_stepsize = 0.3; 
    double golden_ratio = 0.38197; 

    // Calculate the initial coords 
    mat A = mol.coords;
    mat B = zeros<mat>(mol.numAtoms, 3);
    mat C = zeros<mat>(mol.numAtoms, 3);
    mat D = zeros<mat>(mol.numAtoms, 3);
    mat opt_geom = zeros<mat>(mol.numAtoms, 3);

    // Calculate distances between points 
    double AB, BD, AD; 
    double golden_tol = 1e-7; 

    // Calculate the forces 
    mat forces = calculate_central_difference(mol, h); 
    mat unit_Force = zeros<mat>(mol.numAtoms, 3);
    cout << "The starting forces are: " << endl; 
    cout << forces << endl;
    cout << "The norm of the force is: " << norm(forces, 2) << endl;

    // Calculate the energy
    double energy = calculate_energy(mol);
    cout << "The starting energy is: " << energy << endl;

    // Calculate the new coordinates
    int count = 0;
    cout << endl << "Starting the optimization" << endl << endl;

    while (norm(forces, "fro") > tol && count < 5e2){
        unit_Force = forces/norm(forces, 2);
        
        // update B 
        int bracket_count = 0; 
        B = A + search_stepsize * unit_Force;
        while (calculate_energy(B, mol.sigma, mol.epsilon) > calculate_energy(A, mol.sigma, mol.epsilon)){
            bracket_count += 1; 
            search_stepsize /= 2; 
            B = A + search_stepsize * unit_Force;
            cout << "Bracketing for line search, finding B point with stepsize " << search_stepsize << endl;
            if (bracket_count > 100){
                cout << "Finding point B failed, please adjust the search step size" << endl; 
                goto step_forward; 
            }
        }

        cout << "Bracketing for line search, B point found" << endl;
        
        // update D
        bracket_count = 0;
        D = B + search_stepsize * unit_Force;
        while(calculate_energy(D, mol.sigma, mol.epsilon) < calculate_energy(B, mol.sigma, mol.epsilon)){
            bracket_count += 1;
            search_stepsize *= 1.2; 
            D = B + search_stepsize * unit_Force;
            cout << "Bracketing for line search, finding D point with stepsize " << search_stepsize << endl;
            if (bracket_count > 100){
                cout << "Finding point D failed, please adjust the search step size" << endl;
                goto step_forward;
            }
        }

        cout << "Bracketing for line search, D point found" << endl;

        // If we fail to find B or Dm it is possible there's no local minima along the 
        // negative gradient direction
        step_forward:
        if (bracket_count > 100){
            cout << "Bracketing failed, use normal steepest descent instead" << endl; 
            search_stepsize = 0.3; 
            mat new_point = A + search_stepsize * unit_Force;
            while(calculate_energy(new_point, mol.sigma, mol.epsilon) > calculate_energy(A, mol.sigma, mol.epsilon)){
                search_stepsize /= 2; 
            }
            A = new_point;
            cout << "Current energy is: " << calculate_energy(A, mol.sigma, mol.epsilon) << endl;
            cout << "new point" << endl; 
            opt_geom = A;
            goto update_force; 
        }

        // Golden search 
        AB = norm(B-A, "fro");
        BD = norm(D-B, "fro");
        AD = norm(D-A, "fro");
        if(AB < BD){
            C = D + golden_ratio*(A-D);
        } else {
            C = B; 
            B = A + golden_ratio*(D-A);
        }

        while(AD > golden_tol){
            if (calculate_energy(B, mol.sigma, mol.epsilon) > calculate_energy(C, mol.sigma, mol.epsilon)){
                A = B; 
                B = C; 
            } else {
                D = C; 
            }

            AB = norm(B-A, "fro");
            BD = norm(D-B, "fro");
            AD = norm(D-A, "fro");
            if (AB < BD){
                C = D + golden_ratio*(A-D);
            } else {
                C = B; 
                B = A + golden_ratio*(D-A);
            }
        }
        if (calculate_energy(B, mol.sigma, mol.epsilon) > calculate_energy (C, mol.sigma, mol.epsilon)){
            cout << "The new point is: " << endl;
            cout << C << endl;
            cout << "The current energy is: " << calculate_energy(C, mol.sigma, mol.epsilon) << endl;
            opt_geom = C;
        } else {
            cout << "The new point is: " << endl;
            cout << B << endl;
            cout << "The current energy is: " << calculate_energy(B, mol.sigma, mol.epsilon) << endl;
            opt_geom = B;
        }

        // Update the forces in the new position 
        update_force:
        forces = calculate_total_force(opt_geom, mol.sigma, mol.epsilon);
        cout << "The force is" << endl;
        cout << forces << endl;
        cout << "The norm of the force is: " << norm(forces, "fro") << endl;
        count += 1;
    }

    // set opt_geom to the final geometry
    mol.coords = opt_geom;
    cout << "Total iterations: " << count << endl; 
    cout << "The final energy is: " << calculate_energy(opt_geom, mol.sigma, mol.epsilon) << endl;

}

// Function to calculate the hessian using finite difference of the gradient 
mat calculate_hessian(Molecule mol, double h){
    mat hessian = zeros<mat>(3*mol.numAtoms, 3*mol.numAtoms);
    int numCoordinates = 3*mol.numAtoms;

    //finite difference of the gradient
    mat gradient = calculate_total_force(mol);

    // calculate the hessian matrix
    for (int i = 0; i < numCoordinates; i++){
        for(int j= i; j < numCoordinates; j++){
            mat save_coords_plus = mol.coords;
            mat save_coords_minus = mol.coords;

            // perturb the coordinates for finite differences 
            save_coords_plus(i/3, i%3) += h;
            save_coords_plus(j/3, j%3) += h;

            save_coords_minus(i/3, i%3) -= h;
            save_coords_minus(j/3, j%3) -= h;

            // central finite difference formula
            mat gradient_plus = calculate_total_force(save_coords_plus, mol.sigma, mol.epsilon);
            mat gradient_minus = calculate_total_force(save_coords_minus, mol.sigma, mol.epsilon);

            double second_derivative = - (gradient_plus(i/3, i%3) - gradient_minus(i/3, i%3))/(2*h); 

            // set the hessian element
            hessian(i,j) = second_derivative;

            // set the hessian element
            hessian(j,i) = second_derivative;
        }
    }
    return hessian;
}



// Function to calculate the hessian element
double calculate_hessian_element(Molecule& mol, double h, int i, int j){
    mat coords = mol.coords;
    double energy = calculate_energy(mol);
    double original_value = coords(i,j);

    // Perturb the coordinate and calculate the energies 
    coords(i,j) += h;
    double energy_plus = calculate_energy(coords, mol.sigma, mol.epsilon);

    coords(i,j) -= 2*h;
    double energy_minus = calculate_energy(coords, mol.sigma, mol.epsilon);

    coords(i,j) = original_value;

    return (energy_plus - 2*energy + energy_minus)/(h*h);
}

// Function to calculate the hessian mixed element 
double calculate_mixed_hessian_element(Molecule& mol, double h, int i, int j, int k, int l){
    mat coords = mol.coords;
    
    // Perturb the first coordinate and calculate the energies
    coords(i,j) += h; 
    double energy_plus1 = calculate_energy(coords, mol.sigma, mol.epsilon);

    coords(i,j) -= 2*h;
    double energy_minus1 = calculate_energy(coords, mol.sigma, mol.epsilon);

    // Perturb the second coordinate and calculate the energies
    coords(k,l) += h;
    double energy_plus2 = calculate_energy(coords, mol.sigma, mol.epsilon);

    coords(k,l) -= 2*h;
    double energy_minus2 = calculate_energy(coords, mol.sigma, mol.epsilon);

    // calculate the mixed partial derivate using finite difference 
    double hessian_element = (energy_plus1 - energy_minus1 - energy_plus2 + energy_minus2)/(4*h*h);

    // reset the coordinates
    coords(i,j) += h;
    coords(k,l) += h;

    return hessian_element;
}

// Function to calculate the hessian matrix
// mat calculate_hessian(Molecule mol, double h){
//     mat hessian = zeros<mat>(3*mol.numAtoms, 3*mol.numAtoms);

//     // Hij = E(+)  - 2 E0 + E(-)/h2 
//     mat coords = mol.coords;

//     // Calculate the energy at the initial point
//     double energy = calculate_energy(coords, mol.sigma, mol.epsilon);

//     for(int atomIndex = 0; atomIndex < mol.numAtoms; atomIndex++){
//         for(int coordIndex = 0; coordIndex < 3; coordIndex++){
//             hessian(atomIndex * 3 + coordIndex, atomIndex * 3 + coordIndex) = calculate_hessian_element(mol, h, atomIndex, coordIndex);

//             // Calculate the off diagonal elements
//             for(int atomIndex2 = atomIndex + 1; atomIndex2 < mol.numAtoms; atomIndex2++){
//                 for(int coordIndex2 = 0; coordIndex2 < 3; coordIndex2++){
//                     hessian(atomIndex * 3 + coordIndex, atomIndex2 * 3 + coordIndex2) = calculate_mixed_hessian_element(mol, h, atomIndex, coordIndex, atomIndex2, coordIndex2);
//                     hessian(atomIndex2 * 3 + coordIndex2, atomIndex * 3 + coordIndex) = hessian(atomIndex * 3 + coordIndex, atomIndex2 * 3 + coordIndex2);
//                 }
//             }
//         }
//     }
//     return hessian;
// }


// mat calculate_hessian(Molecule mol, double h){
//     mat hessian = zeros<mat>(3*mol.numAtoms, 3*mol.numAtoms);

//     int numAtoms = mol.numAtoms;
//     int coordinates = 3*numAtoms;

//     mat coords = mol.coords;

//     // Calculate the energy at the initial point
//     double energy = calculate_energy(coords, mol.sigma, mol.epsilon);

//     for (int i = 0; i < coordinates; i++){
//         for(int j= i; j < coordinates; j++){
//             mat save_coords_plus = coords;
//             mat save_coords_minus = coords;

//             // perturb the coordinates for finite differences 
//             save_coords_plus(i/3, i%3) += h;
//             save_coords_plus(j/3, j%3) += h;

//             save_coords_minus(i/3, i%3) -= h;
//             save_coords_minus(j/3, j%3) -= h;

//             // central finite difference formula
//             double energy_plus = calculate_energy(save_coords_plus, mol.sigma, mol.epsilon);
//             double energy_minus = calculate_energy(save_coords_minus, mol.sigma, mol.epsilon);

//             double second_derivative = (energy_plus - 2*energy + energy_minus)/(h*h);

//             // set the hessian element
//             hessian(i,j) = second_derivative;
//             hessian(j,i) = second_derivative;
//         }
//     }

//     return hessian;
// }

// Function to generate the mass weighted hessian matrix
mat generate_mass_weighted_hessian(mat hessian, Molecule& mol){

    mat mass_matrix = zeros<mat>(3*mol.numAtoms, 3*mol.numAtoms);

    // generate the mass matrix
    for (int i = 0; i < mol.numAtoms; i++){
        for (int j = 0; j < 3; j++){
            mass_matrix(i*3 + j, i*3 + j) = mol.masses(i);
        }
    }

    cout << mass_matrix << endl;

    // generate the inverse of the mass matrix
    mat inverse_mass_matrix = inv(mass_matrix);

    // generate the square root of the inverse of the mass matrix 
    mat sqrt_inverse_mass_matrix = sqrt(inverse_mass_matrix);

    // transorm the hessian matrix
    mat mass_weighted_hessian = sqrt_inverse_mass_matrix * hessian * sqrt_inverse_mass_matrix;

    return mass_weighted_hessian;
}

// Function to calculate the vibrational frequencies
vec calculate_vibrational_frequencies(mat mass_weighted_hessian){
    vec eigenvalues = eig_sym(mass_weighted_hessian);
    vec frequencies = sqrt(eigenvalues);
    return frequencies;
}