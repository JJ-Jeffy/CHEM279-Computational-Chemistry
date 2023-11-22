// Problem Set 3 
// University of California, Berkeley 

// Created by: Jeffy Jeffy 
// date: 10/16/23

// Description: This c++ file contains the implementation for the class molecule along with structs that define shells and atoms

#include "molecule_test.h"
#include <fstream> 

// Constructor for the molecule class 
Molecule::Molecule(string filename){
    ifstream infile(filename);
    if(!infile){
        cout << "File does not exist" << endl;
        exit(1);
    }

    // map for the elements 
    map<string, int> element = {{"H", 1}, {"C",6}, {"N",7}, {"O",8}, {"F",9}};
    map<string, int> element_valence = {{"H", 1}, {"C",4}, {"N",5}, {"O",6}, {"F",7}};


    infile >> numAtoms >> multiplicity;
    atomList = vector<atom>(numAtoms);
    for(int i = 0; i < numAtoms; i++){
        // initialize the atom vector coord 
        atomList[i].coords = zeros<vec>(3);
        infile >> atomList[i].atomName >> atomList[i].coords(0) >> atomList[i].coords(1) >> atomList[i].coords(2);
        atomList[i].atomicNumber = element[atomList[i].atomName];
        atomList[i].valence = element_valence[atomList[i].atomName];
    }
    

    infile.close();

    // count the number of electrons 
    numElectrons = countElectrons();

    // build basis functions for all the atoms 
    numBasisSets = 0; 
    for(int i = 0; i< numAtoms; i++){
        atomList[i].basisFunctionList = buildBasisFunctions(atomList[i]);
        // calculate the normalization constant
        calcNormalizationAll(atomList[i].basisFunctionList);
        numBasisSets += atomList[i].basisFunctionList.size();
        for(int j = 0; j < atomList[i].basisFunctionList.size(); j++){
            basisFunctionListAll.push_back(atomList[i].basisFunctionList[j]);
        }
    }

}

// Function: Count the number of electrons in the molecule 
int Molecule::countElectrons(){
    int count = 0; 
    for(int i = 0; i < numAtoms; i++){
        if (atomList[i].atomicNumber == 1){
            count += 1;
        }
        else if (atomList[i].atomicNumber == 6){
            count += 4;
        }
        else if (atomList[i].atomicNumber == 7){
            count += 5;
        }
        else if (atomList[i].atomicNumber == 8){
            count += 6;
        }
        else if (atomList[i].atomicNumber == 9){
            count += 7;
        }
    }

    // multiplicity 
    if (multiplicity == 0){
        p = count/2;
        q = count - p; 
    }
    if (multiplicity == 1){
        count -= 1;
        p = count/2;
        q = count - p; 
    }
    if (multiplicity == 3){
        p = 4; 
        q = 1; 
    }



    return count/2;
}

// Build basis functions 
vector<basisFunction> Molecule::buildBasisFunctions(atom &atom){
    vector<basisFunction> basisFunctionList; 

    int atomNum = atom.atomicNumber;
    if(atomNum == 1){
        struct basisFunction s = {"1s", atom.coords, s_lmn, H_exp, H_coeff, zeros<vec>(3)};
        basisFunctionList.push_back(s);
    } else if(atomNum == 6){
        struct basisFunction twos = {"2s", atom.coords, s_lmn, C_exp, C_s_coeff, zeros<vec>(3)};
        struct basisFunction px = {"2px", atom.coords, px_lmn, C_exp, C_p_coeff, zeros<vec>(3)};
        struct basisFunction py = {"2py", atom.coords, py_lmn, C_exp, C_p_coeff, zeros<vec>(3)};
        struct basisFunction pz = {"2pz", atom.coords, pz_lmn, C_exp, C_p_coeff, zeros<vec>(3)};
        basisFunctionList.push_back(twos);
        basisFunctionList.push_back(px);
        basisFunctionList.push_back(py);
        basisFunctionList.push_back(pz);
    } else if(atomNum == 7){
        struct basisFunction twos = {"2s", atom.coords, s_lmn, N_exp, C_s_coeff, zeros<vec>(3)};
        struct basisFunction px = {"2px", atom.coords, px_lmn, N_exp, C_p_coeff, zeros<vec>(3)};
        struct basisFunction py = {"2py", atom.coords, py_lmn, N_exp, C_p_coeff, zeros<vec>(3)};
        struct basisFunction pz = {"2pz", atom.coords, pz_lmn, N_exp, C_p_coeff, zeros<vec>(3)};
        basisFunctionList.push_back(twos);
        basisFunctionList.push_back(px);
        basisFunctionList.push_back(py);
        basisFunctionList.push_back(pz);
    } else if(atomNum == 8){
        struct basisFunction twos = {"2s", atom.coords, s_lmn, O_exp, C_s_coeff, zeros<vec>(3)};
        struct basisFunction px = {"2px", atom.coords, px_lmn, O_exp, C_p_coeff, zeros<vec>(3)};
        struct basisFunction py = {"2py", atom.coords, py_lmn, O_exp, C_p_coeff, zeros<vec>(3)};
        struct basisFunction pz = {"2pz", atom.coords, pz_lmn, O_exp, C_p_coeff, zeros<vec>(3)};
        basisFunctionList.push_back(twos);
        basisFunctionList.push_back(px);
        basisFunctionList.push_back(py);
        basisFunctionList.push_back(pz);
    } else if(atomNum == 9){
        struct basisFunction twos = {"2s", atom.coords, s_lmn, F_exp, C_s_coeff, zeros<vec>(3)};
        struct basisFunction px = {"2px", atom.coords, px_lmn, F_exp, C_p_coeff, zeros<vec>(3)};
        struct basisFunction py = {"2py", atom.coords, py_lmn, F_exp, C_p_coeff, zeros<vec>(3)};
        struct basisFunction pz = {"2pz", atom.coords, pz_lmn, F_exp, C_p_coeff, zeros<vec>(3)};
        basisFunctionList.push_back(twos);
        basisFunctionList.push_back(px);
        basisFunctionList.push_back(py);
        basisFunctionList.push_back(pz);
    }

    return basisFunctionList;
}

// Function that prints the geometry information
void Molecule::print_geom_info(){
    cout << "The molecule has " << numAtoms << " atoms" << endl;
    cout << "The molecule has " << numElectrons << " electrons" << endl;
    cout << "The molecule has " << p << " alpha electrons" << endl;
    cout << "The molecule has " << q << " beta electrons" << endl;
    for(int i = 0; i < numAtoms; i++){
        cout << "~~~~~~~~~~~" << endl;
        cout << "The atom is " << atomList[i].atomName << " and the coordinates are: " << endl; 
        cout << atomList[i].coords << endl;
        cout << "The atomic number is " << atomList[i].atomicNumber << endl;
        cout << "The number of basis sets for this atom is " << atomList[i].basisFunctionList.size() << endl;
        for (int j = 0; j < atomList[i].basisFunctionList.size(); j++){
            cout << "The basis function is " << atomList[i].basisFunctionList[j].name << endl;
            cout << "The quantum numbers are: " << endl;
            cout << atomList[i].basisFunctionList[j].quantumNumbers << endl;
            cout << "The exponents are: " << endl;
            cout << atomList[i].basisFunctionList[j].exponents << endl;
            cout << "The coefficients are: " << endl;
            cout << atomList[i].basisFunctionList[j].coefficients << endl;
            cout << "The normalization constants are: " << endl;
            cout << atomList[i].basisFunctionList[j].normalization << endl;
        }
        cout << "~~~~~~~~~~~" << endl;
    }
}

// Function: Calculate the Normalization Constants for a given basis function 
// Input: A basis function 
// Output: The normalization constant for that basis function
void Molecule::calculateNormalization(basisFunction &basisFunction){
    // Calculate the normalization constant for the basis function
    double norm1 = overlapIntegral(basisFunction.center, basisFunction.center, basisFunction.exponents(0), basisFunction.exponents(0), basisFunction.quantumNumbers, basisFunction.quantumNumbers);
    double norm2 = overlapIntegral(basisFunction.center, basisFunction.center, basisFunction.exponents(1), basisFunction.exponents(1), basisFunction.quantumNumbers, basisFunction.quantumNumbers);
    double norm3 = overlapIntegral(basisFunction.center, basisFunction.center, basisFunction.exponents(2), basisFunction.exponents(2), basisFunction.quantumNumbers, basisFunction.quantumNumbers);

    norm1 = pow(norm1, -0.5);
    norm2 = pow(norm2, -0.5);
    norm3 = pow(norm3, -0.5);

    basisFunction.normalization = {norm1, norm2, norm3};
}

// Function: Calculate the Normalization Constants for all of the basis functions
// Input: A vector of basis functions
// Output: The normalization constants for all of the basis functions
void Molecule::calcNormalizationAll(vector<basisFunction>& basisFunctionList){
    for(int i = 0; i < basisFunctionList.size(); i++){
        calculateNormalization(basisFunctionList[i]);
    }
}

//Function: Calculate the contracted overlap integral for two primitive Gaussian shells
//Input: Two basis functions 
//Output: The contracted overlap integral for two primitive Gaussian shells
double Molecule::contractedOverlapIntegral(basisFunction basisFunction1, basisFunction basisFunction2){
    double overlap = 0; 

    // loop over all the exponent combinations 
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            overlap += basisFunction1.coefficients(i)*basisFunction2.coefficients(j)*basisFunction1.normalization(i)*basisFunction2.normalization(j)*
                overlapIntegral(basisFunction1.center, basisFunction2.center, basisFunction1.exponents(i), basisFunction2.exponents(j), basisFunction1.quantumNumbers, basisFunction2.quantumNumbers);
        }
    }

    return overlap;
}

// Function: Calculate the overlap integral for two primitive Gaussian shells
mat Molecule::overlapMatrix(){
    mat overlap = arma::zeros<mat>(numBasisSets, numBasisSets);
    for(int i = 0; i < numBasisSets; i++){
        for(int j = 0; j < numBasisSets; j++){
            overlap(i,j) = contractedOverlapIntegral(basisFunctionListAll[i], basisFunctionListAll[j]);
        }
    }
    return overlap;
}

// Function: Calculate the gamma value for two given basis functions
// Input: Two basis functions
// Output: The gamma value for the two basis functions
double Molecule::calcGamma(basisFunction basisFunction1, basisFunction basisFunction2){
    // get the lmn values 
    vec lmn1 = basisFunction1.quantumNumbers;
    vec lmn2 = basisFunction2.quantumNumbers;

    // make sure only the s orbitals are being used
    if(lmn1(0) != 0 || lmn1(1) != 0 || lmn1(2) != 0){
        cout << "Error: Only s orbitals are allowed" << endl;
    }

    // extract other information from the basis sets
    vec da=basisFunction1.coefficients % basisFunction1.normalization;
    vec db=basisFunction2.coefficients % basisFunction2.normalization;
    vec alphaa = basisFunction1.exponents;
    vec alphab = basisFunction2.exponents;
    vec Ra = basisFunction1.center;
    vec Rb = basisFunction2.center;
    int len = basisFunction1.exponents.size();

    double sum = 0; 
    for (int k1=0; k1<len; k1++){
        for(int k2=0; k2<len; k2++){
            double sigmaA = 1.0/(alphaa(k1) + alphaa(k2));

            for(int j1=0; j1<len; j1++){
                for(int j2=0; j2<len; j2++){
                    double sigmaB = 1.0/(alphab(j1) + alphab(j2));

                    double I2e = I2e_pG(Ra, Rb, sigmaA, sigmaB); 

                    sum += da(k1)*da(k2)*db(j1)*db(j2)*I2e;
                }
            }
        }
    }

    return sum*27.2114;
}

// Function: create the gamma matrix 
mat Molecule::calcGammaMatrix(){
    mat gamma = arma::zeros<mat>(numAtoms, numAtoms);
    for(int i = 0; i < numAtoms; i++){
        for(int j = 0; j < numAtoms; j++){
            gamma(i,j) = calcGamma(atomList[i].basisFunctionList[0], atomList[j].basisFunctionList[0]);
        }
    }
    return gamma;
}

// Function: Calculate the hamiltonian core matrix 
mat Molecule::calcHCoreMatrix(){
    mat hcore = arma::zeros<mat>(numBasisSets, numBasisSets);
    int mu = 0;
    double para_val = 0; 

    for (int A = 0; A < numAtoms; A++){
        string name_A = atomList[A].atomName;
        int z_val_A = atomList[A].valence;
        int numBasis_A = atomList[A].basisFunctionList.size();
        double gammaAA = calcGamma(atomList[A].basisFunctionList[0], atomList[A].basisFunctionList[0]);
        // cout << "gammaAA: " << gammaAA << endl;

        // loop over all the basis sets in atomA
        for(int  i = 0; i < numBasis_A; i++){
            string name_ao = atomList[A].basisFunctionList[i].name;
            // find the related parameter 
            vec para = atom_para[name_A]; 
            if (name_ao == "1s" || name_ao == "2s"){
                para_val = para(0);
            } else if (name_ao == "2px" || name_ao == "2py" || name_ao == "2pz"){
                para_val = para(1);
            }
            hcore(mu, mu) =  - para_val - (z_val_A - 0.5)*gammaAA;

            int nu = 0;

            // loop over all the atoms in the molecule 
            for(int B = 0; B < numAtoms; B++){
                string name_B = atomList[B].atomName;
                int z_val_B = atomList[B].valence;
                int numBasis_B = atomList[B].basisFunctionList.size();
                double gammaAB = calcGamma(atomList[A].basisFunctionList[0], atomList[B].basisFunctionList[0]);

                // add the ZB and gammaAB terms
                if(A != B){
                    hcore(mu, mu) -= (z_val_B*gammaAB); 
                }

                // loop over all the basis sets in atomB 
                for (int j = 0; j < numBasis_B; j++){
                    if(mu != nu){
                        int beta_A = atom_para[name_A](2);
                        int beta_B = atom_para[name_B](2);
                        hcore(mu,nu) = - 0.5*(beta_A+beta_B)*contractedOverlapIntegral(atomList[A].basisFunctionList[i], atomList[B].basisFunctionList[j]);
                    }
                    nu += 1;
                }
            }
            mu += 1;
        }
    }
    return hcore;
}

//Function: Update the density matrix 
mat Molecule::updateDensityMatrix(mat coeff, string spin){
    mat matrix = zeros<mat>(numBasisSets, numBasisSets);
    // p(mu,vu) = sum over p electrons of C(p,mu)*C(p,vu)
    if(spin == "alpha"){
        for(int mu = 0; mu < numBasisSets; mu++){
            for(int vu = 0; vu < numBasisSets; vu++){
                for(int i = 0; i < p; i++){
                    matrix(mu, vu) += coeff(mu, i)*coeff(vu, i);
                }
            }
        }
    // p(mu,vu) = sum over q electrons of C(q,mu)*C(q,vu)
    } else if(spin == "beta"){
        for(int mu = 0; mu < numBasisSets; mu++){
            for(int vu = 0; vu < numBasisSets; vu++){
                for(int i = 0; i < q; i++){
                    matrix(mu, vu) += coeff(mu, i)*coeff(vu, i);
                }
            }
        }
    }
    return matrix;
}

// Function: Calculate the total density for a certain atom using the P_alpha and P_beta matrices 
vec Molecule::totalDensity(mat P_alpha, mat P_beta){
    vec totalDensity = zeros<vec>(numAtoms); 

    // add alpha and beta density together
    mat newMatrix = P_alpha + P_beta;
    int index = 0; 

    for (int i = 0; i < numAtoms; i++){
        int numBasis_A = atomList[i].basisFunctionList.size();
        for(int j = 0; j < numBasis_A; j++){
            totalDensity(i) += newMatrix(index, index);
            index += 1;
        }
    }
    
   return totalDensity;
}

// Function: Calculate the fock matric for alpha electrons
mat Molecule::calcFockMatrix(mat density_matrix, vec totalDensityVec){
    mat fock = zeros<mat>(numBasisSets, numBasisSets);
    int mu = 0; 
    double para_val = 0;

    for (int A = 0; A < numAtoms; A++){
        string name_A = atomList[A].atomName;
        int z_val_A = atomList[A].valence;
        int numBasis_A = atomList[A].basisFunctionList.size();
        double gammaAA = calcGamma(atomList[A].basisFunctionList[0], atomList[A].basisFunctionList[0]);
        double density = totalDensityVec(A);


        // loop over all the basis sets in atomA
        for(int  i = 0; i < numBasis_A; i++){
            string name_ao = atomList[A].basisFunctionList[i].name;
            // find the related parameter 
            vec para = atom_para[name_A]; 
            if (name_ao == "1s" || name_ao == "2s"){
                para_val = para(0);
            } else if (name_ao == "2px" || name_ao == "2py" || name_ao == "2pz"){
                para_val = para(1);
            }

            fock(mu, mu) = - para_val + ((density - z_val_A) - (density_matrix(mu, mu) - 0.5))*gammaAA;

            int nu = 0;

            // loop over all the atoms in the molecule 
            for(int B = 0; B < numAtoms; B++){
                string name_B = atomList[B].atomName;
                int z_val_B = atomList[B].valence;
                int numBasis_B = atomList[B].basisFunctionList.size();
                double gammaAB = calcGamma(atomList[A].basisFunctionList[0], atomList[B].basisFunctionList[0]);
                double density = totalDensityVec(B);

                // add the ZB and gammaAB terms
                if(A != B){
                    fock(mu, mu) += (density - z_val_B)*gammaAB;
                }

                // loop over all the basis sets in atomB 
                for (int j = 0; j < numBasis_B; j++){
                    if(mu != nu){
                        int beta_A = atom_para[name_A](2);
                        int beta_B = atom_para[name_B](2);
                        fock(mu,nu) = - 0.5*(beta_A+beta_B)*contractedOverlapIntegral(atomList[A].basisFunctionList[i], atomList[B].basisFunctionList[j])
                            - (density_matrix(mu, nu) * gammaAB);
                    }
                    nu += 1;
                }
            }
            mu += 1;
        }
    }
    return fock;
}

// Function: Calculate nuclear repulsion energy for the molecule
double Molecule::calcNuclearRepulsionEnergy(){
    double energy = 0; 
    for(int i = 0; i < numAtoms; i++){
        for(int j = i+1; j < numAtoms; j++){
            double distance = norm(atomList[i].coords - atomList[j].coords);
            energy += (atomList[i].valence*atomList[j].valence)/distance;
        }
    }
    return energy;
}

// Function: Calculate electronic repulsion 
double Molecule::calcElectronicEnergy(mat density_matrix, mat Fock_matrix){
    double energy = 0; 
    mat hcore = calcHCoreMatrix(); 
    for(int mu = 0; mu < numBasisSets; mu++){
        for(int nu = 0; nu < numBasisSets; nu++){
            energy += density_matrix(mu, nu)*(hcore(mu, nu) + Fock_matrix(mu, nu));
        }
    }
    return energy*0.5;
}

// Function: Calculate the classic SCF algorithm 
void Molecule::runSCF(double tolerance, int maxIterations){
    // step 1: Guess the density matrices 
    mat P_alpha = zeros<mat>(numBasisSets, numBasisSets);
    mat P_beta = zeros<mat>(numBasisSets, numBasisSets);
    vec totalDensityvec = zeros<vec>(numAtoms);

    bool converged = false;
    int iteration = 0;

    cout << "SCF Cycle Begins " << endl;
    cout << " -------------------------------" << endl;

    // step 6: Check if the density matrices have converged
    while(!converged && iteration < maxIterations){

        cout << "Iteration: " << iteration << endl;

        // step 2: Calculate the Fock matrix 
        mat F_alpha = calcFockMatrix(P_alpha, totalDensityvec);
        mat F_beta = calcFockMatrix(P_beta, totalDensityvec);

        cout << "F_alpha: " << endl;
        cout << F_alpha << endl;
        cout << "F_beta: " << endl;
        cout << F_beta << endl;

        // step 3: diagonalize the Fock matrix and obtain the eigenvalues and eigenvectors
        vec eigenvalues_alpha;
        mat eigenvectors_alpha;
        eig_sym(eigenvalues_alpha, eigenvectors_alpha, F_alpha);

        vec eigenvalues_beta;
        mat eigenvectors_beta;
        eig_sym(eigenvalues_beta, eigenvectors_beta, F_beta);

        cout << "Solving the eigenvalue problem : " << iteration << endl;
        

        // save the old density matrices
        mat P_alpha_old = P_alpha;
        mat P_beta_old = P_beta;

        cout << "Ea" << endl; 
        cout << eigenvalues_alpha << endl;
        cout << "Eb" << endl;
        cout << eigenvalues_beta << endl;

        cout << "Ca" << endl;
        cout << eigenvectors_alpha << endl;
        cout << "Cb" << endl;
        cout << eigenvectors_beta << endl;

        
        // step 4: Calculate the new density matrices
        P_alpha = updateDensityMatrix(eigenvectors_alpha, "alpha");
        cout << "P_alpha: " << endl;
        cout << P_alpha << endl;

        P_beta = updateDensityMatrix(eigenvectors_beta, "beta");
        cout << "P_beta: " << endl; 
        cout << P_beta << endl;

        totalDensityvec = totalDensity(P_alpha, P_beta);
        cout << "Total Density: " << endl;
        cout << totalDensityvec << endl;

        // step 5: check for convergence 
        if(abs(P_alpha - P_alpha_old).max() < tolerance && abs(P_beta - P_beta_old).max() < tolerance){
            converged = true;
            if (converged){
                cout << "SCF has converged" << endl;
                cout << " ---------------------------- " << endl;
                cout << "The nuclear repulsion energy is: " << endl;
                double nre = calcNuclearRepulsionEnergy()*27.211;
                cout << nre << " eV" << endl;
                double ee_alpha = calcElectronicEnergy(P_alpha, F_alpha);
                double ee_beta = calcElectronicEnergy(P_beta, F_beta);
                cout << "The electronic energy is: " << endl;
                cout << (ee_alpha + ee_beta) << " eV" << endl;
                cout << "The total energy is: " << endl;
                totalEnergy = (nre + ee_alpha + ee_beta);
                cout << totalEnergy << " eV" << endl;
            } else {
                cout << "SCF has not converged" << endl;
            }
        }
        iteration += 1;
    }

}

    


