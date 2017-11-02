/* Object oriented Hamiltonian of the C-C structure.
   Copyright (C) 2015-2017 Free Software Foundation, Inc.
   Contributed by Eduardo Ramos Fernandez <eduradical951@gmail.com>
This file is part of TBMD.
TBMD is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 3, or (at your option) any later
version.
TBMD is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.
You should have received a copy of the GNU General Public License
along with TBMD; see the file COPYING.  If not see
<http://www.gnu.org/licenses/> */

#include "Hamiltonian.h"

// Constructors/Destructors

Hamiltonian::Hamiltonian() {
}

Hamiltonian::Hamiltonian(uint N) {
    H_ = Matrix(4 * N, 4 * N);
    n_ = 4 * N;
    Eigenvectors_ = Matrix(4 * N, 4 * N);
    eigenvalues_ = Vector(4 * N);
}

Hamiltonian::~Hamiltonian() {
}

void Hamiltonian::print() {
    std::cout << std::endl;
    std::cout.precision(16);
    H_.raw_print();
    std::cout << std::endl;
}


// Methods

void Hamiltonian::build2ParticleHamiltonian(const Particle& p1, const Particle& p2, double distance, Vector3D cosines) {
    uint N = n_ / 4;
    double bl_row = p1.id(), bl_col = p2.id();

    double dx = cosines.x, dy = cosines.y, dz = cosines.z;
    // Off-diagonal constants                                                                                                                                                                                                                   
    double v_sssi = -5.0, v_spsi = 4.7, v_ppsi = 5.5, v_pppi = -1.55;
    double scaling = scale_func(distance);
    v_sssi *= scaling;
    v_spsi *= scaling;
    v_ppsi *= scaling;
    v_pppi *= scaling;
    // Diagonal constants
    double E_s = -2.99, E_p = 3.71; // Diagonal and off-diagonal constants  

    // Calculate one submatrix (interaction of A with B)                                                                                                                                                                                        
    int block_col = bl_row, block_row = bl_col;
    H_(4 * block_col, 4 * block_row) = v_sssi;
    H_(4 * block_col + 1, 4 * block_row + 1) = dx * dx * (v_ppsi - v_pppi) + v_pppi;
    H_(4 * block_col + 2, 4 * block_row + 2) = dy * dy * (v_ppsi - v_pppi) + v_pppi;
    H_(4 * block_col + 3, 4 * block_row + 3) = dz * dz * (v_ppsi - v_pppi) + v_pppi;
    H_(4 * block_col + 1, 4 * block_row + 2) = dx * dy * (v_ppsi - v_pppi);
    H_(4 * block_col + 2, 4 * block_row + 1) = H_(4 * block_col + 1, 4 * block_row + 2);
    H_(4 * block_col + 1, 4 * block_row + 3) = dx * dz * (v_ppsi - v_pppi);
    H_(4 * block_col + 3, 4 * block_row + 1) = H_(4 * block_col + 1, 4 * block_row + 3);
    H_(4 * block_col + 2, 4 * block_row + 3) = dy * dz * (v_ppsi - v_pppi);
    H_(4 * block_col + 3, 4 * block_row + 2) = H_(4 * block_col + 2, 4 * block_row + 3);
    H_(4 * block_col + 1, 4 * block_row) = -dx*v_spsi;
    H_(4 * block_col, 4 * block_row + 1) = -H_(4 * block_col + 1, 4 * block_row);
    H_(4 * block_col + 2, 4 * block_row) = -dy*v_spsi;
    H_(4 * block_col, 4 * block_row + 2) = -H_(4 * block_col + 2, 4 * block_row);
    H_(4 * block_col + 3, 4 * block_row) = -dz*v_spsi;
    H_(4 * block_col, 4 * block_row + 3) = -H_(4 * block_col + 3, 4 * block_row);

    // Calculate the second submatrix (interaction of B with A) to make Hamiltonian Hermitian                                                                                                                                                   

    block_col = bl_col, block_row = bl_row;
    H_(4 * block_col, 4 * block_row) = H_(4 * block_row, 4 * block_col);
    H_(4 * block_col + 1, 4 * block_row + 1) = H_(4 * block_row + 1, 4 * block_col + 1);
    H_(4 * block_col + 2, 4 * block_row + 2) = H_(4 * block_row + 2, 4 * block_col + 2);
    H_(4 * block_col + 3, 4 * block_row + 3) = H_(4 * block_row + 3, 4 * block_col + 3);
    H_(4 * block_col + 1, 4 * block_row + 2) = H_(4 * block_row + 1, 4 * block_col + 2);
    H_(4 * block_col + 2, 4 * block_row + 1) = H_(4 * block_col + 1, 4 * block_row + 2);
    H_(4 * block_col + 1, 4 * block_row + 3) = H_(4 * block_row + 1, 4 * block_col + 3);
    H_(4 * block_col + 3, 4 * block_row + 1) = H_(4 * block_col + 1, 4 * block_row + 3);
    H_(4 * block_col + 2, 4 * block_row + 3) = H_(4 * block_row + 2, 4 * block_col + 3);
    H_(4 * block_col + 3, 4 * block_row + 2) = H_(4 * block_col + 2, 4 * block_row + 3);
    H_(4 * block_col + 1, 4 * block_row) = -H_(4 * block_row + 1, 4 * block_col);
    H_(4 * block_col, 4 * block_row + 1) = -H_(4 * block_col + 1, 4 * block_row);
    H_(4 * block_col + 2, 4 * block_row) = -H_(4 * block_row + 2, 4 * block_col);
    H_(4 * block_col, 4 * block_row + 2) = -H_(4 * block_col + 2, 4 * block_row);
    H_(4 * block_col + 3, 4 * block_row) = -H_(4 * block_row + 3, 4 * block_col);
    H_(4 * block_col, 4 * block_row + 3) = -H_(4 * block_col + 3, 4 * block_row);



    for (uint i = 0; i < N; i++) {
        H_(4 * i, 4 * i) = E_s;
    }

    for (uint i = 0; i < N; i++) {
        H_(4 * i + 1, 4 * i + 1) = E_p;
        H_(4 * i + 2, 4 * i + 2) = E_p;
        H_(4 * i + 3, 4 * i + 3) = E_p;
    }

}

double Hamiltonian::scale_func(double d) {
    const double r_1 = 2.45, r_m = 2.6;

    const double A_s = 2.89933161442; // A=[(r_0)^n]*e^[n*(r_0/r_c)^n_c]	


    const double B_s = -0.0126200997391; // B=-n/(r_c^n_c)


    const double n = 2.0;
    const double n_c = 6.5;
    const double cs_0 = 6.7392620074314E-3;
    const double cs_1 = -8.1885359517898E-2;
    const double cs_2 = 0.1932365259144;
    const double cs_3 = 0.354287433238;
    double scaling;

    if (d > r_m)
        scaling = 0.0;

    else if (d <= r_1)
        scaling = A_s * exp(B_s * pow(d, n_c)) / pow(d, n);

    else
        scaling = cs_3 * (d - r_1)*(d - r_1)*(d - r_1) + cs_2 * (d - r_1)*(d - r_1) + cs_1 * (d - r_1) + cs_0;

    return scaling;
}

void Hamiltonian::build(const ParticleList& particle_list) {
    uint k;
    NearestNeighboursList::const_iterator nn;
    ParticleList::const_iterator p;
    const NearestNeighboursList* nn_list;
    Particle* p_next;
    for (p = particle_list.begin(); p != particle_list.end(); p++) {
        nn_list = &p->nearestNeighboursList();
        nn = nn_list->begin();
        k = 0;
        while ((nn != nn_list->end()) && (k < p->noNearestNeighbours())) {
            p_next = nn->particle_ptr;
            if (p_next->id() > p->id()) {
                build2ParticleHamiltonian(*p, *p_next, nn->distance, nn->cosines);
            }
            k++;
            nn++;
        }
    }
}

void Hamiltonian::calculateCharge(ParticleList& particle_list){
    ParticleList::iterator p;
    double charge=0.0;
    int j,k;
    int N=n_/4;
    for (p = particle_list.begin(); p != particle_list.end(); p++) {
        charge=0.0;
        for(j=0;j<4;j++){
         for (k=0;k<2*N;k++){
          //Eigenvectors_.print();
          charge+=2*Eigenvectors_(4*p->id()+j,k)*Eigenvectors_(4*p->id()+j,k);
         }  
        } p->charge() = charge;
    }
}


void Hamiltonian::computeEigs() {
    arma::eig_sym(eigenvalues_, Eigenvectors_, H_);
}
