/* Repulsive energy of the Hamiltonian.
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

#include "RepulsiveEnergy.h"

// Constructors/Destructors

RepulsiveEnergy::RepulsiveEnergy() {
}

RepulsiveEnergy::~RepulsiveEnergy() {
}

// Methods

// Cosntants for the functions below for Carbon
const double p_0 = -2.5909765118191;
const double p_1 = 0.5721151498619;
const double p_2 = -1.7896349903996E-3;
const double p_3 = 2.3539221516757E-5;
const double p_4 = -1.24251169551587E-7;
const double d_m = 2.6;
const double d_1 = 2.57;
const double Phi_0 = 8.18555;
const double m = 3.30304;
const double m_c = 8.6655;
const double d_c = 2.1052;
const double d_0 = 1.64;
const double c_0 = 2.2504290109E-8;
const double c_1 = -1.4408640561E-6;
const double c_2 = 2.1043303374E-5;
const double c_3 = 6.6024390226E-5;

const double A1 = Phi_0 * pow(d_0, m) * exp(m * pow((d_0/d_c), m_c));
const double A2 = -m * pow((1/d_c), m_c);

double phi(double r) {
    double phi;
    if (r > d_m)
        phi = 0;
    else if (r <= d_1)
     	phi = A1 * pow(r, -m) * exp(A2 * pow(r, m_c));
    else
        phi = c_0 + (c_1 * (r - d_1)) + (c_2 * (r - d_1)*(r - d_1)) + (c_3 * (r - d_1)*(r - d_1)*(r - d_1));
    return phi;
}

double dPhi(double r){
	double value = 0.;
	if (r <= d_1) {
	//	value = - phi(r) * ( ((1/r)*m*m_c*pow((r/d_c),m_c)) + (m/r) );
		value =  (phi(r)/r) * ( (A2*pow(r,m_c)*m_c) - m);
 	} else if (d_1 < r && r <= d_m) {
		value = c_1 + (2.*c_2*(r-d_1)) + (3.*c_3*(r-d_1)*(r-d_1));
	} else { 
		value = 0;
	}
 	return value;
}

double poly_func(double x) {
    return p_0 + (p_1 * x) + (p_2 * x * x) + (p_3 * x * x * x) + (p_4 * x * x * x * x);
}

double dfSumPhi(Particle& p) {
    double Ai = A(p);
    return (p_1 + (2.0 * p_2 * Ai) + (3.0 * p_3 * Ai * Ai) + (4.0 * p_4 * Ai * Ai * Ai));
}

double fSumPhi(Particle& p) {
    double Ai = A(p);
    return (p_0 + (p_1 * Ai) + (p_2 * Ai * Ai) + (p_3 * Ai * Ai * Ai) + (p_4 * Ai * Ai * Ai * Ai));
}

double A(Particle& p) {
    NearestNeighboursList::iterator nn;
    NearestNeighboursList& nn_list = p.nearestNeighboursList();
    nn = nn_list.begin();
    uint k = 0;
    double sum = 0.0;
    while (nn != nn_list.end() && (k < p.noNearestNeighbours())) {
        sum += phi(nn->distance);
        k++;
        nn++;
    }
    return sum;
}

double RepulsiveEnergy::computeEnergy(ParticleList& particle_list) {
    ParticleList::iterator p;
    tot_energy_ = 0;
    for (p = particle_list.begin(); p < particle_list.end(); p++)
        tot_energy_ += fSumPhi(*p);
}

void RepulsiveEnergy::updateParticleNetForce(Particle& p) {
    double dfdAj, dPhidxk, dPhidyk, dPhidzk, factor, dfdAk;
    Force f;
    NearestNeighboursList::iterator nn;
    NearestNeighboursList& nn_list = p.nearestNeighboursList();
    nn = nn_list.begin();
    uint k = 0;
    f.x = 0;
    f.y = 0;
    f.z = 0;
    dfdAk = dfSumPhi(p);
    while (nn != nn_list.end() && (k < p.noNearestNeighbours())) {
        dfdAj = dfSumPhi(*(nn->particle_ptr));
        factor = (dfdAk + dfdAj);
        dPhidxk = nn->cosines.x * dPhi(nn->distance);

        Particle* p_next;
        p_next = nn->particle_ptr;
        f.x += factor * dPhidxk;

        dPhidyk = nn->cosines.y * dPhi(nn->distance);
        f.y += factor * dPhidyk;
        dPhidzk = nn->cosines.z * dPhi(nn->distance);
        f.z += factor * dPhidzk;
        k++;
        nn++;
    }
    if (fabs(f.x) < TOL) {
        f.x = 0.0;
    }
    if (fabs(f.y) < TOL)
        f.y = 0.0;
    if (fabs(f.z) < TOL)
        f.z = 0.0;

    p.force() += f;
}
