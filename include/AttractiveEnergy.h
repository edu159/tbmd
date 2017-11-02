/* Attractive energy of the Hamiltonian.
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

#ifndef ATTRACTIVEENERGY_H
#define ATTRACTIVEENERGY_H

#include "Energy.h"
#include "Hamiltonian.h"
#include "Particle.h"
#include "common.h"

/**
 * class AtractiveEnergy
 * 
 */

class AttractiveEnergy : public Energy {
private:

    // Private attributes
    Hamiltonian hamiltonian_;
    Matrix rho_matrix_;

public:

    // Getters

    const Hamiltonian& hamiltonian() {
        return hamiltonian_;
    }

    // Setters

    void hamiltonian(const Hamiltonian& hamiltonian) {
        hamiltonian_ = hamiltonian;
    }

public:

    // Constructors/Destructors
    AttractiveEnergy();
    AttractiveEnergy(uint N);
    virtual ~AttractiveEnergy();


    // Methods
    void permuteCosines(double &dx, double&dy, double&dz, Vector3D cosines, int per);
    double calcOffD(int dir, double s, double ds_dr, double distance, Vector3D cosines, int per);
    double calcDiag(int dir, double s, double ds_dr, double distance, Vector3D cosines, int per);
    double calcTopR(int dir, double s, double ds_dr, double distance, Vector3D cosines, int per);
    double scale_func_copy(double d);
    double scale_func_derivative(double d, double scaling);    

    void calcDeltaH(arma::mat& DeltaH, int dir, double distance, Vector3D cosines);
    void calc1D_Ham(int ncoord, const Particle& p1, const Particle& p2, Hamiltonian& H_delta, double delta_x);
    Force compute2ParticleForce(const Particle& p1, const Particle& p2, double distance, Vector3D cosines);
    virtual double computeEnergy(ParticleList& particle_list);
    virtual void updateParticleNetForce(Particle& p);
    void update(ParticleList& particle_list);
    void GetDensityMatrix(const Matrix& EV);

};

void GetdH(int Ne, Particle p, int c, Matrix& dH);
double GetTraceMod(int Ne, Matrix rho, Matrix dH);

#endif // ATTRACTIVEENERGY_H
