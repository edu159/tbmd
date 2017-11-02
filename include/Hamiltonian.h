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


#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "common.h"
#include "Particle.h"

/**
 * class Hamiltonian
 * 
 */

class Hamiltonian {
private:

    // Private attributes
    Matrix H_;
    Matrix Eigenvectors_;
    Vector eigenvalues_;

public:
    uint n_;

    // Constructors/Destructors
    Hamiltonian();
    Hamiltonian(uint N);
    virtual ~Hamiltonian();

    // Public methods

    const Vector& eigenValues() const {
        return eigenvalues_;
    }

    const Matrix& eigenVectors() const {
        return Eigenvectors_;
    }
    void build2ParticleHamiltonian(const Particle& p1, const Particle& p2, double distance, Vector3D cosines);
    void build(const ParticleList& particle_list);
    void computeEigs();
    void print();
    void calculateCharge(ParticleList& particle_list);

    // Operators
    double& operator()(unsigned row, unsigned col);
    double operator()(unsigned row, unsigned col) const;

private:
    // Private methods
    double scale_func(double d);
};

inline
double& Hamiltonian::operator()(unsigned row, unsigned col) {
    return H_(row, col);
}

inline
double Hamiltonian::operator()(unsigned row, unsigned col) const {
    return H_(row, col);
}

#endif // HAMILTONIAN_H
