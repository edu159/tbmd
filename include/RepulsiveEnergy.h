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

#ifndef REPULSIVEENERGY_H
#define REPULSIVEENERGY_H

#include "Energy.h"
#include "common.h"
#include "Particle.h"

/**
 * class RepulsiveEnergy
 * 
 */
class RepulsiveEnergy : public Energy {
public:

    // Constructors/Destructors
    RepulsiveEnergy();
    virtual ~RepulsiveEnergy();

    // Public methods
    double computeEnergy(ParticleList& particle_list);
    void updateParticleNetForce(Particle& p);


private:


};
// Private methods
double poly_func(double x);
double phi(double r);
double A(Particle& p);
double dfSumPhi(Particle& p);
double dPhi(double r);

#endif // REPULSIVEENERGY_H
