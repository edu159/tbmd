/* Base class for Hamiltonian energies.
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


#ifndef ENERGY_H
#define ENERGY_H

#include "common.h"
#include "Particle.h"


/**
 * class Energy
 * 
 */

class ParticleSystem;

class Energy {
protected:
    // Protected attributes
    double tot_energy_;

public:
    void energy(double energy);

    const double energy() {
        return tot_energy_;
    }

public:

    // Constructors/Destructors
    Energy();
    virtual ~Energy();

    // Public methods
    virtual double computeEnergy(ParticleList& particle_list) = 0;
    virtual void updateParticleNetForce(Particle& p) = 0;

};

#endif // ENERGY_H
