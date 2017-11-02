/* Carbon particle object. 
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

#include "Particle.h"
#include "common.h"
// Constructors/Destructors

Particle::Particle() {
    vel_.x = 0.0;
    vel_.y = 0.0;
    vel_.z = 0.0;
    pos_.x = 0.0;
    pos_.y = 0.0;
    pos_.z = 0.0;
    dis_.x = 0.0;
    dis_.y = 0.0;
    dis_.z = 0.0;
    force_.x = 0.0;
    force_.y = 0.0;
    force_.z = 0.0;
    charge_ = 0.0;
    no_NN_ = 0;
    //TODO:CHange this
    NN_list_ = NearestNeighboursList(MAX_NNEIGHBOURS);
    N_list_ = NeighbourList(MAX_NEIGHBOURS);
}

Particle::~Particle() {
}

// Methods
void Particle::print() const {
    std::cout << "Particle " << id_ << ": " << std::endl;
    std::cout << "Force:" << force_.x << " "<<  force_.y << " " << force_.z << std::endl;
    std::cout << "Velocity:" << vel_.x << " "<<  vel_.y << " " << vel_.z << std::endl;
    std::cout << "Position:" << pos_.x << " "<<  pos_.y << " " << pos_.z << std::endl;
    std::cout << "Number of NN : " << no_NN_ << std::endl;
    std::cout << std::endl;
}

