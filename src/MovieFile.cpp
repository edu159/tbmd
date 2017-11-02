/* Movie file generator as a XYZ file.
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


#include "MovieFile.h"

MovieFile::MovieFile() {}
MovieFile::MovieFile(std::string fname): OutputFile(fname) {}
MovieFile::~MovieFile(){}

void MovieFile::writeFrame(ParticleSystem* particle_system) {
	ParticleList::iterator p;
	fd_ << particle_system->noParticles() << std::endl;
	fd_ << std::endl;
	ParticleList&  particle_list = particle_system->particleList();
	for (p = particle_list.begin(); p != particle_list.end(); p++) {
		fd_ << "C " << p->position().x << " " << p->position().y << " " << p->position().z << std::endl;
	}
	fd_ << std::endl;
}
