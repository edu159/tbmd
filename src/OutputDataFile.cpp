/* Output file writer of info about the MD simulation.
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

#include "OutputDataFile.h"

OutputDataFile::OutputDataFile() {}
OutputDataFile::OutputDataFile(std::string fname): OutputFile(fname) {
	setModes();
	mode_ = BASIC;
	nmodes_fname_ = "NormalModes"; 

}

OutputDataFile::~OutputDataFile(){}

void OutputDataFile::mode(std::string mode) {
	mode_ = mode_list_[mode];
}

bool OutputDataFile::isValidMode(std::string mode) {
	if (mode_list_.find(mode) != mode_list_.end()) 
			return true;
	return false;
}

void OutputDataFile::setModes() {
	mode_list_["BASIC"] = BASIC;
	mode_list_["DETAILED"] = DETAILED;
	mode_list_["DEBUG"] = DEBUG;
	mode_list_["LAST_STEP"] = LAST_STEP;
}

void OutputDataFile::writeStepData(ParticleSystem* particle_system, uint step) {
	switch(mode_) {
	case BASIC:
		basicMode(particle_system, step);
		break;
	case DETAILED:
		detailedMode(particle_system, step);
		break;
	case LAST_STEP:
		detailedMode(particle_system, step);
		break;
	case DEBUG:
		debugMode(particle_system, step);
		break;
	}
}

void OutputDataFile::writeNormalModes(ParticleSystem* particle_system) {
	uint N = particle_system->noParticles();
	std::ostringstream fname;
	fname << nmodes_fname_;
	std::ofstream fd_nm;
	fd_nm.open(fname.str().c_str());
	Vector* normal_modes;
	Matrix* normal_modes_v;
	normal_modes = particle_system->getNormalModes();
	normal_modes_v = particle_system->getNormalModesV();
	fd_nm << std::setiosflags(std::ios::fixed)
				<< std::setprecision(5)
				<< std::setw(12);
	for (uint i = 0; i < (3 * N); i++) {
		fd_nm << std::setw(12) << std::left << normal_modes->at(i) << " ";
		for (uint j = 0; j < (3 * N); j++) {
			fd_nm << std::setw(12) << std::left << normal_modes_v->at(j, i) << " ";
		}
		fd_nm << std::endl;
	}
	fd_nm.close();
}


void OutputDataFile::basicMode(ParticleSystem* particle_system, uint step) {
	fd_ << step << " " << particle_system->energy() << " " << particle_system->getTemperature() << std::endl; 
}

void OutputDataFile::detailedMode(ParticleSystem* particle_system, uint step) {
		fd_ << step << " " << particle_system->kineticEnergy() << " " << particle_system->bondEnergy() 
			  << " " << particle_system->kineticEnergy() << " " << particle_system->repEnergy() 
			  << " " << particle_system->energy() << " " << particle_system->getTemperature() << std::endl;
}
void OutputDataFile::debugMode(ParticleSystem* particle_system, uint step) {
}
