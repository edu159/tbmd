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

#ifndef OUTPUT_DATA_FILE_H
#define OUTPUT_DATA_FILE_H

#include "File.h"
#include "Particle.h"
#include "ParticleSystem.h"
#include "common.h"

typedef std::map<std::string, OutputModeType> ModeListType;

class OutputDataFile: public OutputFile {
private:
	OutputModeType mode_;
	ModeListType mode_list_;
	std::string nmodes_fname_;
	
public:
	OutputDataFile();
	OutputDataFile(std::string fname);
	OutputModeType mode(){return mode_;};
	void mode(OutputModeType mode) {mode_ = mode;}
	virtual ~OutputDataFile();
	void mode(std::string mode);
	bool isValidMode(std::string mode);
	void writeStepData(ParticleSystem* particle_system, uint step);
	void writeNormalModes(ParticleSystem* particle_system);
	void nmodesFileName(std::string fname) {nmodes_fname_ = fname;};
private:
	void setModes(); 
	void basicMode(ParticleSystem* particle_system, uint step);
	void detailedMode(ParticleSystem* particle_system, uint step);
	void debugMode(ParticleSystem* particle_system, uint step);
};

#endif
