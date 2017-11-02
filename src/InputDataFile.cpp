/* Input data file for an MD simulation.
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


#include "InputDataFile.h"

// Constructors
InputDataFile::InputDataFile(){}
InputDataFile::~InputDataFile(){}
InputDataFile::InputDataFile(std::string fname):InputFile(fname){}

// Methods
void InputDataFile::configureSystem(ParticleSystem* p_system) {
	std::string data_line;
	TokenList tokens;
	uint line_no = 0;
	Particle p;
	bool head_line = false;
	uint no_particles = 0;
	uint N;
	while (readLine(data_line)) {
		line_no++;
		if (data_line != "") {
			tokenizeString(data_line, " \t", true, tokens);
			if (tokens[0].at(0) != '%') {
				if (!head_line) {
					if (check_integer_value(tokens[0])) {
						delete p_system;
						N = atoi(tokens[0].c_str());
						p_system = new ParticleSystem(N);
						head_line = true;
					}
					else
						error(DATA_FILE_ERR_MSG, PARSE_ERR_MSG, "The size of the system has to be a valid integer", line_no);
				}
				else {
					if (tokens.size() == 4) {
						if (check_double_value(tokens[1]))
							p.position().x = atof(tokens[1].c_str());
						else
							error(DATA_FILE_ERR_MSG, PARSE_ERR_MSG, "Invalid double value integer", line_no);
						if (check_double_value(tokens[2]))
							p.position().y = atof(tokens[2].c_str());
						else
							error(DATA_FILE_ERR_MSG, PARSE_ERR_MSG, "Invalid double value integer", line_no);
						if (check_double_value(tokens[3]))
							p.position().z = atof(tokens[3].c_str());
						else
							error(DATA_FILE_ERR_MSG, PARSE_ERR_MSG, "Invalid double value integer", line_no);
						p_system->addParticle(p);
						no_particles++;
					}
					else
							error(DATA_FILE_ERR_MSG, PARSE_ERR_MSG, "Ivalid format. Line format is C x y z", line_no);

				}
			}	 
		}
	}
	if (!head_line)
		error(DATA_FILE_ERR_MSG, PARSE_ERR_MSG, "The size of the system was not specified");
	if (no_particles != N) 
		error(DATA_FILE_ERR_MSG, PARSE_ERR_MSG, "The size of the system and the number of particles specified does not match");
	p_system->checkErrors();
}
