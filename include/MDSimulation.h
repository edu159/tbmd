/* Class which implements an MD simulation.
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


#ifndef MDSIMULATION_H
#define MDSIMULATION_H

#include "common.h"
#include "ParticleSystem.h"
#include "File.h"
#include "ConfigFile.h"
#include "InputDataFile.h"
#include "OutputDataFile.h"
#include "MovieFile.h"

/**
 * class MDSimulation
 * 
 */

class MDSimulation {
public:

private:

    // Private attributes
    ParticleSystem* particle_system_;
    double time_step_;
		MovieFile* movie_file_;
		InputDataFile* input_data_file_;
		ConfigFile* config_file_;
		OutputDataFile* output_data_file_;
    bool compute_normal_modes_;
		bool create_movie_file_;
		uint number_of_steps_;
		uint movie_frame_rate_;
		uint output_data_rate_;
		bool quiet_;
		

public:

    // Setters

    void timeStep(double time_step) {
        time_step_ = time_step;
    }

    // Getters

    double timeStep() {
        return time_step_;
    }

public:

    // Constructors/Destructors
    MDSimulation();
    MDSimulation(std::string input_fname, std::string config_fname, std::string output_name);
    virtual ~MDSimulation();
		void configureSimulation(const ParamList& param_list);

    // Public methods
    void step(gsl_rng * r);
    void run();
		void quietMode(bool mode) {quiet_ = mode;}
};

#endif // MDSIMULATION_H
