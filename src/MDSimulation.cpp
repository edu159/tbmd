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


#include "MDSimulation.h"
#include "Hamiltonian.h"
#include <time.h>
// Constructors/Destructors

MDSimulation::MDSimulation() {
}

MDSimulation::MDSimulation(std::string data_fname, std::string config_fname, std::string output_fname) {
	ParamList param_list;
	movie_frame_rate_ = DEFAULT_FRAME_RATE;
	output_data_rate_ = DEFAULT_DATA_RATE;
	create_movie_file_ = false;
	quiet_ = false;
	time_step_ = DEFAULT_TIME_STEP;
	particle_system_ = new ParticleSystem();
	input_data_file_ = new InputDataFile(data_fname);
	output_data_file_ = new OutputDataFile(output_fname);
	config_file_ = new ConfigFile(config_fname);
	movie_file_ = new MovieFile();
	input_data_file_->configureSystem(particle_system_);
	param_list = config_file_->getParams();
	configureSimulation(param_list);
 	particle_system_->updateNeighbourList(); 
  particle_system_->updateParticlesForce();
 	particle_system_->setThermalVelocities();
}


MDSimulation::~MDSimulation() {
	delete config_file_;
	delete input_data_file_;
	delete particle_system_;
	delete output_data_file_;
	delete movie_file_;
}

// Methods

void MDSimulation::step(gsl_rng * r) {
	particle_system_->moveParticles(time_step_, r);
  particle_system_->getTemperature();
}

void MDSimulation::run() {
    uint s, i, actual_no_bars;
    s = 0;
    double t1 = time(NULL);
    double t2;
    //Initialise random number generator:
    const gsl_rng_type * T; //gsl_rng_type holds static information about each type of generator
    gsl_rng * r; //gsl_rng describes an instance of a generator created from a given gsl_rng_type
    T = gsl_rng_ranlxs2; //picks the random number generator. Performance: 769 k doubles/sec, mt19937
    r = gsl_rng_alloc(T); // Returns a pointer to a newly-created instance of a random number generator of type T.
    gsl_rng_set(r, clock()); // Seeds the instance r with the given seed (here with the value of clock()).
		uint no_bars = number_of_steps_ / 10;
    //End initialisation
		actual_no_bars = 0;
		bool expected_time = true;
		std::cout << std::endl;
		std::cout << "PROGRESS: "<< std::endl;
    while (s < number_of_steps_) {
			if (create_movie_file_ && ((s % movie_frame_rate_) == 0))
				movie_file_->writeFrame(particle_system_);
			step(r);
			if (!quiet_) {
				if ((s % no_bars) == 0) {
					actual_no_bars++;
					std::cout << actual_no_bars << "0% --- complete" << std::endl;
				}
			}
			s++;
			if ((output_data_file_->mode() != LAST_STEP) && (s % output_data_rate_) == 0 && (s != number_of_steps_))
				output_data_file_->writeStepData(particle_system_, s);	

    }
		
		// Always print the las step
		output_data_file_->writeStepData(particle_system_, s);	
		// Print elapsed time
		t2 = time(NULL);
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
		std::cout << "*Simulation time = " << difftime(t2, t1) << " s." << std::endl;
		// Normal Modes
    if (compute_normal_modes_) {
				std::cout << "*Computing normal modes...";
        particle_system_->computeNormalModes();
				output_data_file_->writeNormalModes(particle_system_);
				std::cout << "Done." << std::endl;
    }
		std::cout << std::endl;
		if (compute_normal_modes_)
			std::cout << "NormalModes-C" << particle_system_->noParticles() << "... created."<< std::endl;
		if (create_movie_file_)
			std::cout << "movie.xyz" <<"... created."<< std::endl;
		std::cout << output_data_file_->filename() << "... created."<< std::endl;
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
}

void MDSimulation::configureSimulation(const ParamList& param_list) { 
	ParamList::const_iterator param;
	// Required parameters
	bool number_of_steps_bool = false;
	bool time_step_bool = false;
	bool cell_dimx_bool = false, cell_dimy_bool = false, cell_dimz_bool = false;
	CellDimensions cell_dims;
	for (param = param_list.begin(); param != param_list.end(); param++) {
		if (param->key == "TEMPERATURE")
			particle_system_->targetTemperature(param->double_value);
		else if (param->key == "NUMBER_OF_STEPS") {
			number_of_steps_ = param->integer_value;
			number_of_steps_bool = true;
		}
		else if (param->key == "CREATE_MOVIE_FILE")
			create_movie_file_ = param->boolean_value;
		else if (param->key == "THERMOSTAT"){
			if (particle_system_->isValidThermostat(param->string_value)){
				particle_system_->thermostat(param->string_value);
			}
			else
				error(CONFIG_FILE_ERR_MSG, PARSE_ERR_MSG , "Please provide a valid thermostat [VV_DAMPED | OVRVO]");
		}
		else if (param->key == "NORMAL_MODES")
			compute_normal_modes_= param->boolean_value;
		else if (param->key == "OUTPUT_FILE_MODE") {
			if (output_data_file_->isValidMode(param->string_value))
				output_data_file_->mode(param->string_value);
			else
				error(CONFIG_FILE_ERR_MSG, PARSE_ERR_MSG , "Please provide a valid mode of output [BASIC | DETAILED | LAST_STEP | DEBUG]");
		}
		else if (param->key == "MD_GAMMA")
			particle_system_->gamma(param->double_value);
		else if (param->key == "TIME_STEP"){
			time_step_ = param->double_value;
			time_step_bool = true;
		}
		else if (param->key == "OUTPUT_DATA_RATE")
			output_data_rate_ = param->integer_value;
		else if (param->key == "MOVIE_FRAME_RATE")
			movie_frame_rate_ = param->integer_value;
		else if (param->key == "CELL_DIMENSION_X") {
			std::ostringstream sstream;
			sstream << 3*RL;
			std::string min_size_str = sstream.str();
			if (param->double_value < MIN_CELL_SIZE)
				error(CONFIG_FILE_ERR_MSG, "Bad value", "Cell dimensions has to be lower than 3*RL = " + min_size_str + ".");
			cell_dimx_bool = true;
			cell_dims.x = param->double_value;
		}
		else if (param->key == "CELL_DIMENSION_Y") {
			std::ostringstream sstream;
			sstream << 3*RL;
			std::string min_size_str = sstream.str();
			if (param->double_value < MIN_CELL_SIZE)
				error(CONFIG_FILE_ERR_MSG, "Bad value", "Cell dimensions has to be lower than 3*RL = " + min_size_str + ".");
			cell_dimy_bool = true;
			cell_dims.y = param->double_value;
		}
		else if (param->key == "CELL_DIMENSION_Z") {
			std::ostringstream sstream;
			sstream << 3*RL;
			std::string min_size_str = sstream.str();
			if (param->double_value < MIN_CELL_SIZE)
				error(CONFIG_FILE_ERR_MSG, "Bad value", "Cell dimensions has to be lower than 3*RL = " + min_size_str + ".");
			cell_dimz_bool = true;
			cell_dims.z = param->double_value;

		}
		else if (param->key == "PERIODIC_BC")
			particle_system_->pbc(param->boolean_value);
		else if (param->key == "MOVIE_FILE_NAME") {
			delete movie_file_;
			movie_file_ = new MovieFile(param->string_value);
		}
		else if (param->key == "NMODES_FILE_NAME") {
			output_data_file_->nmodesFileName(param->string_value);
		}
	}
	if (!number_of_steps_bool)
		error(CONFIG_FILE_ERR_MSG, "Missing parameter", "You must specify the number of steps to perform.");
	if (!time_step_bool)
		error(CONFIG_FILE_ERR_MSG, "Missing parameter", "You must specify the time step.");

	if (particle_system_->pbc())
		if (cell_dimx_bool && cell_dimy_bool && cell_dimz_bool)
			particle_system_->cellDimensions(cell_dims);
}
