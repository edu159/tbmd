/* System of C particles. Normal modes calculation.
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

#ifndef PARTICLESYSTEM_H
#define PARTICLESYSTEM_H

#include "common.h"
#include "Particle.h"
#include "RepulsiveEnergy.h"
#include "AttractiveEnergy.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/**
 * class ParticleSystem
 * 
 */
typedef Vector3D CellDimensions;
typedef std::map<std::string, ThermostatType> ThermostatListType;
class ParticleSystem {
private:
    // Private attributes
		ThermostatListType thermostat_list_;
    uint no_particles_;
    ParticleList particle_list_;
    AttractiveEnergy attract_energy_;
    RepulsiveEnergy rep_energy_;
		ThermostatType thermostat_;
    double target_temperature_;
    double kinetic_energy_;
    double maxDisplacement1_; //Largest displacement of any atom in the system.
    double maxDisplacement2_; //Second largest displacement of any atom in the system.
    Matrix D_;
    Vector* NormalModes_;
    Matrix* nm_E_Vectors_;
		double gamma_;
    CellDimensions cell_dims_;
    bool pbc_;

public:

    // Constructors/Destructors
    ParticleSystem();
    ParticleSystem(uint N);
    ~ParticleSystem();

    // Setters

    void noParticles(uint no_particles) {
        no_particles_ = no_particles;
    }

    void particleList(ParticleList p_list) {
        particle_list_ = p_list;
    }

    void attractiveEnergy(AttractiveEnergy attract_energy) {
        attract_energy_ = attract_energy;
    }

    void repulsiveEnergy(RepulsiveEnergy rep_energy) {
        rep_energy_ = rep_energy;
    }

    // Getters

    uint noParticles() {
        return no_particles_;
    }

    ParticleList& particleList() {
        return particle_list_;
    }

    const AttractiveEnergy& attractiveEnergy() {
        return attract_energy_;
    }

    const RepulsiveEnergy& repulsiveEnergy() {
        return rep_energy_;
    }

    const Particle& particle(uint id) {
        return particle_list_[id];
    }

public:

    // Public methods
    void updateNeighbourList();
    void clearNeighbourList();
    void updateNN();
    void clearNN();
    void clearParticlesForce();
    void resetDisplacements();
    void updateParticlesForce();
    void moveParticleInPeriodicBC(Particle& p, Displacement & dPos);
    void moveParticles(double dt, gsl_rng * r);
    void setThermalVelocities();
    void addParticle(Particle p);
    void printParticles();
    double energy();
		double kineticEnergy() {return kinetic_energy_;}
		double bondEnergy() {return attract_energy_.energy();}
		double repEnergy() {return rep_energy_.energy();}
		void updateRepEnergy();
		void updateBondEnergy();
    double targetTemperature() {return target_temperature_;}
		void targetTemperature(double target_temperature)  {target_temperature_ = KB * target_temperature;}
		double gamma() {return gamma_;}
		void gamma(double gamma) {gamma_ = gamma;}
		double getTemperature();
    void computeNormalModes();
		Vector* getNormalModes();
		Matrix* getNormalModesV();
		void thermostat(std::string thermostat) {thermostat_ = thermostat_list_[thermostat];}
		bool isValidThermostat(std::string);
		void setThermostats();
		bool pbc() {return pbc_;}
		void pbc(bool pbc) { pbc_ = pbc;}
		void cellDimensions(const CellDimensions& cell_dims) {cell_dims_ = cell_dims;}
		void checkErrors();

private:
    void computeKineticEnergy();
    void updateD_Row(double& coord, Particle& p, int ncoord);

};

#endif // PARTICLESYSTEM_H
