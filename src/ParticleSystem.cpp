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

#include "ParticleSystem.h"
// Constructors/Destructors

ParticleSystem::ParticleSystem() {
	setThermostats();
}

ParticleSystem::ParticleSystem(uint N) {
	setThermostats();
	no_particles_ = 0;
	gamma_ = DEFAULT_GAMMA;
	thermostat_ = DEFAULT_THERMOSTAT;
	target_temperature_ = 0.0;
	particle_list_ = ParticleList(N);
	attract_energy_ = AttractiveEnergy(N);
	rep_energy_ = RepulsiveEnergy();
	maxDisplacement1_ = 0;
	maxDisplacement2_ = 0;
	D_ = arma::zeros(3 * N, 3 * N);
	nm_E_Vectors_ = new Matrix();
	NormalModes_ = new Vector();
}

ParticleSystem::~ParticleSystem() {
	delete nm_E_Vectors_;
	delete NormalModes_;
}

// Methods

void ParticleSystem::updateRepEnergy() {
	rep_energy_.computeEnergy(particle_list_);
}

void ParticleSystem::updateBondEnergy() {
	attract_energy_.update(particle_list_);
	attract_energy_.computeEnergy(particle_list_);
}

void ParticleSystem::updateNeighbourList() {
	int u, v, w;
	double DSQ; // squared distance between the two atoms in Angstroms
	arma::vec dx = arma::zeros<arma::vec>(3);
	arma::vec dy = arma::zeros<arma::vec>(3);
	arma::vec dz = arma::zeros<arma::vec>(3);

	Neighbour n;
	ParticleList::iterator p, p_next; // pointer 
	clearNeighbourList(); // clear number of neighbours
	for (p = particle_list_.begin(); p != particle_list_.end(); p++) { // choose one atom to focus on (atom a)
		for (p_next = particle_list_.begin(); p_next != particle_list_.end(); p_next++) { // loop over all atoms b, s.t. b<a 
			if (p->id() < p_next->id()) { // making sure b<a
				if (!pbc_) {
					if (p->distance(*p_next) <= RL) { // check if within LIST radius
						// Add Neighbour n to both particles (to avoid calculating the same thing twice)
						n.particle_ptr = &(*p_next); // get address for pointer to the atom
						p->addNeighbour(n);
						n.particle_ptr = &(*p);
						p_next->addNeighbour(n);
					}
				} 
				else {
					dx(0) = cell_dims_.x - p_next->position().x + p->position().x;
					dx(1) = p_next->position().x - p->position().x;
					dx(2) = cell_dims_.x + p_next->position().x - p->position().x;

					dy(0) = cell_dims_.y - p_next->position().y + p->position().y;
					dy(1) = p_next->position().y - p->position().y;
					dy(2) = cell_dims_.y + p_next->position().y - p->position().y;

					dz(0) = cell_dims_.z - p_next->position().z + p->position().z;
					dz(1) = p_next->position().z - p->position().z;
					dz(2) = cell_dims_.z + p_next->position().z - p->position().z;

					for (u = 0; u < 3; u++) {
						for (v = 0; v < 3; v++) {
							for (w = 0; w < 3; w++) {
								DSQ = dx(u) * dx(u) + dy(v) * dy(v) + dz(w) * dz(w);
								if (DSQ < RL * RL) {
									n.particle_ptr = &(*p_next); // get address for pointer to the atom
									p->addNeighbour(n);
									n.particle_ptr = &(*p);
									p_next->addNeighbour(n);
								}
							}
						}
					}
				}
			}
		}
	}
	updateNN();
}

void ParticleSystem::updateNN() { //TODO: Should be entirely replaced by updateNeighbourList()
	int u, v, w;
	double DSQ; // squared distance between the two atoms in Angstroms
	arma::vec dx = arma::zeros<arma::vec>(3);
	arma::vec dy = arma::zeros<arma::vec>(3);
	arma::vec dz = arma::zeros<arma::vec>(3);

	NearestNeighbour nn; // structure of type NearestNeighbour
	ParticleList::iterator p; // pointer
	Particle* p_next;
	NeighbourList::const_iterator n;
	clearNN(); // clear number of neighbours
	uint k;
	for (p = particle_list_.begin(); p != particle_list_.end(); p++) { // choose one atom to focus on (atom a)      
		NeighbourList& Neighbour_List = p->neighbourList();
		n = Neighbour_List.begin();
		k = 0;
		while (n != Neighbour_List.end()&&(k < p->noNeighbours())) { //TODO: Include argument that checks if Number_of_Neighbours is okay
			p_next = n->particle_ptr;
			if (p->id() < p_next->id()) { // making sure b<a
				nn.distance = p->distance(*p_next); // gives distance between p and p_next

				if (!pbc_) {
					nn.distance = p->distance(*p_next); // gives distance between p and p_next

					if (nn.distance <= RC) { // check if within cut-off radius
						nn.particle_ptr = &(*p_next); // get address for pointer to the atom b
						nn.cosines = (p_next->position() - p->position()) / nn.distance;
						p->addNearestNeighbour(nn);
						nn.particle_ptr = &(*p);
						nn.cosines = -nn.cosines;
						p_next->addNearestNeighbour(nn);
					}
				} else {
					dx(0) = -(cell_dims_.x - p_next->position().x + p->position().x);
					dx(1) = p_next->position().x - p->position().x;
					dx(2) = cell_dims_.x + p_next->position().x - p->position().x;

					dy(0) = -(cell_dims_.y - p_next->position().y + p->position().y);
					dy(1) = p_next->position().y - p->position().y;
					dy(2) = cell_dims_.y + p_next->position().y - p->position().y;

					dz(0) = -(cell_dims_.z - p_next->position().z + p->position().z);
					dz(1) = p_next->position().z - p->position().z;
					dz(2) = cell_dims_.z + p_next->position().z - p->position().z;
					for (u = 0; u < 3; u++) {
						for (v = 0; v < 3; v++) {
							for (w = 0; w < 3; w++) {
								DSQ = dx(u) * dx(u) + dy(v) * dy(v) + dz(w) * dz(w);
								if (DSQ <= RC * RC) {
									nn.distance = sqrt(DSQ);
									nn.particle_ptr = &(*p_next); // get address for pointer to the atom b
									nn.cosines.x = dx(u) / nn.distance; //TODO: is the sign of dx important?
									nn.cosines.y = dy(v) / nn.distance;
									nn.cosines.z = dz(w) / nn.distance;
									p->addNearestNeighbour(nn);
									nn.particle_ptr = &(*p);
									nn.cosines = -nn.cosines;
									p_next->addNearestNeighbour(nn);
								}
							}
						}
					}
				}
			}
			n++;
			k++;
		}
	}
}

void ParticleSystem::clearNN() {
	ParticleList::iterator p;
	for (p = particle_list_.begin(); p != particle_list_.end(); p++)
		p->clearNNList();
}

void ParticleSystem::clearNeighbourList() {
	ParticleList::iterator p;
	for (p = particle_list_.begin(); p != particle_list_.end(); p++)
		p->clearNeighbourList();
}

void ParticleSystem::clearParticlesForce() {
	ParticleList::iterator p;
	for (p = particle_list_.begin(); p != particle_list_.end(); p++)
		p->setForceByCoords(0.0, 0.0, 0.0);
}

void ParticleSystem::updateParticlesForce() {
	attract_energy_.update(particle_list_);
	ParticleList::iterator p;
	clearParticlesForce();
	for (p = particle_list_.begin(); p != particle_list_.end(); p++) {
		attract_energy_.updateParticleNetForce(*p);
		rep_energy_.updateParticleNetForce(*p);
	}
}

void ParticleSystem::moveParticleInPeriodicBC(Particle& p, Displacement & dPos) {
	Position NewPosition;
	NewPosition = p.position() + dPos;

	if (NewPosition.x < 0) {
		NewPosition.x = cell_dims_.x + (p.position().x + dPos.x);
	} else if (NewPosition.x > cell_dims_.x) {
		NewPosition.x = p.position().x + dPos.x - cell_dims_.x;
	} else if (NewPosition.y < 0) {
		NewPosition.y = cell_dims_.y + (p.position().y + dPos.y);
	} else if (NewPosition.y > cell_dims_.y) {
		NewPosition.y = p.position().y + dPos.y - cell_dims_.y;
	} else if (NewPosition.z < 0) {
		NewPosition.z = cell_dims_.z + (p.position().z + dPos.z);
	} else if (NewPosition.z > cell_dims_.z) {
		NewPosition.z = p.position().z + dPos.z - cell_dims_.z;
	}

	p.position() = NewPosition;

}

void ParticleSystem::moveParticles(double dt, gsl_rng * r) {
	ParticleList::iterator p;
	if (thermostat_ == VV_DAMPED) {
		Displacement dPos;
		double MagnitudeOfNetDisplacement;

		for (p = particle_list_.begin(); p != particle_list_.end(); p++) {

			dPos = (p->velocity() + MD_COEFF * (-gamma_ * C_MASS * p->velocity() + p->force()) * dt) * dt;
			p->displacement() += dPos;
			if (!pbc_) {
				p->position() += dPos;
			} else {
				moveParticleInPeriodicBC((*p), dPos);
			}
			dPos = p->displacement();

			MagnitudeOfNetDisplacement = norm(dPos.x, dPos.y, dPos.z);
			if (MagnitudeOfNetDisplacement > maxDisplacement1_) {
				maxDisplacement2_ = maxDisplacement1_;
				maxDisplacement1_ = MagnitudeOfNetDisplacement;
			} else if (MagnitudeOfNetDisplacement > maxDisplacement2_) {
				maxDisplacement2_ = MagnitudeOfNetDisplacement;
			}
		}

		updateNN();

		for (p = particle_list_.begin(); p != particle_list_.end(); p++) {
			p->velocity() += MD_COEFF * (-gamma_ * C_MASS * p->velocity() + p->force()) * dt;
		}

		attract_energy_.update(particle_list_);
		clearParticlesForce();

		for (p = particle_list_.begin(); p != particle_list_.end(); p++) {
			attract_energy_.updateParticleNetForce(*p);
			rep_energy_.updateParticleNetForce(*p);
			p->velocity() += MD_COEFF * (-gamma_ * C_MASS * p->velocity() + p->force()) * dt;
		}


	} 
	else if (thermostat_ == OVRVO) {
		double a, b, roota, factor1, factor2, MagnitudeOfNetDisplacement;
		Velocity v1q, v2q, v3q;
		Displacement dis;

		a = exp(-gamma_ * dt), b = 1.0;
		roota = sqrt(a);
		factor1 = sqrt((1 - a) * target_temperature_ / C_MASS);
		factor2 = b * dt / 2.0;

		for (p = particle_list_.begin(); p != particle_list_.end(); p++) {
			v1q = roota * p->velocity();
			v1q.x += factor1 * gsl_ran_gaussian(r, 1);
			v1q.y += factor1 * gsl_ran_gaussian(r, 1);
			v1q.z += factor1 * gsl_ran_gaussian(r, 1);
			Force f = p->force();
			v2q = v1q + factor2 * f / C_MASS;
			p->velocity2() = v2q;
			dis = 2 * factor2*v2q;
			p->displacement() += dis; //p->setDisplacementByCoords(p->displacement().x+dis.x, p->displacement().y+dis.y, p->displacement().z+dis.z);
			if (!pbc_) {
				p->position() += dis;
			} else {
				moveParticleInPeriodicBC((*p), dis);
			}
			dis = p->displacement(); //get the net displacement
			MagnitudeOfNetDisplacement = norm(dis.x, dis.y, dis.z);

			if (MagnitudeOfNetDisplacement > maxDisplacement1_) {
				maxDisplacement2_ = maxDisplacement1_;
				maxDisplacement1_ = MagnitudeOfNetDisplacement;
			} else if (MagnitudeOfNetDisplacement > maxDisplacement2_) {
				maxDisplacement2_ = MagnitudeOfNetDisplacement;
			}
		}

		updateNN();
		updateParticlesForce();

		for (p = particle_list_.begin(); p != particle_list_.end(); p++) {

			v3q = roota * (p->velocity2() + factor2 * p->force() / C_MASS);
			v3q.x += factor1 * gsl_ran_gaussian(r, 1);
			v3q.y += factor1 * gsl_ran_gaussian(r, 1);
			v3q.z += factor1 * gsl_ran_gaussian(r, 1);
			p->velocity() = v3q;
		}
	}

	//Check if neighbour lists have to be updated:
	if ((maxDisplacement1_ + maxDisplacement2_) >= (RL - RC)) { //update neighbour list
		ParticleList::iterator p;
		for (p = particle_list_.begin(); p != particle_list_.end(); p++) { //TODO: maybe it's easier/faster to define one null vector and write p->displacement()=NULL_VECTOR
			p->displacement().x = 0;
			p->displacement().y = 0;
			p->displacement().z = 0;
		}
		maxDisplacement1_ = 0;
		maxDisplacement2_ = 0;

		updateNeighbourList();
	}
}

void ParticleSystem::addParticle(Particle p) {
	p.id(no_particles_);
	particle_list_[no_particles_] = p;
	no_particles_++;
}

void ParticleSystem::printParticles() {
	ParticleList::iterator p;
	for (p = particle_list_.begin(); p != particle_list_.end(); p++) {
		p->print();
	}
}

void ParticleSystem::setThermalVelocities() {

	//Initialise random number generator:
	const gsl_rng_type * T; //gsl_rng_type holds static information about each type of generator
	gsl_rng * r; //gsl_rng describes an instance of a generator created from a given gsl_rng_type
	T = gsl_rng_ranlxs2; //picks the random number generator. Performance: 769 k doubles/sec, mt19937
	r = gsl_rng_alloc(T); // Returns a pointer to a newly-created instance of a random number generator of type T.
	gsl_rng_set(r, clock()); // Seeds the instance r with the given seed (here with the value of clock()).
	//End initialisation


	ParticleList::iterator p;
	double sigma = sqrt((3 * PI - 8) / (PI)) * sqrt(target_temperature_ / C_MASS), RandomNumber, Vx = 0, Vy = 0, Vz = 0;

	for (p = particle_list_.begin(); p != particle_list_.end(); p++) {

		RandomNumber = gsl_ran_gaussian(r, sigma);
		p->velocity().x = RandomNumber;
		Vx += RandomNumber;

		RandomNumber = gsl_ran_gaussian(r, sigma);
		p->velocity().y = RandomNumber;
		Vy += RandomNumber;

		RandomNumber = gsl_ran_gaussian(r, sigma);
		p->velocity().z = RandomNumber;
		Vz += RandomNumber;
	}
	gsl_rng_free(r);


	for (p = particle_list_.begin(); p != particle_list_.end(); p++) {
		p->velocity().x -= Vx / no_particles_;
		p->velocity().y -= Vy / no_particles_;
		p->velocity().z -= Vz / no_particles_;
	}
	getTemperature();
}

void ParticleSystem::computeKineticEnergy() {
	ParticleList::iterator p;

	//Get the velcoity of the center of mass
	double dVx, dVy, dVz;
	Velocity V_com, V_p; // com = center of mass
	V_com.x = 0., V_com.y = 0., V_com.z = 0.;
	for (p = particle_list_.begin(); p != particle_list_.end(); p++) {
		V_com += p->velocity();
	}
	kinetic_energy_ = 0.;
	for (p = particle_list_.begin(); p != particle_list_.end(); p++) {
		V_p = p->velocity();
		dVx = V_p.x - V_com.x / no_particles_;
		dVy = V_p.y - V_com.y / no_particles_;
		dVz = V_p.z - V_com.z / no_particles_;
		kinetic_energy_ += 0.5 * C_MASS * ((dVx * dVx)+(dVy * dVy)+(dVz * dVz)); //Note that the kinetic energy is taken relative to the center of mass.
	}
}

double ParticleSystem::energy() {
	attract_energy_.computeEnergy(particle_list_);
	rep_energy_.computeEnergy(particle_list_);
	computeKineticEnergy();
	return rep_energy_.energy() + attract_energy_.energy() + kinetic_energy_;
}

double ParticleSystem::getTemperature() {
	computeKineticEnergy();
	double T;
	if (target_temperature_ == 0.0)
		T = 0.0;
	else
		T = 2. / (3. * (no_particles_ - 1) * KB) * kinetic_energy_;
	return T;
}

void ParticleSystem::computeNormalModes() {
	ParticleList::iterator p;
	delete NormalModes_;
	NormalModes_ = new Vector(3 * no_particles_);
	NormalModes_->zeros();
	delete nm_E_Vectors_;
	nm_E_Vectors_ = new Matrix(3 * no_particles_, 3 * no_particles_);
	for (p = particle_list_.begin(); p != particle_list_.end(); p++) {
		updateD_Row(p->position().x, *p, 0);
		updateD_Row(p->position().y, *p, 1);
		updateD_Row(p->position().z, *p, 2);
	}
	D_ = D_ / (-12.0 * C_MASS * DELTA);

	arma::eig_sym(*NormalModes_, *nm_E_Vectors_, D_);
}

void ParticleSystem::updateD_Row(double& coord, Particle& p, int ncoord) {
	double r = coord;

	coord += DELTA;
	updateNN();
	attract_energy_.update(particle_list_);
	clearParticlesForce();
	int k = 0;
	ParticleList::iterator p_next;
	for (p_next = particle_list_.begin(); p_next != particle_list_.end(); p_next++) {
		attract_energy_.updateParticleNetForce(*p_next);
		rep_energy_.updateParticleNetForce(*p_next);
		D_(3 * p.id() + ncoord, 3 * p_next->id() + 0) += 8.0 * p_next->force().x;
		D_(3 * p.id() + ncoord, 3 * p_next->id() + 1) += 8.0 * p_next->force().y;
		D_(3 * p.id() + ncoord, 3 * p_next->id() + 2) += 8.0 * p_next->force().z;
		k++;
	}

	coord = r;
	coord -= DELTA;
	updateNN();
	attract_energy_.update(particle_list_);
	clearParticlesForce();
	k = 0;
	for (p_next = particle_list_.begin(); p_next != particle_list_.end(); p_next++) {
		attract_energy_.updateParticleNetForce(*p_next);
		rep_energy_.updateParticleNetForce(*p_next);
		D_(3 * p.id() + ncoord, 3 * p_next->id() + 0) -= 8.0 * p_next->force().x;
		D_(3 * p.id() + ncoord, 3 * p_next->id() + 1) -= 8.0 * p_next->force().y;
		D_(3 * p.id() + ncoord, 3 * p_next->id() + 2) -= 8.0 * p_next->force().z;
		k++;
	}


	coord = r;
	coord -= 2.0 * DELTA;
	updateNN();
	attract_energy_.update(particle_list_);
	clearParticlesForce();
	k = 0;
	for (p_next = particle_list_.begin(); p_next != particle_list_.end(); p_next++) {
		attract_energy_.updateParticleNetForce(*p_next);
		rep_energy_.updateParticleNetForce(*p_next);
		D_(3 * p.id() + ncoord, 3 * p_next->id() + 0) += p_next->force().x;
		D_(3 * p.id() + ncoord, 3 * p_next->id() + 1) += p_next->force().y;
		D_(3 * p.id() + ncoord, 3 * p_next->id() + 2) += p_next->force().z;
		k++;
	}


	coord = r;
	coord += 2.0 * DELTA;
	updateNN();
	attract_energy_.update(particle_list_);
	clearParticlesForce();
	k = 0;
	for (p_next = particle_list_.begin(); p_next != particle_list_.end(); p_next++) {
		attract_energy_.updateParticleNetForce(*p_next);
		rep_energy_.updateParticleNetForce(*p_next);
		D_(3 * p.id() + ncoord, 3 * p_next->id() + 0) -= p_next->force().x;
		D_(3 * p.id() + ncoord, 3 * p_next->id() + 1) -= p_next->force().y;
		D_(3 * p.id() + ncoord, 3 * p_next->id() + 2) -= p_next->force().z;
		k++;
	}

	coord = r;
}

Vector* ParticleSystem::getNormalModes() {
	return NormalModes_;
}

Matrix* ParticleSystem::getNormalModesV() {
	return nm_E_Vectors_;
}

bool ParticleSystem::isValidThermostat(std::string thermostat) {
	ThermostatListType::iterator p;
	for (p = thermostat_list_.begin(); p != thermostat_list_.end(); p++)
		if (thermostat_list_.find(thermostat) != thermostat_list_.end())
			return true;
	return false;
}

void ParticleSystem::setThermostats() {
	thermostat_list_["VV_DAMPED"] = VV_DAMPED;
	thermostat_list_["OVRVO"] = OVRVO;
}

void ParticleSystem::checkErrors() {
	ParticleList::iterator p1, p2;
	for (p1 = particle_list_.begin(); p1 != particle_list_.end(); p1++)
		for (p2 = particle_list_.begin(); p2 != particle_list_.end(); p2++) {
			if ((p1->id() != p2->id()) && norm(p1->position().x - p2->position().x, p1->position().y - p2->position().y, p1->position().z - p2->position().z) == 0.0)
				error("","Geometric error","Two carbon atoms cannot lay at the same point. Check the input data file.");
		}
}

