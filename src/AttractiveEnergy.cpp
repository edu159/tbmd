/* Attractive energy of the Hamiltonian.
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

#include "AttractiveEnergy.h"

using namespace std;
using namespace arma;

// Constructors/Destructors

AttractiveEnergy::AttractiveEnergy() {
}

AttractiveEnergy::AttractiveEnergy(uint N) {
	hamiltonian_ = Hamiltonian(N);
	rho_matrix_ = Matrix(4 * N, 4 * N);
}

AttractiveEnergy::~AttractiveEnergy() {

}

// Methods

void AttractiveEnergy::permuteCosines(double &dx, double&dy, double&dz, Vector3D cosines, int per) {
	if (per == 0) {
		dx = cosines.x; dy = cosines.y; dz = cosines.z;
	} 
	else if (per == 1) {
		dx = cosines.y; dy = cosines.z; dz = cosines.x;
	} 
	else if (per == 2) {
		dx = cosines.z; dy = cosines.x; dz = cosines.y;
	} 
	else 
		std::cout << "Wrong per, AttractiveEnergy.cpp" << std::endl;
}


double AttractiveEnergy::calcTopR(int dir, double s, double ds_dr, double distance, Vector3D cosines, int per) {
	double element;
	double dx, dy, dz;
	permuteCosines(dx, dy, dz, cosines, per);

	if (per == dir) {
		element = Vspsi*(ds_dr*dx*dx+s/distance-s*dx*dx/distance);
	} else  {
		if ( (per==0 && dir==1) || (per==1 && dir==2) || (per==2 && dir==0) ) {
			element = Vspsi*dx*dy*(ds_dr-s/distance);
		}
		else if ( (per==0 && dir==2) || (per==1 && dir==0) || (per==2 && dir==1) ) {
			element = Vspsi*dx*dz*(ds_dr-s/distance);
		} else { std::cout << "Wrong per or dir, AttractiveEnergy.cpp" << std::endl; }
	}
	return element;
}

double AttractiveEnergy::calcOffD(int dir, double s, double ds_dr, double distance, Vector3D cosines, int per) {
	double element;
	double dx, dy, dz;
	double W = Vppsi-Vpppi;
	permuteCosines(dx, dy, dz, cosines, per);

	if (per==dir) {
		element = W*dy*(ds_dr*dx*dx+s/distance-2*dx*dx*s/distance);
	} else {
		if ( (per==0 && dir==1) || (per==1 && dir==2) || (per==2 && dir==0)) {
			element = W*dx*(ds_dr*dy*dy+s/distance-2*dy*dy*s/distance);
		} else if ( (per==0 && dir==2) || (per==1 && dir==0) || (per==2 && dir==1)) {
			element = W*dx*dy*dz*(ds_dr-2*s/distance);
		} else { std::cout << "Wrong per or dir, AttractiveEnergy.cpp" << std::endl; }
	}
	return element;
}

double AttractiveEnergy::calcDiag(int dir, double s, double ds_dr, double distance, Vector3D cosines, int per) {
	double element;
	double dx, dy, dz;
	double W = Vppsi-Vpppi;
	permuteCosines(dx, dy, dz, cosines, per);

	if (per == dir) {
		element = ds_dr*dx*(dx*dx*W+Vpppi)+2*s*W*dx*(1-dx*dx)/distance;
	} else {
		if ( (per==0 && dir==1) || (per==1 && dir==2) || (per==2 && dir==0) ) {
			element = ds_dr*(dx*dx*dy*W+dy*Vpppi) - 2*s*W*dx*dx*dy/distance;
		} else if ( (per==0 && dir==2) || (per==1 && dir==0) || (per==2 && dir==1) ) {
			element = ds_dr*(dx*dx*dz*W+dz*Vpppi) - 2*s*W*dx*dx*dz/distance;
		} else { std::cout << "Wrong dir, AttractiveEnergy.cpp" << std::endl;}
	}
	return element;
}

void AttractiveEnergy::calcDeltaH(arma::mat& DeltaH, int dir, double distance, Vector3D cosines) {
	double s = scale_func_copy(distance);
	double ds_dr = scale_func_derivative(distance, s);
	double cos;
	if (dir==0) {cos = cosines.x;}
	else if (dir==1) {cos = cosines.y;}
	else if (dir==2) {cos=cosines.z;}

	DeltaH(0, 4) = Vsssi * ds_dr * cos;

	for (int i=0; i<3; i++) {
		DeltaH(0,5+i) = calcTopR(dir, s, ds_dr, distance, cosines, i);
		DeltaH(i+1,4) = -DeltaH(0,5+i);
		DeltaH(i+1,5+i) = calcDiag(dir, s, ds_dr, distance, cosines, i);
	}
	DeltaH(1,6) = calcOffD(dir, s, ds_dr, distance, cosines, 0);
	DeltaH(2,5) = DeltaH(1,6);
	DeltaH(2,7) = calcOffD(dir, s, ds_dr, distance, cosines, 1);
	DeltaH(3,6) = DeltaH(2,7);
	DeltaH(1,7) = calcOffD(dir, s, ds_dr, distance, cosines, 2);
	DeltaH(3,5) = DeltaH(1,7);

	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			DeltaH(i+4,j) = DeltaH(i,4+j);
			if ((i==0) || (j==0)){
				DeltaH(i+4,j) *= -1.0;
			}
		}
	}
	DeltaH(4,0) = DeltaH(0,4);
}

double AttractiveEnergy::scale_func_copy(double d) {

	double scaling;

	const double r_1 = 2.45, r_m = 2.6;

	const double A_s = 2.89933161442;

	const double cs_0 = 6.7392620074314E-3;
	const double cs_1 = -8.1885359517898E-2;
	const double cs_2 = 0.1932365259144;
	const double cs_3 = 0.354287433238;

	if (d > r_m)
		scaling = 0.0;

	else if (d <= r_1)
		scaling = A_s * exp(Bs * pow(d, Fnc)) / pow(d, Fn);

	else
		scaling = cs_3 * (d - r_1)*(d - r_1)*(d - r_1) + cs_2 * (d - r_1)*(d - r_1) + cs_1 * (d - r_1) + cs_0;

	return scaling;
}

double AttractiveEnergy::scale_func_derivative(double d, double scaling) {

	double ds_dr;

	const double r_1 = 2.45, r_m = 2.6;

	const double A_s = 2.89933161442; // A=[(r_0)^n]*e^[n*(r_0/r_c)^n_c]	

	const double cs_0 = 6.7392620074314E-3;
	const double cs_1 = -8.1885359517898E-2;
	const double cs_2 = 0.1932365259144;
	const double cs_3 = 0.354287433238;

	if (d > r_m)
		ds_dr = 0.0;

	else if (d <= r_1)
		ds_dr = (-Fn+Fnc*Bs*pow(d, Fnc))*scaling/d;

	else
		ds_dr = 3.0 * cs_3 *(d - r_1)*(d - r_1) + 2.0* cs_2 *(d - r_1) + cs_1;

	return ds_dr;
}


Force AttractiveEnergy::compute2ParticleForce(const Particle& p1, const Particle& p2, double distance, Vector3D cosines) {
	const Matrix& Ev = hamiltonian_.eigenVectors();
	uint N = hamiltonian_.n_ / 4.0;
	Force f;
	f.x = 0.0;
	f.y = 0.0;
	f.z = 0.0;
	uint i, j, k;

	arma::mat DeltaH = zeros(8,8);
	calcDeltaH(DeltaH, 0, distance, cosines); 
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			for (k = 0; k < (2 * N); k++) {
				f.x += 4.0 * Ev(4 * p1.id() + i, k) * (DeltaH(i,j+4))   * Ev(4 * p2.id() + j, k);
			}
		}
	}

	DeltaH.zeros();	
	calcDeltaH(DeltaH,1, distance, cosines); 
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			for (k = 0; k < (2 * N); k++) {
				f.y += 4.0 * Ev(4 * p1.id() + i, k) * (DeltaH(i,j+4))   * Ev(4 * p2.id() + j, k);
			}
		}
	}

	DeltaH.zeros();	
	calcDeltaH(DeltaH,2, distance, cosines); 
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			for (k = 0; k < (2 * N); k++) {
				f.z += 4.0 * Ev(4 * p1.id() + i, k) * (DeltaH(i,j+4))   * Ev(4 * p2.id() + j, k);
			}
		}
	}

	if (fabs(f.x) < TOL)
		f.x = 0.0;
	if (fabs(f.y) < TOL)
		f.y = 0.0;
	if (fabs(f.z) < TOL)
		f.z = 0.0;
	return f;
}

double AttractiveEnergy::computeEnergy(ParticleList& particle_list) {
	const Vector& Evals = hamiltonian_.eigenValues();
	tot_energy_ = 0;
	for (int i = 0; i < 2 * particle_list.size(); i++)
		tot_energy_ += Evals(i);
	tot_energy_ *= 2;
}

void AttractiveEnergy::updateParticleNetForce(Particle& p) {
	uint k = 0;
	NearestNeighboursList::const_iterator nn;
	NearestNeighboursList & nn_list = p.nearestNeighboursList();
	Particle* p_next;
	nn = nn_list.begin();
	Force f;

	while (nn != nn_list.end() && (k < p.noNearestNeighbours())) {
		p_next = nn->particle_ptr;
		if (p_next->id() > p.id()) {
			f = compute2ParticleForce(p, *p_next, nn->distance, nn->cosines);
			p.force() += f;
			p_next->force() -= f;
		}
		k++;
		nn++;
	}
}

void AttractiveEnergy::GetDensityMatrix(const Matrix& EV) {
	uint M, i, j, s;
	uint Ne = hamiltonian_.n_;
	M = Ne / 2;
	rho_matrix_.zeros();
	for (i = 0; i < Ne; i++) {
		for (j = 0; j < Ne; j++) {
			for (s = 0; s < M; s++) rho_matrix_(i, j) += (EV(i, s) * EV(j, s));
		}
	}
}

void AttractiveEnergy::update(ParticleList& particle_list) {
	hamiltonian_.build(particle_list);
	hamiltonian_.computeEigs();
	hamiltonian_.calculateCharge(particle_list);
	const Matrix& EV = hamiltonian_.eigenVectors();
	GetDensityMatrix(EV);
	rho_matrix_ *= -1;
}
