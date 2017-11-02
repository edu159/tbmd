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

#ifndef PARTICLE_H
#define PARTICLE_H


#include "common.h"


class Particle;

struct NearestNeighbour { //TODO: Check if this is still needed once NeighbourLists have been introduced completely
    double distance;
    Vector3D cosines;
    Particle* particle_ptr;
};

struct Neighbour { //TODO: Check if this is still needed once NeighbourLists have been introduced completely
    Particle* particle_ptr;
};

typedef std::vector<NearestNeighbour> NearestNeighboursList; //TODO: See above
typedef std::vector<Neighbour> NeighbourList;
typedef std::vector<Particle> ParticleList;

/**
 * class Particle
 * 
 */

class Particle {
private:
    uint id_;
    Position pos_;
    Velocity vel_;
    Velocity vel2_;
    Force force_;
    double charge_;
    NearestNeighboursList NN_list_; //TODO: Check if this is necessary
    NeighbourList N_list_;
    uint no_NN_;
    uint Number_of_Neighbours_;
    Displacement dis_;


public:

    // Setters - Functions that may be used to set the value of the particle's properties

    inline void position(const Position& pos) {
        pos_ = pos;
    }

    inline void displacement(const Displacement& dis) {
        dis_ = dis;
    }

    inline void velocity(const Velocity& vel) {
        vel_ = vel;
    }

    inline void velocity2(const Velocity& vel) {
        vel2_ = vel;
    }

    inline void force(const Force& force) {
        force_ = force;
    }

    inline void charge (double charge) {
	    charge_ = charge;
    }

    inline void nearestNeighboursList(const NearestNeighboursList& NN_list) {
        NN_list_ = NN_list;
    }//TODO: Check if required

    inline void neighbourList(const NeighbourList& N_list) {
        N_list_ = N_list;
    }

    inline void noNearestNeighbours(uint no_NN) {
        no_NN_ = no_NN;
    }

    inline void noNeighbours(uint Number_Of_Neighbours) {
        Number_of_Neighbours_ = Number_Of_Neighbours;
    }

    inline void id(uint id) {
        id_ = id;
    }

    inline void setPositionByCoords(double x, double y, double z) {
        pos_.x = x;
        pos_.y = y;
        pos_.z = z;
    }

    inline void setDisplacementByCoords(double x, double y, double z) {
        dis_.x = x;
        dis_.y = y;
        dis_.z = z;
    }

    inline void setForceByCoords(double x, double y, double z) {
        force_.x = x;
        force_.y = y;
        force_.z = z;
    }

    inline void setVelocityByCoords(double x, double y, double z) {
        vel_.x = x;
        vel_.y = y;
        vel_.z = z;
    }

    inline void setVelocity2ByCoords(double x, double y, double z) {
        vel2_.x = x;
        vel2_.y = y;
        vel2_.z = z;
    }

    // Geters

    inline void getPositionByCoords(double& x, double& y, double& z) const {
        x = pos_.x;
        y = pos_.y;
        z = pos_.z;
    }

    inline void getDisplacementByCoords(double& x, double& y, double& z) const {
        x = dis_.x;
        y = dis_.y;
        z = dis_.z;
    }

    inline void getForceByCoords(double& fx, double& fy, double& fz) const {
        fx = force_.x;
        fy = force_.y;
        fz = force_.z;
    }

    inline void getVelocityByCoords(double& vx, double& vy, double& vz) const {
        vx = vel_.x;
        vy = vel_.y;
        vz = vel_.z;
    }

    inline void getVelocity2ByCoords(double& vx, double& vy, double& vz) const {
        vx = vel2_.x;
        vy = vel2_.y;
        vz = vel2_.z;
    }

    inline Velocity& velocity() {
        return vel_;
    }

    inline Velocity& velocity2() {
        return vel2_;
    }

    inline Position& position() {
        return pos_;
    }

    inline const Position& position()const {
        return pos_;
    }

    inline Displacement& displacement() {
        return dis_;
    }

    inline Force& force() {
        return force_;
    }
	
    inline double& charge() {
	    return charge_;
    }

    const NearestNeighboursList& nearestNeighboursList() const {
        return NN_list_;
    }

    NearestNeighboursList& nearestNeighboursList() {
        return NN_list_;
    }

    const NeighbourList& neighbourList() const {
        return N_list_;
    }

    NeighbourList& neighbourList() {
        return N_list_;
    }

    inline uint noNearestNeighbours() const {
        return no_NN_;
    }

    inline uint noNeighbours() const {
        return Number_of_Neighbours_;
    }

    inline uint id() const {
        return id_;
    }

public:

    // Constructors/Destructors
    Particle();
    virtual ~Particle();

    // Public methods

    inline double distance(Particle p) {
        return norm(pos_.x - p.position().x, pos_.y - p.position().y, pos_.z - p.position().z);
    }

    void addNearestNeighbour(const NearestNeighbour& nn) {
        NN_list_[no_NN_] = nn;
        no_NN_++;
    }

    void addNeighbour(const Neighbour& n) {
        N_list_[Number_of_Neighbours_] = n;
        Number_of_Neighbours_++;
    }

    void clearNNList() {
        no_NN_ = 0;
    }

    void clearNeighbourList() {
        Number_of_Neighbours_ = 0;
    }

    void print() const;

};

#endif // PARTICLE_H
