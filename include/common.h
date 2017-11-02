/* Inlined functions, constants and structures used in various files.
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

#ifndef COMMON_H
#define COMMON_H

// Common includes
#include <cmath>
#include <armadillo>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>


// Common error messages
const std::string PARSE_ERR_MSG = "Parse error";

// Common functions

inline bool check_double_value(std::string value) {
    if ((std::string::npos == value.find_first_not_of("-.1234567890")) &&
            (std::count(value.begin(), value.end(), '.') <= 1))
        return true;
    return false;
}

inline bool check_integer_value(std::string value) {
    if (std::string::npos == value.find_first_not_of("1234567890"))
        return true;
    return false;
}

inline bool check_boolean_value(std::string value) {
    if ((value == "true") | (value == "false"))
        return true;
    return false;
}

inline void error(std::string error_location, std::string error_type, std::string error_msg, uint line_no = 0) {
    std::cerr << "[*]" << error_location << ":" << std::endl;
    std::cerr << "  " << error_type << ": " << error_msg;
    if (line_no != 0)
        std::cout << " -- line " << line_no;
    std::cout << std::endl;
    exit(EXIT_FAILURE);
}

uint tokenizeString(const std::string& i_source,
        const std::string& i_seperators,
        bool i_discard_empty_tokens,
        std::vector<std::string>& o_tokens);

enum OutputModeType {BASIC, DETAILED, DEBUG, LAST_STEP};
enum ThermostatType {VV_DAMPED, OVRVO};

// Constants
const double C_MASS = 1243.71;
const double MD_COEFF = 1.0 / (2.0 * C_MASS);
const double KB = 8.6173e-5;
const double PI = 3.14159265359;
const double RC = 2.60;
const double RL = RC + 1;
const double TOL = 1e-15;
const double DEFAULT_GAMMA = 0.1;
const uint DEFAULT_FRAME_RATE = 10;
const uint DEFAULT_DATA_RATE = 10;
const double DEFAULT_TARGET_TEMPERATURE = 0.0;
const double DELTA = 1e-4;
const double Fn = 2.0; //n
const double Fnc = 6.5; // nc
const double Bs = -0.0126200997391; // Bs=-n/(r_c^n_c)
const double Vspsi = 4.7;
const double Vsssi = -5.0;
const double Vppsi = 5.5;
const double Vpppi = -1.55;
const std::string PATH = "/home/tsmuser/Dropbox/graphs4gpp/"; // output.dat and movie.xyz will be here
const ThermostatType DEFAULT_THERMOSTAT = OVRVO; //should later include Anderson, Langevin and maybe Nose-Hoover 
const double MIN_CELL_SIZE = 3*RL;
const double DEFAULT_TIME_STEP = 0.1;
const uint MAX_NEIGHBOURS = 50;
const uint MAX_NNEIGHBOURS = 40;

// Type definitions
typedef unsigned int uint;
typedef arma::mat Matrix;
typedef arma::vec Vector;
typedef double (*fun1D_ptr)(double);


// Common functions

inline double norm(double x, double y, double z) {
    return sqrt(x * x + y * y + z * z);
}

inline double finite_difference(double x, double h, fun1D_ptr func) {
    return (-func(x + (2.0 * h)) + (8.0 * func(x + h)) - (8.0 * func(x - h)) + func(x - (2.0 * h))) / (12.0 * h);
}

struct Vector3D {
    double x;
    double y;
    double z;

    Vector3D& operator+=(const Vector3D& rhs) {
        this->x += rhs.x;
        this->y += rhs.y;
        this->z += rhs.z;
        return *this;
    }

    Vector3D& operator-=(const Vector3D& rhs) {
        this->x -= rhs.x;
        this->y -= rhs.y;
        this->z -= rhs.z;
        return *this;
    }

    Vector3D& operator-() {
        this->x = -this->x;
        this->y = -this->y;
        this->z = -this->z;
        return *this;
    }

    Vector3D& operator=(const Vector3D& rhs) {
        this->x = rhs.x;
        this->y = rhs.y;
        this->z = rhs.z;
        return *this;
    }

    Vector3D& operator=(const double rhs) {
        this->x = rhs;
        this->y = rhs;
        this->z = rhs;
        return *this;
    }

    Vector3D& operator*=(const double rhs) {
        this->x *= rhs;
        this->y *= rhs;
        this->z *= rhs;
        return *this;
    }

    Vector3D& operator+=(const double rhs) {
        this->x += rhs;
        this->y += rhs;
        this->z += rhs;
        return *this;
    }

    Vector3D& operator/=(const double rhs) {
        this->x /= rhs;
        this->y /= rhs;
        this->z /= rhs;
        return *this;
    }
};

inline Vector3D operator+(Vector3D lhs, const Vector3D& rhs) {
    lhs += rhs;
    return lhs;
}

inline Vector3D operator-(Vector3D lhs, const Vector3D& rhs) {
    lhs -= rhs;
    return lhs;
}

inline Vector3D operator*(Vector3D lhs, const double rhs) {
    lhs *= rhs;
    return lhs;
}

inline Vector3D operator*(const double lhs, Vector3D rhs) {
    rhs *= lhs;
    return rhs;
}

inline Vector3D operator/(Vector3D lhs, const double rhs) {
    lhs /= rhs;
    return lhs;
}

inline Vector3D operator/(const double lhs, Vector3D rhs) {
    rhs /= lhs;
    return rhs;
}
typedef Vector3D Force;
typedef Vector3D Velocity;
typedef Vector3D Position;
typedef Vector3D Displacement;
#endif //COMMON_H
