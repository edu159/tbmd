/* Configuration file parser.
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

#include "ConfigFile.h"

ConfigFile::ConfigFile() {
}

ConfigFile::~ConfigFile() {
    close();
}

ConfigFile::ConfigFile(std::string fname) : InputFile(fname) {
    ;
    setOptions();
}

ParamList ConfigFile::getParams() {
    std::string param_line;
    TokenList tokens;
    uint line_no = 0;
    while (readLine(param_line)) {
        line_no++;
        if (param_line != "") {
            tokenizeString(param_line, " \t", true, tokens);
            parseParamLine(tokens, line_no);
        }
    }
    return params_;
}

void ConfigFile::parseParamLine(const TokenList tokens, uint line_no) {
    if (tokens[0].at(0) != '%') {
        if (tokens.size() == 3) {
            if (valid_options_.find(tokens[0]) == valid_options_.end())
                error(CONFIG_FILE_ERR_MSG, PARSE_ERR_MSG, "Token not recognised", line_no);
            if (tokens[1] == "=")
                addParam(tokens[0], tokens[2], valid_options_[tokens[0]], line_no);
            else
                error(CONFIG_FILE_ERR_MSG, PARSE_ERR_MSG, "Line format is PARAM = VALUE", line_no);

        } else
            error(CONFIG_FILE_ERR_MSG, PARSE_ERR_MSG, "Line format is PARAM = VALUE", line_no);
    }
}

void ConfigFile::setOptions() {
    valid_options_["TEMPERATURE"] = DOUBLE;
    valid_options_["CREATE_MOVIE_FILE"] = BOOLEAN;
    valid_options_["NUMBER_OF_STEPS"] = INTEGER;
    valid_options_["THERMOSTAT"] = STRING;
    valid_options_["NORMAL_MODES"] = BOOLEAN;
    valid_options_["OUTPUT_FILE_MODE"] = STRING;
    valid_options_["TIME_STEP"] = INTEGER;
    valid_options_["MD_GAMMA"] = DOUBLE;
    valid_options_["TIME_STEP"] = DOUBLE;
    valid_options_["MOVIE_FRAME_RATE"] = INTEGER;
    valid_options_["OUTPUT_DATA_RATE"] = INTEGER;
    valid_options_["CELL_DIMENSION_X"] = DOUBLE;
    valid_options_["CELL_DIMENSION_Y"] = DOUBLE;
    valid_options_["CELL_DIMENSION_Z"] = DOUBLE;
    valid_options_["PERIODIC_BC"] = BOOLEAN;
    valid_options_["MOVIE_FILE_NAME"] = STRING;
    valid_options_["NMODES_FILE_NAME"] = STRING;
}


void ConfigFile::addParam(const std::string param_name, const std::string param_value, ParamType param_type, uint line_no) {
    Param p;
    p.key = param_name;
    switch (param_type) {
        case DOUBLE:
            if (check_double_value(param_value))
                p.double_value = atof(param_value.c_str());
            else
                error(CONFIG_FILE_ERR_MSG, PARSE_ERR_MSG, "Double value not well formatted", line_no);
            break;
        case INTEGER:
            if (check_integer_value(param_value))
                p.integer_value = atoi(param_value.c_str());
            else
                error(CONFIG_FILE_ERR_MSG, PARSE_ERR_MSG, "Integer value not well formatted", line_no);
            break;
        case STRING:
            p.string_value = param_value;
            break;
        case BOOLEAN:
            if (check_boolean_value(param_value))
                (param_value == "true") ? (p.boolean_value = true) : (p.boolean_value = false);
            else
                error(CONFIG_FILE_ERR_MSG, PARSE_ERR_MSG, "Boolean value can only be true or false", line_no);
            break;
    }
    params_.push_back(p);
}
