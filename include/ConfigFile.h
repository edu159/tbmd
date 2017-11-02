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


#ifndef CONFIG_FILE_H
#define CONFIG_FILE_H
#include "common.h"
#include "File.h"

struct Param {
    std::string key;
    double double_value;
    std::string string_value;
    uint integer_value;
    bool boolean_value;
};

enum ParamType {
    INTEGER, DOUBLE, STRING, BOOLEAN
};

// Errors
const std::string CONFIG_FILE_ERR_MSG = "Configuration file";

typedef std::vector<Param> ParamList;
typedef std::vector<std::string> TokenList;
typedef std::map<std::string, ParamType> OptionList;

class ConfigFile : public InputFile {
private:
    ParamList params_;
    OptionList valid_options_;

private:
    void parseParamLine(const TokenList token_list, uint line_no);
    void setOptions();
    void addParam(const std::string param_name, const std::string param_value, ParamType param_type, uint line_no);
public:
    ConfigFile();
    virtual ~ConfigFile();
    ConfigFile(std::string fname);
    ParamList getParams();
};


#endif 
