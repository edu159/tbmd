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


#ifndef INPUT_DATA_FILE_H
#define INPUT_DATA_FILE_H
#include "common.h"
#include "File.h"
#include "ParticleSystem.h"

// Errors
const std::string DATA_FILE_ERR_MSG = "Data file";

typedef std::vector<std::string> TokenList;

class InputDataFile : public InputFile {
    //void parseDataLine(const TokenList token_list, uint line_no);
public:
    InputDataFile();
    virtual ~InputDataFile();
    InputDataFile(std::string fname);
    void configureSystem(ParticleSystem* p_system);

};

#endif 
