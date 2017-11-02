/* Object oriented implementation of a File.
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


#ifndef FILE_H
#define FILE_H
#include <fstream>
#include <string>
#include <common.h>

class File {
protected:
    std::string fname_;
    std::fstream fd_;
public:
    File();
    virtual ~File();
    virtual void open() = 0;
    void close();

    std::string filename() {
        return fname_;
    };

    void filename(std::string fname) {
        fname_ = fname;
    };
};

class InputFile : public File {
public:
    InputFile();
    InputFile(std::string fname);
    virtual ~InputFile();
    std::istream& readLine(std::string& line);
    virtual void open();
};

class OutputFile : public File {
public:
    OutputFile();
    OutputFile(std::string fname);
    virtual ~OutputFile();
    void writeLine(std::string line);
    virtual void open();
};

#endif
