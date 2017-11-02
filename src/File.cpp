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


#include "File.h"

// File Methods
File::File() {}
File::~File() {}

void File::close() {
	fd_.close();
}


// InputFile Methods
InputFile::InputFile() {}
InputFile::~InputFile() {}

InputFile::InputFile(std::string fname) {
	fname_ = fname;
	open();
}

void InputFile::open() {
	fd_.open(fname_.c_str(), std::ios_base::in);
	if (!fd_.good())
		error("", "Input error", "File '" + fname_ + "' does not exist.");
}

std::istream& InputFile::readLine(std::string& line) {
	return std::getline(fd_, line);
}

// OutputFile Methods
OutputFile::OutputFile() {}
OutputFile::~OutputFile() {}

OutputFile::OutputFile(std::string fname) {
	fname_ = fname;
	open();
}

void OutputFile::open() {
	fd_.open(fname_.c_str(), std::ios_base::out);
	if (!fd_.good())
		error("", "Input error", "File " + fname_ + " does not exist.");
}

void OutputFile::writeLine(std::string line) {
	fd_ << line;
}

