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

#include "common.h"

uint tokenizeString(const std::string& i_source,
    		    const std::string& i_seperators,
    		    bool i_discard_empty_tokens,
    		    std::vector<std::string>& o_tokens) 
{
    size_t prev_pos = 0;
    size_t pos = 0;
    size_t number_of_tokens = 0;
    o_tokens.clear();
    pos = i_source.find_first_of(i_seperators, pos);
    while (pos != std::string::npos)
    {
    	std::string token = i_source.substr(prev_pos, pos - prev_pos);
    	if (!i_discard_empty_tokens || token != "")
    	{
    		o_tokens.push_back(i_source.substr(prev_pos, pos - prev_pos));
    		number_of_tokens++;
    	}

    	pos++;
    	prev_pos = pos;
    	pos = i_source.find_first_of(i_seperators, pos);
    }

    if (prev_pos < i_source.length())
    {
    	o_tokens.push_back(i_source.substr(prev_pos));
    	number_of_tokens++;
    }

    return number_of_tokens;
}


