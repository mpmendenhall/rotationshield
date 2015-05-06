/* 
 * SMExcept.hh, part of the RotationShield program
 * Copyright (c) 2007-2014 Michael P. Mendenhall
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/// \file SMExcept.hh \brief exception handling class
#ifndef SMEXCEPT_HH
/// make sure this file is only included once
#define SMEXCEPT_HH

#include <exception>
#include "QFile.hh"

/// exception class for error handling
class SMExcept: public std::exception, public Stringmap {
public:
    /// constructor
    SMExcept(std::string tp);
    /// destructor
    ~SMExcept() throw() {}
    /// display error
    virtual const char* what() const throw();
    /// string for holding error message
    mutable std::string msg;
};

#endif
