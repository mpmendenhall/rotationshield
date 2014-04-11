/* 
 * strutils.hh, part of the RotationShield program
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

/// \file strutils.hh \brief String manipulation utilities
#ifndef STRUTILS_HH
/// Make sure this header is only loaded once
#define STRUTILS_HH

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>

/// convert a double to a string
std::string dtos(double d);
/// convert an int to a string
std::string itos(int i);

/// convert a vector of doubles to a string list
std::string vtos(const double* st, const double* en, std::string sep = ",");
/// convert a vector of doubles to a string list
std::string vtos(const std::vector<double>& ds,std::string sep = ",");
/// convert a vector of floats to a string list
std::string vtos(const float* st, const float* en, std::string sep = ",");
/// convert a vector of doubles to a string list
std::string vtos(const std::vector<float>& ds,std::string sep = ",");
/// convert a vector of ints to a string list
std::string vtos(const int* st, const int* en, std::string sep = ",");
/// convert a vector of ints to a string list
std::string vtos(const std::vector<int>& ds,std::string sep = ",") ;

/// convert a char to a string
std::string ctos(char c);
/// convert a string to lowercase
std::string lower(std::string s);
/// convert a string to uppercase
std::string upper(std::string s);
/// replace all of one character in a string with another
std::string replace(std::string s, char o, char n);
/// check whether string a begins with string b
bool startsWith(const std::string& a, const std::string& b);
/// split a string into substrings on given split characters
std::vector<std::string> split(const std::string& s, const std::string splitchars = " \t\r\n");
/// join a list of strings into a single string
std::string join(const std::vector<std::string>& ss, const std::string& sep = " ");
/// strip junk chars off start and end of string
std::string strip(const std::string& s, const std::string stripchars = " \t\r\n");
/// split a string into a vector of doubles
std::vector<double> sToDoubles(const std::string& s, const std::string splitchars = ", \t\r\n");
/// split a string into a vector of floats
std::vector<float> sToFloats(const std::string& s, const std::string splitchars = ", \t\r\n");
/// split a string into a vector of ints
std::vector<int> sToInts(const std::string& s, const std::string splitchars = ", \t\r\n");
/// read in an array from a file
std::vector< std::vector<float> > readArray(std::ifstream& fin, unsigned int minitems = 1, const std::string splitchars = ", \t\r\n");

#endif
