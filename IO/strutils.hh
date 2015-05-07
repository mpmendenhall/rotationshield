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
using std::string;
#include <vector>
using std::vector;
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>

/// convert a double to a string
string dtos(double d);
/// convert an int to a string
string itos(int i);

/// convert a vector of doubles to a string list
string vtos(const double* st, const double* en, string sep = ",");
/// convert a vector of doubles to a string list
string vtos(const vector<double>& ds,string sep = ",");
/// convert a vector of floats to a string list
string vtos(const float* st, const float* en, string sep = ",");
/// convert a vector of doubles to a string list
string vtos(const vector<float>& ds,string sep = ",");
/// convert a vector of ints to a string list
string vtos(const int* st, const int* en, string sep = ",");
/// convert a vector of ints to a string list
string vtos(const vector<int>& ds,string sep = ",") ;

/// convert a char to a string
string ctos(char c);
/// convert a string to lowercase
string lower(string s);
/// convert a string to uppercase
string upper(string s);
/// replace all of one character in a string with another
string replace(string s, char o, char n);
/// check whether string a begins with string b
bool startsWith(const string& a, const string& b);
/// split a string into substrings on given split characters
vector<string> split(const string& s, const string splitchars = " \t\r\n");
/// join a list of strings into a single string
string join(const vector<string>& ss, const string& sep = " ");
/// strip junk chars off start and end of string
string strip(const string& s, const string stripchars = " \t\r\n");
/// split a string into a vector of doubles
vector<double> sToDoubles(const string& s, const string splitchars = ", \t\r\n");
/// split a string into a vector of floats
vector<float> sToFloats(const string& s, const string splitchars = ", \t\r\n");
/// split a string into a vector of ints
vector<int> sToInts(const string& s, const string splitchars = ", \t\r\n");
/// read in an array from a file
vector< vector<float> > readArray(std::ifstream& fin, unsigned int minitems = 1, const string splitchars = ", \t\r\n");

#endif
