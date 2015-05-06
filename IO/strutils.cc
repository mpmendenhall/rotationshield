/* 
 * strutils.cc, part of the RotationShield program
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

#include "strutils.hh"
#include <stdlib.h>

std::string dtos(double d) {
    char c[16];
    sprintf(c,"%g",d);
    return std::string(c);
}

std::string itos(int i) {
    char c[32];
    sprintf(c,"%i",i);
    return std::string(c);    
}

std::string vtos(const double* st, const double* en, std::string sep) {
    std::string s = "";
    if(st==en)
        return s;
    s = dtos(*st);
    for(const double* it = st+1; it != en; it++)
        s += sep+dtos(*it);
    return s;
}

std::string vtos(const std::vector<double>& ds,std::string sep) { return vtos(&*ds.begin(),&*ds.end(),sep); }

std::string vtos(const float* st, const float* en, std::string sep) {
    std::string s = "";
    if(st==en)
        return s;
    s = dtos(*st);
    for(const float* it = st+1; it != en; it++)
        s += sep+dtos(*it);
    return s;
}

std::string vtos(const std::vector<float>& ds,std::string sep) { return vtos(&*ds.begin(),&*ds.end(),sep); }

std::string vtos(const int* st, const int* en, std::string sep) {
    std::string s = "";
    if(st==en)
        return s;
    s = itos(*st);
    for(const int* it = st+1; it != en; it++)
        s += sep+itos(*it);
    return s;
}

std::string vtos(const std::vector<int>& ds,std::string sep) { return vtos(&*ds.begin(),&*ds.end(),sep); }

std::string ctos(char c) {
    char ch[3];
    sprintf(ch,"%c",c);
    return std::string(ch);        
}

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), (int(*)(int))tolower);
    return s;
}

std::string upper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), (int(*)(int))toupper);
    return s;
}

std::string replace(std::string s, char o, char n) {
    std::string::size_type found = s.find_first_of(o);
    while( found != std::string::npos ) {
        s[found] = n;
        found = s.find_first_of(o,found+1);
    }
    return s;
}

bool startsWith(const std::string& a, const std::string& b) { return a.substr(0,b.size()) == b; }

std::vector<std::string> split(const std::string& s, const std::string splitchars) {
    std::vector<std::string> v;
    size_t p = 0;
    while(p<s.size()) {
        size_t wstart = s.find_first_not_of(splitchars,p);
        if(wstart == std::string::npos)
            break;
        p = s.find_first_of(splitchars,wstart);
        if(p == std::string::npos)
            p = s.size();
        v.push_back(s.substr(wstart,p-wstart));
    }
    return v;
}

std::string join(const std::vector<std::string>& ss, const std::string& sep) {
    std::string s = "";
    if(!ss.size())
        return s;
    s = ss[0];
    for(std::vector<std::string>::const_iterator it = ss.begin()+1; it < ss.end(); it++)
        s += sep + *it;
    return s;
}

std::string strip(const std::string& s, const std::string stripchars) {
    size_t wstart = s.find_first_not_of(stripchars);
    if(wstart == std::string::npos)
        return "";
    size_t wend = s.find_last_not_of(stripchars);
    return s.substr(wstart,wend-wstart+1);
}

std::vector<double> sToDoubles(const std::string& s, const std::string splitchars) {
    std::vector<double> v;
    std::vector<std::string> words = split(s,splitchars);
    for(unsigned int i=0; i<words.size(); i++)
        v.push_back(atof(words[i].c_str()));
    return v;
}

std::vector<float> sToFloats(const std::string& s, const std::string splitchars) {
    std::vector<float> v;
    std::vector<std::string> words = split(s,splitchars);
    for(unsigned int i=0; i<words.size(); i++)
        v.push_back(atof(words[i].c_str()));
    return v;
}

std::vector<int> sToInts(const std::string& s, const std::string splitchars) {
    std::vector<int> v;
    std::vector<std::string> words = split(s,splitchars);
    for(unsigned int i=0; i<words.size(); i++)
        v.push_back(atoi(words[i].c_str()));
    return v;
}

std::vector< std::vector<float> > readArray(std::ifstream& fin, unsigned int minitems, const std::string splitchars) {
    std::vector< std::vector<float> > a;
    std::string s;
    while (fin.good()) {
        std::getline(fin,s);
        std::vector<float> v = sToFloats(s,splitchars);
        if(v.size() >= minitems)
            a.push_back(v);
    }
    return a;
}
