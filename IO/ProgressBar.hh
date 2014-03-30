/* 
 * ProgressBar.hh, part of the RotationShield program
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

#ifndef PROGRESSBAR_HH
#define PROGRESSBAR_HH 1

#include <stdio.h>
#include <string>

/// class for printing a progress bar to stdout
class ProgressBar {
public:

	/// constructor
	ProgressBar(unsigned int nt, unsigned int nm=1, bool v=true, const std::string& label="");
	
	/// destructor
	~ProgressBar() { if(verbose) printf("* Done.\n"); }
	
	/// update status at i items completed
	void update(unsigned int i);
	
protected:

	const unsigned int ntotal;	//< total number of items to completion
	const unsigned int nsteps;	//< number of steps to mark = ntotal/nmod
	const unsigned int nmod;	//< number of items per marked step
	
	unsigned int c;				//< number of items completed
	unsigned int s;				//< steps displayed
	const bool verbose;			//< whether to display the progress bar
};

#endif
