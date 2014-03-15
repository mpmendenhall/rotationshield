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
