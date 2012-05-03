#ifndef PROGRESSBAR_HH
#define PROGRESSBAR_HH 1

#include <stdio.h>
#include <string>

class ProgressBar {
public:
	ProgressBar(unsigned int nt, unsigned int nm, bool v, const std::string& label=""): ntotal(nt), nsteps(nt/nm), nmod(nm), c(0), s(0), verbose(v) {
		if(verbose) {
			printf("%s+",label.c_str()); 
			for(unsigned int i=0; i<nsteps; i++)
				printf("-");
			printf("\n|");
			fflush(stdout);
		}
	}
	
	~ProgressBar() { if(verbose) printf("* Done.\n"); }
	
	void update(unsigned int i) {
		if(i<=c) return;
		c = i;
		
		for(; s<c/nmod; s++) {
			if(verbose) {
				printf("*");
				fflush(stdout);
			}
		}
		
	}
	
protected:
	const unsigned int ntotal;
	const unsigned int nsteps;
	const unsigned int nmod;
	
	unsigned int c;
	unsigned int s;
	const bool verbose;
};

#endif
