#include "ProgressBar.hh"

ProgressBar::ProgressBar(unsigned int nt, unsigned int nm, bool v, const std::string& label): ntotal(nt), nsteps(nt/nm), nmod(nm), c(0), s(0), verbose(v) {
	if(verbose) {
		printf("%s+",label.c_str()); 
		for(unsigned int i=0; i<nsteps; i++)
			printf("-");
		printf("\n|");
		fflush(stdout);
	}
}
	
void ProgressBar::update(unsigned int i) {
	if(i<=c) return;
	c = i;
	
	for(; s<c/nmod; s++) {
		if(verbose) {
			printf("*");
			fflush(stdout);
		}
	}
}
