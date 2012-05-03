#include "Permutation.hh"

Permutation::Permutation(unsigned int n): data(std::vector<unsigned int>(n)) {
	for(unsigned int i=0; i<n; i++)
		data[i] = i;
}

Permutation& Permutation::nshuffle(int n) {
	for(unsigned int i=0; i<size()/n; i++) 
		for(int j=0; j<n; j++)
			data[j*size()/n+i] = i*n+j;
	return *this;
}

void Permutation::swap(unsigned int a, unsigned int b) {
	std::swap(data[a],data[b]);
}

Permutation Permutation::inverse() const {
	Permutation inv = Permutation(size());
	for(unsigned int i=0; i<size(); i++) inv[data[i]] = i;
	return inv;
}

Permutation& Permutation::invert() {
	std::vector<unsigned int> newdat = std::vector<unsigned int>(size());
	for(unsigned int i=0; i<size(); i++) newdat[data[i]] = i;
	data = newdat;
	return *this;
}

const Permutation Permutation::operator*(const Permutation& p) const {
	assert(p.size() == size());
	Permutation o = Permutation(size());
	for(unsigned int i=0; i<size(); i++) o[i] = p[data[i]];
	return o;
}

int main() {
	return 0;
}