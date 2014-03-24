#include "CMatrix.hh"

#include <stdlib.h>
#include <iomanip>
#include <algorithm>

std::vector<cmatrix_fft> CMatrix::ffters = std::vector<cmatrix_fft>();

// Constructor and Destructor ------______------______-------_______------______------

CMatrix::CMatrix(unsigned int ncyc):
ncycles(ncyc), data(), kdata(), has_realspace(false), has_kspace(false) { }

void CMatrix::writeToFile(std::ostream& o) const {
	o.write((char*)&ncycles,sizeof(int));
	o.write((char*)&has_realspace,sizeof(bool));
	o.write((char*)&has_kspace,sizeof(bool));
	if(has_realspace)
		o.write((char*)&data.front(),sizeof(double)*ncycles);
	if(has_kspace)
		o.write((char*)&kdata.front(),sizeof( complex<double> )*(ncycles/2+1));
}

CMatrix CMatrix::readFromFile(std::istream& s) {
	int ncyc;
	bool hrd, hkd;
	s.read((char*)&ncyc,sizeof(int));
	s.read((char*)&hrd,sizeof(bool));
	s.read((char*)&hkd,sizeof(bool));
	CMatrix foo = CMatrix(ncyc);
	if(hrd) {
		foo.alloc_data();
		s.read((char*)&foo.data.front(),sizeof(double)*ncyc);
	}
	if(hkd) {
		foo.alloc_kdata();
		s.read((char*)&foo.kdata.front(),sizeof( complex<double> )*(ncyc/2+1));
	}
	return foo;
}

// Special Matrices ------______------______-------_______------______------

CMatrix CMatrix::identity(unsigned int nc) {
	CMatrix m = CMatrix(nc);
	m[0] = 1.0;
	return m;
}

void CMatrix::zero() {
	has_kspace = true;
	has_realspace = true;
	data = std::vector<double>(ncycles);
	kdata = std::vector< complex<double> >(ncycles/2+1);
}

CMatrix CMatrix::ramp(unsigned int nc, double r0) {
	CMatrix m = CMatrix(nc);
	for(unsigned int i=0; i<nc; i++)
		m[i] = i+r0;
	return m;
}

CMatrix CMatrix::random(unsigned int nc) {
	CMatrix m = CMatrix(nc);
	for(unsigned int i=0; i<nc; i++)
		m[i] = double(rand())/double(RAND_MAX);
	return m;
}

// Data Providers ------______------______-------_______------______------

void CMatrix::alloc_data() const {
	data = std::vector<double>(ncycles);
	has_realspace = true;
}

void CMatrix::alloc_kdata() const {
	kdata = std::vector< complex<double> >(ncycles);
	has_kspace = true;
}


void CMatrix::calculateRealData() const {
	cmatrix_fft& f = get_ffter();
	for(unsigned int i=0; i<ncycles/2+1; i++)
		f.kspace[i] = kdata[i];
	fftw_execute(f.reverseplan);
	alloc_data();
	for(unsigned int n=0; n<ncycles; n++) data[n] = double(f.realspace[n])/double(ncycles);	
}


double CMatrix::operator[](unsigned int i) const {
	assert(i<ncycles);
	if(has_realspace)
		return data[i];
	if(!has_kspace) {
		alloc_data();
		return data[i];
	}
	calculateRealData();
	return data[i];
}

double& CMatrix::operator[](unsigned int i) {
	assert(i<ncycles);
	if(has_realspace)
		return data[i];
	if(!has_kspace) {
		alloc_data();
		return data[i];
	}
	calculateRealData();
	has_kspace = false; // messing with real space data invalidates kspace data
	return data[i];	
}


void CMatrix::calculateKData() const {
	cmatrix_fft& f = get_ffter();
	for(unsigned int n=0; n<ncycles; n++) f.realspace[n] = data[n];
	fftw_execute(f.forwardplan);
	alloc_kdata();
	kdata.assign(f.kspace,f.kspace+ncycles/2+1);
	has_kspace = true;
}

std::vector< complex<double> >& CMatrix::getKData() {
	if(has_kspace)
		return kdata;
	if(!has_realspace) {
		alloc_kdata();
		return kdata;
	}
	calculateKData();
	has_realspace = false;
	return kdata;
}

const std::vector< complex<double> >& CMatrix::getKData() const {
	if(has_kspace)
		return kdata;
	if(!has_realspace) {
		alloc_kdata();
		return kdata;
	}
	calculateKData();
	return kdata;
}

std::vector<double>& CMatrix::getRealData() {
	if(has_realspace)
		return data;
	if(!has_kspace) {
		alloc_data();
		return data;
	}
	calculateRealData();
	has_kspace = false;
	return data;
}

const std::vector<double>& CMatrix::getRealData() const {
	if(has_realspace)
		return data;
	if(!has_kspace) {
		alloc_data();
		return data;
	}
	calculateRealData();
	return data;
}

// Matrix properties ------______------______-------_______------______------

double CMatrix::norm_L2() const {
	const std::vector< complex<double> >& v = getKData();
	std::vector<double> vn;
	for(std::vector< complex<double> >::const_iterator it = v.begin(); it < v.end(); it++)
		vn.push_back(abs(*it));
	return *std::max_element(vn.begin(),vn.end());
}

double CMatrix::det() const {
	const std::vector< complex<double> >& v = getKData();
	double d = v.begin()->real();
	for(unsigned int i=1; i<v.size(); i++)
		d *= norm(v[i]);
	if(!(ncycles%2)) d /= v.back().real();
	return d;
}

double CMatrix::trace() const {
	if(has_realspace) return ncycles*data[0];
	else if(has_kspace) {
		double s = kdata.begin()->real();
		for(std::vector< complex<double> >::const_iterator it = kdata.begin()+1; it < kdata.end(); it++)
			s += it->real()*2;
		if(!(ncycles%2)) s -= kdata.back().real();
		return s;
	} else return 0;
}

// "GUI" ------______------______-------_______------______------


void CMatrix::printRow(int r) const {
	printf("| ");
	for(unsigned int i=0; i<ncycles; i++) printf("%.3g ",(*this)[(i+(ncycles-r))%ncycles]);
	printf("|");
}

void CMatrix::display() const {
	std::cout << "CMatrix " << ncycles << " " << has_realspace << " " << has_kspace << std::endl;
	for(int r=0; r<ncycles; r++) {
		printRow(r);
		printf("\n");
	}
	std::cout << "CMatrix " << ncycles << " " << has_realspace << " " << has_kspace << std::endl;
}

void CMatrix::displayK() const {
	std::cout << "{ ";
	for(int i=0; i<ncycles/2+1; i++) std::cout << getKData()[i] << " ";
	std::cout << "}" << std::endl;
}

std::ostream& operator<<(std::ostream& o, const CMatrix& m) {
	for(unsigned int r=0; r<m.ncycles; r++) {
		o << "| ";
		for(int c=0; c<m.ncycles; c++)
			o << m[(c+(m.ncycles-r))%m.ncycles] << " ";
		o << "|\n";
	}
	return o;
}


// Matrix Ops ------______------______-------_______------______------


CMatrix& CMatrix::operator+=(const CMatrix& m) {

	assert(ncycles == m.ncycles);
	bool kdata_ok, realdata_ok;
	kdata_ok = realdata_ok = false;
	
	if(!has_realspace && !has_kspace) {
		*this = m;
		return *this;
	}
	
	if(has_kspace && (m.has_kspace || !has_realspace)) {
		const std::vector< complex<double> >& kd = m.getKData();
		for(unsigned int i=0; i<ncycles/2+1; i++)
			kdata[i] += kd[i];
		kdata_ok = true;
		has_kspace = true;
	}

	if(has_realspace && (m.has_realspace || !has_kspace)) {
		for(unsigned int i=0; i<ncycles; i++)
			data[i] += m[i];
		realdata_ok = true;
	}


	
	has_realspace = realdata_ok;
	has_kspace = kdata_ok;
	
	return *this;
}

const CMatrix CMatrix::operator+(const CMatrix& m) const {
	CMatrix r = *this;
	r += m;
	return r;
}

CMatrix& CMatrix::operator-=(const CMatrix& m) {
	assert(ncycles == m.ncycles);
	bool kdata_ok, realdata_ok;
	kdata_ok = realdata_ok = false;
	
	if(!has_realspace && !has_kspace) {
		*this = -m;
		return *this;
	}
	
	if(has_kspace && (m.has_kspace || !has_realspace)) {
		const std::vector< complex<double> >& kd = m.getKData();
		for(unsigned int i=0; i<ncycles/2+1; i++)
			kdata[i] -= kd[i];
		kdata_ok = true;
		has_kspace = true;
	}

	if(has_realspace && (m.has_realspace || !has_kspace)) {
		for(unsigned int i=0; i<ncycles; i++)
			data[i] -= m[i];
		realdata_ok = true;
	}
	
	
	has_realspace = realdata_ok;
	has_kspace = kdata_ok;
	
	return *this;
}

const CMatrix CMatrix::operator-(const CMatrix& m) const {
	CMatrix r = *this;
	r -= m;
	return r;
}

CMatrix& CMatrix::operator*=(double c) {
	
	if(has_realspace)
		for(unsigned int i=0; i<ncycles; i++) data[i] *= c;
	
	if(has_kspace)
		for(unsigned int i=0; i<ncycles/2+1; i++)
			kdata[i] *= c;
	
	return *this;
}

const CMatrix CMatrix::operator*(double c) const {
	CMatrix r = *this;
	r *= c;
	return r;
}

const CMatrix CMatrix::operator-() const {
	CMatrix r = *this;
	r *= -1.0;
	return r;
}

CMatrix& CMatrix::operator*=(const CMatrix& m) {
	getKData();
	const std::vector< complex<double> >& mkd = m.getKData();
	has_realspace = false;
	for(unsigned int i=0; i<ncycles/2+1; i++)
		kdata[i] *= mkd[i];
	return *this;
}

const CMatrix CMatrix::operator*(const CMatrix& m) const {
	CMatrix r = *this;
	r *= m;
	return r;
}


const VarVec<double> CMatrix::operator*(const VarVec<double>& v) const {
	
	const std::vector< complex<double> >& kd = getKData();
	
	cmatrix_fft& f = get_ffter();
	f.realspace[0] = v[0];
	for(unsigned int i=1; i<ncycles; i++) f.realspace[i] = v[ncycles - i];
	fftw_execute(f.forwardplan);
	
	for(unsigned int i=0; i<ncycles/2+1; i++)
		f.kspace[i] *= kd[i];
	
	fftw_execute(f.reverseplan);
	
	VarVec<double> out = VarVec<double>(ncycles);
	out[0] = f.realspace[0];
	for(unsigned int i=1; i<ncycles; i++) out[i] = f.realspace[ncycles-i];
	out /= double(ncycles);
	return out;
}

CMatrix& CMatrix::invert() {
	std::vector< complex<double> >& kd = getKData();
	for(unsigned int i=0; i<ncycles/2+1; i++)
		kd[i] = 1./kd[i];
	return *this;
}

const CMatrix CMatrix::inverse() const {
	CMatrix I = *this;
	I.invert();
	return I;
}

CMatrix CMatrix::cmatrixFromColumn(int ncyc, double* coldat) {
	CMatrix m = CMatrix(ncyc);
	m[0] = coldat[0];
	for(int n=1; n<ncyc; n++) m[n] = coldat[ncyc-n];
	return m;
}

const CMatrix CMatrix::transpose() const {
	CMatrix m = CMatrix(ncycles);
	std::vector<double>& md = m.getRealData(); 
	md[0] = (*this)[0];
	for(int n=1; n<ncycles; n++) md[n] = (*this)[ncycles-n];
	return m;
}

// Inverter handler ------______------______-------_______------______------

void CMatrix::add_inverter() const {
	cmatrix_fft NI;
	ffters.push_back(NI);
	cmatrix_fft& newInverter = ffters.back();
	
	newInverter.ncyc = ncycles;
	newInverter.realspace = new double[ncycles];
	newInverter.kspace = new  complex<double> [ncycles/2+1];
	
	FILE* fin = fopen("fftw_wisdom","r");
	if(fin) {
		fftw_import_wisdom_from_file(fin);
		fclose(fin);
	}
	
	newInverter.forwardplan = fftw_plan_dft_r2c_1d(ncycles,
												   newInverter.realspace,
												   (fftw_complex*)newInverter.kspace,
												   FFTW_EXHAUSTIVE);
	newInverter.reverseplan = fftw_plan_dft_c2r_1d(ncycles,
												   (fftw_complex*)newInverter.kspace,
												   newInverter.realspace,
												   FFTW_EXHAUSTIVE);
	
	FILE* fout = fopen("fftw_wisdom","w");
	if(fout) {
		fftw_export_wisdom_to_file(fout);
		fclose(fout);
	}
}

cmatrix_fft& CMatrix::get_ffter() const {
	std::vector<cmatrix_fft>::iterator it;
	for(it = ffters.begin(); it != ffters.end(); it++)
		if(it->ncyc == ncycles) return *it;
	add_inverter();
	return ffters[ffters.size()-1];
}
