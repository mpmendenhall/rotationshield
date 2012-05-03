/// \file "ScanRange.hh" \brief Class for scanning a variable over an interval
#ifndef SCANRANGE_HH
/// Make sure this header is only loaded once
#define SCANRANGE_HH 1

/// Class for scanning over an interval in a specified number of uniform steps
class ScanRange {
public:
	/// Constructor
	/**	\param start first point in scan
	 \param end last point in scan
	 \param npts number of points in scan */
	ScanRange(double start, double end, int npts): nmax(npts), s(start), e(end), n(0) { }
	
	/// Destructor (nothing to do)
	~ScanRange() {}
	
	/// Step to next point in interval
	double next() { if(nmax>1) return s + (e-s)*double(n++)/double(nmax - 1); else return 0.5*(s+e)*(++n); }
	/// Get current location
	double current() const { if(nmax>1) return s + (e-s)*double(n-1)/double(nmax - 1); else return 0.5*(s+e); }
	/// Check whether there are more points to go after this one
	bool goOn() const { return n<=nmax; }
	/// Return the current step number #n
	int getN() const { return n; }
	/// Print the status of the scan to stdout
	void printStatus() const { printf("%g [%i/%i]",current(),n,nmax); }
	
	const int nmax; //< total number of steps over interval
	const double s; //< Starting point of interval
	const double e; //< Ending point of interval
private:
	int n; //< Step number in [0,#nmax)
};


#endif
