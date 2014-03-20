#ifndef ANGLES_HH
#define ANGLES_HH 1

#include <vector>
#include <set>
#include <cmath>
#include <iostream>

/// normalize an angle to [t0,t0+2*PI)
float normalizeAngle(float a, float theta0 = -M_PI);

/// utility class for math on angular intervals
struct angular_interval {
	/// constructor
	angular_interval(double t0=0, double t1=0, bool nrm = true);
	/// normalize so start is in [0,2pi)
	void normalize();
	/// add a rotating angle
	void add(double a);
	
	double th0;	//< start angle
	double th1;	//< end angle
};

/// interval comparison for ordering
inline bool operator<(const angular_interval& i1, const angular_interval& i2) { return i1.th0==i2.th0 ? i1.th1 < i2.th1 : i1.th0 < i2.th0; }

/// unions of closed angular intervals
class Angular_Interval_Set {
public:
	/// add an interval to the set
	void add_interval(angular_interval i);
	/// remove an interval from the set
	void subtract_interval(angular_interval i);
	/// add an interval to the set
	void add_interval(double a, double b) { add_interval(angular_interval(a,b)); }
	/// remove an interval from the set
	void subtract_interval(double a, double b) { subtract_interval(angular_interval(a,b)); }
	
	/// get list of intervals in set, optionally merging wrap-around
	std::vector<angular_interval> get_intervals(bool merge_wrap = true) const;

protected:
	std::set<angular_interval> endpts;	//< interval endpoints
};

/// string output representation for angular intervals
std::ostream& operator<<(std::ostream& o, const angular_interval& i);
/// string output representation for Angular_Interval_Set
std::ostream& operator<<(std::ostream& o, const Angular_Interval_Set& S);

#endif
