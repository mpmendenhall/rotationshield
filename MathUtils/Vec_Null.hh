#ifndef VEC_NULL_HH
#define VEC_NULL_HH 1

/// Special 0 case
template<typename T>
class Vec<0,T> {
public:
	
	/// default constructor for zero vector
	Vec() { }
	
	/// print to stdout
	void display(const char* suffix = "\n") const { printf("< >%s",suffix); };
	
	/// dot product with another vector
	T dot(const Vec<0,T> v) const { return T(); }
	/// square magnitude \f$ v \cdot v \f$
	T mag2() const { return T(); }
	/// magnitude \f$ \sqrt{v\cdot v} \f$
	T mag() const { return T(); }
	/// sum of vector elements
	T sum() const { return T(); }
	/// product of vector elements
	T prod() const { return T(); }
	
	/// this vector, normalized to magnitude 1
	Vec<0,T> normalized() const { return Vec<0,T>(); }
	/// project out component parallel to another vector
	Vec<0,T> paraProj(const Vec<0,T> v) const { return Vec<0,T>(); }
	/// project out component orthogonal to another vector
	Vec<0,T> orthoProj(const Vec<0,T> v) const { return Vec<0,T>(); }
	/// angle with another vector
	T angle(const Vec<0,T> v) const { return T(); }
	
	/// mutable element access operator
	T& operator[](unsigned int i) { assert(false); }
	/// immutable element access operator
	T operator[](unsigned int i) const { assert(false); }
	
	/// unary minus operator
	const Vec<0,T> operator-() const { return Vec<0,T>(); }
	
	/// inplace addition
	Vec<0,T>& operator+=(const Vec<0,T>& rhs) { return *this; }
	/// inplace subtraction
	Vec<0,T>& operator-=(const Vec<0,T>& rhs) { return *this; }
	
	/// inplace multiplication
	Vec<0,T>& operator*=(T c) { return *this; }
	/// inplace elementwise multiplication
	Vec<0,T>& operator*=(const Vec<0,T>& other) { return *this; }
	/// inplace division
	Vec<0,T>& operator/=(T c) { return *this; }
	/// inplace elementwise division
	Vec<0,T>& operator/=(const Vec<0,T>& other) { return *this; }
	
	/// addition operator
	const Vec<0,T> operator+(const Vec<0,T>& other) const { return Vec<0,T>(); }
	/// subtraction operator
	const Vec<0,T> operator-(const Vec<0,T>& other) const { return Vec<0,T>(); }
	
	/// multiplication operator
	const Vec<0,T> operator*(T c) const { return Vec<0,T>(); }
	/// elementwise multiplication operator
	const Vec<0,T> operator*(const Vec<0,T>& other) const { return Vec<0,T>(); }
	/// division operator
	const Vec<0,T> operator/(T c) const { return Vec<0,T>(); }
	/// elementwise division operator
	const Vec<0,T> operator/(const Vec<0,T>& other) const { return Vec<0,T>(); }
	
	/// equality operator
	bool operator==(const Vec<0,T>& rhs) const { return true; }
	/// inequality operator
	bool operator!=(const Vec<0,T>& rhs) const { return false; }
	/// dictionary-order comparison operator
	bool operator<(const Vec<0,T>& rhs) const { return false; }
	/// dictionary-order comparison operator
	bool operator<=(const Vec<0,T>& rhs) const { return true; }
	/// dictionary-order comparison operator
	bool operator>(const Vec<0,T>& rhs) const { return false; }
	/// dictionary-order comparison operator
	bool operator>=(const Vec<0,T>& rhs) const { return true; }
};

#endif
