#ifndef BINARYOUTPUTOBJECT_HH
#define BINARYOUTPUTOBJECT_HH 1

#include <iostream>
#include <stdio.h>
#include <string>

/// virtual base class with functions for objects that write themselves to file
class BinaryOutputObject {
public:
	/// constructor
	BinaryOutputObject() {}
	/// destructor
	virtual ~BinaryOutputObject() {}
	
	/// write string to file
	static void writeString(const std::string& s, std::ostream& o) {
		o.write(s.c_str(), s.size());
		o.flush();
	}
	/// check whether file matches string
	static bool checkString(const std::string& s, std::istream& i, bool throwIfMismatch = true) {
		std::string s2(s.size(),'x');
		i.read(&s2[0],s2.size());
		if(throwIfMismatch && s2!=s) {
			std::cout << "Mismatched file string '" << s << "' read as '" << s2 << "'!\n";
			throw;
		}
		return s2==s;
	}
};

#endif
