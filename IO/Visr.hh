#ifndef VISCONTROLLER_HH
#define VISCONTROLLER_HH 1

#include "Vec.hh"
#include "MiscUtils.hh"

#define WITH_OPENGL 1 

#ifndef WITH_OPENGL

namespace vsr {
	
	class Visr {
	public:
		Visr() {}
		~Visr() {}
		void clearWindow() {}
		void resetViewTransformation() {}
		void startRecording(bool newseg = false) {}
		void stopRecording() {}
		void pause() {}
		void setColor(float r, float g, float b) {}
		void setColor(float r, float g, float b, float a) {}
		void line(vec3 s, vec3 e) {}
		void plane(vec3 o, vec3 dx, vec3 dy) {}
		void quad(float* xyz) {}
		void filledquad(float* xyz) {}
		void dot(vec3 p) {}
		
		void startLines() {}
		void vertex(vec3 v) {}
		void endLines() {}
		
		static Visr* W;
	};
	
}

#else

#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <vector>

namespace vsr {
	
	class Visr {
	public:
		Visr();
		~Visr() {}
		void clearWindow();
		void resetViewTransformation();
		void startRecording(bool newseg = false);
		void stopRecording();
		void pause();
		void setColor(float r, float g, float b);
		void setColor(float r, float g, float b, float a);
		void line(vec3 s, vec3 e);
		void plane(vec3 o, vec3 dx, vec3 dy);
		void quad(float* xyz);
		void filledquad(float* xyz);
		void dot(vec3 p);
		
		void startLines();
		void vertex(vec3 v);
		void endLines();
		
		static Visr* W;
	};
}

#endif

#endif
