#ifndef VISCONTROLLER_HH
#define VISCONTROLLER_HH 1

#include "Vec.hh"
#include "MiscUtils.hh"
#include <vector>
#include <deque>

namespace vsr {
	
	void initWindow(const std::string& title = "OpenGL Viewer Window");
	void clearWindow();
	void resetViewTransformation();
	void startRecording(bool newseg = false);
	void stopRecording();
	void pause();
	void setColor(float r, float g, float b, float a = 1);
	void line(vec3 s, vec3 e);
	void plane(vec3 o, vec3 dx, vec3 dy);
	void quad(float* xyz);
	void filledquad(float* xyz);
	void dot(vec3 p);
	
	void startLines();
	void vertex(vec3 v);
	void endLines();
	void doGlutLoop();
}

#endif
