#ifndef VISCONTROLLER_HH
#define VISCONTROLLER_HH 1

#include "Vec.hh"
#include "Typedefs.hh"
#include <vector>
#include <deque>

namespace vsr {

	/// initialize visualization window
	void initWindow(const std::string& title = "OpenGL Viewer Window");
	/// enter main drawing loop
	void doGlutLoop();
	
	/// reset view to default
	void resetViewTransformation();
	/// update view window after changes
	void updateViewWindow();
	/// start recording a series of draw commands; newseg=true to erase all previous series
	void startRecording(bool newseg = false);
	/// stop recording a series of draw commands
	void stopRecording();
	/// wait for [ENTER] to be pushed in vis window
	void pause();
	
	/// clear window to blank screen
	void clearWindow();
	/// set background clear color
	void setClearColor(float r, float g, float b, float a=0);
	/// set color for subsequent draws
	void setColor(float r, float g, float b, float a = 1);
	/// draw specified line
	void line(vec3 s, vec3 e);
	/// draw specified plane, centered at o, +/- x and y vectors
	void plane(vec3 o, vec3 dx, vec3 dy);
	/// draw quadrangle outline
	void quad(float* xyz);
	/// draw filled quadrangle
	void filledquad(float* xyz);
	/// draw dot at location
	void dot(vec3 p);
	
	/// start a polygon/series-of-lines
	void startLines();
	/// next vertex in line series
	void vertex(vec3 v);
	/// end series of lines
	void endLines();
	
	/// set the pause flag, cleared when [ENTER] pressed in vis window
	void set_pause();
	/// get current state of pause flag
	bool get_pause();
	/// set the kill flag to end visualization thread
	void set_kill();
}

/// virtual base class for visualizable objects
class Visualizable {
public:
	/// constructor
	Visualizable() {}
	/// destructor
	virtual ~Visualizable() {}
	

	/// visualize "top level"
	void visualize() const { vsr::startRecording(true); vsr::clearWindow(); _visualize(); vsr::stopRecording(); }
		
	/// visualize without clearing screen
	virtual void _visualize() const = 0;
	
	static bool vis_on;
};


#endif
