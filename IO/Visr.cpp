#include "Visr.hh"
#include <string>

#ifdef WITH_OPENGL
#include <GL/freeglut.h>

namespace vsr {
	
	void initWindow(const std::string& windowTitle = "OpenGL Viewer Window");
	void* doGlutLoop( void *vptr_args );
	void redrawDisplay();
	
	void resetViewTransformation();
	void reshapeWindow(int width, int height);
	void keypress(unsigned char key, int x, int y);
	void specialKeypress(int key, int x, int y);
	void startMouseTracking(int button, int state, int x, int y);
	void mouseTrackingAction(int x, int y);
	
	pthread_mutex_t displayLock;
	GLuint displayList;
	std::vector<GLuint> displaySegs;
	int clickx0,clicky0;
	int modifier;
	float winwidth, winheight;
	float viewrange;
	float xrot,yrot,zrot;
	float xtrans, ytrans;
	float ar;
	bool pause_display;
	
	Visr::Visr() {
		pause_display = false;
		pthread_mutexattr_t displayLockAttr;
		pthread_mutexattr_init(&displayLockAttr);
		pthread_mutexattr_settype(&displayLockAttr, PTHREAD_MUTEX_RECURSIVE);
		pthread_mutex_init(&displayLock,&displayLockAttr);
		initWindow("OpenGL Viewer Window");
		resetViewTransformation();
		startRecording();
		printf("Drawing initial teapot...\n");
		clearWindow();
		glColor3f(0.0, 0.0, 1.0);
		glutWireTeapot(0.5);
		stopRecording();
	}
	
	void Visr::setColor(float r, float g, float b) {
		glColor3f(r,g,b);
	}
	
	void Visr::setColor(float r, float g, float b, float a) {
		glColor4f(r,g,b,a);
	}
	
	void Visr::clearWindow() {
		glClearColor(0.0, 0.0, 0.0, 0.0);
		glClearDepth(100.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	}
	
	void Visr::pause() {
		pause_display = true;
		while(pause_display) usleep(50000);
	}
	
	void Visr::resetViewTransformation() {	
		pthread_mutex_lock(&displayLock);
		printf("Re-setting view...\n");
		viewrange = 1.0;
		xtrans = ytrans = 0;
		glLineWidth(1.5/viewrange);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslatef(0,0,1.0*viewrange);
		glMatrixMode(GL_PROJECTION);	
		glLoadIdentity();
		glOrtho(-viewrange*ar, viewrange*ar, -viewrange, viewrange, -1000,1000);
		pthread_mutex_unlock(&displayLock);
	}
	
	void Visr::startRecording(bool newseg) {
		printf("Recording GL commands..."); fflush(stdout);
		pthread_mutex_lock(&displayLock);
		printf(" startRecording Lock acquired.\n");
		glFlush();
		glFinish();
		if(!newseg) {
			if(displaySegs.size() && glIsList(displaySegs.back()))
				glDeleteLists(displaySegs.back(),1);
			else
				displaySegs.push_back(glGenLists(1));
		} else
			displaySegs.push_back(glGenLists(1));
		glNewList(displaySegs.back(), GL_COMPILE);
	}
	
	void Visr::stopRecording() {
		glEndList();
		glutPostRedisplay();
		glFlush();
		glFinish();
		pthread_mutex_unlock(&displayLock);
		printf("Recording complete.\n");
	}
	
	void Visr::line(vec3 s, vec3 e) {
		glBegin(GL_LINES);
		glVertex3f(s[0],s[1],s[2]);
		glVertex3f(e[0],e[1],e[2]);
		glEnd();
	}
	
	void Visr::startLines() {
		glBegin(GL_LINE_STRIP);
	}
	
	void Visr::vertex(vec3 v) {
		glVertex3f(v[0],v[1],v[2]);
	}
	
	void Visr::endLines() {
		glEnd();
	}
	
	void Visr::plane(vec3 o, vec3 dx, vec3 dy) {
		glBegin(GL_QUADS);
		glVertex3f(o[0]+0.5*dx[0]+0.5*dy[0],o[1]+0.5*dx[1]+0.5*dy[1],o[2]+0.5*dx[2]+0.5*dy[2]);
		glVertex3f(o[0]-0.5*dx[0]+0.5*dy[0],o[1]-0.5*dx[1]+0.5*dy[1],o[2]-0.5*dx[2]+0.5*dy[2]);
		glVertex3f(o[0]-0.5*dx[0]-0.5*dy[0],o[1]-0.5*dx[1]-0.5*dy[1],o[2]-0.5*dx[2]-0.5*dy[2]);
		glVertex3f(o[0]+0.5*dx[0]-0.5*dy[0],o[1]+0.5*dx[1]-0.5*dy[1],o[2]+0.5*dx[2]-0.5*dy[2]);
		glEnd();
	}
	
	void Visr::quad(float* xyz)
	{
		glBegin(GL_LINE_LOOP);
		glVertex3f(xyz[0],xyz[1],xyz[2]);
		glVertex3f(xyz[3],xyz[4],xyz[5]);
		glVertex3f(xyz[9],xyz[10],xyz[11]);
		glVertex3f(xyz[6],xyz[7],xyz[8]);
		glEnd();
	}
	
	void Visr::filledquad(float* xyz)
	{
		glBegin(GL_QUADS);
		glVertex3f(xyz[0],xyz[1],xyz[2]);
		glVertex3f(xyz[3],xyz[4],xyz[5]);
		glVertex3f(xyz[9],xyz[10],xyz[11]);
		glVertex3f(xyz[6],xyz[7],xyz[8]);
		glEnd();
	}
	
	void Visr::dot(vec3 p)
	{
		mdouble l = viewrange/100.0;
		
		glBegin(GL_QUADS);
		glVertex3f(p[0]+l,p[1],p[2]);
		glVertex3f(p[0],p[1]+l,p[2]);
		glVertex3f(p[0]-l,p[1],p[2]);
		glVertex3f(p[0],p[1]-l,p[2]);
		glVertex3f(p[0]+l,p[1],p[2]);
		glVertex3f(p[0],p[1],p[2]+l);
		glVertex3f(p[0]-l,p[1],p[2]);
		glVertex3f(p[0],p[1],p[2]-l);
		glVertex3f(p[0],p[1]+l,p[2]);
		glVertex3f(p[0],p[1],p[2]+l);
		glVertex3f(p[0],p[1]-l,p[2]);
		glVertex3f(p[0],p[1],p[2]-l);
		glEnd();
	}
	
	
	void* doGlutLoop( void *vptr_args ) {
		int c = 0;
		printf("Starting GLUT main loop...\n");
		usleep(1000000);
		while(1) {
			if(pthread_mutex_trylock(&displayLock)) {
				usleep(10000);
				continue;
			}
			glutMainLoopEvent();
			c = (c+1)%50;
			if(!c)
				glutPostRedisplay();
			pthread_mutex_unlock(&displayLock);
			usleep(1000);
		}
		return 0x0;
	}
	
	
	void initWindow(const std::string& windowTitle) {
		
		printf("Initializing OpenGL visualization window...\n");
		
		int a = 1;
		char programname[] = "glviewer";
		char* pnameptr = programname;
		glutInit(&a,&pnameptr);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
		glutInitWindowSize(600, 600);
		glutInitWindowPosition(100, 100);
		glutCreateWindow(windowTitle.c_str());
		
		ar = 1.0;
		
		glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_ONE_MINUS_DST_ALPHA,GL_DST_ALPHA);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_FOG);
		float fadecolor[4] = {0,0,0,1.0};
		glFogfv(GL_FOG_COLOR,fadecolor);
		glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		
		glutDisplayFunc(&redrawDisplay);
		glutMouseFunc(&startMouseTracking);
		glutMotionFunc(&mouseTrackingAction);
		glutReshapeFunc(&reshapeWindow);
		glutKeyboardFunc(&keypress);
		glutSpecialFunc(&specialKeypress);
		//glutIdleFunc(&redrawIfUnlocked);
		
		printf("Launching visualization thread...\n");
		pthread_mutex_unlock(&displayLock);
		pthread_t thread;
		pthread_create( &thread, NULL, &doGlutLoop, NULL );
		
		printf("Window init done.\n");
	}
	
	void reshapeWindow(int width, int height) {	
		
		printf("Resizing view window to %ix%i...",width,height); fflush(stdout);
		while(pthread_mutex_trylock(&displayLock)) {
			printf("."); fflush(stdout);
			usleep(100000);
		}
		printf(" reshapeWindow Lock acquired.\n");
		
		glViewport(0,0,width,height);
		winwidth = width; winheight = height;
		ar = float(width)/float(height);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(-viewrange*ar+xtrans, viewrange*ar+xtrans, -viewrange+ytrans, viewrange+ytrans, 1000, -1000);
		
		glFlush();
		glFinish();
		
		pthread_mutex_unlock(&displayLock);
		
		printf("Resize Done.\n");
	}
	
	
	
	void redrawDisplay() {
		glCallLists(displaySegs.size(),GL_UNSIGNED_INT,&displaySegs.front());
		glutSwapBuffers();
		glFlush();
		glFinish();
	}
	
	
	void keypress(unsigned char key, int x, int y) {
		if(key == 32 || key == 13) // spacebar or return
			pause_display = false;
		//if(key == 27) // escape
		//	resetViewTransformation();
	}
	
	void specialKeypress(int key, int x, int y) {
		
	}
	
	void startMouseTracking(int button, int state, int x, int y) {
		modifier = glutGetModifiers();
		if(state == GLUT_DOWN) {
			clickx0 = x;
			clicky0 = y;
		}
	}
	
	void mouseTrackingAction(int x, int y) {
		
		if(modifier == GLUT_ACTIVE_SHIFT) {
			float s = (1 - 0.005*(x-clickx0));
			if( (viewrange > 1.0e-2 || s > 1.0) && (viewrange < 1.0e3 || s < 1.0) ) viewrange *= s;
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			glOrtho(-viewrange*ar+xtrans, viewrange*ar+xtrans, -viewrange+ytrans, viewrange+ytrans, 1000, -1000);
			glLineWidth(1.5/viewrange);
			
		} else if(modifier == GLUT_ACTIVE_CTRL) {
			xtrans -= ar*2.0*(x-clickx0)*viewrange/winwidth;
			ytrans += 2.0*(y-clicky0)*viewrange/winheight;
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			glOrtho(-viewrange*ar+xtrans, viewrange*ar+xtrans, -viewrange+ytrans, viewrange+ytrans, 1000, -1000);
		} else {
			
			glMatrixMode(GL_MODELVIEW);
			
			float m[16];
			glGetFloatv(GL_MODELVIEW_MATRIX , m);
			
			if(modifier == (GLUT_ACTIVE_CTRL | GLUT_ACTIVE_SHIFT))
				glRotatef( -0.2*(x-clickx0),0,0,1);
			else {
				glRotatef( -0.2*(y-clicky0),m[0],m[4],m[8]);
				glRotatef( -0.2*(x-clickx0),m[1],m[5],m[9]);
			}
		}
		
		clickx0 = x; clicky0 = y;
		
		glFlush();
		glFinish();
	}
}

#endif

namespace vsr {
	Visr* Visr::W = new Visr();
}
