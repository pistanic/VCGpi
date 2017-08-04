///////////////////////////////////////////////////////////////////////////////////
// The Following program generates a cube, sphere, and tetrahedron in a normalized 
// 1024 x 1024 window using the OpenGL/GLUT libraries.
//	 
// Program Commands: 
// 	Quit the program: q
// 	Increase theta value by 5 degrees: R
// 	Increase Phi value by 0.1: Z
//	Decrease theta value by 5 degrees: r 
//	Decrease Phi value by 0.1: z
//	Reset theta and phi to origional valyes at progran start up: i
//
///////////////////////////////////////////////////////////////////////////////////

extern "C"
{
	#include <wiringPi.h>
}

#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <GL/glut.h>

#define BUTTON_PIN 0 // GPIO pin 17 
//////////////////////////////////////////////////////////////////////////////////

// Frame rate update period(ms)
int frameRate = 33;

// Parameter used for adjusting eye position 
float theta = 0.0; 

// Parameter used for adjusting eye position 
float phi = 0.0; 

GLfloat sqrt2 = std::sqrt(2); 

// Paramerter used to specify the eye position
GLfloat eye0 = 5.0;
GLfloat eye1 = 1.0;

//RGBA color of diffuse light 
const GLfloat diffuseLight[] = 
{
	1.0, 1.0, 1.0
};

// Position of light source
const GLfloat lightPosition[] = 
{
	0.25, 0.5, 1.0, 0.0
};

// Convert Degrees to radians 
inline double degToRad(double x)
{
	return x * M_PI / 180.0; 
}

// Draw blue (0.0, 0.0, 0.9) cube with sides of unit length,
// centered at (-1.5, 1.5, 0.0)
void drawCube()
{
	glPushMatrix(); 
	glColor3f(0.0, 0.0, 0.9); // Blue
	glTranslatef(-1.5, 1.5, 0.0);
	glutSolidCube(1.0);
	glPopMatrix(); 

}

// Draw magenta (0.9, 0.0, 0.9) sphere with radius 1/sqrt(2)
// centered at (0.0, 0.0, 0.0)
void drawSphere()
{
	glPushMatrix(); 
	glColor3f(0.9,0.0,0.9); // Magenta 
	glutSolidSphere((1.0/sqrt2), 15, 15);
	glPopMatrix(); 
}

// Draw cyan (0.0, 0.9, 0.9) tetrahedron with sides length sqrt(2)
// centered at (1.5, -1.5, 0.0)
void drawTetrahedron()
{
	// Tetrahedron Vertices 
	static const GLfloat vertices[][3] = 
	{

	  {0.0, 0.0, 0.61237},
	  {0.5, 0.28867, -0.20412},
	  {-0.5, 0.28867, -0.24012}, 
	  {0.0, -0.57735, -0.20412}

	};
	
	// Tetrahedon Faces 
	// triplets are sets of vertex indecies for face
	// Note: specified in CCW order
	static const int faces[][3] = 
	{

	  {0, 1, 2},
	  {0, 2, 3},
	  {0, 3, 1},
	  {3, 2, 1}

	};
	
	static const GLfloat normals[][3] = 
	{

	  {0.0, 0.942809, 0.3333334},
	  {-0.816497, -0.471405, 0.333333},	
	  {0.816497, -0.471405, 0.333333}, 
	  {0.0, 0.0, -1.0}
	 
	};

	// Draw tetrahedron 
	glBegin(GL_TRIANGLES);
	for(int i = 0; i < 4; ++i) // For each face
	{

	  glNormal3f(normals[i][0], normals[i][1], normals[i][2]);
 
	  for(int j = 0; j < 3; ++j) // For each vertex
	  {

	    int v = faces[i][j];
	    // Specify Colors 
	    glColor3f(0.0,0.9,0.9); // Cyan 
	    // Specify position of the vertex
	    glVertex3f(vertices[v][0], vertices[v][1], vertices[v][2]); 
	    
	  }

	}

	glEnd(); 	

}

void display()
{
 	

	//Clear window 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Set the current matric to model view 
	glMatrixMode(GL_MODELVIEW);
	
	// Save the current modelview matrix
	glPushMatrix(); 


	//Set eye position 
	GLfloat eyeX = eye0 * cos(degToRad(theta));
	GLfloat eyeY = eye0 * sin(degToRad(theta));
	GLfloat eyeZ = eye1 * phi; 
	gluLookAt(eyeX, eyeY, eyeZ, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0); 

	// Set the position of light to 0
	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

	// Draw Cube 
	drawCube(); 

	// Draw Sphere 
	drawSphere(); 
	
	// Draw Tetrahedron with light enabled
	glPushMatrix();
	glTranslatef(1.5, -1.5, 0.0);
	glScalef(sqrt2, sqrt2, sqrt2);
	drawTetrahedron(); 
	glPopMatrix(); 

	// Restore old modelview matrix
	glPopMatrix(); 
	
	// Flush graphics output to framebuffer
	glutSwapBuffers(); 
	

}

void reshape(GLint width, GLint height)
{
	//Compute Aspect ration 
	GLfloat aspectRatio = static_cast<GLfloat>(width) / 
	  ((height) ? height : 1.0);

	// Set viewport to entire window
	glViewport(0, 0, width, height);
	
	// Initalize the projection matrix
	// Done to maintain the aspect ratio in case of 
	// non square window
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 

	// Perspective projection 
	gluPerspective(70.0, aspectRatio, 1.0, 100.0);
	
	// Initalize the model view matrix 
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity(); 
	

}

void keyboard(unsigned char key, int x, int y)
{

	switch(key) {

	case 'R':
		theta += 5;
	 	glutPostRedisplay(); 
		break;
	case 'r':
		theta -= 5;
		glutPostRedisplay();
		break;
	case 'Z':
		phi += 0.1;
		glutPostRedisplay(); 
		break;
	case 'z':
		phi -= 0.1;
		glutPostRedisplay(); 
		break; 
	case 'i':
		theta = 0; 
		phi = 0; 
	 	glutPostRedisplay(); 
		break; 

	case 'q':
		exit(0);
		break;;

	}

}

void voltageInput(int value)
{
	if(!digitalRead(BUTTON_PIN))
	{
		theta += 5; 
		glutPostRedisplay(); 
	}
	glutTimerFunc(15,voltageInput,0);
}


int main(int argc, char **argv)
{
	
	wiringPiSetup();
	
	pinMode(BUTTON_PIN, INPUT);
	pullUpDnControl(BUTTON_PIN, PUD_UP); 
	
	const int winWidth = 1024; // normalized window width  
	const int winHeight = 1024; //normalized window height 
	
	// Glut initialization 
	glutInit(&argc, argv);
 
	// Specift the type of display mode
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); 

	// Set the window size
	glutInitWindowSize(winWidth, winHeight);
	
	// Create new window with the same name as
	// executable file 
	glutCreateWindow(argv[0]);

	// Display callback 
	glutDisplayFunc(display);

	// Regester a reshape function
	glutReshapeFunc(reshape);

	// Register a keyboard call back function
	glutKeyboardFunc(keyboard);

	voltageInput(0); 
	 
		// Set the clear color
	glClearColor(0.0, 0.0, 0.0, 0.0);

	// Set Color of diffuse light and enable 
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);  
		
	// Front material to track color 
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	glFrontFace(GL_CCW);

	glEnable(GL_DEPTH_TEST);
	
	glEnable(GL_CULL_FACE);

	glEnable(GL_NORMALIZE);

	// Enable the Glut main loop.
	glutMainLoop();
	
	
	return 0;
}



// Function for diffusion light with position
// (0.5, 0.25, 1.0) and color white (1.0, 1.0, 1.0)
/*
void myLight()
{

  GLfloat light_position[] = {0.5, 0.25, 1.0};

  glLightfv(GL_LIGHT0, GL_POSITION, light_position);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  

}
*/


