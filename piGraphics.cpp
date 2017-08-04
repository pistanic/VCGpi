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


#include <wiringPi.h>
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

void drawCube()
{

	glPushMatrix(); 
	glColor3f(0.9,0.3,0.9);
	glutSolidCube(1.0);
	glPopMatrix(); 

}

void drawHexagon()
{

	GLfloat firstHex[6][3];
	GLfloat secondHex[6][3];
	glBegin(GL_POLYGON);
	for (GLint i = 0; i < 6; i++)
	{
		GLfloat angle = 2 * M_PI / 6 * (i + 0.5);
		GLfloat x_tempVertex = sqrt(2)*cos(angle);
		GLfloat z_tempVertex = sqrt(2)*sin(angle);
		firstHex[i][0] = x_tempVertex;
		firstHex[i][1] = 0;
		firstHex[i][2] = z_tempVertex;
		glVertex3f(x_tempVertex, 0, z_tempVertex);
	}

	glEnd();

	glBegin(GL_POLYGON);

	for (GLint i = 0; i < 6; i++)
	{
		GLfloat angle = 2 * M_PI / 6 * (i + 0.5);
		GLfloat x_tempVertex = sqrt(2)*cos(angle);
		GLfloat z_tempVertex = sqrt(2)*sin(angle);
		secondHex[i][0] = x_tempVertex;
		secondHex[i][1] = sqrt(2);
		secondHex[i][2] = z_tempVertex;
		glVertex3f(x_tempVertex, sqrt(2), z_tempVertex);
	}

	glEnd();

// Begin filling the rectangles between the two hexagons
// drawing begins at ( 0, 0, 0 ) and moves clockwise

	glColor3f(0.9,0.0,0.9);
	glBegin(GL_QUADS);
	
	for (GLint i = 0; i < 5; i++)
	{
		glVertex3f(firstHex[i][0], firstHex[i][1], firstHex[i][2]); // top left
		glVertex3f(secondHex[i][0], secondHex[i][1], secondHex[i][2]); // bottom left
		glVertex3f(secondHex[i+1][0], secondHex[i+1][1], secondHex[i+1][2]); // top right
		glVertex3f(firstHex[i+1][0], firstHex[i+1][1], firstHex[i+1][2]); // bottom right
	}

	// corner case (connecting the end to the beginning)
	glVertex3f(firstHex[5][0], firstHex[5][1], firstHex[5][2]); // top left
	glVertex3f(secondHex[5][0], secondHex[5][1], secondHex[5][2]); // bottom left
	glVertex3f(secondHex[0][0], secondHex[0][1], secondHex[0][2]); // top right
	glVertex3f(firstHex[0][0], firstHex[0][1], firstHex[0][2]); // bottom right

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

	drawCube(); 

	// Draw Tetrahedron with light enabled
	glPushMatrix();
	//drawHexagon(); 
	
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


