///////////////////////////////////////////////////////////////////////////////////
// The Following program generates a Polygon with voltage controls via GPIO ports on 
// the raspberry pi mk II 
// 1024 x 1024 window using the OpenGL/GLUT libraries.
//	 
// Program Commands: 
// 	Quit the program: q
//	Increase theta value by 5 degrees: r 
//	Decrease Phi value by 0.1: z
//	Groound GPIO pin 17 to increase thetae by 5 
///////////////////////////////////////////////////////////////////////////////////


#include <wiringPi.h>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <GL/glut.h>

#define BUTTON_PIN 0 // GPIO pin 17 
//////////////////////////////////////////////////////////////////////////////////

// Back Ground Variable

const double pi2 = M_PI*2;
double mapPoints1[100][1000][2];
double mapPoints2[100][1000][2];

// Frame rate update period(ms)
int frameRate = 33;

// Parameter used for adjusting eye position 
float theta = 0.0; 

// Parameter used for adjusting eye position 
float phi = 0.0; 

GLfloat sqrt2 = std::sqrt(2); 

// Paramerter used to specify the eye position
GLfloat eye_x = 5.0;
GLfloat eye_y = 5.0;
GLfloat eye_z = 1.0;

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

// Function for diffusion light with position
// (0.5, 0.25, 1.0) and color white (1.0, 1.0, 1.0)

void myLight()
{

  GLfloat light_position[] = {0.5, 0.25, 1.0};

  glLightfv(GL_LIGHT0, GL_POSITION, light_position);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  

}


// Draw the x, y, and z axes.
void drawAxes()
{
        const GLfloat length = 100.0;

        glLineWidth(2.0);
        glBegin(GL_LINES);

        // Draw the positive x axis (in bright red).
        glColor3f(1.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(length, 0.0, 0.0);

        // Draw the negative x axis (in dark red).
        glColor3f(0.5, 0.0, 0.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(-length, 0.0, 0.0);

        // Draw the positive y axis (in bright green).
        glColor3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, length, 0.0);

        // Draw the negative y axis (in dark green).
        glColor3f(0.0, 0.5, 0.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, -length, 0.0);

        // Draw the positive z axis (in bright blue).
        glColor3f(0.0, 0.0, 1.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, length);

        // Draw the negative z axis (in dark blue).
        glColor3f(0.0, 0.0, 0.5);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, -length);

        glEnd();
}

void drawHexagon()
{
	

	GLfloat Base_Hex[8][3];
	glBegin(GL_POLYGON);		
	glColor3f(0.9,0.0,0.9);
	glNormal3f(0,0,1);
	for (GLint i = 0; i < 8; ++i)
	{

		GLfloat angle = (M_PI * i ) / 4;
		GLfloat x_tempVertex = cos(angle);
		GLfloat y_tempVertex = sin(angle);
		Base_Hex[i][0] = x_tempVertex;
		Base_Hex[i][1] = y_tempVertex;
		Base_Hex[i][2] = 1;
		glVertex3f(Base_Hex[i][0], Base_Hex[i][1], Base_Hex[i][2]);//Draw Top
	}
	glEnd();

	
	
	glBegin(GL_POLYGON);
	glNormal3f(0,0,-1);
	for(GLint i= 0; i < 8; ++i)
	{
		glVertex3f(Base_Hex[i][1], Base_Hex[i][0], -Base_Hex[i][2]);//Draw Base
	}
	glEnd();

	glLineWidth(3); 
	glBegin(GL_LINES);
	for (GLint i = 0; i < 7; ++i)
	{
		glColor3f((0.0 + 0.1*i),(0.0 + 0.1*i), 0.9);
		glVertex3f(Base_Hex[i][0], Base_Hex[i][1], Base_Hex[i][2]);
		glVertex3f(Base_Hex[i+1][0], Base_Hex[i+1][1], Base_Hex[i+1][2]);	
	}

	glEnd();
	
	GLfloat Side_Hex1[8][3]; // Side Hex1

	glBegin(GL_POLYGON);
	glColor3f(0.9,0.3,0.9);
	glNormal3f(0,1,0);
	for (GLint i = 0; i < 8; ++i)
	{

		GLfloat angle = (M_PI * i) / 4;
		GLfloat x_tempVertex = cos(angle);
		GLfloat z_tempVertex = sin(angle);
		Side_Hex1[i][0] = x_tempVertex; 
		Side_Hex1[i][1] = 1;		
		Side_Hex1[i][2] = z_tempVertex;
		glVertex3f(Side_Hex1[i][2], Side_Hex1[i][1], Side_Hex1[i][0]);
	}

	glEnd();
		
	glBegin(GL_POLYGON);
	glColor3f(0.9,0.3,0.9);
	glNormal3f(0,-1,0);
	for (GLint i = 0; i < 8; ++i)
	{
		glVertex3f(Side_Hex1[i][0], -Side_Hex1[i][1], Side_Hex1[i][2]);
	
	}	
	glEnd();

	GLfloat Side_Hex2[8][3];

	glBegin(GL_POLYGON);
	glColor3f(0.9,0.6,0.9);
	glNormal3f(1,0,0);
	for (GLint i = 0; i < 8; ++i)
	{

		GLfloat angle = (M_PI * i) / 4;
		GLfloat y_tempVertex = cos(angle);
		GLfloat z_tempVertex = sin(angle);
		Side_Hex2[i][0] = 1; 
		Side_Hex2[i][1] = y_tempVertex;		
		Side_Hex2[i][2] = z_tempVertex;
		glVertex3f(Side_Hex2[i][0], Side_Hex2[i][1], Side_Hex2[i][2]);
	}

	glEnd();
		
	glLineWidth(3); 
	glBegin(GL_LINES);
	for (GLint i = 0; i < 7; ++i)
	{
		glColor3f((0.0 + 0.1*i),(0.0 + 0.1*i), 0.9);
		glVertex3f(Side_Hex2[i][0], Side_Hex2[i][1], Side_Hex2[i][2]);
		glVertex3f(Side_Hex2[i+1][0], Side_Hex2[i+1][1], Side_Hex2[i+1][2]);	
	}

	glEnd();

	glBegin(GL_POLYGON);
	glColor3f(0.9,0.6,0.9);
	glNormal3f(-1,0,0);
	for (GLint i = 0; i < 8; ++i)
	{
		glVertex3f(-Side_Hex2[i][0], Side_Hex2[i][2], Side_Hex2[i][1]);
	
	}	
	glEnd();

	glBegin(GL_TRIANGLES);
	glColor3f(0.3,0.3,0.9); 
	glVertex3f(Side_Hex2[1][0], Side_Hex2[1][1], Side_Hex2[1][2]);
	glVertex3f(Side_Hex2[0][0], Side_Hex2[0][1], Side_Hex2[0][2]); 
	glVertex3f(Side_Hex1[1][0], Side_Hex1[1][1], Side_Hex1[1][2]); 
	
	glEnd(); 

	glBegin(GL_TRIANGLES);
	glColor3f(0.3,0.3,0.9); 
	glVertex3f(Side_Hex2[2][0], Side_Hex2[2][1], Side_Hex2[2][2]);
	glVertex3f(Side_Hex2[1][0], Side_Hex2[1][1], Side_Hex2[1][2]); 
	glVertex3f(Base_Hex[1][0], Base_Hex[1][1], Base_Hex[1][2]); 
	
	glEnd(); 

	glBegin(GL_TRIANGLES);
	glColor3f(0.3,0.3,0.9); 
	glVertex3f(Base_Hex[2][0], Base_Hex[2][1], Base_Hex[2][2]);
	glVertex3f(Base_Hex[1][0], Base_Hex[1][1], Base_Hex[1][2]);
	glVertex3f(Side_Hex1[1][0], Side_Hex1[1][1], Side_Hex1[1][2]); 
	
	glEnd(); 


	glBegin(GL_TRIANGLES);
	glColor3f(0.7,0.7,0.9); 
	glVertex3f(Base_Hex[1][0], Base_Hex[1][1], Base_Hex[1][2]);
	glVertex3f(Side_Hex2[1][0], Side_Hex2[1][1], Side_Hex2[1][2]); 
	glVertex3f(Side_Hex1[1][0], Side_Hex1[1][1], Side_Hex1[1][2]); 
	
	glEnd(); 

}

//////////////////////////////////////////////////////////////////////////////////////
// Back Ground Code is based on : example2.c written by Miguel A Sepulveda. 
// Origonal code: http://www.linuxfocus.org/English/January1998/article17.html
/////////////////////////////////////////////////////////////////////////////////////

void NonLinearMap(double *x, double *y)
{
	

	static double K = 1.04295;
	*y += K * sin(*x);
	*x += *y; 
	*x = fmod(*x, pi2);
	if(*x < 0.0) *x += pi2; 
}

void MapCalc(void)
{
	
	const int NumberSteps = 1000;
	const int NumberOrbits = 100; 
	const double Delta_x = pi2/(NumberOrbits-1); 
	int step, orbit; 
	glColor3f(1.0,1.0,1.0);
	
	for (orbit = 0; orbit < NumberOrbits; orbit++)
	{
      		double x, y;
      		y = 3.1415;
      		x = Delta_x * orbit;
      		for (step = 0; step < NumberSteps; step++)
		{
        			NonLinearMap(&x, &y);
				mapPoints1[orbit][step][1] = x*1.2;
				mapPoints1[orbit][step][2] = y-cos(degToRad(3.6*orbit));
      		};
    	};
	
    	for (orbit = 0; orbit < NumberOrbits; orbit++)
	{
     		double x, y;
      		x = 3.1415;
        	y = Delta_x * orbit;

       		for (step = 0; step < NumberSteps; step++)
		{
          		NonLinearMap(&x, &y);
			mapPoints2[orbit][step][1] = x*1.2;
			mapPoints2[orbit][step][2] = y-cos(degToRad(3.6*orbit));
        	};

     	};

}

void MapHandler(void)
{

    	for (GLint i = 0; i < 100; ++i)
	{
	        glBegin(GL_POINTS);

       		for (GLint j = 0; j < 1000; ++j)
		{
          		glVertex3f(0, mapPoints1[i][j][2], mapPoints1[i][j][1]);
        	}

       	 	glEnd();

     	}


    	for (GLint i = 0; i < 100; ++i)
	{
	        glBegin(GL_POINTS);

       		for (GLint j = 0; j < 1000; ++j)
		{
          		glVertex3f(0, mapPoints2[i][j][2], mapPoints2[i][j][1]);
        	};

       	 	glEnd();

     	}


}

/////////////////////////////////////////////////////////////////////////////////////

void display()
{
 	

	//Clear window 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Set the current matric to model view 
	glMatrixMode(GL_MODELVIEW);
	
	// Save the current modelview matrix
	glPushMatrix(); 


	//Set eye position 
	GLfloat eyeX = 5.0;
	GLfloat eyeY = 0.0;
	GLfloat eyeZ = 1.0; 
	gluLookAt(eyeX, eyeY, eyeZ, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0); 

	// Set the position of light to 0
	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

	glPushMatrix();
	drawAxes(); 
	glPopMatrix();	

	// Draw Tetrahedron with light enabled
	glPushMatrix();
	glRotatef(degToRad(45+theta), 0.0, 1.0, 0.0);
	glRotatef(degToRad(phi), 0.0, 0.0, 1.0);
	drawHexagon(); 
	glPopMatrix(); 

	glPushMatrix();
	glTranslatef(-1.0,-3.5,-4.0);
	MapHandler(); 
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
	gluPerspective(70.0, aspectRatio, 1.0, 150.0);
	
	// Initalize the model view matrix 
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity(); 
	

}

void keyboard(unsigned char key, int x, int y)
{

	switch(key) {
	case 'r':	
		theta += 5;
		glutPostRedisplay(); 
		break;

	case 'Z':
		phi += 5;
		glutPostRedisplay(); 
		break;
	case 'z':
		phi -= 5;
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
	MapCalc(); 
	
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



