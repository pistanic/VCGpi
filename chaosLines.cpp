///////////////////////////////////////////////////////////////////////////////////
// The Following program generates a Polygon with voltage controls via GPIO ports on 
// the raspberry pi mk II 
// 1024 x 1024 window using the OpenGL/GLUT libraries.
//	 
// Program Commands: 
// 	Quit the program: q
//	"Walk" : w,s,a,d
//	Subdevide Trapizod: x
//	Randomize topographic points: f
// 
///////////////////////////////////////////////////////////////////////////////////

#include <queue>
#include <iostream>
#include <time.h>
#include <cmath>
#include <cstdlib>
#include <GL/glut.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Subdivision_method_3.h>
#include <SPL/cgalUtil.hpp>
#include <SPL/math.hpp>
#include <SPL/Arcball.hpp>

#define BUTTON_PIN 0 // GPIO pin 17

using SPL::norm;
using SPL::radToDeg;
using SPL::degToRad;
using SPL::normalize;
using SPL::angleBetweenVectors;

////////////////////////////////////////////////////////////////////////////////
// Types
////////////////////////////////////////////////////////////////////////////////

// Basic types.
typedef double Real;
typedef CGAL::Cartesian<Real> Kernel;
typedef Kernel::Point_3 Point3;
typedef Kernel::Point_2 Point2;
typedef Kernel::Vector_3 Vector3;
typedef CGAL::Bbox_3 Bbox_3;
typedef SPL::Arcball<Kernel> ArcBall;
typedef SPL::Rotation_3<Kernel> Rotation3;

///////////////////////////////////////////////////////////////////////////
// Following containers were provided by Dr. Michael Adams
///////////////////////////////////////////////////////////////////////////

// A vertex type (to be used with Polyhedron_3) that includes some 
// // user-defined information.
template <class Refs, class Traits, class P>
struct MyVertex : public CGAL::HalfedgeDS_vertex_base<Refs, Traits, P>
{
        MyVertex() : CGAL::HalfedgeDS_vertex_base<Refs, Traits, P>() {}
        MyVertex(const Kernel::Point_3& p) :
                CGAL::HalfedgeDS_vertex_base<Refs, Traits, P>(p) {}
        Point3 position; // The limit position of the vertex.
        Vector3 normal; // The normal of the limit surface at the limit position.
};

// An items type (to be used with Polyhedron_3) that makes use of the
// // above vertex type.
struct MyItems : public CGAL::Polyhedron_items_3 {
        template <class Refs, class Traits>
        struct Vertex_wrapper {
                typedef MyVertex<Refs, CGAL::Tag_true, typename Traits::Point_3>
                  Vertex;
        };
};

/////////////////////////////////////////////////////////////////////////////
// End of Provided Code
/////////////////////////////////////////////////////////////////////////////

typedef CGAL::Polyhedron_3<Kernel, MyItems> Polyhedron;


struct Info
{

      //static constexpr Real sphereScale = 0.1;
        const Real sphereScale = 0.05;
        Info() : eyePos(0, 0, 10), sceneCenter(0, 0, 0), eyeUpDir(0, 1, 0),
          scale(1.0), rot(Vector3(0, 0, 1), 0), rotateMethod(0), arcBallRadius(1.0){}
        Polyhedron mesh;  // The polyhedral mesh.
        Bbox_3 boundBox;  // The bounding box for the mesh.
        Real smallEdge;  // The length of the shortest edge in the mesh.
        int viewportWidth;  // The viewport width.
        int viewportHeight;  // The viewport height.
        int rotateMethod; 
	Point3 eyePos;  // The eye position.
        Point3 sceneCenter;  // The center of the scene.
        Vector3 eyeUpDir;  // The eye's up direction.
        Real scale;  // The scaling used in drawing the mesh object.
	Real arcBallRadius;
	Rotation3 rot;
	ArcBall arcBall; 
        bool displayMeshEdge = false;
        int numSubDevision = 0;
        int divLevel = 0;
        std::queue <Polyhedron> prevPoly; //store previous polyhedrons  

};

Info info; 
void polyCalc(Polyhedron &mesh);  
//////////////////////////////////////////////////////////////////////////////////

// Back Ground Variable

const double pi2 = M_PI*2;
double mapPoints1[100][1000][3];
double mapPoints2[100][1000][3];

// Frame rate update period(ms)
int frameRate = 33;

// Parameter used for adjusting eye position 
float theta = 0.0; 

// Parameter used for adjusting eye position 
float phi = 0.0; 

float walk_x = 0.0;

float walk_y = 0.0; 

GLfloat sqrt2 = std::sqrt(2); 

bool dance = false; 

//RGBA color of diffuse light 
const GLfloat diffuseLight[] = 
{
	1.0, 1.0, 1.0
};

// Position of light source
const GLfloat lightPosition[] = 
{
	0.0, 0.0, 0.0, 0.0
};

// Trapoziod Mesh 
const std::string defaultMesh(
        "OFF\n"
        "6 8 0\n"
        "-0.5 -0.5  0\n"
        " 0.5 -0.5  0\n"
        " 0.5  0.5  0\n"
        "-0.5  0.5  0\n"
        " 0  0  1.35\n"
	" 0  0 -1.35\n"
        "3 0 1 4\n"
        "3 1 2 4\n"
        "3 2 3 4\n"
        "3 0 4 3\n"
	"3 2 1 5\n"
	"3 1 0 5\n"
	"3 0 3 5\n"
	"3 3 2 5\n"
);

// Function for diffusion light with position
// (0.5, 0.25, 1.0) and color white (1.0, 1.0, 1.0)

void myLight()
{

  GLfloat light_position[] = {-7.0, -7.0, 7.0};

  glLightfv(GL_LIGHT0, GL_POSITION, light_position);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  

}

// Draw a sphere.
void drawSphere()
{
        glutSolidSphere(1.0, 32, 32);
}

// Draw a rectangular prism with the specified center of gravity, axis,
// axis length, and "radius".
void drawRectPrism(const Point3& center, Real length, Real radius,
  const Vector3& axis)
{
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();

        Vector3 dir = axis / norm(axis);
        Real theta = angleBetweenVectors(dir, Vector3(0.0, 0.0, 1.0));
        glTranslatef(center.x(), center.y(), center.z());
        if (fabs(theta) > 1e-6) 
	{   
	      Vector3 rotAxis = CGAL::cross_product(Vector3(0.0, 0.0, 1.0), dir);
              glRotatef(radToDeg(theta), rotAxis.x(), rotAxis.y(), rotAxis.z());

	}
        glScalef(radius, radius, length);
        glutSolidCube(1.0);

        glPopMatrix();
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
//////////////////////////////////////////////////////////////////////////////////////
// Land scape Code is based on : example2.c written by Miguel A Sepulveda. 
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
	
	for (orbit = 0; orbit < NumberOrbits; orbit++)
	{
      		double x, y, z;
      		y = 3.1415;
      		x = Delta_x * orbit;
		z = cos(degToRad(3.6*orbit));
      		for (step = 0; step < NumberSteps; step++)
		{
        			NonLinearMap(&x, &y);
				mapPoints1[orbit][step][0] = x*1.1;
				mapPoints1[orbit][step][1] = y*1.1;
				mapPoints1[orbit][step][2] = z*1.4; 
      		}
    	}
	
    	for (orbit = 0; orbit < NumberOrbits; orbit++)
	{
     		double x, y, z;
      		x = 3.1415;
        	y = Delta_x * orbit;
		z = cos(degToRad(3.6*orbit)); 
       		for (step = 0; step < NumberSteps; step++)
		{
          		NonLinearMap(&x, &y);
			mapPoints2[orbit][step][0] = x*1.1;
			mapPoints2[orbit][step][1] = y*1.1;
			mapPoints2[orbit][step][2] = z*1.1;
        	}

     	}

}

void MapHandler(void)
{
	GLint rand; 
	
	// Tunnel 
    	for (GLint i = 0; i < 100; ++i)
	{
		GLint isOdd = i%2; 

		if(i < 25) glColor3f(0.0,1.0,1.0);
		else if(i < 50) glColor3f(0.5,1.0,1.0);
		else if(i < 75) glColor3f(1.0,0.5,0.0);
		else if(i < 100) glColor3f(0.5,1.0,0.5);

		glPointSize(2);
	        glBegin(GL_POINTS);

       		for (GLint j = 0; j < 1000; ++j)
		{
			GLfloat zpos  = 0.01; 
		       	if(dance) rand = std::rand() % 10; 
			else zpos = 0;

			zpos = 0.0 - zpos * float(rand);

   			glVertex3f(mapPoints1[i][j][0], 
				mapPoints1[i][j][1],
				mapPoints1[i][j][2]-zpos);
		
		}

       	 	glEnd();		
	
     	}

    	for (GLint i = 0; i < 100; ++i)
	{

		if(i < 25) glColor3f(0.0,1.0,1.0);
		else if(i < 50) glColor3f(0.5,1.0,1.0);
		else if(i < 75) glColor3f(1.0,0.5,0.0);
		else if(i < 100) glColor3f(0.5,1.0,0.5);

	        glBegin(GL_LINES);

       		for (GLint j = 0; j < 20; ++j)
		{

   			glVertex3f(mapPoints2[i][j*50][0], 
				mapPoints2[i][j*50][1],
				mapPoints2[i][j*50][2]);

   			glVertex3f(mapPoints2[i][1+j*50][0], 
				mapPoints2[i][1+j*50][1],
				mapPoints2[i][1+j*50][2]);
			
		}

       	 	glEnd();		
	
     	}


	// Mounds
    	for (GLint i = 0; i < 100; ++i)
	{

	      	GLint isOdd = i%2;

		if(i < 25) glColor3f(0.0,1.0,1.0);
		else if(i < 50) glColor3f(0.5,1.0,1.0);
		else if(i < 75) glColor3f(1.0,0.5,0.0);
		else if(i < 100) glColor3f(0.5,1.0,0.5);
		
		glPointSize(1);		
		glBegin(GL_POINTS);
       		for (GLint j = 0; j < 1000; ++j)
		{
			GLfloat zpos = 0.01; 

       			if(dance) rand = std::rand() % 10; 
			else zpos = 0;

			zpos = 0.0 - zpos * float(rand);
	
   			glVertex3f(mapPoints2[i][j][0], 
				mapPoints2[i][j][1],
				mapPoints2[i][j][2]-zpos);

		}
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
	GLfloat eyeX = -4.0+walk_x;
	GLfloat eyeY = 0.0+walk_y;
	GLfloat eyeZ = -0.25; 
	gluLookAt(eyeX, eyeY, eyeZ, 1.0, 0.0, -0.25, 0.0, 0.0, 1.0); 

	// Set the position of light to 0
	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);


	glPushMatrix();
//	drawAxes(); 
	glPopMatrix();	

	// Draw Tetrahedron with light enabled
	glPushMatrix();
	glRotatef(degToRad(theta), 0.0, 1.0, 0.0);
	glRotatef(degToRad(phi), 0.0, 0.0, 1.0);
	//drawHexagon(); 
	glPopMatrix(); 

	glPushMatrix();
	//GLint rand = -50;
	glColor3f(0.9,0.9,0.9);
	glRotatef(degToRad(theta-15), 0.0, 1.0, 0.0);
	glRotatef(degToRad(phi-120), 0.0, 0.0, 1.0);
	glTranslatef(-3.5,-3.5,0.0);
	MapHandler(); 
	glPopMatrix();

	glEnable(GL_LIGHTING);

	glTranslatef(0.0,0.0,0.0);
	glColor3f(0.9,0.0,0.9);

	// For each vertex in the mesh...
        for (Polyhedron::Vertex_const_iterator vertexIter =
          info.mesh.vertices_begin(); vertexIter != info.mesh.vertices_end();
          ++vertexIter) 
	{
                // Draw a sphere at the vertex position.
                Point3 v = vertexIter->point();
                glPushMatrix();
                glTranslatef(v.x(),v.y(),v.z());
                glScalef(info.sphereScale * info.smallEdge, info.sphereScale *
                  info.smallEdge, info.sphereScale * info.smallEdge);
                drawSphere();
                glPopMatrix();
        }

        // For each edge in the mesh...
        for (Polyhedron::Edge_const_iterator edgeIter = info.mesh.edges_begin();
          edgeIter != info.mesh.edges_end(); ++edgeIter) 
	{
                // Draw a rectangular prism along the extent of the edge.
                Point3 v0 = edgeIter->vertex()->point();
                Point3 v1 = edgeIter->opposite()->vertex()->point();
                Point3 midpoint = CGAL::midpoint(v0, v1);
                Vector3 axis = v1 - v0;
                drawRectPrism(midpoint, norm(axis), 0.25 * info.sphereScale *
                  info.smallEdge, axis);
        }
	// Flush graphics output to framebuffer
	glPopMatrix();
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
	gluPerspective(45.0, aspectRatio, 1.0, 1000.0);
	
	// Initalize the model view matrix 
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity(); 
	

}

void keyboard(unsigned char key, int x, int y)
{

	switch(key) {
	case 'r':	
		theta += 15;
		glutPostRedisplay(); 
		break;
	case 'R':	
		theta -= 15;
		glutPostRedisplay(); 
		break;
	case 'Z':
		phi += 15;
		glutPostRedisplay(); 
		break;
	case 'z':
		phi -= 15;
		glutPostRedisplay(); 
		break;
	case 'w':
		walk_x += 0.05;
		glutPostRedisplay(); 
		break;
	case 's':
		walk_x -= 0.05;
		glutPostRedisplay(); 
		break;
	case 'd':
		walk_y -= 0.05;
		glutPostRedisplay(); 
		break;
	case 'a':
		walk_y += 0.05;
		glutPostRedisplay(); 
		break;
	case 'f':
		if(dance) dance = false; 
		else dance = true; 
		glutPostRedisplay();
		break;  
	case 'S':
		dance = false; 
		glutPostRedisplay(); 
		break;
	case 'i':
		theta = 0; 
		phi = 0; 
		walk_x = 0;
		walk_y = 0;  
		dance = false; 
	 	glutPostRedisplay(); 
		break; 
        case 'x': // Increase the number of sub devision levels by one 
                info.prevPoly.push(info.mesh); // Save Current mesh polyhedron
                info.numSubDevision += 1;
                ++info.divLevel;
		//
		if(info.numSubDevision > 1)
		{
		  info.numSubDevision = 0;
		  info.divLevel = 0;
		  info.mesh = info.prevPoly.front(); 
		}
		//
                CGAL::Subdivision_method_3::CatmullClark_subdivision(info.mesh, info.numSubDevision);
                polyCalc(info.mesh);
                glutPostRedisplay();
                break;
        //case 'X': // Decrease the number of sub devision levels by one
                //if(info.divLevel != 0)
                //{
                 // info.numSubDevision -= 1;
                 //--info.divLevel;
                 // info.mesh = info.prevPoly.top(); // Restore Previous mesh 
                 // info.prevPoly.pop();
                //}
                polyCalc(info.mesh);
                glutPostRedisplay();
                break;
	case 'q':
		exit(0);
		break;;

	}

}

/*
void voltageInput(int value)
{
	if(!digitalRead(BUTTON_PIN))
	{
		theta += 5; 
		glutPostRedisplay(); 
	}
	glutTimerFunc(15,voltageInput,0);
}
*/

int main(int argc, char **argv)
{
	
	MapCalc(); 

	srand(time(NULL));
	
	const int winWidth = 1024; // normalized window width  
	const int winHeight = 1024; //normalized window height 
	std::stringstream inStream(defaultMesh);
        inStream >> info.mesh;

        // Compute the bounding box of the mesh.
        if (info.mesh.size_of_vertices() > 0) {
                Point3 v = info.mesh.vertices_begin()->point();
                info.boundBox = Bbox_3(v.x(), v.y(), v.z(), v.x(), v.y(), v.z());
        } else {
                info.boundBox = Bbox_3(0, 0, 0, 0, 0, 0);
        }
        for (Polyhedron::Vertex_const_iterator vertexIter =
          info.mesh.vertices_begin(); vertexIter != info.mesh.vertices_end();
          ++vertexIter) {
                const Point3& v = vertexIter->point();
                info.boundBox = info.boundBox + Bbox_3(v.x(), v.y(), v.z(), v.x(),
                  v.y(), v.z());
        }

        // Compute the length of the shortest edge in the mesh.
        info.smallEdge = -1.0;
        for (Polyhedron::Edge_const_iterator edgeIter = info.mesh.edges_begin();
          edgeIter != info.mesh.edges_end(); ++edgeIter) {
                Point3 v0 = edgeIter->vertex()->point();
                Point3 v1 = edgeIter->opposite()->vertex()->point();
                Real length = norm(v1 - v0);
                if (info.smallEdge < 0.0 || length < info.smallEdge) {
                        info.smallEdge = length;
                }
        }

	polyCalc(info.mesh);

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

void polyCalc(Polyhedron &mesh)
{
	
	for(auto vertIt = mesh.vertices_begin(); vertIt != mesh.vertices_end();
									 ++vertIt)
	{
	 #ifdef DEBUG
   	   std::cerr << "Here we go\n";
	 #endif   

	//Relevant variable declarations
	  bool edge = false;

	  double Bn = 0; 
	  double An = 0; 
   	  double x, x0, xN;
  	  double y, y0, yN;
	  double z, z0, zN;   

	 
 
	  const int n = vertIt-> vertex_degree();
	  #ifdef DEBUG
   	    std::cerr << "Vertex Degree: "<< n << "\n";
	  #endif   
	
          auto Circu = vertIt-> vertex_begin(); //Half edge circulator  
	 
	  do
	  {
	    if(Circu-> is_border_edge())
	    {
	      edge = true; 
	    }
	    ++Circu;      
	  }while(Circu != vertIt-> vertex_begin());// Circulate over surronding	 
  	  				// half edges to determine if border edge 
	  Bn = (1.0/(double)n)*( (5.0/8.0) - 
	       pow( ( (3.0/8.0)+ ((1.0/4.0)*cos((2.0*M_PI)/(double)n)) ),2.0) );
 
 	  An = 1.0/( (3.0/(8.0*Bn)) + (double)n );
	
///////////////////////////////////////////////////////////////////////////////////
//	 Interior Point Calulation 
//////////////////////////////////////////////////////////////////////////////////
	  if(!edge) 
	  {
	 	x = x0 = xN = 0;
	  	y = y0 = yN = 0;
		z = z0 = zN = 0;   

		Circu = vertIt-> vertex_begin(); // Reinitalize ciruclator 
	  	x0 = vertIt-> point().x(); 
	  	y0 = vertIt-> point().y();	
	  	z0 = vertIt-> point().z();

	  	do
	  	{
		#ifdef DEBUG
		  std::cerr << "\n Alph and Beta values: "<< An <<" "<< Bn;
		  std::cerr << "\nWorking interior x/y/z vertex limit\n";
	      	  std::cerr << "x,y,z: " << x <<" "<< y <<" "<< z <<"\n\n"; 
	        #endif
 
	  	  x += Circu-> opposite()-> vertex()-> point().x()*An; 
	          y += Circu-> opposite()-> vertex()-> point().y()*An;	
	      	  z += Circu-> opposite()-> vertex()-> point().z()*An;
	      	  ++Circu;
 
	    	}while(Circu != vertIt-> vertex_begin());

	    	x += (1.0-An*n)*x0;
	    	y += (1.0-An*n)*y0;
	    	z += (1.0-An*n)*z0;
	    	#ifdef DEBUG
		  std::cerr << "Stored interior x/y/z vertex limit\n";
	      	  std::cerr << "x,y,z: " << x <<" "<< y <<" "<< z <<"\n"; 
	        #endif
		// Store Points to Vertex
		Point3 vlim(x,y,z); 
	    	vertIt-> position = vlim; 

///////////////////////////////////////////////////////////////////////////////////
//	    Normal Vector Calulation 
///////////////////////////////////////////////////////////////////////////////////

	    	double taoi; 
	    	std::vector<double> tao;

	    	Circu = vertIt-> vertex_begin(); //Half edge Circulator 

	    	for(int i = 0; i < n; ++i)//Generate t0 and tao vector 
	    	{
	      	  tao.push_back(cos( (2.0*M_PI*(double)i) / (double)n ));
	      	  taoi = tao[i];
	      	  #ifdef DEBUG
	            std::cerr << "Current Vertex xN " 
			<< Circu-> opposite()-> vertex()-> point().x()<< "\n";
	            std::cerr << "Current Vertex yN " 
			<< Circu-> opposite()-> vertex()-> point().y()<< "\n";
                    std::cerr << "Current Vertex zN " 
			<< Circu-> opposite()-> vertex()-> point().z()<< "\n";
	            std::cerr << "tao (" << i << ") =  " << taoi << "\n";
	      	  #endif
  	      	  xN += Circu-> opposite()-> vertex()-> point().x()*taoi; 
	      	  yN += Circu-> opposite()-> vertex()-> point().y()*taoi;	
	      	  zN += Circu-> opposite()-> vertex()-> point().z()*taoi;
	      	  #ifdef DEBUG
	          std::cerr << "xN, yN, zN: "<< xN <<" "<< yN <<" "<< zN <<"\n\n"; 
	      	  #endif 
	      	  ++Circu;   
	    	}
	    	Vector3 t0(xN,yN,zN); //Tangest Vector t0
	    	#ifdef DEBUG
   	          std::cerr << "\n\nVector t0 calculated"<< t0  <<"\n\n";
	    	#endif   

    	    	xN = yN = zN = 0; // Clear normal values 

 	    	Circu = vertIt-> vertex_begin();

	    	for(int i = 0; i < n; ++i)//Generate t1 with previously
	    	{			      // // generated tao vector 
	      	  taoi = tao[n-i-1];
	      	  #ifdef DEBUG
		    std::cerr << "(tao i-n-1) ("<< i <<") " << taoi << "\n";
	            std::cerr << "Current Vertex xN " 
			<< Circu-> opposite()-> vertex()-> point().x()<< "\n";
	            std::cerr << "Current  Vertex yN " 
			<< Circu-> opposite()-> vertex()-> point().y()<< "\n";
                    std::cerr << "Current Vertex zN " 
			<< Circu-> opposite()-> vertex()-> point().z()<< "\n";
  	          #endif
 	       	  xN += Circu-> opposite()-> vertex()-> point().x()*taoi; 
	       	  yN += Circu-> opposite()-> vertex()-> point().y()*taoi;	
	       	  zN += Circu-> opposite()-> vertex()-> point().z()*taoi;
	       	  ++Circu;
	    	  #ifdef DEBUG
	      std::cerr << "xN2,yN2,zN2: "<< xN << " " << yN 
	                << " " << zN <<"\n\n"; 
	    	  #endif   
	    	}	

	    	Vector3 t1(xN,yN,zN); // Tanget Vector t1

	    	#ifdef DEBUG
   	      	  std::cerr << "\n\nVector t1 calulated"<< t1 <<"\n\n";
	    	#endif   

		// Store Points to Vertex
	    	Vector3 aNormal(CGAL::cross_product(t1,t0)); 
	    	Vector3 Normal = SPL::normalize(aNormal);
	   	vertIt-> normal = Normal;  
	   
	  }// End of interior point calulations

///////////////////////////////////////////////////////////////////////////////////
//	 Boundry Point 
///////////////////////////////////////////////////////////////////////////////////
	  else 
	  {
 		x = x0 = xN = 0;
	 	y = y0 = yN = 0;
	  	z = z0 = zN = 0; 	
   
		Circu = vertIt-> vertex_begin(); // Reinitalize ciruclator 	 

		x0 = Circu-> vertex()-> point().x(); 
	  	y0 = Circu-> vertex()-> point().y();	
	 	z0 = Circu-> vertex()-> point().z();

	    	#ifdef DEBUG
	      	  int neb = 0; 
	        #endif

	    	do
	    	{	
		  if(Circu -> is_border_edge())
	      	  {
	      	    #ifdef DEBUG
   	       	      std::cerr << "Neighbour "<< neb <<"\n";
       		      std::cerr << "x for Neighbour: " << neb << " = " 
	                << Circu-> opposite()-> vertex()-> point().x() <<"\n"; 
       		      std::cerr << "y for Neighbour: " << neb << " = " 
	                << Circu-> opposite()-> vertex()-> point().y() <<"\n"; 
	       	      std::cerr << "z for Neighbour: " << neb << " = " 
	                << Circu-> opposite()-> vertex()-> point().z() <<"\n"; 
        	      ++neb; 
		    #endif  
       		    x += Circu-> opposite()-> vertex()-> point().x()*(1.0/6.0); 
	            y += Circu-> opposite()-> vertex()-> point().y()*(1.0/6.0);
		    z += Circu-> opposite()-> vertex()-> point().z()*(1.0/6.0);
	      	  }

	      	  ++Circu; 

		}while(Circu != vertIt-> vertex_begin()); 

		x += (2.0/3.0)*x0;
	    	y += (2.0/3.0)*y0;
	   	z += (2.0/3.0)*z0;

///////////////////////////////////////////////////////////////////////////////////
//		Normal vector calulation 
///////////////////////////////////////////////////////////////////////////////////
	
		Vector3 tNxt,tPrv,tmp;
 	
		Circu = vertIt-> vertex_begin(); 

		double xNxt,yNxt,zNxt;
		double xPrv,yPrv,zPrv;
		double xS,yS,zS; 
		int face = 0;  

		xS = yS = zS = 0; 

		xN = vertIt-> point().x(); 
		yN = vertIt-> point().y(); 
		zN = vertIt-> point().z(); 

		for(int i = 0; i < n; ++i)
		{
		  if(!Circu-> is_border())
		  {
		    ++face; 	
	            xNxt = Circu-> next()-> vertex()-> point().x();
		    yNxt = Circu-> next()-> vertex()-> point().y();
	 	    zNxt = Circu-> next()-> vertex()-> point().z(); 
	            xPrv = Circu-> prev()-> vertex()-> point().x(); 
	  	    yPrv = Circu-> prev()-> vertex()-> point().y();
	  	    zPrv = Circu-> prev()-> vertex()-> point().z(); 
	 	    xNxt -= xN;
		    yNxt -= yN; 
	 	    zNxt -= zN;
		    xPrv -= xN; 
		    yPrv -= yN; 
		    zPrv -= zN;
	    	    #ifdef DEBUG
		      std::cerr << "\n\n NXT DIFFERANCE \n";
	      	      std::cerr << "xN,yN,zN: "
			 << xNxt <<" "<< yNxt <<" "<< zNxt <<"\n"; 
		      std::cerr << "\n\n PRV DIFFERANCE \n";
	      	      std::cerr << "xP,yP,zP: "
			 << xPrv <<" "<< yPrv <<" "<< zPrv <<"\n"; 
		      std::cerr << "\nFace: "<< face <<"\n";
	            #endif

		    tNxt = Vector3(xNxt,yNxt,zNxt);
		    tPrv = Vector3(xPrv,yPrv,zPrv);
		    tmp = CGAL::cross_product(tPrv,tNxt);
		    tmp = SPL::normalize(tmp); 
		    xS += tmp.x();
		    yS += tmp.y();
		    zS += tmp.z();
		    #ifdef DEBUG   
		    std::cerr << "\n\n SUMM POINTS \n";
	      	    std::cerr << "xS,yS,zS: "
		      << xS <<" "<< yS <<" "<< zS <<"\n"; 
	            #endif
		  }
		  ++Circu;
		}
		
		//tmp = SPL::normalize(Vector3(xS,yS,zS));
	        tmp = ((1.0/(double)face)*Vector3(xS,yS,zS)); 
 	        vertIt-> normal = tmp;    	// Stor Limit Point 
	    	Point3 vlim(x,y,z);
	    	vertIt -> position = vlim; 
		

	  } //End of Boundry Calulation 
	
	}

} 



