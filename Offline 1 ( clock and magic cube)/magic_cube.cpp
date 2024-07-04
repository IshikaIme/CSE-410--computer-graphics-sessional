#include <windows.h>  // for MS Windows
#include <GL/glut.h>  // GLUT, include glu.h and gl.h
#include <cmath>
#include<vector>
GLfloat eyex = 4, eyey = 4, eyez = 4;
GLfloat centerx = 0, centery = 0, centerz = 0;
GLfloat upx = 0, upy = 1, upz = 0;
GLfloat vertexAx, vertexBx, vertexCx;
GLfloat vertexAy, vertexBy, vertexCy;
GLfloat vertexAz, vertexBz, vertexCz;
GLfloat scalefactor=1.0;
GLfloat distance=0.0;
GLfloat COGTriangleX, COGTriangleY, COGTriangleZ;

/* Initialize OpenGL Graphics */
void initGL() {
    // Set "clearing" or background color
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);   // Black and opaque
    glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
}



void drawTriangle(int x, int y, int z)
{   
    glPushMatrix();
    float translateDist = (1- scalefactor)/3;
    glTranslatef(translateDist, translateDist, translateDist);

    vertexAx= 0.0; 
    vertexAy= 0.0;
    vertexAz= scalefactor* 1.0;

    vertexBx= 0.0; 
    vertexBy= scalefactor* 1.0;
    vertexBz= 0.0;

    vertexCx= scalefactor* 1.0; 
    vertexCy= 0.0;
    vertexCz= 0.0;
    glBegin(GL_TRIANGLES);  // Begin drawing triangles
    glColor3d(x,y,z) ; // Set color to red

    // Front-facing triangles
    glVertex3f(vertexAx, vertexAy, vertexAz);  // Top vertex
    glVertex3f(vertexBx, vertexBy, vertexBz);  // Bottom-left vertex
    glVertex3f(vertexCx,vertexCy, vertexCz);  // Bottom-right vertex
        
    glEnd();  // End drawing triangles

    glLineWidth(5);
    glColor3f(1.0f, 1.0f, 1.0f);
    // glBegin(GL_LINES);
    
    // glVertex3f(vertexAx, vertexAy, vertexAz);
    // glVertex3f(COGTriangleX, COGTriangleY, COGTriangleZ);

    // glEnd();
    glPopMatrix();
}

void drawCylinder(double height, double radius, int segments) {
    double tempx = radius, tempy = 0;
    double currx, curry;
    glBegin(GL_QUADS);
        for (int i = 1; i <= segments; i++) {
            double theta = i * ( 70.528* 1 * M_PI /180) / segments;
            currx = radius * cos(theta);
            curry = radius * sin(theta);

           // GLfloat c = (2+cos(theta))/3;
           // glColor3f(c,c,c);
           glColor3f(1.0f,0.0f,0.0f);       // Red
            glVertex3f(currx, curry, height/2);
            glVertex3f(currx, curry, -height/2);

            glVertex3f(tempx, tempy, -height/2);
            glVertex3f(tempx, tempy, height/2);

            tempx = currx;
            tempy = curry;
        }
    glEnd();
}


void draw_Octahedron()
{

    //drawing upper pyramid
    drawTriangle(0,1,0);

    glPushMatrix();

    glRotatef(90,0,0,1);
    drawTriangle(1,0,0);
    glRotatef(90,0,0,1);
    drawTriangle(1,1,0);
    glRotatef(90,0,0,1);
    drawTriangle(1,1,1);
        

     glPopMatrix();   

    //drawing lower pyramid

     glPushMatrix();

    glRotatef(180,1,0,0);
    drawTriangle(1,0,1);


    glPushMatrix();

    glRotatef(90,0,0,1);
    drawTriangle(1,1,0);
    glRotatef(90,0,0,1);
    drawTriangle(1,1,1);
    glRotatef(90,0,0,1);
    drawTriangle(1,0,0);
        

    glPopMatrix();  

   glPopMatrix();

}


// generate vertices for +X face only by intersecting 2 circular planes
// (longitudinal and latitudinal) at the given longitude/latitude angles
std::vector<float> buildUnitPositiveX(int subdivision, float centerX, float centerY, float centerZ)
{
    const float DEG2RAD = acos(-1) / 180.0f;

    std::vector<float> vertices;
    float n1[3];        // normal of longitudinal plane rotating along Y-axis
    float n2[3];        // normal of latitudinal plane rotating along Z-axis
    float v[3];         // direction vector intersecting 2 planes, n1 x n2
    float a1;           // longitudinal angle along Y-axis
    float a2;           // latitudinal angle along Z-axis

    // compute the number of vertices per row, 2^n + 1
    int pointsPerRow = (int)pow(2, subdivision) + 1;

    // rotate latitudinal plane from 45 to -45 degrees along Z-axis (top-to-bottom)
    for(unsigned int i = 0; i < pointsPerRow; ++i)
    {
        // normal for latitudinal plane
        // if latitude angle is 0, then normal vector of latitude plane is n2=(0,1,0)
        // therefore, it is rotating (0,1,0) vector by latitude angle a2
        a2 = DEG2RAD * (45.0f - 90.0f * i / (pointsPerRow - 1));
        n2[0] = -sin(a2);
        n2[1] = cos(a2);
        n2[2] = 0;

        // rotate longitudinal plane from -45 to 45 along Y-axis (left-to-right)
        for(unsigned int j = 0; j < pointsPerRow; ++j)
        {
            // normal for longitudinal plane
            // if longitude angle is 0, then normal vector of longitude is n1=(0,0,-1)
            // therefore, it is rotating (0,0,-1) vector by longitude angle a1
            a1 = DEG2RAD * (-45.0f + 90.0f * j / (pointsPerRow - 1));
            n1[0] = -sin(a1);
            n1[1] = 0;
            n1[2] = -cos(a1);

            // find direction vector of intersected line, n1 x n2
            v[0] = n1[1] * n2[2] - n1[2] * n2[1];
            v[1] = n1[2] * n2[0] - n1[0] * n2[2];
            v[2] = n1[0] * n2[1] - n1[1] * n2[0];

            // normalize direction vector
            float scale = distance / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            v[0] *= scale;
            v[1] *= scale;
            v[2] *= scale;

            // v[0]*= distance;
            // v[1]*= distance;
            // v[2]*= distance;

            // add a vertex into array
            vertices.push_back((v[0]+centerX));
            vertices.push_back((v[1]+ centerY));
            vertices.push_back((v[2]+centerZ));
        }
    }

    return vertices;
}

void drawSphereFace() {
    std::vector<float> vertices = buildUnitPositiveX(6,0.0,0.0,0.0);// subdivision 6
    glBegin(GL_POLYGON);
    for (int i = 0; i < vertices.size(); i += 3) {
        GLfloat x = vertices[i] ;
        GLfloat y = vertices[i + 1];
        GLfloat z = vertices[i + 2];
        glVertex3f(x, y, z);
    }
    glEnd();
}

void drawFullSphere()
{
   // glScalef(sclX, sclX, sclX);
    drawSphereFace();
    distance= (1-scalefactor)/ sqrt(3);
    glPushMatrix();
    glRotatef(90, 0, 1, 0);
    glColor3f(1.0f, 0.5f, 0.0f); 
    glTranslatef(scalefactor, 0,0 ); 
    drawSphereFace();
    glPopMatrix();

    glPushMatrix();
    glRotatef(180, 0, 1, 0);
    glColor3f(0.0f, 1.0f, 0.0f);     // Green  
    glTranslatef(scalefactor, 0,0 ); 
    drawSphereFace();
    glPopMatrix();

    glPushMatrix();
    glRotatef(270, 0, 1, 0);
    glColor3f(1.0f, 0.3f, 0.5f);  
    glTranslatef(scalefactor, 0,0 ); 
    drawSphereFace();
    glPopMatrix();

    glPushMatrix();
    glRotatef(90, 0, 0, 1);
    glColor3f(1.0f, 0.5f, 0.7f);  
    glTranslatef(scalefactor, 0,0 ); 
    drawSphereFace();
    glPopMatrix();

    glPushMatrix();
    glRotatef(-90, 0, 0, 1);
    glColor3f(1.0f, 0.7f, 0.2f);  
    glTranslatef(scalefactor, 0,0 ); 
    drawSphereFace();
    glPopMatrix();
}


// Global variables

// bool isAxes = true, isCube = false, isPyramid = false;

/* Draw axes: X in Red, Y in Green and Z in Blue */
void drawAxes() {
    glLineWidth(3);
    glBegin(GL_LINES);
        glColor3f(1,0,0);   // Red
        // X axis
        glVertex3f(0,0,0);
        glVertex3f(1,0,0);

        glColor3f(0,1,0);   // Green
        // Y axis
        glVertex3f(0,0,0);
        glVertex3f(0,1,0);

        glColor3f(0,0,1);   // Blue
        // Z axis
        glVertex3f(0,0,0);
        glVertex3f(0,0,1);
    glEnd();
}

/* Draw a pyramid centered at the origin */
void drawPyramid() {
    glBegin(GL_TRIANGLES);           // Begin drawing the pyramid with 4 triangles
        // Front
        glColor3f(1.0f, 0.0f, 0.0f);     // Red
        glVertex3f( 0.0f, 1.0f, 0.0f);
        glColor3f(0.0f, 1.0f, 0.0f);     // Green
        glVertex3f(-1.0f, -1.0f, 1.0f);
        glColor3f(0.0f, 0.0f, 1.0f);     // Blue
        glVertex3f(1.0f, -1.0f, 1.0f);

        // Right
        glColor3f(1.0f, 0.0f, 0.0f);     // Red
        glVertex3f(0.0f, 1.0f, 0.0f);
        glColor3f(0.0f, 0.0f, 1.0f);     // Blue
        glVertex3f(1.0f, -1.0f, 1.0f);
        glColor3f(0.0f, 1.0f, 0.0f);     // Green
        glVertex3f(1.0f, -1.0f, -1.0f);

        // Back
        glColor3f(1.0f, 0.0f, 0.0f);     // Red
        glVertex3f(0.0f, 1.0f, 0.0f);
        glColor3f(0.0f, 1.0f, 0.0f);     // Green
        glVertex3f(1.0f, -1.0f, -1.0f);
        glColor3f(0.0f, 0.0f, 1.0f);     // Blue
        glVertex3f(-1.0f, -1.0f, -1.0f);

        // Left
        glColor3f(1.0f,0.0f,0.0f);       // Red
        glVertex3f( 0.0f, 1.0f, 0.0f);
        glColor3f(0.0f,0.0f,1.0f);       // Blue
        glVertex3f(-1.0f,-1.0f,-1.0f);
        glColor3f(0.0f,1.0f,0.0f);       // Green
        glVertex3f(-1.0f,-1.0f, 1.0f);
    glEnd();   // Done drawing the pyramid
}



void drawSingleCylinder()
{

    float cylinderHeight = scalefactor * sqrt(2.0);

    float cylinderRadius =   (1-scalefactor)/ sqrt(3);
    float cylinderCenter =0;


    glPushMatrix();
    glRotated( -45,0,1,0);
    glTranslatef(scalefactor/sqrt(2), 0,0 );
    glRotated( -70.5/2,0,0,1);
    
    drawCylinder(cylinderHeight, cylinderRadius , 100);
    glPopMatrix();

    
}

void draw4cylinders()
{
    drawSingleCylinder();
    glRotated( 90,0,1,0);
        drawSingleCylinder();
    glRotated( 90,0,1,0);
        drawSingleCylinder();
    glRotated( 90,0,1,0);
        drawSingleCylinder();
    glRotated( 90,0,1,0);

}

void drawAllCylinders()
{
    draw4cylinders();
    glRotated( 90,1,0,0);
    draw4cylinders();
    glRotated( 90,0,0,1);
    draw4cylinders();

}

/*  Handler for window-repaint event. Call back when the window first appears and
    whenever the window needs to be re-painted. */
void display() {
    // glClear(GL_COLOR_BUFFER_BIT);            // Clear the color buffer (background)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);             // To operate on Model-View matrix
    glLoadIdentity();                       // Reset the model-view matrix

    // default arguments of gluLookAt
    // gluLookAt(0,0,0, 0,0,-100, 0,1,0);
    
    // control viewing (or camera)
    gluLookAt(eyex,eyey,eyez,
              centerx,centery,centerz,
              upx,upy,upz);
    // draw
   // drawAxes();

    //drawSphereFace();
    drawAllCylinders();
    draw_Octahedron();
    drawFullSphere();
    
    
    // if (isAxes) drawAxes();
    // if (isCube) drawCube();
    // if (isPyramid) drawPyramid();

    glutSwapBuffers();  // Render now
}

/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshapeListener(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0) height = 1;                // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset the projection matrix
    /*if (width >= height) {
        // aspect >= 1, set the height from -1 to 1, with larger width
        gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
    } else {
        // aspect < 1, set the width to -1 to 1, with larger height
        gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
    }*/
    // Enable perspective projection with fovy, aspect, zNear and zFar
    gluPerspective(45.0f, aspect, 0.1f, 100.0f);
}




/* Callback handler for normal-key event */
void keyboardListener(unsigned char key, int x, int y) {

    double v = 0.25;
    double rate = 0.1;
    double s;
    double directionx = centerx- eyex;
    double directiony = centery- eyey;
    double directionz = centerz- eyez;
    double length = sqrt(directionx * directionx + directiony * directiony + directionz * directionz);

    directionx /= length;
    directiony /= length;
    directionz /= length;

    double crossX = directiony * 0 - 1 * directionz;
    double crossY = - directionx * 0 + 0 * directionz;
    double crossZ = directionx * 1 - 0 * directiony;
    double crossUnit =  sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ);


    double rightx = crossX /crossUnit;
    double righty = crossY/crossUnit;
    double rightz= crossZ /crossUnit;


    double pgupX = righty * directionz - rightz * directiony;
    double pgupY = - rightx * directionz + rightz * directionx;
    double pgupZ = rightx * directiony - righty * directionx;

    centerx = directionx + eyex;
    centery = directiony + eyey;
    centerz = directionz + eyez;
 
    switch (key) {
    // Control eye (location of the eye)
    // control eyex
    case '1':
        rightx = rightx * cos(rate) + directionx * sin(rate);
        righty = righty * cos(rate) + directiony * sin(rate);
        rightz = rightz * cos(rate) + directionz * sin(rate);
        
        directionx = directionx * cos(rate) - rightx * sin(rate);
        directiony = directiony * cos(rate) - righty * sin(rate);
        directionz = directionz * cos(rate) - rightz * sin(rate);
        break;
    case '2':
        rightx = rightx * cos(-rate) + directionx * sin(-rate);
        righty = righty * cos(-rate) + directiony * sin(-rate);
        rightz = rightz * cos(-rate) + directionz * sin(-rate);
        
        directionx = directionx * cos(-rate) - rightx * sin(-rate);
        directiony = directiony * cos(-rate) - righty * sin(-rate);
        directionz = directionz * cos(-rate) - rightz * sin(-rate);
        break;

    // control eyey
    case '3':
        eyey += v;
        break;
    case '4':
        eyey -= v;
        break;
    // control eyez
    case '5':
        eyez += v;
        break;
    case '6':
        eyez -= v;
        break;

    case 'a':
        eyex += v * (upy*directionz);
        eyez += v * (-directionx*upy);
        s = sqrt(eyex*eyex + eyez*eyez) / (4 * sqrt(2));
        eyex /= s;
        eyez /= s;
        break;
    case 'd':
        eyex += v * (-upy*directionz);
        eyez += v * (directionx*upy);
        s = sqrt(eyex*eyex + eyez*eyez) / (4 * sqrt(2));
        eyex /= s;
        eyez /= s;
        break;
    case 'w':
        eyey += v;
        break;
    case 's':
        eyey -= v;
        break;

    case '.':
        if(scalefactor<1.0)
            scalefactor+=0.1;
        // cylinderHeight+= (sqrt(2.0))/16.0;
        // cylinderRadius-= 1.0 / (16 * sqrt(3));
        // cylinderCenter+= 1.0/32.0;
  
        break;

    case ',':
        if(scalefactor>0.0)
            scalefactor-=0.1;
        // cylinderHeight-= (sqrt(2.0))/16.0;
        // cylinderRadius += 1.0 / (16 * sqrt(3));
        // cylinderCenter-= 1.0/32.0;
     

         break;
    

    case 27:    // ESC key
        exit(0);    // Exit window
        break;
    }
    glutPostRedisplay();    // Post a paint request to activate display()
}

/* Callback handler for special-key event */
void specialKeyListener(int key, int x,int y) {
    double v = 0.25;
    double rate = 0.1;
    double directionx = centerx- eyex;
    double directiony = centery- eyey;
    double directionz = centerz- eyez;
    double length = sqrt(directionx * directionx + directiony * directiony + directionz * directionz);

    directionx /= length;
    directiony /= length;
    directionz /= length;

    double crossX = directiony * 0 - 1 * directionz;
    double crossY = - directionx * 0 + 0 * directionz;
    double crossZ = directionx * 1 - 0 * directiony;
    double crossUnit =  sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ);


    double rightx = crossX /crossUnit;
    double righty = crossY/crossUnit;
    double rightz= crossZ /crossUnit;


    double pgupX = righty * directionz - rightz * directiony;
    double pgupY = - rightx * directionz + rightz * directionx;
    double pgupZ = rightx * directiony - righty * directionx;

    switch (key) {
    case GLUT_KEY_UP:
        eyex = eyex+ directionx *rate;
        eyey = eyey+ directiony *rate;
        eyez = eyez+ directionz *rate;
        break;
    case GLUT_KEY_DOWN:
        eyex = eyex - directionx *rate;
        eyey = eyey - directiony *rate;
        eyez = eyez - directionz *rate;
        break;

    case GLUT_KEY_RIGHT:
        eyex = eyex - rightx *rate;
        eyey = eyey - righty *rate;
        eyez = eyez - rightz *rate;
        break;

    case GLUT_KEY_LEFT:
        eyex = eyex + rightx *rate;
        eyey = eyey + righty *rate;
        eyez = eyez + rightz *rate;
        break;
    case GLUT_KEY_PAGE_UP:
        eyex = eyex + pgupX *rate;
        eyey = eyey + pgupY *rate;
        eyez = eyez + pgupZ *rate;
        break;

    case GLUT_KEY_PAGE_DOWN:
        eyex = eyex - pgupX *rate;
        eyey = eyey - pgupY *rate;
        eyez = eyez - pgupZ *rate;
        break;
    default:
        return;
    }
    glutPostRedisplay();    // Post a paint request to activate display()
}



/* Main function: GLUT runs as a console application starting at main()  */
int main(int argc, char** argv) {


    COGTriangleX= 1.0/3;
    COGTriangleY= 1.0/3;
    COGTriangleZ= 1.0/3;

    glutInit(&argc, argv);                      // Initialize GLUT
    glutInitWindowSize(640, 640);               // Set the window's initial width & height
    glutInitWindowPosition(50, 50);             // Position the window's initial top-left corner
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color
    glutCreateWindow("OpenGL 3D Drawing");      // Create a window with the given title
    glutDisplayFunc(display);                   // Register display callback handler for window re-paint
    glutReshapeFunc(reshapeListener);           // Register callback handler for window re-shape
    glutKeyboardFunc(keyboardListener);         // Register callback handler for normal-key event
    glutSpecialFunc(specialKeyListener);        // Register callback handler for special-key event
    initGL();                                   // Our own OpenGL initialization
    glutMainLoop();                             // Enter the event-processing loop
    return 0;
}
