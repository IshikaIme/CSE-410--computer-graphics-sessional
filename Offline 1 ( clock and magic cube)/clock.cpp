#include <windows.h>  // For MS Windows
#include <GL/glut.h>  // GLUT, includes glu.h and gl.h
#include <cmath>

float minuteAngle = M_PI / 30.0f;
float t = 0.0, theta_t;
float minuteStart = 0.9f, minuteEnd = 1.0f, secondAngle = 0, clockRadius = 45.0f, hourAngle = 0, secAngleHand = 0, minAngleHand = 0, hourAngleHand = 0;
float stepStart = 0.8f, stepEnd = 1.0f;
float clockInnerRadius = 10.0f, pendulumRadius = 5.0f, pendulumHand = 6.0f;
float pendulumLength = 40.0f;  // Length of the pendulum arm

void init()
{
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);  // Set background color to black
    glMatrixMode(GL_PROJECTION);
    glOrtho(-150, 150, -100, 100, -100, 100);
}

void timer(int value)
{
    glutPostRedisplay();
    secAngleHand += 6;
    minAngleHand += 6.0 / 60;
    hourAngleHand += 6.0 / 3600;

    glutTimerFunc(1000, timer, 0);
}

void timerPendulum(int value)
{
    glutPostRedisplay();
    float omega = 2.0 * M_PI / 2.0;  // Angular frequency for 2-second period
    theta_t = 45 * cos(omega * t);   // initial angle 45 deegree
    t += 0.01;

    glutTimerFunc(10, timerPendulum, 0);
}

void MinutesDesign()
{
    glLineWidth(2);
    glEnable(GL_LINE_SMOOTH);
    glBegin(GL_LINES);
    for (int i = 0; i < 60; i++)
    {
        float c = cos(i * minuteAngle), s = sin(i * minuteAngle);
        if (i % 5 == 0) {
            glColor3f(1.0f, 1.0f, 1.0f);
            glVertex2f(clockRadius * minuteStart * c, clockRadius * minuteStart * s);
            glVertex2f(clockRadius * minuteEnd * c, clockRadius * minuteEnd * s);
        }
        else {
            glColor3f(0.0f, 1.0f, 0.0f);
            glVertex2f(clockRadius * stepStart * c, clockRadius * stepStart * s);
            glVertex2f(clockRadius * stepEnd * c, clockRadius * stepEnd * s);
        }
    }
    glEnd();
}

void handsDesign()
{
    // Second hand
    glColor4f(1.0f, 1.0f, 0.0f, 1.0f); // yellow

    glLineWidth(3);
    glPushMatrix();
    glTranslatef(-1, 3, 0);
    glPushMatrix();
    glEnable(GL_LINE_SMOOTH);
    glRotatef(-secAngleHand, 0, 0, 1);
    glPushMatrix();
    glTranslatef(1, -3, 0);

    glBegin(GL_LINES);
    {
        glVertex2i(-1, 3);
        glVertex2i(20, -22);
    }
    glEnd();

    glPopMatrix();
    glPopMatrix();
    glPopMatrix();

    // Minute hand
    glColor4f(1.0f, 0.5f, 0.0f, 1.0f); // orange/brown

    glLineWidth(2);
    glPushMatrix();
    glTranslatef(-1, 3, 0);
    glPushMatrix();
    glEnable(GL_LINE_SMOOTH);
    glRotatef(-minAngleHand, 0, 0, 1);
    glPushMatrix();
    glTranslatef(1, -3, 0);

    glBegin(GL_LINES);
    {
        glVertex2i(-1, 3);
        glVertex2i(10, -15);
    }
    glEnd();

    glPopMatrix();
    glPopMatrix();
    glPopMatrix();

    // Hour hand
    glColor4f(1.0f, 0.0f, 0.0f, 1.0f); // red

    glPushMatrix();
    glTranslatef(-1, 3, 0);
    glPushMatrix();
    glEnable(GL_LINE_SMOOTH);
    glRotatef(-hourAngleHand, 0, 0, 1);
    glPushMatrix();
    glTranslatef(1, -3, 0);

    glBegin(GL_LINES);
    {
        glVertex2i(-1, 3);
        glVertex2i(3, -8);
    }
    glEnd();

    glPopMatrix();
    glPopMatrix();
    glPopMatrix();
}


void CaseDesign()
{
    float x, y, angle;

    glColor4f(0.0f, 1.0f, 1.0f, 0.0f);  // Cyan

    glPointSize(2.5);
    glBegin(GL_POINTS);
    for (angle = 0.0f; angle <= 2.0 * M_PI * 100; angle += 0.1f)
    {
        x = clockRadius * cos(angle);
        y = clockRadius * sin(angle);
        glVertex2f(x, y);
    }
    glEnd();
}

void PendulumDesign()
{
    glColor3f(1.0f, 1.0f, 1.0f);  // White

    glPushMatrix();
    glTranslatef(0, -(pendulumLength + 5), 0);
    glRotatef(theta_t, 0, 0, 1);

    glBegin(GL_LINES);
    {
        glVertex2i(0, 0);
        glVertex2i(0, -pendulumLength);
    }
    glEnd();

    glColor3f(1.0f, 1.0f, 0.0f);  // Yellow
    glTranslatef(0, -pendulumLength, 0);
    glBegin(GL_POLYGON);
    for (float theta = 0; theta < 360; theta += 10)
    {
        float x = pendulumRadius * cos(theta / 180.0 * M_PI);
        float y = pendulumRadius * sin(theta / 180.0 * M_PI);
        glVertex2f(x, y);
    }
    glEnd();
    // glColor3f(1.0f, 0.0f, 0.0f);  // Red
    // glTranslatef(0, -pendulumLength, 0);
    // glBegin(GL_POLYGON);
    // for (float theta = 0; theta < 360; theta += 10)
    // {
    //     float x = pendulumRadius*cos(theta / 180.0 * M_PI);
    //     float y = pendulumRadius * sin(theta / 180.0 * M_PI);
    //     glVertex2f(x, y);
    // }
    // glEnd();


    // glColor3f(1.0f, 0.0f, 0.0f);  // Red
    // glTranslatef(0, -pendulumRadius, 0);
    // glBegin(GL_LINES);
    // {
    //     glVertex2i(0, 0);
    //     glVertex2i(0, -pendulumHand);
    // }
    // glEnd();

    glPopMatrix();
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    CaseDesign();
    MinutesDesign();
    handsDesign();
    PendulumDesign();
    glFlush();
}

int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(1080 / 1.5, 780 / 1.5);
    glutInitWindowPosition(100, 10);
    glutCreateWindow("OpenGL 2D Drawing");
    init();
    glutTimerFunc(0, timer, 0);
    glutTimerFunc(0, timerPendulum, 0);
    glutDisplayFunc(display);
    glutMainLoop();
    return 0;
}
