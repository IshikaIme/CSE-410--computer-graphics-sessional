#include<bits/stdc++.h>
#include "bitmap_image.hpp"
using namespace std;
#define PI acos(-1)
static unsigned long int g_seed = 1;

inline int random()
{
 g_seed = (214013 * g_seed + 2531011);
 return (g_seed >> 16) & 0x7FFF;
}

class Point
{

    public:
    double x,y,z,n;
    Point()
    {
        x = y = z = 0.0;
        n = 1.0;
    }

	Point(double a, double b, double c)
    {
        x=a;
        y=b;
        z=c;
        n= 1.0;
    }

    Point(double a, double b, double c, double n)
    {
        x= a;
        y= b;
        z= c;
        n= n;
    }


    Point(const Point &p)
    {
        x= p.x;
        y= p.y;
        z= p.z;
        n= p.n;

    }

    Point normalize()
    {
        double length = sqrt(x*x+y*y+z*z);
        x= x/length;
        y= y/length;
        z= z/length;
        Point p(x,y,z);
        return p;

    }



    Point operator+ (Point b)
    {

        double sumx, sumy, sumz;
        sumx= x + b.x;
        sumy = y + b.y;
        sumz= z + b.z;
        Point sum( sumx, sumy, sumz);
        return sum;
    }


    Point operator- (Point b)
    {

        double rx, ry, rz;
        rx= x - b.x;
        ry = y - b.y;
        rz= z - b.z;
        Point r( rx, ry, rz);
        return r;
    }

    double operator* (Point a)
    {

        double rx, ry, rz;
        rx= x *  a.x;
        ry = y *  a.y;
        rz= z*  a.z;
        ;
        return rx + ry + rz;
    }


    Point operator* (double a)
    {

        double rx, ry, rz;
        rx= x * a;
        ry = y * a;
        rz= z * a;
        Point r( rx, ry, rz);
        return r;
    }

    Point  X(Point b)  {
        Point temp(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x);
        return temp;
    }


    string print(Point a)  {
        string str="";
        str=to_string(a.x);;
        return str;

    }

    friend ostream& operator<<(ostream& os, const Point& p);
    friend istream& operator>>(istream &in, Point &p);


};

Point crossProduct(Point p1, Point p2) 
{
    double cx = p1.y * p2.z - p1.z * p2.y;
    double cy = p1.z * p2.x - p1.x * p2.z;
    double cz = p1.x * p2.y - p1.y * p2.x;
    return Point(cx, cy, cz);
}

class Triangle
{

    public:
    Point points[3];
    int rgb[3];

    Triangle()
    {
      
        rgb[0] = random()%256;
        rgb[1] = random()%256;
        rgb[2] = random()%256;
    }

    Triangle( Point v1, Point v2, Point v3)
    {
        points[0] = v1;
        points[1] = v2;
        points[2]= v3;
        rgb[0] = random()%256;
        rgb[1] = random()%256;
        rgb[2] = random()%256;
    }

};

class Matrix
{
    public:
    vector<vector<double> > ara;
    int dimension;

    Matrix()
    {
        ara.resize(4,vector<double>(4,0));
        dimension = 4;
    }

    Matrix(int dim)
    {
        ara.resize(dim,vector<double>(dim,0));
        dimension = dim;
    }



    void IdentityMatrix()
    {
        for (int i=0; i<dimension; i++)
        {
            for (int j=0; j<dimension; j++)
                {
                    if(i==j)
                        ara[i][j]=1;
                    else
                        ara[i][j]=0;
                }
        }
    }

    Matrix Multiplication( Matrix m2)
    {
        Matrix result;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                for (int k = 0; k < dimension; ++k) {
                   result.ara[i][j] += ara[i][k] * m2.ara[k][j];
            }
        }
    }

    return result;
    }


    Point Multiplication(Point point)
    {
        Point result;
        result.x = ara[0][0] * point.x + ara[0][1] * point.y + ara[0][2] * point.z + ara[0][3] * point.n;
        result.y = ara[1][0] * point.x + ara[1][1] * point.y + ara[1][2] * point.z + ara[1][3] * point.n;
        result.z = ara[2][0] * point.x + ara[2][1] * point.y + ara[2][2] * point.z + ara[2][3] * point.n;
        result.n = ara[3][0] * point.x + ara[3][1] * point.y + ara[3][2] * point.z + ara[3][3] * point.n;
        result.x = result.x / result.n;
        result.y = result.y / result.n;
        result.z = result.z / result.n;
        result.n = 1;
        return result;
    }

      Point R(Point x, Point a, double theta) 
    {   
        double rad = theta*(PI)/180;
        Point ret = x*cos(rad)+a*(a*x)*(1-cos(rad))+(a.X(x))*sin(rad);
        return ret;
    }

};

 ostream& operator<<(ostream& os, const Point& p)
{
    os << p.x << ' ' << p.y << ' ' << p.z ;
    return os;
}

istream& operator>>(istream &in, Point &p)
{
    in >> p.x >> p.y >> p.z;
    return in;
}

double getValueZ( Triangle tr, double x, double y)
{
    Point normal_vector = crossProduct(tr.points[1]- tr.points[0], tr.points[2]- tr.points[0]);
    double d= - (normal_vector.x * tr.points[0].x + normal_vector.y * tr.points[0].y +normal_vector.z * tr.points[0].z) ;
    double valueZ= (-d - normal_vector.x * x- normal_vector.y * y) / (normal_vector.z);
    return valueZ; 
}


int main()
{
    ifstream fin("scene.txt");
    ofstream foutStage1("stage1.txt");
    ofstream foutStage2("stage2.txt");
    ofstream foutStage3("stage3.txt");
    
    Point eye, look, up;
    double fovX, rr, tt, fovY, aspectRatio, near, far;

    cout << std::fixed << std::setprecision(7);
    foutStage1 << std::fixed << std::setprecision(7);
    foutStage2 << std::fixed << std::setprecision(7);
    foutStage3 << std::fixed << std::setprecision(7);

    fin>> eye >> look >> up;
    //cout<< eye <<" "<<look<<" "<< up;

    fin >> fovY >> aspectRatio >> near >> far;


    //############################################
    //stage 2 starts here
    Matrix viewMatrix;
    viewMatrix.IdentityMatrix();
    Point l = look - eye;
    l.normalize();
    Point r = l.X(up);
    r.normalize();
    Point u = r.X(l);

    viewMatrix.ara[0][0]=r.x;
    viewMatrix.ara[0][1]=r.y;
    viewMatrix.ara[0][2]=r.z;
    viewMatrix.ara[1][0]=u.x;
    viewMatrix.ara[1][1]=u.y;
    viewMatrix.ara[1][2]=u.z;
    viewMatrix.ara[2][0]=-l.x;
    viewMatrix.ara[2][1]=-l.y;
    viewMatrix.ara[2][2]=-l.z;
    
    Matrix trans;
    trans.IdentityMatrix();
    trans.ara[0][3] = -eye.x;
    trans.ara[1][3] = -eye.y;
    trans.ara[2][3] = -eye.z;
    viewMatrix= trans.Multiplication(viewMatrix);
    

    //###################################################


    //########################################################
    // stage 3: projection matrix starts here
    fovX = fovY * aspectRatio;
    tt = near * tan((fovY/2)*(PI/180));
    rr = near * tan((fovX/2)*(PI/180));
    Matrix projectionMatrix;
    projectionMatrix.IdentityMatrix();
    projectionMatrix.ara[0][0] = near/rr;
    projectionMatrix.ara[1][1] = near/tt;
    projectionMatrix.ara[2][2] = -(far+near)/(far-near);
    projectionMatrix.ara[3][2] = -1;
    projectionMatrix.ara[2][3] = -2.0*far*near/(far-near);
    projectionMatrix.ara[3][3] = 0;
    

    Matrix idMat(4);

    idMat.IdentityMatrix();

    stack<Matrix>states;
    states.push(idMat);

    int trianglesNumber = 0;





    //stage 1 starts here 

    while(true)
    {
        string command;
        fin>>command;

        if (command== "triangle")
            {
                trianglesNumber++;
                Point point[3], viewPoint[3], projectionPoint[3];
                fin>> point[0] >> point[1] >> point[2];
                for(int i=0; i< 3; i++)
                {
                    point[i]= states.top().Multiplication(point[i]);
                    viewPoint[i] = viewMatrix.Multiplication(point[i]);
                    projectionPoint[i] = projectionMatrix.Multiplication(viewPoint[i]);

                    foutStage1<< point[i]<<endl;
                    foutStage2<< viewPoint[i]<<endl;
                    foutStage3<< projectionPoint[i]<<endl;


                }
            foutStage1<<endl;
            foutStage2<<endl;
            foutStage3<<endl;
            }

        else if (command== "translate")
            {
                double tx,ty,tz;
                fin>> tx>>ty >>tz;
                Matrix translateMatrix;
                translateMatrix.IdentityMatrix();
                translateMatrix.ara[0][3] = tx;
                translateMatrix.ara[1][3] = ty;
                translateMatrix.ara[2][3] = tz;
                Matrix res = (states.top()).Multiplication(translateMatrix);
                states.pop();
                states.push(res);
                
            }

        else if (command=="scale")
            {
                double sx, sy, sz;
                fin >>sx>>sy>>sz;
                Matrix transformMatrix;
                transformMatrix.IdentityMatrix();
                transformMatrix.ara[0][0] = sx;
                transformMatrix.ara[1][1] = sy;
                transformMatrix.ara[2][2] = sz;
                Matrix res = (states.top()).Multiplication(transformMatrix);
                states.pop();
                states.push(res);
               
            }

        else if (command=="rotate")
            {
                double ax,ay,az;
                double angle;
                fin>> angle >> ax >> ay >> az;
                Point p( ax,ay,az);
                Matrix Rot;
                Rot.IdentityMatrix();
                Point i(1,0,0);
                Point j(0,1,0);
                Point k(0,0,1);
                
                p.normalize();
                

                Point c1= Rot.R(i, p, angle);
                Point c2= Rot.R(j, p, angle);
                Point c3= Rot.R(k, p, angle);
            
                Rot.ara[0][0] = c1.x;
                Rot.ara[1][0] = c1.y;
                Rot.ara[2][0] = c1.z;

                Rot.ara[0][1] = c2.x;
                Rot.ara[1][1] = c2.y;
                Rot.ara[2][1] = c2.z;

                Rot.ara[0][2] = c3.x;
                Rot.ara[1][2] = c3.y;
                Rot.ara[2][2] = c3.z;

                Matrix mat = states.top().Multiplication(Rot);
                states.pop();
                states.push(mat);
                
            }

        else if (command=="push")
            {
                states.push(states.top());
                
            }

        else if (command=="pop")
            {
                if(states.empty())
                    cout<<"Stack is empty, cant pop error";
                states.pop();
                
            }
        else if (command=="end")
                break;

        else
            {
                cout<<"Invalid command"<<endl;
            }

    }

    fin.close();
    foutStage1.close();
    foutStage2.close();
    foutStage3.close();

    
    //########################stage 4 ######################################
    fin.open("config.txt");
    ofstream fout;
    fout.open("z_buffer.txt");

    int screenWidth, screenHeight;
    fin>>screenWidth>>screenHeight;

    double minX,maxX, minY, maxY;

    //calculating the left limit, right limit, top limit, bottom limit, dx, dy
    

    double dx = 2.0/screenWidth;
    double dy = 2.0/screenHeight;

    double topY= 1.0 - dy/2;
    double bottomY = (-1.0) + dy/2;
    double leftX = (-1.0) + dx/2;
    double rightX = 1.0  - dx/2;

    double zMax= 1.0;
    double zMin= -1.0;
   

    vector<vector<double>> zBuffer(screenHeight,vector<double>(screenWidth,zMax));;
    bitmap_image image(screenWidth, screenHeight);
    for(int i=0;i<screenWidth;i++)
    {
        for(int j=0;j<screenHeight;j++) image.set_pixel(i,j,0,0,0);
        
    }


    fin.close();

    fin.open("stage3.txt");
    
    for(int t=0; t<trianglesNumber;t++)
    {
        Point p1,p2,p3;
        fin>>p1>>p2>>p3;

        Triangle triangle;
        triangle.points[0] = p1;
        triangle.points[1] = p2;
        triangle.points[2] = p3;

       
        // calculate min max values of triangles
        double minX = min(p1.x,min(p2.x,p3.x));
        double maxX = max(p1.x,max(p2.x,p3.x));
        double minY = min(p1.y,min(p2.y,p3.y));
        double maxY = max(p1.y,max(p2.y,p3.y));

        //clipping 
        minX = max ( minX, leftX) ;
        maxX= min (maxX, rightX);

        minY= max( minY , bottomY);
        maxY= min ( maxY, topY);

        int bottomScanLineY = round((topY- minY)/dy);
        int topScanLineY = round ((topY - maxY)/ dy);

        int scanY;
        // scanning from top to bottom
        for(scanY =topScanLineY; scanY <= bottomScanLineY ; scanY++)
        {  
            double disY = topY- scanY * dy;
            vector<double> intersectionPointX(2), intersectionPointZ(2);
            int index=0;;
            int iterate;
            for(iterate=0; iterate<3; iterate++)
            {
                int next ;
                next = (iterate+1) %3;
 
                if (triangle.points[iterate].y == triangle.points[next].y) 
                {
                    //if one side of the triangle lies on the scanLine then ignore
                    continue;
                }

                //cout<<"Hii"<<endl;

                //checking if any point intersects with the scanline in triangle (can be highest 2 points)
                double minYAtScanLine = min ( triangle.points[iterate].y , triangle.points[next].y);
                double maxYAtScanLine = max ( triangle.points[iterate].y , triangle.points[next].y);

                if(disY >= minYAtScanLine && disY <= maxYAtScanLine)
                {
                    // x= x1 + (x1-x2) * ((y-y1)/(y1-y2))
                    intersectionPointX[index] =triangle.points[iterate].x + (triangle.points[iterate].x - triangle.points[next].x) *((disY - triangle.points[iterate].y ) /(triangle.points[iterate].y - triangle.points[next].y));
                    index++;

                }
            }

            vector<double> newIntersectionPointX(2);
            newIntersectionPointX = intersectionPointX;

                //cliping on x 
            for (int i=0; i<2; i++ )
            {
                if(minX> intersectionPointX[i] )
                    intersectionPointX[i] = minX;
                else if (maxX< intersectionPointX[i] )
                    intersectionPointX[i] = maxX;

            }
                

            intersectionPointZ[0]= intersectionPointZ[1] - (intersectionPointZ[1] - intersectionPointZ[0]) * ((newIntersectionPointX[1] - intersectionPointX[0])/ (newIntersectionPointX[1] - newIntersectionPointX[0]) );
            intersectionPointZ[1]= intersectionPointZ[1] - (intersectionPointZ[1] - intersectionPointZ[0]) * ((newIntersectionPointX[1] - intersectionPointX[1])/ (newIntersectionPointX[1] - newIntersectionPointX[0]) );

            // a is left point and b is right point
            double xa= intersectionPointX[0];
            double za= intersectionPointZ[0];
            
            double xb= intersectionPointX[1];
            double zb= intersectionPointZ[1];

            if(intersectionPointX[0]>= intersectionPointX[1])
            {
                swap(xa,xb);
                swap(za,zb);
                swap(intersectionPointX[0], intersectionPointX[1]);
            }

            int startX = round((xa-leftX)/dx);
            int endX = round((xb-leftX)/dx);
            
            for(int scanX =startX; scanX<= endX; scanX++)
            {
                double currX= leftX + scanX *dx;
                double currY = disY;
                double valuez= getValueZ(triangle, currX, currY);
                
                if(valuez>= zMin && valuez < zBuffer[scanY][scanX])
                {
                    zBuffer[scanY][scanX] = valuez;
                    image.set_pixel(scanX,scanY,triangle.rgb[0],triangle.rgb[1],triangle.rgb[2]);
                }
            }
        }
    }
    

     for (int i = 0; i < screenHeight; i++) 
    {
        for (int j = 0; j < screenWidth; j++)
        {
            if (zBuffer[i][j] < zMax) fout << setprecision(6) << fixed << zBuffer[i][j] << "\t";
        }
        fout << endl;
    }
    fin.close();
    fout.close();
    

    image.save_image("out.bmp");

    // a. Free image memory
    image.clear();  // This clears the memory used by the image

    // b. Free z-buffer memory
    for (int i = 0; i < zBuffer.size(); ++i)
    {
        zBuffer[i].clear();
    }
    zBuffer.clear();  // This clears the memory used by the z-buffer

    return 0;



}





