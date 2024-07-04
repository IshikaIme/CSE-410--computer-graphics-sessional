#include<bits/stdc++.h>
using namespace std;
#define PI acos(-1)

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



class Triangle
{

    public:
    Point ver[3];
    Triangle( Point v1, Point v2, Point v3)
    {
        ver[0] = v1;
        ver[1] = v2;
        ver[2]= v3;
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

  //  int trianglesNumber = 0;

    //stage 1 starts here 

    while(true)
    {
        string command;
        fin>>command;

        if (command== "triangle")
            {
               // trianglesNumber++;
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





}





