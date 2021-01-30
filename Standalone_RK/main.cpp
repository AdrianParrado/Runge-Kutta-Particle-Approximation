#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
using namespace std;

struct ArrStruct
{
    float particles[3][2];
};
//float particleInit(float x_, float y_, float z_);

ArrStruct rungeKutta(ArrStruct var, float h, float b1, float b2, float b3, float u_, float v_, float w_);
ArrStruct getParticles();
void printParticles(ArrStruct var, float T_);

int main()
{
//Parameters
    //Time-Step
    float T = 0;
    float dt = 1;
    float T_max;
    cin >> T_max;
    stringstream ss;
    string T_;

    //RK Coefficients
    float a_21 = 0.333;
    float a_31 = 2;
    float a_32 = -1;

    float b_1 = 0;
    float b_2 = 0.75;
    float b_3 = 0.25;

    //Flow Velocity
    float u = 1;
    float v = 0.5;
    float w = 0;

    /*srand(time(0));
    float x = (float)(rand())/(float)(rand());
    float y = (float)(rand())/(float)(rand());
    float z = (float)(rand())/(float)(rand());*/

    //cout <<"At T=: " << T << ", the particle position is: "<< "x= " << x << ", y= " << y << ", z= " << z << endl;


    //RK Loop
    ArrStruct particles;
    ArrStruct newParticles;
    particles = getParticles();
    printParticles(particles, T);
    newParticles = particles;
    while(T < T_max)
    {
        T += dt;

        newParticles = rungeKutta(newParticles, dt, b_1, b_2, b_3, u, v, w);
        printParticles(newParticles, T);
        /*for(int i = 0; i < 3; ++i)
        {
            for(int j = 0; j < 10; ++j)
            {
                //cout << "particles[" << i << "][" << j << "]: " << particles[i][j] << endl;
                switch(i)
                {
                    case 0:
                        particles[i][j] = rungeKutta(particles[i][j], dt, a_21, a_31, a_32, b_1, b_2, b_3, u);
                        break;
                    case 1:
                        particles[i][j] = rungeKutta(particles[i][j], dt, a_21, a_31, a_32, b_1, b_2, b_3, v);
                        break;
                    case 2:
                        particles[i][j] = rungeKutta(particles[i][j], dt, a_21, a_31, a_32, b_1, b_2, b_3, w);
                        break;
                }
                cout << "new particles[" << i << "][" << j << "]: " << particles[i][j] << endl;
                //cout << "particles[" << i << "][" << j << "]: " << particles[i][j] << endl;
            }
        }*/


        /*x = rungeKutta(x, dt, a_21, a_31, a_32, b_1, b_2, b_3, u);
        y = rungeKutta(y, dt, a_21, a_31, a_32, b_1, b_2, b_3, v);
        z = rungeKutta(z, dt, a_21, a_31, a_32, b_1, b_2, b_3, w);*/

        //cout << endl << "At T= " << T << ", the particle position is: ";
        //cout << "x= " << x << ", y= " << y << ", z= " << z << endl;
    }
    return 0;
}

/*float particleInit(float x_, float y_, float z_)
{
    srand(time(0));

    x_ = (float)(rand())/(float)(rand());
    y_ = (float)(rand())/(float)(rand());
    z_ = (float)(rand())/(float)(rand());
}
*/
ArrStruct getParticles()
{
    //Particle Initialisation
    srand(time(0));

    ArrStruct var;

    for(int i = 0; i < 3; ++i)
    {
        for(int j = 0; j < 2; ++j)
        {
             var.particles[i][j] = (float)rand() / (float)rand();
        }
    }

    return var;
}

void printParticles(ArrStruct var, float T_)
{
    //Create, fill and close
        stringstream ss;
        ss << T_;
        string Ts;
        ss >> Ts;

        cout << "Creating File "  << Ts << endl;
        ofstream fileOut;
        fileOut.open("T= "+Ts+".txt");
        fileOut << "x,y,z" << endl;

        for(int j = 0; j < 2; ++j)
        {
            fileOut << var.particles[0][j] << "," << var.particles[1][j] << "," << var.particles[2][j] << endl;
        }
        fileOut.close();

        cout << "File Created " << T_ << endl;
}

ArrStruct rungeKutta(ArrStruct var, float h, float b1, float b2, float b3, float u_, float v_, float w_)
{

    //RK Substeps
     /* float Y1 = y_;
        float Y2 = y_ + h*(a21*u_);
        float Y3 = y_ + h*(a31*u_ + a32*u_);
     */
    //Calculate new value of y
        //y_ = y_ + h*(b1*u_ + b2*u_ + b3*u_);

    for(int i = 0; i < 3; ++i)
        {
            for(int j = 0; j < 2; ++j)
            {
                //cout << "particles[" << i << "][" << j << "]: " << particles[i][j] << endl;
                switch(i)
                {
                    case 0:
                        var.particles[i][j] += h*(b1*u_ + b2*u_ +b3*u_);
                        //particles[i][j] = rungeKutta(particles[i][j], dt, a_21, a_31, a_32, b_1, b_2, b_3, u);
                        break;
                    case 1:
                        var.particles[i][j] += h*(b1*v_ + b2*v_ +b3*v_);
                        //particles[i][j] = rungeKutta(particles[i][j], dt, a_21, a_31, a_32, b_1, b_2, b_3, v);
                        break;
                    case 2:
                        var.particles[i][j] += h*(b1*w_ + b2*w_ +b3*w_);
                        //particles[i][j] = rungeKutta(particles[i][j], dt, a_21, a_31, a_32, b_1, b_2, b_3, w);
                        break;
                }
                //cout << "new particles[" << i << "][" << j << "]: " << particles[i][j] << endl;
                //cout << "particles[" << i << "][" << j << "]: " << particles[i][j] << endl;
            }
        }

    return var;
}
