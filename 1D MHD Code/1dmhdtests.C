//#include<iostream>
//#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<array>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<sstream>
#include "1dmhdtests.h"
//------------READ PARAMETERS FROM SETTINGS FILE------------------------
//----------------------------------------------------------------------
void ReadParams(int &test, int &limiter, int &ncells, double &x0, double &x1, double &tstop, double & C, double &gamma, double & omega){
    
    std::string temp;
    std::ifstream file;
    
    //Change name of settings file
    file.open("inputs.txt");
    
    std::getline(file, temp);
    std::istringstream input0(temp);
    input0 >> test;
    
    std::getline(file, temp);
    std::istringstream input1(temp);
    input1 >> limiter;
    
    std::getline(file, temp);
    std::istringstream input2(temp);
    input2 >> ncells;
    
    std::getline(file, temp);
    std::istringstream input3(temp);
    input3 >> x0;
    
    std::getline(file, temp);
    std::istringstream input4(temp);
    input4 >> x1;
    
    std::getline(file, temp);
    std::istringstream input5(temp);
    input5 >> tstop;
    
    std::getline(file, temp);
    std::istringstream input6(temp);
    input6 >> C;
    
    std::getline(file, temp);
    std::istringstream input7(temp);
    input7 >> gamma;
    
    std::getline(file, temp);
    std::istringstream input8(temp);
    input8 >> omega;
    
}
//--------------------------TESTS---------------------------------------
//----------------------------------------------------------------------
void test_modifiedSod(double &x0,double &dx, std::vector<std::array<double,8> >& U_P){
    double rho_L, vx_L, vy_L, vz_L, p_L, Bx_L, By_L, Bz_L;
    double rho_R, vx_R, vy_R, vz_R, p_R, Bx_R, By_R, Bz_R;
    
    // Left initial conditions
    rho_L = 1.0;
    vx_L = 0.75;
    vy_L = 0.0;
    vz_L = 0.0;
    p_L = 1.0;
    Bx_L = 0.0;
    By_L = 0.0;
    Bz_L = 0.0;
    
    // Right initial conditions
    rho_R = 0.125;
    vx_R = 0.0;
    vy_R = 0.0;
    vz_R = 0.0;
    p_R = 0.1;
    Bx_R = 0.0;
    By_R = 0.0;
    Bz_R = 0.0;
    
    for(int i=0; i<int(U_P.size()); ++i){
        double x = x0 + (i-1.5)*dx;
            
        if (x<=0.3){
            U_P[i][0]=rho_L;
            U_P[i][1]=vx_L;
            U_P[i][2]=vy_L;
            U_P[i][3]=vz_L;
            U_P[i][4]=p_L;
            U_P[i][5]=Bx_L;
            U_P[i][6]=By_L;
            U_P[i][7]=Bz_L;
        }
        
        else{
            U_P[i][0]=rho_R;
            U_P[i][1]=vx_R;
            U_P[i][2]=vy_R;
            U_P[i][3]=vz_R;
            U_P[i][4]=p_R;
            U_P[i][5]=Bx_R;
            U_P[i][6]=By_R;
            U_P[i][7]=Bz_R;
        }
    }
}

void test_BrioWu(double &x0,double &dx,std::vector<std::array<double,8> >& U_P){

    double rho_L, vx_L, vy_L, vz_L, p_L, Bx_L, By_L, Bz_L;
    double rho_R, vx_R, vy_R, vz_R, p_R, Bx_R, By_R, Bz_R;

    // Left initial conditions
    rho_L = 1.0;
    vx_L = 0.0;
    vy_L = 0.0;
    vz_L = 0.0;
    p_L = 1.0;
    Bx_L = 0.75;
    By_L = 1.0;
    Bz_L = 0.0;
    
    // Right initial conditions
    rho_R = 0.125;
    vx_R = 0.0;
    vy_R = 0.0;
    vz_R = 0.0;
    p_R = 0.1;
    Bx_R = 0.75;
    By_R = -1.0;
    Bz_R = 0.0;
    
    for(int i=0; i<int(U_P.size()); ++i){
        double x = x0 + (i-1.5)*dx;
            
        if (x<=400.0){
            U_P[i][0]=rho_L;
            U_P[i][1]=vx_L;
            U_P[i][2]=vy_L;
            U_P[i][3]=vz_L;
            U_P[i][4]=p_L;
            U_P[i][5]=Bx_L;
            U_P[i][6]=By_L;
            U_P[i][7]=Bz_L;
        }
        
        else{
            U_P[i][0]=rho_R;
            U_P[i][1]=vx_R;
            U_P[i][2]=vy_R;
            U_P[i][3]=vz_R;
            U_P[i][4]=p_R;
            U_P[i][5]=Bx_R;
            U_P[i][6]=By_R;
            U_P[i][7]=Bz_R;
        }
    }
}

void test_DaiWoodward(double &x0,double &dx, std::vector<std::array<double,8> >& U_P){
    double rho_L, vx_L, vy_L, vz_L, p_L, Bx_L, By_L, Bz_L;
    double rho_R, vx_R, vy_R, vz_R, p_R, Bx_R, By_R, Bz_R;
    
    // Left initial conditions
    rho_L = 1.0;
    vx_L = 10.0;
    vy_L = 0.0;
    vz_L = 0.0;
    p_L = 20.0;
    Bx_L = 5.0/sqrt(4.0*M_PI);
    By_L = 5.0/sqrt(4.0*M_PI);
    Bz_L = 0.0/sqrt(4.0*M_PI);
    
    // Right initial conditions
    rho_R = 1.0;
    vx_R = -10.0;
    vy_R = 0.0;
    vz_R = 0.0;
    p_R = 1.0;
    Bx_R = 5.0/sqrt(4.0*M_PI);
    By_R = 5.0/sqrt(4.0*M_PI);
    Bz_R = 0.0/sqrt(4.0*M_PI);
    
    for(int i=0; i<int(U_P.size()); ++i){
        double x = x0 + (i-1.5)*dx;
            
        if (x<=0.5){
            U_P[i][0]=rho_L;
            U_P[i][1]=vx_L;
            U_P[i][2]=vy_L;
            U_P[i][3]=vz_L;
            U_P[i][4]=p_L;
            U_P[i][5]=Bx_L;
            U_P[i][6]=By_L;
            U_P[i][7]=Bz_L;
        }
        
        else{
            U_P[i][0]=rho_R;
            U_P[i][1]=vx_R;
            U_P[i][2]=vy_R;
            U_P[i][3]=vz_R;
            U_P[i][4]=p_R;
            U_P[i][5]=Bx_R;
            U_P[i][6]=By_R;
            U_P[i][7]=Bz_R;
        }
    }
}

void test_BrioWuSwitched(double &x0,double &dx, std::vector<std::array<double,8> >& U_P){
    double rho_L, vx_L, vy_L, vz_L, p_L, Bx_L, By_L, Bz_L;
    double rho_R, vx_R, vy_R, vz_R, p_R, Bx_R, By_R, Bz_R;

    // Left initial conditions
    rho_L = 1.0;
    vx_L = 0.0;
    vy_L = 0.0;
    vz_L = 0.0;
    p_L = 1.0;
    Bx_L = 0.75;
    By_L = 0.0;
    Bz_L = 1.0;
    
    // Right initial conditions
    rho_R = 0.125;
    vx_R = 0.0;
    vy_R = 0.0;
    vz_R = 0.0;
    p_R = 0.1;
    Bx_R = 0.75;
    By_R = 0.0;
    Bz_R = -1.0;
    
    for(int i=0; i<int(U_P.size()); ++i){
        double x = x0 + (i-1.5)*dx;
            
        if (x<=400.0){
            U_P[i][0]=rho_L;
            U_P[i][1]=vx_L;
            U_P[i][2]=vy_L;
            U_P[i][3]=vz_L;
            U_P[i][4]=p_L;
            U_P[i][5]=Bx_L;
            U_P[i][6]=By_L;
            U_P[i][7]=Bz_L;
        }
        
        else{
            U_P[i][0]=rho_R;
            U_P[i][1]=vx_R;
            U_P[i][2]=vy_R;
            U_P[i][3]=vz_R;
            U_P[i][4]=p_R;
            U_P[i][5]=Bx_R;
            U_P[i][6]=By_R;
            U_P[i][7]=Bz_R;
        }
    }
}
void test_DaiWoodwardSwitched(double &x0,double &dx, std::vector<std::array<double,8> >& U_P){
    double rho_L, vx_L, vy_L, vz_L, p_L, Bx_L, By_L, Bz_L;
    double rho_R, vx_R, vy_R, vz_R, p_R, Bx_R, By_R, Bz_R;
    
    // Left initial conditions
    rho_L = 1.0;
    vx_L = 10.0;
    vy_L = 0.0;
    vz_L = 0.0;
    p_L = 20.0;
    Bx_L = 5.0/sqrt(4.0*M_PI);
    By_L = 0.0/sqrt(4.0*M_PI);
    Bz_L = 5.0/sqrt(4.0*M_PI);
    
    // Right initial conditions
    rho_R = 1.0;
    vx_R = -10.0;
    vy_R = 0.0;
    vz_R = 0.0;
    p_R = 1.0;
    Bx_R = 5.0/sqrt(4.0*M_PI);
    By_R = 0.0/sqrt(4.0*M_PI);
    Bz_R = 5.0/sqrt(4.0*M_PI);
    
    for(int i=0; i<int(U_P.size()); ++i){
        double x = x0 + (i-1.5)*dx;
            
        if (x<=0.5){
            U_P[i][0]=rho_L;
            U_P[i][1]=vx_L;
            U_P[i][2]=vy_L;
            U_P[i][3]=vz_L;
            U_P[i][4]=p_L;
            U_P[i][5]=Bx_L;
            U_P[i][6]=By_L;
            U_P[i][7]=Bz_L;
        }
        
        else{
            U_P[i][0]=rho_R;
            U_P[i][1]=vx_R;
            U_P[i][2]=vy_R;
            U_P[i][3]=vz_R;
            U_P[i][4]=p_R;
            U_P[i][5]=Bx_R;
            U_P[i][6]=By_R;
            U_P[i][7]=Bz_R;
        }
    }
}
    
void test_ToroTest1(double &x0,double &dx,std::vector<std::array<double,8> >& U_P){
    double rho_L, vx_L, vy_L, vz_L, p_L, Bx_L, By_L, Bz_L;
    double rho_R, vx_R, vy_R, vz_R, p_R, Bx_R, By_R, Bz_R;
    
    // Left initial conditions
    rho_L = 1.0;
    vx_L = 0.0;
    vy_L = 0.0;
    vz_L = 0.0;
    p_L = 1.0;
    Bx_L = 0.0;
    By_L = 0.0;
    Bz_L = 0.0;
    
    // Right initial conditions
    rho_R = 0.125;
    vx_R = 0.0;
    vy_R = 0.0;
    vz_R = 0.0;
    p_R = 0.1;
    Bx_R = 0.0;
    By_R = 0.0;
    Bz_R = 0.0;
    
    
    for(int i=0; i<int(U_P.size()); ++i){
        double x = x0 + (i-1.5)*dx;
            
        if (x<=0.5){
            U_P[i][0]=rho_L;
            U_P[i][1]=vx_L;
            U_P[i][2]=vy_L;
            U_P[i][3]=vz_L;
            U_P[i][4]=p_L;
            U_P[i][5]=Bx_L;
            U_P[i][6]=By_L;
            U_P[i][7]=Bz_L;
        }
        
        else{
            U_P[i][0]=rho_R;
            U_P[i][1]=vx_R;
            U_P[i][2]=vy_R;
            U_P[i][3]=vz_R;
            U_P[i][4]=p_R;
            U_P[i][5]=Bx_R;
            U_P[i][6]=By_R;
            U_P[i][7]=Bz_R;
        }
    }
}

void test_Li25D(double &x0,double &dx,std::vector<std::array<double,8> >& U_P){
    double rho_L, vx_L, vy_L, vz_L, p_L, Bx_L, By_L, Bz_L;
    double rho_R, vx_R, vy_R, vz_R, p_R, Bx_R, By_R, Bz_R;
    
    // Left initial conditions
    rho_L = 1.08;
    vx_L = 1.2;
    vy_L = 0.01;
    vz_L = 0.5;
    p_L = 0.95;
    Bx_L = 2.0/sqrt(4.0*M_PI);
    By_L = 3.6/sqrt(4.0*M_PI);
    Bz_L = 2.0/sqrt(2.0*M_PI);
    
    // Right initial conditions
    rho_R = 1.0;
    vx_R = 0.0;
    vy_R = 0.0;
    vz_R = 0.0;
    p_R = 1.0;
    Bx_R = 2.0/sqrt(4.0*M_PI);
    By_R = 4.0/sqrt(4.0*M_PI);
    Bz_R = 2.0/sqrt(4.0*M_PI);
    
    
    for(int i=0; i<int(U_P.size()); ++i){
        double x = x0 + (i-1.5)*dx;
            
        if (x<=0.5){
            U_P[i][0]=rho_L;
            U_P[i][1]=vx_L;
            U_P[i][2]=vy_L;
            U_P[i][3]=vz_L;
            U_P[i][4]=p_L;
            U_P[i][5]=Bx_L;
            U_P[i][6]=By_L;
            U_P[i][7]=Bz_L;
        }
        
        else{
            U_P[i][0]=rho_R;
            U_P[i][1]=vx_R;
            U_P[i][2]=vy_R;
            U_P[i][3]=vz_R;
            U_P[i][4]=p_R;
            U_P[i][5]=Bx_R;
            U_P[i][6]=By_R;
            U_P[i][7]=Bz_R;
        }
    }
}
