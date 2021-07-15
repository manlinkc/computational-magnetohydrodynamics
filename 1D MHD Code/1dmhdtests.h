#include<iostream> 
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<array>
#include<algorithm>

void ReadParams(int &test, int &limiter, int &ncells, double &x0, double &x1, double &tstop, double & C, double &gamma, double & omega);

void test_modifiedSod(double &x0,double &dx, std::vector<std::array<double,8> >& U_P);
void test_BrioWu(double &x0,double &dx, std::vector<std::array<double,8> >& U_P);
void test_DaiWoodward(double &x0,double &dx, std::vector<std::array<double,8> >& U_P);
void test_BrioWuSwitched(double &x0,double &dx,  std::vector<std::array<double,8> >& U_P);
void test_DaiWoodwardSwitched(double &x0,double &dx,  std::vector<std::array<double,8> >& U_P);
void test_ToroTest1(double &x0,double &dx, std::vector<std::array<double,8> >& U_P);
void test_Li25D(double &x0,double &dx,std::vector<std::array<double,8> >& U_P);
