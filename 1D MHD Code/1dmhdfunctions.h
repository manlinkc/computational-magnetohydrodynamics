#include<iostream> 
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<array>
#include<algorithm>

// Declare functions used in main code

double computeTimeStep(std::vector<std::array<double,8> >& U_P, double &C, double dx, double gamma, double ncells);

// Declare functions used to convert between conservative and primitive variables
double magnitude(const double a, const double b, const double c);
std::array<double,8> Primitive_to_Conservative(std::array<double,8> & U_Pi, const double gamma);
std::array<double,8> Conservative_to_Primitive(std::array<double,8> & U_Ci, const double &gamma);

// Declare functions for the limiters
double Minbeelimiter (const int &i, const int &j, const std::vector<std::array<double,8> > & U_C);
double Superbeelimiter (const int &i, const int &j, const std::vector<std::array<double,8> > & U_C);
double VanLeerlimiter (const int &i, const int &j, const std::vector<std::array<double,8> > & U_C);
double VanAlbadalimiter (const int &i, const int &j, const std::vector<std::array<double,8> > & U_C);

// Declare functions use to compute the flux
std::array<double,8> f_mhd(std::array<double,8> & U_Ci, const double & gamma);

// Declare functions used in main code
std::array<double,8> getHLLCflux(std::array<double,8> & UC_leftstate, std::array<double,8> & UC_rightstate, const double & gamma);
double getfastwavespeed(std::array<double,8>& UP_state,const double &gamma);

// Declare functions for boundary conditions
void Transmissiveboundaries(std::vector<std::array<double, 8> > &U, const int &ncells);
