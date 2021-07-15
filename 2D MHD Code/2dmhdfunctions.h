#include<iostream> 
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<array>
#include<algorithm>


// Declare functions for setup
//void initialise(std::string& test, std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain);
void ReadParams(int &testnum, int &limiternum, int &bconditions, int &xCells, int &yCells, double &x0, double &x1, double &y0, double &y1, double &tstop, double &C, double &gamma, double &omega);


double computeTimeStep(const std::vector < std::vector < std::array<double, 8> > >& UP, const double &C, const double dx, const double dy, const double gamma);

std::array<double,8> Primitive_to_Conservative(std::array<double,8> & UP_ij, const double gamma);
std::array<double,8> Conservative_to_Primitive(std::array<double,8> & UC_ij, const double gamma);

double magnitude(const double a, const double b, const double c);


// Declare functions for boundary conditions
void transmissive_bc(std::vector< std::vector <std::array<double,8> > > &UP, int &xCells, int &yCells);
void periodic_bc(std::vector< std::vector <std::array<double,8> > > &UP, int &xCells, int &yCells);
void KH_bc(std::vector< std::vector <std::array<double,8> > > &UP, int &xCells, int &yCells);


// Declare functions for getting limiters
double computeR_x(const int i, const int j, const int var, const std::vector < std::vector < std::array<double, 8> > > &UC);
double computeR_y(const int i, const int j, const int var, const std::vector < std::vector < std::array<double, 8> > > &UC);
double Minbeelimiter (double r);
double Superbeelimiter (double r);
double VanLeerlimiter (double r);
double VanAlbadalimiter (double r);

//Declare functions for each test
void ToroTest1_x(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain);
void BrioWu_x(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain);
void BrioWu_y(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain);
void DaiWoodward_x(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain);
void OrzagTang(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain);
void KelvinHelmholtz(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain);
void ToroTest1_y(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain);
void ToroTest1_nonalign(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain);
void CylindricalExplosion(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain);
void BrioWu_nonalign(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain);

// Declare functions for getting flux
std::array<double,8> f_euler_xdirection(std::array<double,8> & U_Ci, const double & gamma);
std::array<double,8> f_euler_ydirection(std::array<double,8> & U_Ci, const double & gamma);
std::array<double,8> getFORCEflux_xdirection(std::array<double,8> & U_Ci, std::array<double,8> & U_Ciplus1, const double & dx, double &dt, const double &gamma);
std::array<double,8> getFORCEflux_ydirection(std::array<double,8> & U_Ci, std::array<double,8> & U_Ciplus1, const double & dx, double &dt, const double &gamma);
std::array<double,8> getHLLCflux_xdirection(std::array<double,8> & UC_leftstate, std::array<double,8> & UC_rightstate, const double & gamma);
std::array<double,8> getHLLCflux_ydirection(std::array<double,8> & UC_leftstate, std::array<double,8> & UC_rightstate, const double & gamma);
double getfastwavespeed(std::array<double,8>& UP_state,const double &gamma);
double getfastwavespeed2(std::array<double,8>& UP_state,const double &gamma);
