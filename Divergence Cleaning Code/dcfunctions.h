#include<iostream> 
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<array>
#include<algorithm>

double getBxtilde(std::array<double, 9> & UCleft, std::array<double, 9> & UCright, double &ch);
std::array<double,2> getDCflux_x(std::array<double,9> & UCleft, std::array<double,9> & UCright, double &ch);
double getBytilde(std::array<double, 9> & UCleft, std::array<double, 9> & UCright, double &ch);
std::array<double,2> getDCflux_y(std::array<double,9> & UCleft, std::array<double,9> & UCright, double &ch);

// Declare functions for setup
void initialise(std::string& test, std::vector< std::vector <std::array<double,9> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain);
double computeTimeStep(const std::vector < std::vector < std::array<double, 9> > >& UP, const double &C, const double dx, const double dy, const double gamma);

std::array<double,9> Primitive_to_Conservative(std::array<double,9> & UP_ij, const double gamma);
std::array<double,9> Conservative_to_Primitive(std::array<double,9> & UC_ij, const double gamma);

double magnitude(const double a, const double b, const double c);

// Declare functions for getting flux
std::array<double,9> f_euler_xdirection(std::array<double,9> & U_Ci, const double & gamma);
std::array<double,9> f_euler_ydirection(std::array<double,9> & U_Ci, const double & gamma);
std::array<double,9> getFORCEflux_xdirection(std::array<double,9> & U_Ci, std::array<double,9> & U_Ciplus1, const double & dx, double &dt, const double &gamma, std::array<double,2> & BPsiflux_x);
std::array<double,9> getFORCEflux_ydirection(std::array<double,9> & U_Ci, std::array<double,9> & U_Ciplus1, const double & dx, double &dt, const double &gamma, std::array<double,2> & BPsiflux_y);
std::array<double,9> getHLLCflux_xdirection(std::array<double,9> & UC_leftstate, std::array<double,9> & UC_rightstate, const double & gamma, std::array<double,2> &BPsi_x);
std::array<double,9> getHLLCflux_ydirection(std::array<double,9> & UC_leftstate, std::array<double,9> & UC_rightstate, const double & gamma, std::array<double,2> &BPsi_y);
double getfastwavespeed(std::array<double,9>& UP_state,const double &gamma);
double getfastwavespeed2(std::array<double,9>& UP_state,const double &gamma);

// Declare functions for boundary conditions
void KH_bc(std::vector< std::vector <std::array<double,9> > > &UP, int &xCells, int &yCells);
void periodic_bc(std::vector< std::vector <std::array<double,9> > > &UP, int &xCells, int &yCells);
void transmissive_bc(std::vector< std::vector <std::array<double,9> > > &UP, int &xCells, int &yCells);

void OrzagTang(std::vector< std::vector <std::array<double,9> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain);
void KelvinHelmholtz(std::vector< std::vector <std::array<double,9> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain);
void ReadParams(int &testnum, int &limiternum, int &bconditions, int &xCells, int &yCells, double &x0, double &x1, double &y0, double &y1, double &tstop, double &C, double &gamma, double &omega);
//Declare functions for getting limiters
double computeR_x(const int i, const int j, const int var, const std::vector < std::vector < std::array<double, 9> > > &UC);
double computeR_y(const int i, const int j, const int var, const std::vector < std::vector < std::array<double, 9> > > &UC);
double Minbeelimiter (double r);
double Superbeelimiter (double r);
double VanLeerlimiter (double r);
double VanAlbadalimiter (double r);
