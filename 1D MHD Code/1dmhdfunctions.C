#include<iostream> 
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<array>
#include<algorithm>
#include "1dmhdfunctions.h"

//--------------------------COMPUTE TIME STEP---------------------------
//----------------------------------------------------------------------
double computeTimeStep(std::vector<std::array<double,8> >& U_P, double &C, double dx, double gamma, double ncells){
	// Function Description: Compute dt for MHD Equations
	
	double dt;
	
	//Define a vector for a 
	std::vector<double>a;
	// Make a the same size as the number of cells + ghost cells in our grid
	a.resize(ncells+4);
	
	// Compute values for a but not including the boundary condition
	for (int i = 2; i < int(a.size())-2; i++)
	{
		// Get primitive variables
		double rho = U_P[i][0];
		double vx = U_P[i][1];
		double vy = U_P[i][2];
		double vz = U_P[i][3];
		double p = U_P[i][4];
        double Bx = U_P[i][5];
		double By = U_P[i][6];
		double Bz = U_P[i][7];
		
		//Get |v| and B.B
		double modv = magnitude(vx,vy,vz);
		double Bsquared = pow(magnitude(Bx,By,Bz),2.0);
		//Get c_s and c_a 
		double cs_squared = (gamma*p)/rho;
		double ca_squared =  Bsquared/rho;
		// Get c_s**2 + c_a**2 
		double cterms = cs_squared + ca_squared;
		
		// Fast wave speed
		double cf = sqrt( 0.5*(cterms + sqrt(cterms*cterms - ((4.0*cs_squared*Bx*Bx)/rho)))); 
		
		// a = |v| + fast wave speed
		a[i] = modv + cf; 
	} 
	
	// Get the maximum value of a to set the time setep
	double amax = 0.0;
	
	// Loop over u but not including the boundary condition
	for (int i = 2; i < int(a.size())-2; i++)
	{
		if (fabs(a[i]) > amax)
		{
			amax = fabs(a[i]);
			}
	} 
	dt = (C*dx)/amax;
	
	return dt;
}

//--------------------------MAGNITUDE-----------------------------------
//----------------------------------------------------------------------
double magnitude(const double a, const double b, const double c){
	//Function Description: Compute magnitude of vector containing three elements a, b, c
    
    double result = sqrt(a*a + b*b + c*c);
    return result;
    
}
//-----------------------PRIMITIVE TO CONSERVATIVE----------------------
//----------------------------------------------------------------------

std::array<double,8> Primitive_to_Conservative(std::array<double,8> & U_Pi, const double gamma){
	//Function description: Convert primitive variables to conservative variables

		std::array<double,8> U_Ci;
		
		double rho = U_Pi[0];
		double vx = U_Pi[1];
		double vy = U_Pi[2];
		double vz = U_Pi[3];
		double p = U_Pi[4];
        double Bx = U_Pi[5];
		double By = U_Pi[6];
		double Bz = U_Pi[7];
		
		double vsquared = pow(magnitude(vx,vy,vz),2.0);
		double Bsquared = pow(magnitude(Bx,By,Bz),2.0);
		
		//Density remains unchanged
		U_Ci[0] = rho;
		
		//Change v to momentum by multiplying rho and v
		U_Ci[1] = rho*vx;
		U_Ci[2] = rho*vy;
		U_Ci[3] = rho*vz;
		
		//Change pressure to Energy
		U_Ci[4] = (p/(gamma-1.0)) + 0.5*rho*vsquared + 0.5*Bsquared;
		
		//Magnetic field remains unchanged
        U_Ci[5] = Bx;
		U_Ci[6] = By;
		U_Ci[7] = Bz;
	
		return U_Ci;
}
//--------------------------CONSERVATIVE TO PRIMATIVE-------------------
//----------------------------------------------------------------------

std::array<double,8> Conservative_to_Primitive(std::array<double,8> & U_Ci, const double &gamma){
		//Function description: Convert conservative variables to primitive variables
	
		std::array<double,8> U_Pi;
		
		double rho = U_Ci[0];
		double momx = U_Ci[1];
		double momy = U_Ci[2];
		double momz = U_Ci[3];
		double U = U_Ci[4];
        double Bx = U_Ci[5];
		double By = U_Ci[6];
		double Bz = U_Ci[7];
		
		double momsquared = pow(magnitude(momx,momy,momz),2.0);
		double Bsquared = pow(magnitude(Bx,By,Bz),2.0);
		
		//Density remains unchanged
		U_Pi[0]=rho;
		
		//Change momentum to velocity by dividing momentum with rho
		U_Pi[1] = momx/rho;
		U_Pi[2] = momy/rho;
		U_Pi[3] = momz/rho;
		
		//Change Energy to Pressure
		U_Pi[4] = (gamma-1.0)*(U - (0.5*Bsquared) - (0.5*(momsquared/rho)));

		//Magnetic field remains unchanged
        U_Pi[5] = Bx;
		U_Pi[6] = By;
		U_Pi[7] = Bz;
		
		return U_Pi;
}

//--------------------------CONSERVATIVE FLUX---------------------------
//----------------------------------------------------------------------
std::array<double,8> f_mhd(std::array<double,8> & U_Ci, const double & gamma){
	//Function Description: Computes mhd flux using conservative variables
		
		//Define array for returning output 
		std::array<double, 8> flux_i;
		
		double rho = U_Ci[0];
		double momx = U_Ci[1];
		double momy = U_Ci[2];
		double momz = U_Ci[3];
		double U = U_Ci[4];
        double Bx = U_Ci[5];
		double By = U_Ci[6];
		double Bz = U_Ci[7];
		
		double momsquared = pow(magnitude(momx,momy,momz),2.0);
		double Bsquared = pow(magnitude(Bx,By,Bz),2.0);
		double p = (gamma-1.0)*(U - (0.5*Bsquared) - (0.5*(momsquared/rho)));
		double vdotB = (momx*Bx + momy*By + momz*Bz)/rho;
		
		//double rhosquared = rho*rho;
		double momxsquared = momx*momx;
		//double momcubed = mom*mom*mom;
		
		flux_i[0] = momx;
		
		flux_i[1] = (momxsquared/rho) + p + (0.5*Bsquared) - (Bx*Bx);
		
		flux_i[2] = ((momx*momy)/rho) - Bx*By;
		
		flux_i[3] = ((momx*momz)/rho) - Bx*Bz;
		
		flux_i[4] = (U + p + 0.5*Bsquared)*(momx/rho) - vdotB*Bx;
    
        flux_i[5] = 0.0;
		
		flux_i[6] = (By*momx - Bx*momy)/rho;
		
		flux_i[7] = (Bz*momx - Bx*momz)/rho;
		
		return flux_i;
	}

//--------------------------MINBEE LIMITER-----------------------------
//----------------------------------------------------------------------
double Minbeelimiter (const int &i, const int &j, const std::vector<std::array<double,8> > & U_C){
	
	double minbee_psi=0.0;
	double r_num = U_C[i][j] - U_C[i-1][j];
	double r_denom = U_C[i+1][j] - U_C[i][j];
	
	//Case when r becomes singular 
    if (r_denom == 0){
        minbee_psi=0;
        }
    else if (r_denom!=0){
		// Construct r 
		double r;
		r = r_num/r_denom;
		
		if (r <= 0){
			minbee_psi=0.0;
		}
		else if (0<r && r<=1){ 
			minbee_psi = r;
		}
		else if (r>1){
			double epsilonR = 2.0/(1.0 + r);
			minbee_psi = std::min(1.0, epsilonR);
        }
	}
	return minbee_psi;
}

//--------------------------SUPERBEE LIMITER-----------------------------
//-----------------------------------------------------------------------
double Superbeelimiter (const int &i, const int &j, const std::vector<std::array<double,8> > & U_C){
	
	double superbee_psi=0.0;
	double r_num = U_C[i][j] - U_C[i-1][j];
	double r_denom = U_C[i+1][j] - U_C[i][j];

	//Case when r becomes singular 
    if (r_denom == 0){
        superbee_psi=0;
        }
     else if (r_denom!=0){
		// Construct r 
		double r;
		r = r_num/r_denom;
		
		if (r <= 0){
			superbee_psi=0.0;
		}
		else if (0<r && r<=0.5){ 
			superbee_psi = 2.0*r;
		}
		else if (0.5<r && r<=1){ 
			superbee_psi = 1.0;
		}
		else if (r>1){
			double epsilonR = 2.0/(1.0 + r);
			superbee_psi = std::min(r, std::min(epsilonR, 2.0));
		}	
	}
	
	return superbee_psi;
}

//--------------------------VAN-LEER LIMITER-----------------------------
//-----------------------------------------------------------------------
double VanLeerlimiter (const int &i, const int &j, const std::vector<std::array<double,8> > & U_C){
	
	double vanleer_psi=0.0;
	double r_num = U_C[i][j] - U_C[i-1][j];
	double r_denom = U_C[i+1][j] - U_C[i][j];
	
	//Case when r becomes singular 
    if (r_denom == 0){
        vanleer_psi=0;
        }
    
    else if (r_denom!=0){
		// Construct r 
		double r;
		r = r_num/r_denom;
		
		if (r <= 0){
			vanleer_psi=0;
		}
		else if (r>0){
			double epsilonL = (2.0*r)/(1.0 + r);
			double epsilonR = 2.0/(1.0 + r);
			
			vanleer_psi = std::min(epsilonL, epsilonR);
		}	
	}
	return vanleer_psi;
}
//--------------------------VAN-ALBADA LIMITER-----------------------------
//-----------------------------------------------------------------------
double VanAlbadalimiter (const int &i, const int &j, const std::vector<std::array<double,8> > & U_C){
	
	double vanalbada_psi=0.0;
	double r_num = U_C[i][j] - U_C[i-1][j];
	double r_denom = U_C[i+1][j] - U_C[i][j];
	
	//Case when r becomes singular 
    if (r_denom == 0){
        vanalbada_psi=0;
        }
    else if (r_denom!=0){
		// Construct r 
		double r;
		r = r_num/r_denom;
		
		if (r <= 0){
			vanalbada_psi=0;
		}
		else if (r>0){
			double rquantity = (r*(1.0+r))/(1.0 + r*r);
			double epsilonR = 2.0/(1.0 + r);
			
			vanalbada_psi = std::min(rquantity, epsilonR);
		}	
	}
	
	return vanalbada_psi;
}


//-------------------------------------HLLC FLUX---------------------------------------------
//-------------------------------------------------------------------------------------------

std::array<double,8> getHLLCflux(std::array<double,8> & UC_leftstate, std::array<double,8> & UC_rightstate, const double & gamma){
    // Function Description: Calculates the HLLC flux

    // Left state (variables)
    std::array<double,8> UP_leftstate;
    UP_leftstate = Conservative_to_Primitive(UC_leftstate,gamma);
    
    double rho_L = UC_leftstate[0];
    //double momx_L = UC_leftstate[1];
    //double momy_L = UC_leftstate[2];
    //double momz_L = UC_leftstate[3];
    double U_L = UC_leftstate[4];
    double Bx_L = UC_leftstate[5];
    double By_L = UC_leftstate[6];
    double Bz_L = UC_leftstate[7];
    double vx_L = UP_leftstate[1];
    double vy_L = UP_leftstate[2];
    double vz_L = UP_leftstate[3];
    double p_L = UP_leftstate[4];

    // Right state (conservative variables)
    std::array<double,8> UP_rightstate;
    UP_rightstate = Conservative_to_Primitive(UC_rightstate,gamma);
    
    double rho_R = UC_rightstate[0];
    //double momx_R = UC_rightstate[1];
    //double momy_R = UC_rightstate[2];
    //double momz_R = UC_rightstate[3];
    double U_R = UC_rightstate[4];
    double Bx_R = UC_rightstate[5];
    double By_R = UC_rightstate[6];
    double Bz_R = UC_rightstate[7];
    double vx_R = UP_rightstate[1];
    double vy_R = UP_rightstate[2];
    double vz_R = UP_rightstate[3];
    double p_R = UP_rightstate[4];
    
    // Get minimum velocity x-direction
    double min_vx = std::min(vx_L, vx_R);
    double max_vx = std::max(vx_L, vx_R);
    double cf_L = getfastwavespeed(UP_leftstate,gamma);
    double cf_R = getfastwavespeed(UP_rightstate,gamma);
    double max_cf = std::max(cf_L, cf_R);
    
    // Wave Speeds SL, SR, SM
    double SL = min_vx - max_cf;
    double SR = max_vx + max_cf;
    
    double p_TL = p_L + 0.5*(Bx_L*Bx_L + By_L*By_L + Bz_L*Bz_L);
    double p_TR = p_R + 0.5*(Bx_R*Bx_R + By_R*By_R + Bz_R*Bz_R);

    double Sstar = (rho_R * vx_R * (SR-vx_R) - rho_L * vx_L * (SL-vx_L) + p_TL - p_TR - Bx_L*Bx_L + Bx_R*Bx_R )/(rho_R*(SR - vx_R) - rho_L*(SL - vx_L));

    // Compute MHD conservative flux left and right state
    std::array<double, 8> f_L;
    std::array<double, 8> f_R;
    f_L = f_mhd(UC_leftstate, gamma);
    f_R = f_mhd(UC_rightstate, gamma);

    // Compute intermediate state in HLL Li_MHDHLLC paper eqn (2)
    std::array<double, 8> UHLL;
    for (int i=0; i<8; i++){
        UHLL[i] = (SR * UC_rightstate[i] - SL * UC_leftstate[i] - f_R[i] + f_L[i])/(SR - SL);
    }
    
    double p_Tstar = rho_L*(SL - vx_L)*(Sstar - vx_L) + p_TL;
    
    // Compute intermediate states in HLLC
    std::array<double,8> UHLLC_L;
    std::array<double,8> UHLLC_R;
    double HLLvdotB = (UHLL[1]/UHLL[0])*UHLL[5] + (UHLL[2]/UHLL[0])*UHLL[6] + (UHLL[3]/UHLL[0])*UHLL[7];

    UHLLC_L[0] = rho_L*((SL - vx_L)/(SL - Sstar)); //rho_star_L
    UHLLC_L[1] = UHLLC_L[0]*Sstar; //momx_star_L
    UHLLC_L[2] = rho_L*vy_L*((SL - vx_L)/(SL - Sstar)) - ((UHLL[5]*UHLL[6] - Bx_L*By_L)/(SL - Sstar)); //momy_star_L
    UHLLC_L[3] = rho_L*vz_L*((SL - vx_L)/(SL - Sstar)) - ((UHLL[5]*UHLL[7] - Bx_L*Bz_L)/(SL - Sstar)); //momz_star_L
    double vdotB_L = vx_L*Bx_L + vy_L*By_L + vz_L*Bz_L;
    UHLLC_L[4] = U_L*((SL - vx_L)/(SL - Sstar)) + (p_Tstar*Sstar - p_TL*vx_L - (UHLL[5]*HLLvdotB - Bx_L*vdotB_L))/(SL - Sstar);
    UHLLC_L[5] = UHLL[5]; // Bx_star_L
    UHLLC_L[6] = UHLL[6]; // By_star_L
    UHLLC_L[7] = UHLL[7]; // Bz_star_L
    
    // HLLC right star state
    UHLLC_R[0] = rho_R*((SR - vx_R)/(SR - Sstar)); //rho_star_R
    UHLLC_R[1] = UHLLC_R[0]*Sstar; //momx_star_R
    UHLLC_R[2] = rho_R*vy_R*((SR - vx_R)/(SR - Sstar)) - ((UHLL[5]*UHLL[6] - Bx_R*By_R)/(SR - Sstar)); //momy_star_R
    UHLLC_R[3] = rho_R*vz_R*((SR - vx_R)/(SR - Sstar)) - ((UHLL[5]*UHLL[7] - Bx_R*Bz_R)/(SR - Sstar)); //momz_star_R
    double vdotB_R = vx_R*Bx_R + vy_R*By_R + vz_R*Bz_R;
    UHLLC_R[4] = U_R*((SR - vx_R)/(SR - Sstar)) + (p_Tstar*Sstar - p_TR*vx_R - (UHLL[5]*HLLvdotB - Bx_R*vdotB_R))/(SR - Sstar); //E_star_R
    UHLLC_R[5] = UHLL[5]; // Bx_star_R
    UHLLC_R[6] = UHLL[6]; // By_star_R
    UHLLC_R[7] = UHLL[7]; // Bz_star_R
    
    // Define fluxes
    std::array<double,8> fHLLC_L;
    std::array<double,8> fHLLC_R;
    for (int i=0; i <8; i++){
        fHLLC_L[i] = f_L[i] + SL*(UHLLC_L[i] - UC_leftstate[i]);
        fHLLC_R[i] = f_R[i] + SR*(UHLLC_R[i] - UC_rightstate[i]);
    }
    
    std::array<double,8> finalflux;
    
    // Return flux conditions
    if (0 <= SL) {
        //std::cout << "f_L" << std::endl;
        finalflux = f_L;
    }
    else if (SL < 0 && 0 <= Sstar) {
        //std::cout << "fHLLC_L" << std::endl;
        finalflux = fHLLC_L;
    }
    else if (Sstar < 0 && 0 <= SR) {
        //std::cout << "fHLLC_R" << std::endl;
        finalflux = fHLLC_R;
    }
    else if (SR < 0) {
        //std::cout << "f_R" << std::endl;
        finalflux = f_R;
    }

    return finalflux;
}

//------------------------------GET FAST WAVE SPEED-------------------------------------------------
//--------------------------------------------------------------------------------------------------

double getfastwavespeed(std::array<double,8>& UP_state,const double &gamma){
    
        // Get primitive variables
        double rho = UP_state[0];
        //double vx = U_P[i][1];
        //double vy = U_P[i][2];
        //double vz = U_P[i][3];
        double p = UP_state[4];
        double Bx = UP_state[5];
        double By = UP_state[6];
        double Bz = UP_state[7];
        
        double Bsquared = pow(magnitude(Bx,By,Bz),2.0);
        //Get c_s and c_a
        double cs_squared = (gamma*p)/rho;
        double ca_squared =  Bsquared/rho;
        // Get c_s**2 + c_a**2
        double cterms = cs_squared + ca_squared;
        
        // Fast speed
        double cf = sqrt( 0.5*(cterms + sqrt(cterms*cterms - ((4.0*cs_squared*Bx*Bx)/rho))));
    
    return cf;
}

void Transmissiveboundaries(std::vector<std::array<double, 8> > &U, const int &ncells){
    //Apply transmissive boundary conditions,the first real cell is at i=2
    for (int j=0; j<8; j++){
        U[0][j] = U[2][j];
        U[1][j] = U[2][j];
        U[ncells+2][j] = U[ncells+1][j];
        U[ncells+3][j] = U[ncells+1][j];
    }
}



