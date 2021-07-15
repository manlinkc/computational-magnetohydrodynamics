#include<iostream> 
#include<fstream>
#include<sstream>
#include<cmath>
#include<vector>
#include<string>
#include<array>
#include<algorithm>
#include "2dmhdfunctions.h"

//--------------------------COMPUTE R------------------------------------
//-----------------------------------------------------------------------
double computeR_x(const int i, const int j, const int var, const std::vector < std::vector < std::array<double, 8> > > &UC){
    double r_num, r_denom;
    double r = 0.0;
    
    r_num = UC[i][j][var] - UC[i-1][j][var];
    r_denom = UC[i+1][j][var] - UC[i][j][var];
    
    if (r_denom == 0){
        r = 0.0;
    }
    else if (r_denom!=0){
        r = r_num/r_denom;
    }
    
    return r;
}

double computeR_y(const int i, const int j, const int var, const std::vector < std::vector < std::array<double, 8> > > &UC){
    double r_num, r_denom;
    double r = 0.0;
    
    r_num = UC[i][j][var] - UC[i][j-1][var];
    r_denom = UC[i][j+1][var] - UC[i][j][var];
    
    if (r_denom == 0){
        r = 0.0;
    }
    else if (r_denom!=0){
        r = r_num/r_denom;
    }
    
    return r;
}

//--------------------------MINBEE LIMITER-------------------------------
//-----------------------------------------------------------------------
double Minbeelimiter (double r){

    double minbee_psi = 0.0;

    if (r <= 0){
        minbee_psi = 0.0;
    }
    else if (0<r && r<=1){
        minbee_psi = r;
    }
    else if (r>1){
        double epsilonR = 2.0/(1.0 + r);
        minbee_psi = std::min(1.0, epsilonR);
    }

    return minbee_psi;
}
//--------------------------SUPERBEE LIMITER-----------------------------
//-----------------------------------------------------------------------
double Superbeelimiter (double r){
    
    double superbee_psi=0.0;
    
    if (r <= 0){
        superbee_psi=0;
    }
    else if (0<r && r<=0.5){
        superbee_psi = 2*r;
    }
    else if (0.5<r && r<=1){
        superbee_psi = 1;
    }
    else if (r>1){
        double epsilonR = 2.0/(1.0 + r);
        superbee_psi = std::min(r, std::min(epsilonR, 2.0));
    }
    
    return superbee_psi;
}

//--------------------------VAN-LEER LIMITER-----------------------------
//-----------------------------------------------------------------------
double VanLeerlimiter (double r){
        
    double vanleer_psi=0.0;
    
    if (r <= 0){
        vanleer_psi=0;
    }
    else if (r>0){
        double epsilonL = (2.0*r)/(1.0 + r);
        double epsilonR = 2.0/(1.0 + r);
        
        vanleer_psi = std::min(epsilonL, epsilonR);
    }
    
    return vanleer_psi;
}

//--------------------------VAN-ALBADA LIMITER-----------------------------
//-------------------------------------------------------------------------
double VanAlbadalimiter (double r){
    
    double vanalbada_psi=0.0;
    
    if (r <= 0){
        vanalbada_psi=0;
    }
    else if (r>0){
        double rquantity = (r*(1.0+r))/(1.0 + r*r);
        double epsilonR = 2.0/(1.0 + r);
        
        vanalbada_psi = std::min(rquantity, epsilonR);
    }

    return vanalbada_psi;
}
//---------------------COMPUTE TIME STEP--------------------------------
//----------------------------------------------------------------------
double computeTimeStep(const std::vector < std::vector < std::array<double, 8> > >& UP, const double &C, const double dx, const double dy, const double gamma){
	// Function Description: Compute dt for Euler's Equation in 2D
	
	double dt = 0.0;
	double a = 0.0;
	double amax = 0.0;
    double a_x, a_y;
    
	// Compute values for a 
	for (int i=0; i < int(UP.size()); i++){
		for (int j=0; j<int(UP[1].size()); j++){
            
            // Get primitive variables
            double rho = UP[i][j][0];
            double vx = UP[i][j][1];
            double vy = UP[i][j][2];
            double vz = UP[i][j][3];
            double p  = UP[i][j][4];
            double Bx = UP[i][j][5];
            double By = UP[i][j][6];
            double Bz = UP[i][j][7];
            
            //Get v.v and B.B
            double modv = magnitude(vx,vy,vz);
            double Bsquared = pow(magnitude(Bx,By,Bz),2.0);
            
            //Get c_s and c_a
            double cs_squared = (gamma*p)/rho;
            double ca_squared =  Bsquared/rho;
            
            // Get c_s**2 + c_a**2
            double cterms = cs_squared + ca_squared;
            // Fast wave in x direction
            double cf_x = sqrt( 0.5*(cterms + sqrt(cterms*cterms - ((4.0*cs_squared*Bx*Bx)/rho))));
            double cf_y = sqrt( 0.5*(cterms + sqrt(cterms*cterms - ((4.0*cs_squared*By*By)/rho))));
            
            // a = |v| + fast wave speed
            //a[i] = vdotv + cf;
            a_x = modv + cf_x;
            a_y = modv + cf_y;
            
            a = std::max(a_x, a_y);
			
			// Replace amax if a is larger
			if (fabs(a) > amax){
				amax = a;
				}				
		}
	}
	
	// need to do min of dx and dy here
	dt = (C*std::min(dx,dy))/amax;
	
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

std::array<double,8> Primitive_to_Conservative(std::array<double,8> & UP_ij, const double gamma){
	//Function description: Turn primitve variables to conservative variables
    
        std::array<double,8> UC_ij;
            
        double rho = UP_ij[0];
        double vx = UP_ij[1];
        double vy = UP_ij[2];
        double vz = UP_ij[3];
        double p = UP_ij[4];
        double Bx = UP_ij[5];
        double By = UP_ij[6];
        double Bz = UP_ij[7];
    
        double vsquared = pow(magnitude(vx,vy,vz),2.0);
        double Bsquared = pow(magnitude(Bx,By,Bz),2.0);
		
        //Density remains unchanged
		UC_ij[0] = rho;
		
        //Change v to momentum by multiplying rho and v
        UC_ij[1] = rho*vx;
        UC_ij[2] = rho*vy;
        UC_ij[3] = rho*vz;
        
        //Change pressure to Energy
        UC_ij[4] = (p/(gamma-1.0)) + 0.5*rho*vsquared + 0.5*Bsquared;
		
        //Magnetic field remains unchanged
        UC_ij[5] = Bx;
        UC_ij[6] = By;
        UC_ij[7] = Bz;
	
		return UC_ij;
}

//--------------------------CONSERVATIVE TO PRIMATIVE-------------------
//----------------------------------------------------------------------

std::array<double,8> Conservative_to_Primitive(std::array<double,8> & UC_ij,const double gamma){
		//Function description: Turn conservative variables to primitive variables
		
		std::array<double,8> UP_ij;
		
        double rho = UC_ij[0];
        double momx = UC_ij[1];
        double momy = UC_ij[2];
        double momz = UC_ij[3];
        double U = UC_ij[4];
        double Bx = UC_ij[5];
        double By = UC_ij[6];
        double Bz = UC_ij[7];
    
        double momsquared = pow(magnitude(momx,momy,momz),2.0);
        double Bsquared = pow(magnitude(Bx,By,Bz),2.0);
		
        //Density remains unchanged
		UP_ij[0]=rho;
		
        //Change momentum to velocity by dividing momentum with rho
        UP_ij[1] = momx/rho;
        UP_ij[2] = momy/rho;
        UP_ij[3] = momz/rho;
		
        //Change Energy to Pressure
        UP_ij[4] = (gamma-1.0)*(U - (0.5*Bsquared) - (0.5*(momsquared/rho)));
    
        //Magnetic field remains unchanged
        UP_ij[5] = Bx;
        UP_ij[6] = By;
        UP_ij[7] = Bz;

		return UP_ij;
}

//--------------------------CONSERVATIVE FLUX---------------------------
//------------------------------X - DIRECTION --------------------------
std::array<double,8> f_euler_xdirection(std::array<double,8> & U_Ci, const double & gamma){
	//Function Description:
	//Computes flux only in x direction
		
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
        double momxsquared = momx*momx;
            
		//Calculate flux
        flux_i[0] = momx;
        
        flux_i[1] = (momxsquared/rho) + p + (0.5*Bsquared) - (Bx*Bx);
        
        flux_i[2] = ((momx*momy)/rho) - Bx*By;
        
        flux_i[3] = ((momx*momz)/rho) - Bx*Bz;
        
        flux_i[4] = (U + p + 0.5*Bsquared)*(momx/rho) - vdotB*Bx;
        
        flux_i[5] = 0;
        
        flux_i[6] = (By*momx - Bx*momy)/rho;
        
        flux_i[7] = (Bz*momx - Bx*momz)/rho;

		
		return flux_i;
	}

//--------------------------CONSERVATIVE FLUX---------------------------
//------------------------------Y - DIRECTION --------------------------
std::array<double,8> f_euler_ydirection(std::array<double,8> & U_Ci, const double & gamma){
	//Function Description:
	//Computes flux only in y direction
		
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
        double momysquared = momy*momy;
		
        //Calculate flux
        flux_i[0] = momy;
        
        flux_i[1] = ((momy*momx)/rho) - By*Bx;
        
        flux_i[2] = (momysquared/rho) + p + (0.5*Bsquared) - (By*By);
        
        flux_i[3] = ((momy*momz)/rho) - By*Bz;
        
        flux_i[4] = (U + p + 0.5*Bsquared)*(momy/rho) - vdotB*By;
        
        flux_i[5] = (Bx*momy - By*momx)/rho;
        
        flux_i[6] = 0;
        
        flux_i[7] = (Bz*momy - By*momz)/rho;
    
		return flux_i;
	}

//--------------------------GET FORCE FLUX------------------------------
//----------------------------X-DIRECTION-------------------------------
	
std::array<double,8> getFORCEflux_xdirection(std::array<double,8> & U_Ci, std::array<double,8> & U_Ciplus1, const double & dx, double &dt, const double &gamma){
	// Function Description: Calculates the force flux
	
	std::array<double, 8> LF_flux;
	std::array<double, 8> R_flux;
	std::array<double, 8> u_iph_nph;
	std::array<double, 8> FORCE_flux;
	
	std::array<double, 8> flux_uiplusone = f_euler_xdirection(U_Ciplus1, gamma);
	std::array<double, 8> flux_ui = f_euler_xdirection(U_Ci, gamma);

	// Compute Lax-Friedrich Flux
    for (int var = 0; var < 8; var++){
        LF_flux[var] = 0.5*(dx/dt)*(U_Ci[var]-U_Ciplus1[var]) + 0.5*(flux_uiplusone[var]+flux_ui[var]);
    }
    
	// Compute term u_iplushalf_nplushalf to use in Richtmyer's flux
    for (int var = 0; var < 8; var++){
        u_iph_nph[var] = 0.5*(U_Ci[var] + U_Ciplus1[var]) - 0.5*(dt/dx)*(flux_uiplusone[var]-flux_ui[var]);
    }
    
	// Compute Richtmyer's flux
	R_flux= f_euler_xdirection(u_iph_nph, gamma);
	
	// Compute FORCE flux
    for (int var =0; var <8; var ++){
        FORCE_flux[var] = 0.5*(LF_flux[var] + R_flux[var]);
    }

	return FORCE_flux;

}

//--------------------------GET FORCE FLUX------------------------------
//----------------------------Y-DIRECTION-------------------------------
	
std::array<double,8> getFORCEflux_ydirection(std::array<double,8> & U_Ci, std::array<double,8> & U_Ciplus1, const double & dx, double &dt, const double &gamma){
	// Function Description: Calculates the force flux
	
	std::array<double, 8> LF_flux;
	std::array<double, 8> R_flux;
	std::array<double, 8> u_iph_nph;
	std::array<double, 8> FORCE_flux;
	
	std::array<double, 8> flux_uiplusone = f_euler_ydirection(U_Ciplus1, gamma);
	std::array<double, 8> flux_ui = f_euler_ydirection(U_Ci, gamma);

	// Compute Lax-Friedrich Flux
    for (int var=0; var <8; var++){
        LF_flux[var] = 0.5*(dx/dt)*(U_Ci[var]-U_Ciplus1[var]) + 0.5*(flux_uiplusone[var]+flux_ui[var]);
    }
    
	// Compute term u_iplushalf_nplushalf to use in Richtmyer's flux
    for (int var=0; var <8; var++){
        u_iph_nph[var] = 0.5*(U_Ci[var] + U_Ciplus1[var]) - 0.5*(dt/dx)*(flux_uiplusone[var]-flux_ui[var]);
    }

	// Compute Richtmyer's flux
	R_flux= f_euler_ydirection(u_iph_nph, gamma);
	
	// Compute FORCE flux
    for (int var =0; var <8; var++){
        FORCE_flux[var] = 0.5*(LF_flux[var] + R_flux[var]);
    }

	return FORCE_flux;

}
//-------------------------------------HLLC FLUX X-DIRECTION ---------------------------------------------
//--------------------------------------------------------------------------------------------------------

std::array<double,8> getHLLCflux_xdirection(std::array<double,8> & UC_leftstate, std::array<double,8> & UC_rightstate, const double & gamma){
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
    f_L = f_euler_xdirection(UC_leftstate, gamma);
    f_R = f_euler_xdirection(UC_rightstate, gamma);

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

//-------------------------------------HLLC FLUX Y-DIRECTION ---------------------------------------------
//--------------------------------------------------------------------------------------------------------
std::array<double,8> getHLLCflux_ydirection(std::array<double,8> & UC_leftstate, std::array<double,8> & UC_rightstate, const double & gamma){
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
    
    // Get minimum velocity y-direction
    double min_vy = std::min(vy_L, vy_R);
    double max_vy = std::max(vy_L, vy_R);
    double cf_L = getfastwavespeed2(UP_leftstate,gamma);
    double cf_R = getfastwavespeed2(UP_rightstate,gamma);
    double max_cf = std::max(cf_L, cf_R);
    
    // Wave Speeds SL, SR, SM
    double SL = min_vy - max_cf;
    double SR = max_vy + max_cf;
    
    double p_TL = p_L + 0.5*(Bx_L*Bx_L + By_L*By_L + Bz_L*Bz_L);
    double p_TR = p_R + 0.5*(Bx_R*Bx_R + By_R*By_R + Bz_R*Bz_R);

    double Sstar = (rho_R * vy_R * (SR-vy_R) - rho_L * vy_L * (SL-vy_L) + p_TL - p_TR - By_L*By_L + By_R*By_R )/(rho_R*(SR - vy_R) - rho_L*(SL - vy_L));

    // Compute MHD conservative flux left and right state
    std::array<double, 8> f_L;
    std::array<double, 8> f_R;
    f_L = f_euler_ydirection(UC_leftstate, gamma);
    f_R = f_euler_ydirection(UC_rightstate, gamma);

    // Compute intermediate state in HLL Li_MHDHLLC paper eqn (2)
    std::array<double, 8> UHLL;
    for (int i=0; i<8; i++){
        UHLL[i] = (SR * UC_rightstate[i] - SL * UC_leftstate[i] - f_R[i] + f_L[i])/(SR - SL);
    }
    
    double p_Tstar = rho_L*(SL - vy_L)*(Sstar - vy_L) + p_TL;
    
    // Compute intermediate states in HLLC
    std::array<double,8> UHLLC_L;
    std::array<double,8> UHLLC_R;
    double HLLvdotB = (UHLL[1]/UHLL[0])*UHLL[5] + (UHLL[2]/UHLL[0])*UHLL[6] + (UHLL[3]/UHLL[0])*UHLL[7];

    UHLLC_L[0] = rho_L*((SL - vy_L)/(SL - Sstar)); //rho_star_L
    UHLLC_L[1] = rho_L*vx_L*((SL - vy_L)/(SL - Sstar)) - ((UHLL[6]*UHLL[5] - By_L*Bx_L)/(SL - Sstar)); //momx_star_L
    UHLLC_L[2] = UHLLC_L[0]*Sstar; //momy_star_L
    UHLLC_L[3] = rho_L*vz_L*((SL - vy_L)/(SL - Sstar)) - ((UHLL[6]*UHLL[7] - By_L*Bz_L)/(SL - Sstar)); //momz_star_L
    double vdotB_L = vx_L*Bx_L + vy_L*By_L + vz_L*Bz_L;
    UHLLC_L[4] = U_L*((SL - vy_L)/(SL - Sstar)) + (p_Tstar*Sstar - p_TL*vy_L - (UHLL[6]*HLLvdotB - By_L*vdotB_L))/(SL - Sstar);
    UHLLC_L[5] = UHLL[5]; // Bx_star_L
    UHLLC_L[6] = UHLL[6]; // By_star_L
    UHLLC_L[7] = UHLL[7]; // Bz_star_L
    
    // HLLC right star state
    UHLLC_R[0] = rho_R*((SR - vy_R)/(SR - Sstar)); //rho_star_R
    UHLLC_R[1] = rho_R*vx_R*((SR - vy_R)/(SR - Sstar)) - ((UHLL[6]*UHLL[5] - By_R*Bx_R)/(SR - Sstar)); //momx_star_R
    UHLLC_R[2] = UHLLC_R[0]*Sstar; //momy_star_R
    UHLLC_R[3] = rho_R*vz_R*((SR - vy_R)/(SR - Sstar)) - ((UHLL[6]*UHLL[7] - By_R*Bz_R)/(SR - Sstar)); //momz_star_R
    double vdotB_R = vx_R*Bx_R + vy_R*By_R + vz_R*Bz_R;
    UHLLC_R[4] = U_R*((SR - vy_R)/(SR - Sstar)) + (p_Tstar*Sstar - p_TR*vy_R - (UHLL[6]*HLLvdotB - By_R*vdotB_R))/(SR - Sstar); //E_star_R
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

double getfastwavespeed2(std::array<double,8>& UP_state,const double &gamma){
    
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
        double cf = sqrt( 0.5*(cterms + sqrt(cterms*cterms - ((4.0*cs_squared*By*By)/rho))));
    
    return cf;
}

//--------------------------BOUNDARY CONDITIONS--------------------------
//-----------------------------------------------------------------------
void transmissive_bc(std::vector< std::vector <std::array<double,8> > > &UP, int &xCells, int &yCells){
    // Note: Psi has same boundary conditions as rho
    
    //Apply transmissive boundary conditions to left and right boundaries
    for (int j=0; j<int(UP[1].size()); j++){
        for (int k=0; k<8; k++){
            UP[0][j][k] = UP[2][j][k];
            UP[1][j][k] = UP[2][j][k];
            UP[xCells+2][j][k] = UP[xCells+1][j][k];
            UP[xCells+3][j][k] = UP[xCells+1][j][k];
            }
    }
    
    //Apply transmissive boundary conditions, to top and bottom boundaries
    for (int i=0; i<int(UP.size()); i++){
        for (int k=0; k<8; k++){
            UP[i][0][k] = UP[i][2][k];
            UP[i][1][k] = UP[i][2][k];
            UP[i][yCells+2][k] = UP[i][yCells+1][k];
            UP[i][yCells+3][k] = UP[i][yCells+1][k];
            }
    }
}

void periodic_bc(std::vector< std::vector <std::array<double,8> > > &UP, int &xCells, int &yCells){
    //Note: Psi has same boundary conditions as rho
    
    //Apply periodic boundary conditions, to left and right boundaries
    for (int j=0; j<int(UP[1].size()); j++){
        for (int k=0; k<8; k++){
            UP[0][j][k] = UP[xCells][j][k];
            UP[1][j][k] = UP[xCells + 1][j][k];
            UP[xCells + 2][j][k] = UP[2][j][k];
            UP[xCells + 3][j][k] = UP[3][j][k];
            }
    }
    
    //Apply periodic boundary conditions, to top and bottom boundaries
    for (int i=0; i<int(UP.size()); i++){
        for (int k=0; k<8; k++){
            UP[i][0][k] = UP[i][yCells][k];
            UP[i][1][k] = UP[i][yCells + 1][k];
            UP[i][yCells + 2][k] = UP[i][2][k];
            UP[i][yCells + 3][k] = UP[i][3][k];
            }
    }
}
void KH_bc(std::vector< std::vector <std::array<double,8> > > &UP, int &xCells, int &yCells){
    //Note: Psi has same boundary conditions as rho
    
    //Apply periodic boundary conditions, to left and right boundaries
    for (int j=0; j<int(UP[1].size()); j++){
        for (int k=0; k<8; k++){
            UP[0][j][k] = UP[xCells][j][k];
            UP[1][j][k] = UP[xCells + 1][j][k];
            UP[xCells + 2][j][k] = UP[2][j][k];
            UP[xCells + 3][j][k] = UP[3][j][k];
            }
    }

    //Apply reflexive boundary conditions, to top and bottom boundaries
    for(int i = 0; i < int(UP.size()); i++){
        // ghost vx = - vx
        // ghost Bx = - Bx
          UP[i][0][2] = (-1.0)*UP[i][3][2];
          UP[i][0][6] = (-1.0)*UP[i][3][6];
          
          UP[i][1][2] = (-1.0)*UP[i][2][2];
          UP[i][1][6] = (-1.0)*UP[i][2][6];
          
          
          UP[i][yCells + 3][2] = (-1.0)*UP[i][yCells][2];
          UP[i][yCells + 3][6] = (-1.0)*UP[i][yCells][6];
           
          UP[i][yCells + 2][2] = (-1.0)*UP[i][yCells + 1][2];
          UP[i][yCells + 2][6] = (-1.0)*UP[i][yCells + 2][6];
          
          for(int var = 0; var < 9; var++)
          {
              // ghost variable = variable apart from vx and Bx
              if (var!=2 && var!=6)
              {
                  UP[i][0][var] = UP[i][3][var];
                  UP[i][1][var] = UP[i][2][var];
                  UP[i][yCells + 3][var] = UP[i][yCells][var];
                  UP[i][yCells + 2][var] = UP[i][yCells + 1][var];
              }
            }
        }
}

void ToroTest1_x(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain){
    
    for (int i=0; i<int(xydomain.size()); i++)
    {
        for (int j=0; j<int(xydomain[1].size()); j++)
        {
            if (xydomain[i][j][0] <= 0.5) //if x<=0.5
            {
                UP[i][j][0] = 1.0; //rho
                UP[i][j][1] = 0.0; //vx
                UP[i][j][2] = 0.0; //vy
                UP[i][j][3] = 0.0; //vz
                UP[i][j][4] = 1.0; //p
                UP[i][j][5] = 0.0; //Bx
                UP[i][j][6] = 0.0; //By
                UP[i][j][7] = 0.0; //Bz
            }
            else if (xydomain[i][j][0] > 0.5) //if x>0.5
            {
                UP[i][j][0] = 0.125; //rho
                UP[i][j][1] = 0.0; //vx
                UP[i][j][2] = 0.0; //vy
                UP[i][j][3] = 0.0; //vz
                UP[i][j][4] = 0.1; //p
                UP[i][j][5] = 0.0; //Bx
                UP[i][j][6] = 0.0; //By
                UP[i][j][7] = 0.0; //Bz
            }
        }
    }
}




void BrioWu_x(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain){
    
    for (int i=0; i<int(xydomain.size()); i++)
    {
        for (int j=0; j<int(xydomain[1].size()); j++)
        {
            if (xydomain[i][j][0] <= 400) //if x<=400
            {
                UP[i][j][0] = 1.0; //rho
                UP[i][j][1] = 0.0; //vx
                UP[i][j][2] = 0.0; //vy
                UP[i][j][3] = 0.0; //vz
                UP[i][j][4] = 1.0; //p
                UP[i][j][5] = 0.75; //Bx
                UP[i][j][6] = 1.0; //By
                UP[i][j][7] = 0.0; //Bz
            }
            else if (xydomain[i][j][0] > 400) //if x>400
            {
                UP[i][j][0] = 0.125; //rho
                UP[i][j][1] = 0.0; //vx
                UP[i][j][2] = 0.0; //vy
                UP[i][j][3] = 0.0; //vz
                UP[i][j][4] = 0.1; //p
                UP[i][j][5] = 0.75; //Bx
                UP[i][j][6] = -1.0; //By
                UP[i][j][7] = 0.0; //Bz
            }
        }
    }
}

void BrioWu_y(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain){
    
    for (int i=0; i<int(xydomain.size()); i++)
    {
        for (int j=0; j<int(xydomain[1].size()); j++)
        {
            if (xydomain[i][j][1] <= 400) //if x<=400
            {
                UP[i][j][0] = 1.0; //rho
                UP[i][j][1] = 0.0; //vx
                UP[i][j][2] = 0.0; //vy
                UP[i][j][3] = 0.0; //vz
                UP[i][j][4] = 1.0; //p
                UP[i][j][5] = 1.0; //Bx
                UP[i][j][6] = 0.75; //By
                UP[i][j][7] = 0.0; //Bz
            }
            else if (xydomain[i][j][1] > 400) //if x>400
            {
                UP[i][j][0] = 0.125; //rho
                UP[i][j][1] = 0.0; //vx
                UP[i][j][2] = 0.0; //vy
                UP[i][j][3] = 0.0; //vz
                UP[i][j][4] = 0.1; //p
                UP[i][j][5] = -1.0; //Bx
                UP[i][j][6] = 0.75; //By
                UP[i][j][7] = 0.0; //Bz
            }
        }
    }
}

void DaiWoodward_x(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain){
    
    for (int i=0; i<int(xydomain.size()); i++)
    {
        for (int j=0; j<int(xydomain[1].size()); j++)
        {
            if (xydomain[i][j][0] <= 0.5) //if x<=0.5
            {
                UP[i][j][0] = 1.0; //rho
                UP[i][j][1] = 10.0; //vx
                UP[i][j][2] = 0.0; //vy
                UP[i][j][3] = 0.0; //vz
                UP[i][j][4] = 20.0; //p
                UP[i][j][5] = 5.0/sqrt(4.0*M_PI); //Bx
                UP[i][j][6] = 0.0/sqrt(4.0*M_PI); //By
                UP[i][j][7] = 5.0/sqrt(4.0*M_PI); //Bz
            }
            else if (xydomain[i][j][0] > 0.5) //if x>0.5
            {
                UP[i][j][0] = 1.0; //rho
                UP[i][j][1] = -10.0; //vx
                UP[i][j][2] = 0.0; //vy
                UP[i][j][3] = 0.0; //vz
                UP[i][j][4] = 1.0; //p
                UP[i][j][5] = 5.0/sqrt(4.0*M_PI); //Bx
                UP[i][j][6] = 0.0/sqrt(4.0*M_PI); //By
                UP[i][j][7] = 5.0/sqrt(4.0*M_PI); //Bz
            }
        }
    }
}

void OrzagTang(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain){
   
    double gamma = 5.0/3.0;
    
    for (int i=0; i<int(xydomain.size()); i++)
    {
        for (int j=0; j<int(xydomain[1].size()); j++)
        {
            double x = xydomain[i][j][0];
            double y = xydomain[i][j][1];
            
            UP[i][j][0] = gamma*gamma;             //rho
            UP[i][j][1] = -sin(2.0*M_PI*y);     //vx
            UP[i][j][2] = sin(2.0*M_PI*x) ;     //vy
            UP[i][j][3] = 0.0;                     //vz
            UP[i][j][4] = gamma;                 //p
            UP[i][j][5] = -sin(2.0*M_PI*y);     //Bx
            UP[i][j][6] = sin(4.0*M_PI*x);         //By
            UP[i][j][7] = 0.0;                     //Bz
            
        }
    }
}

void KelvinHelmholtz(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain){
    
    double sigma = 0.1;
    double x0 = 0.0;
    double y0 = -1.0;
    double x1 = 1.0;
    double y1 = 1.0;
    double dx = (x1-x0)/256.0;
    double dy = (y1-y0)/512.0;
    
    for (int i=0; i<int(xydomain.size()); i++)
    {
        for (int j=0; j<int(xydomain[1].size()); j++)
        {
            //double x = xydomain[i][j][0];
            //double y = xydomain[i][j][1];
            double x = x0 + (i-1.5)*dx;
            double y = y0 + (j-1.5)*dy;
            
            UP[i][j][0] = 1.0;                     //rho
            UP[i][j][1] = (0.5)*(tanh(20.0*y)); //vx
            UP[i][j][2] = (0.01)*sin(2.0*M_PI*x)*exp((-1.0)*((y*y)/(sigma*sigma)));    //vy
            //UP[i][j][2] = 0.0;
            UP[i][j][3] = 0.0;                     //vz
            UP[i][j][4] = 3.0/5.0;                 //p
            UP[i][j][5] = (0.1)*cos(M_PI/3.0);     //Bx
            UP[i][j][6] = 0.0;                    //By
            UP[i][j][7] = (0.1)*sin(M_PI/3.0);  //Bz
            
        }
    }
}

void ToroTest1_y(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain){
    
    for (int i=0; i<int(xydomain.size()); i++)
    {
        for (int j=0; j<int(xydomain[1].size()); j++)
        {
            if (xydomain[i][j][1] <= 0.5) //if y<=0.5
            {
                UP[i][j][0] = 1.0; //rho
                UP[i][j][1] = 0.0; //vx
                UP[i][j][2] = 0.0; //vy
                UP[i][j][3] = 0.0; //vz
                UP[i][j][4] = 1.0; //p
                UP[i][j][5] = 0.0; //Bx
                UP[i][j][6] = 0.0; //By
                UP[i][j][7] = 0.0; //Bz
            }
            else if (xydomain[i][j][1] > 0.5) //if y>0.5
            {
                UP[i][j][0] = 0.125; //rho
                UP[i][j][1] = 0.0; //vx
                UP[i][j][2] = 0.0; //vy
                UP[i][j][3] = 0.0; //vz
                UP[i][j][4] = 0.1; //p
                UP[i][j][5] = 0.0; //Bx
                UP[i][j][6] = 0.0; //By
                UP[i][j][7] = 0.0; //Bz
            }
        }
    }
}

void ToroTest1_nonalign(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain){
    
    for (int i=0; i<int(xydomain.size()); i++)
    {
        for (int j=0; j<int(xydomain[1].size()); j++)
        {
            if (xydomain[i][j][0] < (1.0-xydomain[i][j][1])) //if x<1-y
            {
                UP[i][j][0] = 1.0; //rho
                UP[i][j][1] = 0.0; //vx
                UP[i][j][2] = 0.0; //vy
                UP[i][j][3] = 0.0; //vz
                UP[i][j][4] = 1.0; //p
                UP[i][j][5] = 0.0; //Bx
                UP[i][j][6] = 0.0; //By
                UP[i][j][7] = 0.0; //Bz

            }
            else if ((1.0-xydomain[i][j][1]) <= xydomain[i][j][0])  //if x=>1-y
            {
                UP[i][j][0] = 0.125; //rho
                UP[i][j][1] = 0.0; //vx
                UP[i][j][2] = 0.0; //vy
                UP[i][j][3] = 0.0; //vz
                UP[i][j][4] = 0.1; //p
                UP[i][j][5] = 0.0; //Bx
                UP[i][j][6] = 0.0; //By
                UP[i][j][7] = 0.0; //Bz

            }
        }
    }
}

void CylindricalExplosion(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain){
    double R=0.4;
    for (int i=0; i<int(xydomain.size()); i++)
    {
        for (int j=0; j<int(xydomain[1].size()); j++)
        {
            double x = xydomain[i][j][0];
            double y = xydomain[i][j][1];
            
            if ( ((x-1.0)*(x-1.0) + (y-1.0)*(y-1.0)) <= R*R ) //if (x-1)**2 + (y-1)**2 <= R*R
            {
                UP[i][j][0] = 1.0; //rho
                UP[i][j][1] = 0.0; //vx
                UP[i][j][2] = 0.0; //vy
                UP[i][j][3] = 0.0; //vz
                UP[i][j][4] = 1.0; //p
                UP[i][j][5] = 0.0; //Bx
                UP[i][j][6] = 0.0; //By
                UP[i][j][7] = 0.0; //Bz

            }
            else if ( ((x-1.0)*(x-1.0) + (y-1.0)*(y-1.0)) > R*R) //if (x-1)**2 + (y-1)**2 > R*R
            {
                UP[i][j][0] = 0.125; //rho
                UP[i][j][1] = 0.0; //vx
                UP[i][j][2] = 0.0; //vy
                UP[i][j][3] = 0.0; //vz
                UP[i][j][4] = 0.1; //p
                UP[i][j][5] = 0.0; //Bx
                UP[i][j][6] = 0.0; //By
                UP[i][j][7] = 0.0; //Bzs
            }
        }
    }
}

void BrioWu_nonalign(std::vector< std::vector <std::array<double,8> > > &UP, std::vector< std::vector <std::array<double,2> > > &xydomain){
    
    for (int i=0; i<int(xydomain.size()); i++)
    {
        for (int j=0; j<int(xydomain[1].size()); j++)
        {
            if (xydomain[i][j][1] < (400 -  xydomain[i][j][0])) //if x<=400
            {
                UP[i][j][0] = 1.0; //rho
                UP[i][j][1] = 0.0; //vx
                UP[i][j][2] = 0.0; //vy
                UP[i][j][3] = 0.0; //vz
                UP[i][j][4] = 1.0; //p
                UP[i][j][5] = 0.75; //Bx
                UP[i][j][6] = 1.0; //By
                UP[i][j][7] = 0.0; //Bz
            }
            else //if x>400
            {
                UP[i][j][0] = 0.125; //rho
                UP[i][j][1] = 0.0; //vx
                UP[i][j][2] = 0.0; //vy
                UP[i][j][3] = 0.0; //vz
                UP[i][j][4] = 0.1; //p
                UP[i][j][5] = 0.75; //Bx
                UP[i][j][6] = -1.0; //By
                UP[i][j][7] = 0.0; //Bz
            }
        }
    }
}



void ReadParams(int &testnum, int &limiternum, int &bconditions, int &xCells, int &yCells, double &x0, double &x1, double &y0, double &y1, double &tstop, double &C, double &gamma, double &omega){
    std::string temp;
    
    std::ifstream file;
    //Change name of settings file
    file.open("inputs.txt");
    
    std::getline(file, temp);
    std::istringstream input0(temp);
    input0 >> testnum;
    
    std::getline(file, temp);
    std::istringstream input1(temp);
    input1 >> limiternum;
    
    std::getline(file, temp);
    std::istringstream input2(temp);
    input2 >> bconditions;
    
    std::getline(file, temp);
    std::istringstream input3(temp);
    input3 >> xCells;
    
    std::getline(file, temp);
    std::istringstream input4(temp);
    input4 >> yCells;
    
    std::getline(file, temp);
    std::istringstream input5(temp);
    input5 >> x0;
    
    std::getline(file, temp);
    std::istringstream input6(temp);
    input6 >> x1;
    
    std::getline(file, temp);
    std::istringstream input7(temp);
    input7 >> y0;
    
    std::getline(file, temp);
    std::istringstream input8(temp);
    input8 >> y1;
    
    std::getline(file, temp);
    std::istringstream input9(temp);
    input9 >> tstop;
    
    std::getline(file, temp);
    std::istringstream input10(temp);
    input10 >> C;
    
    std::getline(file, temp);
    std::istringstream input11(temp);
    input11 >> gamma;
    
    std::getline(file, temp);
    std::istringstream input12(temp);
    input12 >> omega;
}

