#include<iostream> 
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<array>
#include<algorithm>
#include<functional>
#include<iostream>
#include<fstream>
#include "1dmhdfunctions.h"
#include "1dmhdtests.h"

/*
 File Description: 1D MHD HLLC approximate solver using MUSCL-Hancock with slope limiters
 Date: Monday 25th January 2021
 */

int main(void){
	
	//---------------------INITIALIZE PARAMETERS------------------------
	//------------------------------------------------------------------
    std::string testName, limiterName;
    int test, limiter, ncells;
    double x0, x1, tstop, C, gamma, omega, dx;
    double tstart = 0.0;
    double dt;
    
    // Function to read in parameters
    ReadParams(test, limiter, ncells, x0, x1, tstop, C, gamma, omega);
    
    dx = (x1 - x0) / ncells;
    
    std::function<void(double&, double&,std::vector<std::array<double,8> >&)> initialiseFn;
    switch(test){
        case 1:
            initialiseFn = test_ToroTest1;
            testName = "TEST1";
            break;
        case 2:
            initialiseFn = test_modifiedSod;
            testName = "TEST2";
            break;
        case 3:
            initialiseFn = test_BrioWu;
            testName = "TEST3";
            break;
        case 4:
            initialiseFn = test_BrioWuSwitched;
            testName = "TEST4";
            break;
        case 5:
            initialiseFn = test_DaiWoodward;
            testName = "TEST5";
            break;
        case 6:
            initialiseFn = test_DaiWoodwardSwitched;
            testName = "TEST6";
            break;
        case 7:
            initialiseFn = test_Li25D;
            testName = "TEST7";
            break;
    }
    
    std::function<double(const int&,const int&,const std::vector<std::array<double,8> > &)> limiterFn;
    switch(limiter){
        case 1:
            limiterFn = Minbeelimiter;
            limiterName = "Minbee";
            break;
        case 2:
            limiterFn = Superbeelimiter;
            limiterName = "Superbee";
            break;
        case 3:
            limiterFn = VanLeerlimiter;
            limiterName = "Van-Leer";
            break;
        case 4:
            limiterFn = VanAlbadalimiter;
            limiterName = "Vanalbada";
            break;
    }

    std::cout << "***SUMMARY***" <<std::endl;
    std::cout << "Test Name = " << testName << std::endl;
    std::cout << "Limiter Name = " << limiterName << std::endl;
    std::cout << "ncells =  " << ncells << std::endl;
    std::cout << "Domain Min = " << x0 << std::endl;
    std::cout << "Domain Max = " << x1 << std::endl;
    std::cout << "tStop = " << tstop << std::endl;
    std::cout << "C = " << C << std::endl;
    std::cout << "omega = " << omega << std::endl;
    std::cout << "gamma = " << gamma << std::endl;
    std::cout << "***SUMMARY***" <<std::endl;
    
    double Coriginal = C;
    double Ctemp = 0.2*Coriginal;
    
    
	//---------------------INITIAL CONDITIONS --------------------------
	//------------------------------------------------------------------
	std::vector<std::array<double, 8> > U_P (ncells + 4); //Primitive Variables (2 dummy cells either side)

	// Set up the initial data
    initialiseFn(x0,dx,U_P);
    
    std::string resolution = std::to_string(ncells);
	// Output the initial data for plotting
	std::ofstream output_initial(testName+"_"+resolution+"_initial.txt");
	
	// Don't include extra ghost cell into initial data output
	for (int i = 2; i < int(U_P.size())-2; i++)
	{
			double x = x0 + (i-1.5)*dx;		
			output_initial << x << " " << U_P[i][0] << " " << U_P[i][1] << " " << U_P[i][2];
			output_initial << " " << U_P[i][3] << " " << U_P[i][4] << " " << U_P[i][5] << " " << U_P[i][6] << " " << U_P[i][7] << std::endl;
	}
	output_initial.close();
	
	//---------------------INITIALIZE VECTORS---------------------------
	//------------------------------------------------------------------
	std::vector<std::array<double, 8> > U_C(ncells + 4); //Conserved variables
	std::vector<std::array<double, 8> > U_Ciplus1(ncells + 4); // Updated Conserved variables
    std::vector<std::array<double, 8> > Uleft(ncells + 4); //Left Data Construction
    std::vector<std::array<double, 8> > Uright(ncells + 4); //Right Data Construction
    std::vector<std::array<double, 8> > Uleft_nph(ncells + 4); //Boundary states UL
    std::vector<std::array<double, 8> > Uright_nph(ncells + 4); //Boundary states UR
    std::vector<std::array<double, 8> > HLLCflux(ncells + 1); //Vector for HLLC flux
    
	//---------------------TIME MARCHING--------------------------------
	//------------------------------------------------------------------
    int itercount = 0;
	
    // Find the solution to time t = tStop by looping over discrete times
	double t = tstart;
	do{
        
        //Small CFL number  for first 10 iteration to remove start up error at discontinuities
        if(itercount <= 10){
            C = Ctemp;
        }
        else{
            C = Coriginal;
        }
        
        itercount += 1;
        
        Transmissiveboundaries(U_P,ncells);
		
        //Compute time step dt and output to terminal every 10 iterations
		dt = computeTimeStep(U_P, C, dx, gamma, ncells);
        t = t + dt;
        
        if(itercount%10 == 0){
            std::cout <<"dt = "<< dt << std::endl;
            std::cout <<"t = "<< t << std::endl;
        }

		// Manually reduce dt if this overshoots tstop
		if( t>tstop ){
			// Go back to tstop and reduce dt
			dt = dt - abs(t-tstop);
			t = tstop;
			std::cout << "Adjust final t to " << t << std::endl;
		}
		
		//Convert primitive variables to conservative variables
		for (int i=0; i<ncells+4; i++){
			U_C[i]= Primitive_to_Conservative(U_P[i], gamma);
		}
        
        //Data reconstruction to get boundary extrapolated values
        for(int i=1; i< int(U_C.size())-1; i++){
            for(int j=0; j<8; j++){
                Uleft[i][j] = U_C[i][j] - 0.5 * limiterFn(i,j,U_C) * (0.5 *(1.0 + omega)*(U_C[i][j] - U_C[i-1][j]) + 0.5 * (1.0 - omega)*(U_C[i+1][j] - U_C[i][j]));
                Uright[i][j] = U_C[i][j] + 0.5 * limiterFn(i,j,U_C) * (0.5 *(1.0 + omega)*(U_C[i][j] - U_C[i-1][j]) + 0.5 * (1.0 - omega)*(U_C[i+1][j] - U_C[i][j]));
            }
        }
        
		// Get the half time step evolution
		for (int i=0; i<int(Uleft.size()); i++){
			
			//Get mhd flux of ulineleft and ulineright
			std::array<double, 8> f_Uleft;
			std::array<double, 8> f_Uright;
			
			f_Uleft = f_mhd(Uleft[i], gamma);
			f_Uright = f_mhd(Uright[i], gamma);
			
            // Note this is same formula as Page 501 Toro (14.3) but with negative sign
			for(int j=0; j<8; j++){
				Uleft_nph[i][j] = Uleft[i][j] - 0.5*(dt/dx)*(f_Uright[j] - f_Uleft[j]);
				Uright_nph[i][j] = Uright[i][j] - 0.5*(dt/dx)*(f_Uright[j] - f_Uleft[j]);
			}
		}

		//Compute HLLC flux
        for (int i=0; i < int(HLLCflux.size()); i++){
            HLLCflux[i] = getHLLCflux(Uright_nph[i+1], Uleft_nph[i+2],gamma);
        }
        
		// Update the data using forward difference
        for (int i=2; i<int(U_C.size())-2; i++){
            for(int j=0; j<8; j++){
                U_Ciplus1[i][j] = U_C[i][j] - dt/dx * (HLLCflux[i-1][j] - HLLCflux[i-2][j]);
            }
        }
        
        // Apply boundary condition again here to UCiplus1s
        Transmissiveboundaries(U_Ciplus1,ncells);
        
        // Update conservative variables
        U_C = U_Ciplus1;
    
		//Convert conservative variables to primitive variables
		for (int i=0; i<ncells+4; i++){
			U_P[i]= Conservative_to_Primitive(U_C[i], gamma);
		}
		
    }while (t < tstop);
		
	//--------------------------OUTPUT FINAL DATA-------------------
	//--------------------------------------------------------------
	// Output the final data for plotting
	std::ofstream output_final(testName+"_"+resolution + "_final.txt");
	
	// Don't include extra ghost cell into final data output	
	for (int i = 2; i < int(U_P.size())-2; i++)
	{
			double x = x0 + (i-1.5)*dx;		
			output_final << x << " " << U_P[i][0] << " " << U_P[i][1] << " " << U_P[i][2];
			output_final << " " << U_P[i][3] << " " << U_P[i][4] << " " << U_P[i][5] << " " << U_P[i][6] << " " << U_P[i][7] << " " << U_P[i][4]/(U_P[i][0]*(gamma-1.0)) << std::endl;
	}
	
	//Close output file
	output_final.close();
	
return 0;
}


