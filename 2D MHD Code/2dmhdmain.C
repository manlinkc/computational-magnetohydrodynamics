#include<iostream> 
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<array>
#include<algorithm>
#include<functional>
#include<iostream>
#include<sstream>
#include<fstream>
#include "2dmhdfunctions.h"

/*
 File Description: 2D MHD
 */

/*
 * TEST1 = Toro's 1st test t=0.25
 * TEST2 = Brio and Wu test t = 80
 * TEST3 = Brio and Wu x-y rotated t = 80
 * TEST4 = Dai and Woodward
 * TEST5 = Orzag Tang
 * TEST6 = Kelvin Helmholtz
 */

int main(void){
	//---------------------INITIALIZE PARAMETERS------------------------
	//------------------------------------------------------------------
    std::string test, limiter;
    int testnum, limiternum, bconditions, xCells, yCells;
    double x0, x1, y0, y1, tstop, C, gamma, omega;
    double tstart = 0.0;
    double dt;
    
    ReadParams(testnum, limiternum, bconditions, xCells, yCells, x0, x1, y0, y1, tstop, C, gamma, omega);

    
    const double dx = (x1 - x0)/ xCells;
    const double dy = (y1 - y0)/ yCells;
    
    std::function<void(std::vector<std::vector<std::array<double, 8> > >&,std::vector<std::vector<std::array<double, 2> > >&)> initialiseFn;
    switch(testnum){
        case 1:
            initialiseFn = ToroTest1_x;
            test = "TEST1";
            break;
        case 2:
            initialiseFn = BrioWu_x;
            test = "TEST2";
            break;
        case 3:
            initialiseFn = DaiWoodward_x;
            test = "TEST3";
            break;
        case 4:
            initialiseFn = OrzagTang;
            test = "TEST4";
            break;
        case 5:
            initialiseFn = KelvinHelmholtz;
            test = "TEST5";
            break;
        case 6:
            initialiseFn = CylindricalExplosion;
            test = "TEST6";
            break;
        case 7:
            initialiseFn = ToroTest1_y;
            test = "TEST7";
            break;
        case 8:
            initialiseFn = BrioWu_y;
            test = "TEST8";
            break;
        case 9:
            initialiseFn = ToroTest1_nonalign;
            test = "TEST9";
            break;
        case 10:
            initialiseFn = BrioWu_nonalign;
            test = "TEST10";
            break;
    }
    
    std::function<double(double &r)> limiterFn;
    switch(limiternum){
        case 1:
            limiterFn = Minbeelimiter;
            limiter = "Minbee";
            break;
        case 2:
            limiterFn = Superbeelimiter;
            limiter = "Superbee";
            break;
        case 3:
            limiterFn = VanLeerlimiter;
            limiter = "Vanleer";
            break;
        case 4:
            limiterFn = VanAlbadalimiter;
            limiter = "Vanalbada";
            break;
    }
    
    std::function<void(std::vector<std::vector<std::array<double, 8> > >&,int &, int &)> boundaryconditionsFn;
    switch(bconditions){
        case 1:
            boundaryconditionsFn = transmissive_bc;
            break;
        case 2:
            boundaryconditionsFn = periodic_bc;
            break;
        case 3:
            boundaryconditionsFn = KH_bc;
            break;
    }
    std::cout << "*** SUMMARY ***" << std::endl;
    std::cout << "Test = " << test << std::endl;
    std::cout << "Limiter = " << limiter << std::endl;
    std::cout << "Boundary Conditions = " << bconditions << std::endl;
    std::cout << "Number of xcells = " << xCells << std::endl;
    std::cout << "Number of ycells = " << yCells << std::endl;
    std::cout << "x0 = " << x0 << std::endl;
    std::cout << "x1 = " << x1 << std::endl;
    std::cout << "y0 = " << y0 << std::endl;
    std::cout << "y1 = " << y1 << std::endl;
    std::cout << "Final Time = " << tstop << std::endl;
    std::cout << "CFL = " << C << std::endl;
    std::cout << "Gamma = " << gamma << std::endl;
    std::cout << "*** SUMMARY ***" << std::endl;

    double Coriginal = C;
    double Ctemp = 0.2*Coriginal;
    
	//---------------------INITIAL CONDITIONS --------------------------
	//------------------------------------------------------------------
	// Define a vector to store primitive variables
	// Add two ghost cells either side for boundary conditions
	std::vector < std::vector < std::array<double, 8> > > UP;
	UP.resize(xCells + 4, std::vector< std::array<double,8 > >(yCells + 4) );
	
	//Construct domain
	std::vector < std::vector < std::array<double,2> > > xydomain;
	xydomain.resize(xCells +4, std::vector< std::array<double,2> >(yCells + 4) );
	
	for (int i=0; i<int(xydomain.size()); i++)
	{
		for (int j=0; j<int(xydomain.size()); j++)
		{
			xydomain[i][j][0] = x0 + (i-1.5)*dx;
			xydomain[i][j][1] = y0 + (j-1.5)*dy;
		}	
	}
    
	initialiseFn(UP,xydomain);

	std::ofstream output_initial(test+"_"+limiter+"_initial.txt");

	// Don't include extra ghost cell into initial data output
	for (int i = 2; i < int(UP.size()-2); i++){
			double x = x0 + (i-1.5)*dx;	
			for (int j = 2; j < int(UP[1].size()-2); j++){
				 double y = y0 + (j-1.5)*dy;
				 
                output_initial << x << " " << y << " " << UP[i][j][0] << " " << UP[i][j][1] << " " << UP[i][j][2] << " " << UP[i][j][3] << " " << UP[i][j][4] << " " << UP[i][j][5] << "  " << UP[i][j][6] << " " << UP[i][j][7] << std::endl;
			}
			//Add extra blank line
			output_initial << " " << std::endl;
	}
	
	//Close output file
	output_initial.close();

	
	//---------------------INITIALIZE VECTORS---------------------------
	//------------------------------------------------------------------
	// Define a vector to store conservative variables
	// Needs to be the same size as UP
	std::vector < std::vector < std::array<double, 8> > > UC;
	UC.resize(xCells + 4, std::vector< std::array<double,8 > >(yCells + 4) );
	
	// Define a vector to store the UC at time i+1 flux values
	std::vector < std::vector < std::array<double, 8> > > UCiplus1;
	UCiplus1.resize(xCells + 4, std::vector< std::array<double,8 > >(yCells + 4) );
	
	// Define a vector to store the UCline after x update
	std::vector < std::vector < std::array<double, 8> > > UCline;
	UCline.resize(xCells + 4, std::vector< std::array<double,8 > >(yCells + 4) );
	
	// Define a vector to store SLIC flux variables for f flux
	std::vector < std::vector < std::array<double, 8> > > SLICflux_f;
	SLICflux_f.resize(xCells + 1, std::vector< std::array<double,8 > >(yCells + 4) );
	
	// Define a vector to store SLIC flux variables for g flux
	std::vector < std::vector < std::array<double, 8> > > SLICflux_g;
	SLICflux_g.resize(xCells + 4, std::vector< std::array<double,8 > >(yCells + 1) );
	
    // Define a vector to store SLIC flux variables for f flux
    std::vector < std::vector < std::array<double, 8> > > HLLCflux_f;
    HLLCflux_f.resize(xCells + 1, std::vector< std::array<double,8 > >(yCells + 4) );
    
    // Define a vector to store SLIC flux variables for g flux
    std::vector < std::vector < std::array<double, 8> > > HLLCflux_g;
    HLLCflux_g.resize(xCells + 4, std::vector< std::array<double,8 > >(yCells + 1) );
    
	// Define a vector to store Uleft quantity
	std::vector < std::vector < std::array<double, 8> > > Uleft;
	Uleft.resize(xCells + 4, std::vector< std::array<double,8 > >(yCells + 4) );
	
	// Define a vector to store Ulineleft quantity
	std::vector < std::vector < std::array<double, 8> > > Ulineleft;
	Ulineleft.resize(xCells + 4, std::vector< std::array<double,8 > >(yCells + 4) );
	
	// Define a vector to store Uright quantity
	std::vector < std::vector < std::array<double, 8> > > Uright;
	Uright.resize(xCells + 4, std::vector< std::array<double,8 > >(yCells + 4) );
	
	// Define a vector to store Ulineright quantity
	std::vector < std::vector < std::array<double, 8> > > Ulineright;
	Ulineright.resize(xCells + 4, std::vector< std::array<double,8 > >(yCells + 4) );
	
	// Define a vector to store Uleft nplushalf quantity
	std::vector < std::vector < std::array<double, 8> > > Uleft_nph;
	Uleft_nph.resize(xCells + 4, std::vector< std::array<double,8 > >(yCells + 4) );
	
	// Define a vector to store Ulineleft nplushalf quantity
	std::vector < std::vector < std::array<double, 8> > > Ulineleft_nph;
	Ulineleft_nph.resize(xCells + 4, std::vector< std::array<double,8 > >(yCells + 4) );
	
	// Define a vector to store Uright nplushalf quantity
	std::vector < std::vector < std::array<double, 8> > > Uright_nph;
	Uright_nph.resize(xCells + 4, std::vector< std::array<double,8 > >(yCells + 4) );
	
	// Define a vector to store Ulineright nplushalf quantity
	std::vector < std::vector < std::array<double, 8> > > Ulineright_nph;
	Ulineright_nph.resize(xCells + 4, std::vector< std::array<double,8 > >(yCells + 4) );
	
	//---------------------TIME MARCHING--------------------------------
	//------------------------------------------------------------------
    int itercount = 0;
    
	// Now find the solution to time t = tStop by looping over 
	// discrete times
	double t = tstart;
	do{
        boundaryconditionsFn(UP,xCells, yCells);
		//----------------COMPUTE TIME STEP----------------------------
		//-------------------------------------------------------------
        //Small CFL number  for first 10 iteration to remove start up error at discontinuities
        if(itercount <= 10){
            C = Ctemp;
        }
        else{
            C = Coriginal;
        }
        
        itercount += 1;
        
        dt = computeTimeStep(UP, C, dx, dy, gamma);
        t = t + dt;
    
        if(itercount%10 == 0){
            std::cout <<"dt = "<< dt << std::endl;
            std::cout <<"t = "<< t << std::endl;
        }
		
		// May want to manually reduce dt is this would overshoot tstop
		if( t>tstop ){
			// Go back to tstop and reduce dt
			dt = dt - abs(t-tstop);
			t = tstop;
			std::cout << "Adjust final t to " << t << std::endl;
		}
		
		//-----------------BOUNDARY CONDITIONS-------------------------
		//-------------------------------------------------------------
		
		//Convert primitive variables to conservative variables
		for (int i=0; i<int(UP.size()); i++){
			for (int j=0; j<int(UP[1].size()); j++){
				UC[i][j]= Primitive_to_Conservative(UP[i][j], gamma);
			}
		}
		
		
		//-----------------UPDATE X DIRECTION USING SLIC---------------------
		//-------------------------------------------------------------------
        for(int i=1; i< int(UC.size())-1; i++){
            for(int j=0; j< int(UC[1].size()); j++){
                for(int var=0; var<8; var++){

                        double r = computeR_x(i,j,var,UC);
                        double Delta = limiterFn(r)*(0.5*(1.0 + omega)*(UC[i][j][var] - UC[i-1][j][var]) + 0.5*(1.0 - omega)*(UC[i+1][j][var] - UC[i][j][var]));
                        
                        Uleft[i][j][var] = UC[i][j][var] - 0.5*Delta;
                        Uright[i][j][var] = UC[i][j][var]+ 0.5*Delta;
                }
            }
        }
        
		// Get the half time step evolution
		for (int i=0; i<int(Uleft.size()); i++){
			for (int j=0; j<int(Uleft[1].size()); j++){
				
				//Get euler flux of ulineleft and ulineright
				std::array<double, 8> f_Uleft;
				std::array<double, 8> f_Uright;
				
				f_Uleft = f_euler_xdirection(Uleft[i][j], gamma);
				f_Uright = f_euler_xdirection(Uright[i][j], gamma);
				
				// Get u left n plus a half and u right n plus a half
				for (int var=0; var<8; var++){
					Uleft_nph[i][j][var] = Uleft[i][j][var] - 0.5*(dt/dx)*(f_Uright[var] - f_Uleft[var]);
					Uright_nph[i][j][var] = Uright[i][j][var] - 0.5*(dt/dx)*(f_Uright[var] - f_Uleft[var]);
				}
			}
		}
        
        for (int i=0; i< int(HLLCflux_f.size()); i++){
            for (int j=0; j< int(HLLCflux_f[1].size()); j++){
                HLLCflux_f[i][j] = getHLLCflux_xdirection(Uright_nph[i+1][j], Uleft_nph[i+2][j],gamma);
            }
        }
        // Update equation for x direction
        // Don't include boundary conditions
        for (int i=2; i<int(UC.size())-2; i++){
            for (int j=0; j<int(UC[1].size()); j++){
                for (int var=0; var<8; var++){
                    UCline[i][j][var] = UC[i][j][var] - dt/dx * (HLLCflux_f[i-1][j][var] - HLLCflux_f[i-2][j][var]);
                }
            }
        }
         
        /*
        //Get the fluxes across spatial domain, uses conservative flux as defined in lectures
        //Only compute flux on real domain not ghost cells
        for (int i=0; i< int(SLICflux_f.size()); i++){
            for (int j=0; j< int(SLICflux_f[1].size()); j++){
                SLICflux_f[i][j] = getFORCEflux_xdirection(Uright_nph[i+1][j], Uleft_nph[i+2][j],dx,dt,gamma);
            }
        }
         
		// Update equation for x direction 
		// Don't include boundary conditions
		for (int i=2; i<int(UC.size())-2; i++){
			for (int j=0; j<int(UC[1].size()); j++){
                for (int var=0; var<8; var++){
                    UCline[i][j][var] = UC[i][j][var] - dt/dx * (SLICflux_f[i-1][j][var] - SLICflux_f[i-2][j][var]);
                }
            }
		}
		*/
        
        boundaryconditionsFn(UCline,xCells, yCells);
		//-----------------UPDATE Y DIRECTION USING SLIC---------------------
		//-------------------------------------------------------------------
        // Get ulineleft
        for(int j=1; j< int(UC[1].size())-1; j++){
            for(int i=0; i< int(UC.size()); i++){
                for(int var=0; var<8; var++){
                    double r = computeR_y(i,j,var,UCline);
                    double Delta = limiterFn(r)*(0.5*(1.0 + omega)*(UCline[i][j][var] - UCline[i][j-1][var]) + 0.5*(1.0 - omega)*(UCline[i][j+1][var] - UCline[i][j][var]));
                    
                    Ulineleft[i][j][var] = UCline[i][j][var] - 0.5*Delta;
                    Ulineright[i][j][var] = UCline[i][j][var] + 0.5*Delta;
                }
            }
        }
        
		// Get the half time step evolution
		for (int i=0; i<int(Ulineleft.size()); i++){
			for (int j=0; j<int(Ulineleft[1].size()); j++){
				
				//Get euler flux of ulineleft and ulineright
				std::array<double, 8> f_Ulineleft;
				std::array<double, 8> f_Ulineright;
				
				//Now need to use the conservative flux in the y direction
				f_Ulineleft = f_euler_ydirection(Ulineleft[i][j], gamma);
				f_Ulineright = f_euler_ydirection(Ulineright[i][j], gamma);
				
				// Get u left n plus a half and u right n plus a half
				for (int var=0; var<8; var++){
					Ulineleft_nph[i][j][var] = Ulineleft[i][j][var] - 0.5*(dt/dy)*(f_Ulineright[var] - f_Ulineleft[var]);
					Ulineright_nph[i][j][var] = Ulineright[i][j][var] - 0.5*(dt/dy)*(f_Ulineright[var] - f_Ulineleft[var]);
				}
			}
		}
		/*
		//Get the fluxes across spatial domain, uses conservative flux as defined in lectures
		//Only compute flux on real domain not ghost cells
		for (int i=0; i< int(SLICflux_g.size()); i++){
			for (int j=0; j< int(SLICflux_g[1].size()); j++){
				SLICflux_g[i][j] = getFORCEflux_ydirection(Ulineright_nph[i][j+1], Ulineleft_nph[i][j+2],dy,dt,gamma);
			}
		}
    
		// Update equation for y direction 
		// Don't include boundary conditions
		for (int j=2; j<int(UCiplus1[1].size())-2; j++){
			for (int i=0; i<int(UCiplus1.size()); i++){
                for (int var=0; var<8; var++){
                    UCiplus1[i][j][var] = UCline[i][j][var] - dt/dy * (SLICflux_g[i][j-1][var] - SLICflux_g[i][j-2][var]);
                }
          
            }
		}
		*/
        
        
        //Get the fluxes across spatial domain, uses conservative flux as defined in lectures
        //Only compute flux on real domain not ghost cells
        for (int i=0; i< int(HLLCflux_g.size()); i++){
            for (int j=0; j< int(HLLCflux_g[1].size()); j++){
                HLLCflux_g[i][j] = getHLLCflux_ydirection(Ulineright_nph[i][j+1], Ulineleft_nph[i][j+2],gamma);
            }
        }
        // Update equation for y direction
        // Don't include boundary conditions
        for (int j=2; j<int(UCiplus1[1].size())-2; j++){
            for (int i=0; i<int(UCiplus1.size()); i++){
                for (int var=0; var<8; var++){
                    UCiplus1[i][j][var] = UCline[i][j][var] - dt/dy * (HLLCflux_g[i][j-1][var] - HLLCflux_g[i][j-2][var]);
                }
          
            }
        }
        

        boundaryconditionsFn(UCiplus1,xCells, yCells);
		//-----------------------UPDATE U's FOR NEXT ITERATION----------------
		//--------------------------------------------------------------------
		
		// Update the conservative variables
		// Only need to update on the real domain 
		for (int i=0; i<int(UC.size()); i++){
			for (int j=0; j<int(UC[1].size()); j++){
                for (int var=0; var<8; var++){
                    UC[i][j][var] = UCiplus1[i][j][var];
                }
            }
		}
		
		//Convert conservative variables to primitive variables
		for (int i=0; i<int(UC.size()); i++){
			for (int j=0; j<int(UC[1].size()); j++){
				UP[i][j]= Conservative_to_Primitive(UC[i][j], gamma);
			}
		}
	
	}while (t < tstop);
	
	
	//--------------------------OUTPUT FINAL DATA-------------------
	//--------------------------------------------------------------
	// Output the initial data for plotting later
	std::ofstream output_final(test+"_"+limiter+"_final.txt");
	

	// Don't include extra ghost cell into initial data output
	for (int i = 2; i < int(UP.size()-2); i++){
			double x = x0 + (i-1.5)*dx;	
			for (int j = 2; j < int(UP[1].size()-2); j++){
				 double y =y0 + (j-1.5)*dy;
				 double ratio = sqrt(UP[i][j][5]*UP[i][j][5]+UP[i][j][6]*UP[i][j][6])/UP[i][j][7];
                 double internalenergy = UP[i][j][4]/(UP[i][j][0]*(gamma - 1.0));
				 output_final << x << " " << y << " " << UP[i][j][0] << " " << UP[i][j][1] << " " << UP[i][j][2] << " " << UP[i][j][3] << " " << UP[i][j][4] << " " << UP[i][j][5] << " " << UP[i][j][6] << " " << UP[i][j][7] << " " << ratio << " " << internalenergy << std::endl;
			}
			//Add extra blank line
			output_final << " " << std::endl;
	}

	//Close output file
	output_final.close();
    
    double L2error = 0.0;
    double maxnabla = 0.0;
    // Loop over all real cells
    for (int i = 2; i < int(UP.size()-2); i++){
            for (int j = 2; j < int(UP[1].size()-2); j++){
                double nablaB2 = ((UP[i+1][j][5] - UP[i-1][j][5])/(2.0*dx)) + ((UP[i][j+1][6] - UP[i][j-1][6])/(2.0*dy));
                maxnabla = std::max(maxnabla, fabs(nablaB2));
                L2error = L2error + fabs(nablaB2);
            }
    }
    std::cout << "L1error at t=final is " << L2error*dx*dy <<std::endl;
    std::cout << "Max div B = " << maxnabla <<std::endl;
	return 0;
}

