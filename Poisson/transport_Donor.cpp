#include<iostream>
#include <stdlib.h>
#include<string>
#include<math.h>
#include <fstream>
#include <iomanip>

using namespace std;
void advectDensity(double *rho_old,double *rho_new, double *Vx, double *Vy, int Nx, int Ny, double dt, double dl);

int main(int argc, char *argv[])
{
	int nx = pow(2,3), ny = pow(2,3);
	int Nx = nx+4,     Ny = ny+4;
	double *rho_old = new double[(Nx)*(Ny)];
	double *rho_new = new double[(Nx)*(Ny)];
	double *Vx = new double[(Nx)*(Ny)];
	double *Vy = new double[(Nx)*(Ny)];
	double *temp;
	int steps=2;
	double dl    = 1./Nx;
	double dt    = dl;
	//cout<<dl<<endl;
	//cout<<dt<<endl;
	for (int i=2; i<Nx-2; i++)
    {
		for (int j=2; j<Ny-2; j++)	
		{	
			if (j<Nx/2. )
				rho_old[i*Ny+j] = 10.;
			else
				rho_old[i*Ny+j] = 0.;
			
			Vx[i*Ny+j]  = 1.;
			Vy[i*Ny+j]  = 1.;
		} 
		
    }
    	//Implement boundary conditions in Ghost Zones-ooooo	
	//Y direction
	for (int i=2; i<Nx-2; i++)
	{
		rho_old[i*Ny+0]    = rho_old[i*Ny+Ny-4]; //bottom row 
		rho_old[i*Ny+1]    = rho_old[i*Ny+Ny-3];
		rho_old[i*Ny+Ny-2] = rho_old[i*Ny+2];	 //top row
		rho_old[i*Ny+Ny-1] = rho_old[i*Ny+3];
	}
	//Y Direction

	for (int j=2; j<Ny-2; j++)
	{
		rho_old[0*Ny+j]      = rho_old[(Ny-4)*Ny+j]; //left column 
		rho_old[1*Ny+j]      = rho_old[(Ny-3)*Ny+j];
		rho_old[(Nx-2)*Ny+j] = rho_old[2*Ny+j];	//right column
		rho_old[(Nx-1)*Ny+j] = rho_old[3*Ny+j];
	}
    for (int i=0; i<steps; i++)
    {cout<<"Advecting"<<endl;
    	advectDensity(rho_old,rho_new, Vx, Vy, Nx, Ny,  dt, dl);
    	temp    = rho_old;
    	rho_old = rho_new;
    	rho_new=temp;
    	
    }
    
    ofstream myfile ("test.txt");
    for (int j=Nx-1; j>-1; j--)
    {
		for (int i=0; i<Ny; i++)
		{	myfile <<setprecision(2)<<rho_old[i*Ny+j]<<" ";
		} 
		myfile <<"\n";
    }
    myfile.close();
}

void advectDensity(double *rho_old,double *rho_new, double *Vx, double *Vy, int Nx, int Ny, double dt, double dl)
{
	double *rho_partial 	= new double[Nx*Ny];
	double *F          		= new double[(Nx)*(Ny)];
	double rho_interp   	= 0.;

	//Calculate x-fluxes everywhere
	for (int i=2; i<Nx-2; i++)
	{cout<<"i: "<<i<<endl;
		for (int j=2; j<Ny-2; j++)
		{	
		//Calculate interpolated density
			if (Vx[i*Ny+j] > 0.)
			{
				rho_interp=rho_old[(i-1)*Ny+j];
			}
			else
			{
				rho_interp=rho_old[i*Ny+j];
			}
			cout<<rho_interp<<" ";
			//Calculate F_1
			F[i*Ny+j] = rho_interp * Vx[i*Ny+j] ;
		
		} 		
	}//end calculating fluxes

	//Partial update density back into rho_partial
	for (int i=2; i<Nx-3; i++)
	{
		for (int j=2; j<Ny-2; j++)
		{
			rho_partial[i*Ny+j] = rho_old[i*Ny+j] - dt/dl * (F[(i+1)*Ny+j] - F[i*Ny+j]);
		}
	}
	// Boundary Fluxes
	for (int j=2; j<Nx-2; j++)
	{int i=Ny-3;
			rho_partial[i*Ny+j] = rho_old[i*Ny+j] - dt/dl * (F[2*Ny+j] - F[i*Ny+j]);
		
	}
	
//Implement boundary conditions in Ghost Zones-ooooo	
	//Y direction
	for (int i=2; i<Nx-2; i++)
	{
		rho_partial[i*Ny+0]    = rho_partial[i*Ny+Ny-4]; //bottom row 
		rho_partial[i*Ny+1]    = rho_partial[i*Ny+Ny-3];
		rho_partial[i*Ny+Ny-2] = rho_partial[i*Ny+2];	 //top row
		rho_partial[i*Ny+Ny-1] = rho_partial[i*Ny+3];
	}
	//Y Direction

	for (int j=2; j<Ny-2; j++)
	{
		rho_partial[0*Ny+j]      = rho_partial[(Ny-4)*Ny+j]; //left column 
		rho_partial[1*Ny+j]      = rho_partial[(Ny-3)*Ny+j];
		rho_partial[(Nx-2)*Ny+j] = rho_partial[2*Ny+j];	//right column
		rho_partial[(Nx-1)*Ny+j] = rho_partial[3*Ny+j];
	}
	//Calculate y-fluxes everywhere
	for (int i=2; i<Nx-2; i++)
	{
		for (int j=2; j<Ny-2; j++)
		{	
		//Calculate interpolated density
			if (Vy[i*Ny+j] >= 0.)
			{
				rho_interp=rho_partial[i*Ny+j-1];
			}
			else
			{	
				rho_interp=rho_partial[i*Ny+j];
			}
			//Calculate F_1
			F[i*Ny+j] = rho_interp * Vy[i*Ny+j];
		
		} 		
	}//end calculating fluxes
	ofstream myfile ("test1.txt");
    for (int j=Nx-3; j>1; j--)
    {
		for (int i=2; i<Ny-2; i++)
		{	myfile <<setprecision(2)<<F[i*Ny+j]<<" ";
		} 
		myfile <<"\n";
    }
    myfile.close();

	//Update density
	for (int i=2; i<Nx-2; i++)
	{
		for (int j=2; j<Ny-3; j++)
		{
			rho_new[i*Ny+j] = rho_partial[i*Ny+j] - dt/dl * (F[i*Ny+j+1] - F[i*Ny+j]);
		}
	}// Boundary Fluxes
	for (int i=2; i<Nx-2; i++)
	{int j=Ny-3;
			rho_new[i*Ny+j] = rho_partial[i*Ny+j] - dt/dl * (F[i*Ny+2] - F[i*Ny+j]);
		
	}
	
	ofstream myfile2 ("test2.txt");
    for (int j=Nx-3; j>1; j--)
    {
		for (int i=2; i<Ny-2; i++)
		{	myfile2 <<setprecision(2)<<rho_new[i*Ny+j]<<" ";
		} 
		myfile2 <<"\n";
    }
    myfile2.close();
	delete [] F;
	//Implement boundary conditions in Ghost Zones-ooooo	
	//Y direction
	for (int i=2; i<Nx-2; i++)
	{
		rho_new[i*Ny+0]    = rho_new[i*Ny+Ny-4]; //bottom row 
		rho_new[i*Ny+1]    = rho_new[i*Ny+Ny-3];
		rho_new[i*Ny+Ny-2] = rho_new[i*Ny+2];	 //top row
		rho_new[i*Ny+Ny-1] = rho_new[i*Ny+3];
	}
	//Y Direction

	for (int j=2; j<Ny-2; j++)
	{
		rho_new[0*Ny+j]      = rho_new[(Ny-4)*Ny+j]; //left column 
		rho_new[1*Ny+j]      = rho_new[(Ny-3)*Ny+j];
		rho_new[(Nx-2)*Ny+j] = rho_new[2*Ny+j];	//right column
		rho_new[(Nx-1)*Ny+j] = rho_new[3*Ny+j];
	}
}
