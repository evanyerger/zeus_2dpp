#include<iostream>
#include <stdlib.h>
#include<string>
#include<math.h>
#include <fstream>
#include <iomanip>

using namespace std;
void advectDensity(double *rho_old,double *rho_new, double *Vx, double *Vy, int Nx, int Ny, double dt, double dl);
void updatePeriodicBC(double *arr,int Nx, int Ny);
double vanLeer(double *arr,double V,int i,int j,int Nx, int Ny,double dt,double dl, int direction);
double vanLeer1(double V,double ip1j,double im1j,double im2j,double ij,double ijm1,double ijm2,double ijp1,int Nx, int Ny,double dt,double dl,int direction);
int main(int argc, char *argv[])
{
	int nx = pow(2,9), ny = pow(2,9);
	int Nx = nx+4,     Ny = ny+4;
	double *rho_old = new double[(Nx)*(Ny)];
	double *rho_new = new double[(Nx)*(Ny)];
	double *Vx = new double[(Nx)*(Ny)];
	double *Vy = new double[(Nx)*(Ny)];
	double *temp;
	int steps=0;
	double dl    = 1./Nx;
	double dt    = dl;
	//cout<<dl<<endl;
	//cout<<dt<<endl;
	if (argc==2)
		sscanf(argv[1],"%d",&steps);
	
	for (int i=2; i<Nx-2; i++)
    {
		for (int j=2; j<Ny-2; j++)	
		{	
		/*
			if (j<Nx/2.  and i<Nx/2.)
				rho_old[i*Ny+j] = i/10.;
			else
				rho_old[i*Ny+j] = 0.;
		*/rho_old[i*Ny+j] = 10*exp(-((i-Nx/2.)*(i-Nx/2.)+(j-Nx/2.)*(j-Nx/2.))/400);
			Vx[i*Ny+j]  = 1.;
			Vy[i*Ny+j]  = 1.;
		} 
		
    }
	updatePeriodicBC(rho_old,Nx,Ny);
	updatePeriodicBC(Vy,Nx,Ny);
	updatePeriodicBC(Vx,Nx,Ny);
    for (int i=0; i<steps; i++)
    {
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
        ofstream myfile3 ("output.txt");
    for (int i=0; i<Nx; i++)
    {
		for (int j=0; j<Ny; j++)
		{	myfile3 <<setprecision(2)<<rho_old[i*Ny+j]<<" ";
		} 
		myfile3 <<"\n";
    }
    myfile3.close();
}

void advectDensity(double *rho_old,double *rho_new, double *Vx, double *Vy, int Nx, int Ny, double dt, double dl)
{	int view=0;
	double *rho_partial 	= new double[Nx*Ny];
	double *F          		= new double[(Nx)*(Ny)];
	double rho_interp   	= 0.;
	//Calculate x-fluxes everywhere
	for (int i=2; i<Nx-2; i++)
	{
		for (int j=2; j<Ny-2; j++)
		{	
		//Calculate interpolated density
		
			//rho_interp = vanLeer(rho_old,Vx[i*Ny+j],i,j,Nx,Ny,dt,dl,0);
			rho_interp=vanLeer1(Vx[i*Ny+j],rho_old[(i+1)*Ny+j] ,rho_old[(i-1)*Ny+j] ,rho_old[(i-2)*Ny+j] ,rho_old[i*Ny+j] ,rho_old[i*Ny+j-1] ,rho_old[i*Ny+j-2] ,rho_old[i*Ny+j+1] ,Nx,Ny,dt,dl,0);

			//Calculate F_1
			F[i*Ny+j]  = rho_interp * Vx[i*Ny+j] ;
		
		} 		
	}//end calculating fluxes	  
	updatePeriodicBC(F,Nx,Ny);
	if (view==1){
		ofstream myfile4 ("Fx.txt");
    for (int j=Nx-1; j>-1; j--)
    {
		for (int i=0; i<Ny; i++)
		{	myfile4 <<setprecision(2)<<F[i*Ny+j]<<" ";
		} 
		myfile4 <<"\n";
    }
    myfile4.close();
    }
	//Partial update density back into rho_partial
	for (int i=2; i<Nx-2; i++)
	{
		for (int j=2; j<Ny-2; j++)
		{
			rho_partial[i*Ny+j] = rho_old[i*Ny+j] - dt/dl * (F[(i+1)*Ny+j] - F[i*Ny+j]);
		}
	}

	updatePeriodicBC(rho_partial,Nx,Ny);
	if (view==1){
	ofstream myfile3 ("rho_partial.txt");
    for (int j=Nx-1; j>-1; j--)
    {
		for (int i=0; i<Ny; i++)
		{	myfile3 <<setprecision(2)<<rho_partial[i*Ny+j]<<" ";
		} 
		myfile3 <<"\n";
    }
    myfile3.close();
    }
	//Calculate y-fluxes everywhere
	for (int i=2; i<Nx-2; i++)
	{
		for (int j=2; j<Ny-2; j++)
		{	
		//Calculate interpolated density
			
			//rho_interp=vanLeer(rho_partial,Vy[i*Ny+j],i,j,Nx,Ny,dt,dl,1);
			rho_interp=vanLeer1(Vy[i*Ny+j],rho_partial[(i+1)*Ny+j] ,rho_partial[(i-1)*Ny+j] ,rho_partial[(i-2)*Ny+j] ,rho_partial[i*Ny+j] ,rho_partial[i*Ny+j-1] ,rho_partial[i*Ny+j-2] ,rho_partial[i*Ny+j+1],Nx,Ny,dt,dl,1);
			//Calculate F_1
			F[i*Ny+j] = rho_interp * Vy[i*Ny+j];
		
		} 		
	}//end calculating fluxes
	updatePeriodicBC(F,Nx,Ny);
	if (view==1){
	  ofstream myfile1 ("Fy.txt");
    for (int j=Nx-1; j>-1; j--)
    {
		for (int i=0; i<Ny; i++)
		{	myfile1 <<setprecision(2)<<F[i*Ny+j]<<" ";
		} 
		myfile1 <<"\n";
    }
    myfile1.close();
	}
	//Update density
	for (int i=2; i<Nx-2; i++)
	{
		for (int j=2; j<Ny-2; j++)
		{
			rho_new[i*Ny+j] = rho_partial[i*Ny+j] - dt/dl * (F[i*Ny+j+1] - F[i*Ny+j]);
		}
	}
	if (view==1){
	ofstream myfile2 ("rho_new.txt");
    for (int j=Nx-1; j>-1; j--)
    {
		for (int i=0; i<Ny; i++)
		{	myfile2 <<setprecision(2)<<rho_new[i*Ny+j]<<" ";
		} 
		myfile2 <<"\n";
    }
    myfile2.close();
    }
	delete [] F;
	updatePeriodicBC(rho_new,Nx,Ny);
}

void updatePeriodicBC(double *arr,int Nx, int Ny)
{//Implement boundary conditions in Ghost Zones-ooooo	
	//Y direction
	for (int i=2; i<Nx-2; i++)
	{
		arr[i*Ny+0]    = arr[i*Ny+Ny-4]; //bottom row 
		arr[i*Ny+1]    = arr[i*Ny+Ny-3];
		arr[i*Ny+Ny-2] = arr[i*Ny+2];	 //top row
		arr[i*Ny+Ny-1] = arr[i*Ny+3];
	}
	//X Direction

	for (int j=2; j<Ny-2; j++)
	{
		arr[0*Ny+j]      = arr[(Ny-4)*Ny+j]; //left column 
		arr[1*Ny+j]      = arr[(Ny-3)*Ny+j];
		arr[(Nx-2)*Ny+j] = arr[2*Ny+j];	//right column
		arr[(Nx-1)*Ny+j] = arr[3*Ny+j];
	}
}
double vanLeer1(double V,double ip1j,double im1j,double im2j,double ij,double ijm1,double ijm2,double ijp1,int Nx, int Ny,double dt,double dl,int direction)
{
	double rho_interp   	= 0.;
	double d_rho			= 0.;
	double delta_rho_plus  	= 0.;
	double delta_rho_minus 	= 0.;
		//Calculate interpolated density
			if (V > 0.)
			{
				//Calculate delta_rho's
				if (direction==0)
				{	
					delta_rho_plus  = ij     - im1j;
					delta_rho_minus = im1j   - im2j;
				}
				else
				{
					delta_rho_plus  = ij  - ijm1;
					delta_rho_minus = ijm1- ijm2;
				}
				if (delta_rho_plus * delta_rho_minus > 0)
					d_rho = 2. * delta_rho_plus * delta_rho_minus / (delta_rho_plus + delta_rho_minus);
				else
					d_rho = 0.;
				
				//Calculate rho_interp
				if (direction==0)
					rho_interp = im1j+ (1.-V*dt/dl)*d_rho/2.;
				else
					rho_interp = ijm1   + (1.-V*dt/dl)*d_rho/2.;

			}
			else
			{
				//Calculate delta_rho's
				if (direction==0)
				{
					delta_rho_plus  = ip1j   - ij;
					delta_rho_minus = ij    - im1j;
				}
				else
				{
					delta_rho_plus  = ijp1   - ij;
					delta_rho_minus = ij     - ijm1;
				}
				if (delta_rho_plus * delta_rho_minus > 0)
					d_rho = 2. * delta_rho_plus * delta_rho_minus / (delta_rho_plus + delta_rho_minus);
				else
					d_rho = 0.;
				
				//Calculate rho_interp
				rho_interp = ij - (1.+V*dt/dl)*d_rho/2.;
			}
	return rho_interp;
}

double vanLeer(double *arr,double V,int i,int j,int Nx, int Ny,double dt,double dl,int direction)
{
	double rho_interp   	= 0.;
	double d_rho			= 0.;
	double delta_rho_plus  	= 0.;
	double delta_rho_minus 	= 0.;
		//Calculate interpolated density
			if (V > 0.)
			{
				//Calculate delta_rho's
				if (direction==0)
				{	
					delta_rho_plus  = arr[i*Ny+j]     - arr[(i-1)*Ny+j];
					delta_rho_minus = arr[(i-1)*Ny+j] - arr[(i-2)*Ny+j];
				}
				else
				{
					delta_rho_plus  = arr[i*Ny+j]   - arr[i*Ny+j-1];
					delta_rho_minus = arr[i*Ny+j-1] - arr[i*Ny+j-2];
				}
				if (delta_rho_plus * delta_rho_minus > 0)
					d_rho = 2. * delta_rho_plus * delta_rho_minus / (delta_rho_plus + delta_rho_minus);
				else
					d_rho = 0.;
				
				//Calculate rho_interp
				if (direction==0)
					rho_interp = arr[(i-1)*Ny+j] + (1.-V*dt/dl)*d_rho/2.;
				else
					rho_interp = arr[i*Ny+j-1]   + (1.-V*dt/dl)*d_rho/2.;

			}
			else
			{
				//Calculate delta_rho's
				if (direction==0)
				{
					delta_rho_plus  = arr[(i+1)*Ny+j]   - arr[i*Ny+j];
					delta_rho_minus = arr[i*Ny+j]       - arr[(i-1)*Ny+j];
				}
				else
				{
					delta_rho_plus  = arr[i*Ny+j+1]   - arr[i*Ny+j];
					delta_rho_minus = arr[i*Ny+j]     - arr[i*Ny+j-1];
				}
				if (delta_rho_plus * delta_rho_minus > 0)
					d_rho = 2. * delta_rho_plus * delta_rho_minus / (delta_rho_plus + delta_rho_minus);
				else
					d_rho = 0.;
				
				//Calculate rho_interp
				rho_interp = arr[i*Ny+j] - (1.+V*dt/dl)*d_rho/2.;
			}
	return rho_interp;
}



