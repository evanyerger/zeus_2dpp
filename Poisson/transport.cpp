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
int main(int argc, char *argv[])
{
	int nx = pow(2,10), ny = pow(2,10);
	int Nx = nx+4,     Ny = ny+4;
	double *rho_old = new double[(Nx)*(Ny)];
	double *rho_new = new double[(Nx)*(Ny)];
	double *Vx = new double[(Nx)*(Ny)];
	double *Vy = new double[(Nx)*(Ny)];
	double *temp;
	int steps=100;
	double dl    = 1./Nx;
	double dt    = dl;
	//cout<<dl<<endl;
	//cout<<dt<<endl;
	for (int i=2; i<Nx-2; i++)
    {
		for (int j=2; j<Ny-2; j++)	
		{	
		/*
			if (j<Nx/2. and i<Nx/2.)
				rho_old[i*Ny+j] = 10.;
			else
				rho_old[i*Ny+j] = 0.;
		*/rho_old[i*Ny+j] = exp(-((i-Nx/2.)*(i-Nx/2.)+(j-Nx/2.)*(j-Nx/2.))/4000);
			Vx[i*Ny+j]  = 1.;
			Vy[i*Ny+j]  = 1.;
		} 
		
    }
	updatePeriodicBC(rho_old,Nx,Ny);
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
{
	double *rho_partial 	= new double[Nx*Ny];
	double *F          		= new double[(Nx)*(Ny)];
	double rho_interp   	= 0.;
	//Calculate x-fluxes everywhere
	for (int i=2; i<Nx-2; i++)
	{
		for (int j=2; j<Ny-2; j++)
		{	/*
		//Calculate interpolated density
			if (Vx[i*Ny+j] > 0.)
			{
				//Calculate delta_rho's
				delta_rho_plus  = rho_old[i*Ny+j]   - rho_old[(i-1)*Ny+j];
				delta_rho_minus = rho_old[(i-1)*Ny+j] - rho_old[(i-2)*Ny+j];
				if (delta_rho_plus * delta_rho_minus > 0)
					d_rho = 2. * delta_rho_plus * delta_rho_minus / (delta_rho_plus + delta_rho_minus);
				else
					d_rho = 0.;
				
				//Calculate d_rho
				rho_interp = rho_old[(i-1)*Ny+j] + (1.-Vx[i*Ny+j]*dt/dl)*d_rho/2.;

			}
			else
			{
				//Calculate delta_rho's
				delta_rho_plus  = rho_old[(i+1)*Ny+j]   - rho_old[i*Ny+j];
				delta_rho_minus = rho_old[i*Ny+j] - rho_old[(i-1)*Ny+j];
				if (delta_rho_plus * delta_rho_minus > 0)
					d_rho = 2. * delta_rho_plus * delta_rho_minus / (delta_rho_plus + delta_rho_minus);
				else
					d_rho = 0.;
				
				//Calculate d_rho
				rho_interp = rho_old[i*Ny+j] - (1.+Vx[i*Ny+j]*dt/dl)*d_rho/2.;
			}*/
			rho_interp = vanLeer(rho_old,Vx[i*Ny+j],i,j,Nx,Ny,dt,dl,0);
			//Calculate F_1
			F[i*Ny+j]  = rho_interp * Vx[i*Ny+j] ;
		
		} 		
	}//end calculating fluxes
	updatePeriodicBC(F,Nx,Ny);
	//Partial update density back into rho_partial
	for (int i=2; i<Nx-2; i++)
	{
		for (int j=2; j<Ny-2; j++)
		{
			rho_partial[i*Ny+j] = rho_old[i*Ny+j] - dt/dl * (F[(i+1)*Ny+j] - F[i*Ny+j]);
		}
	}/*
	// Boundary Fluxes
	for (int j=2; j<Nx-2; j++)
	{int i=Ny-3;
			rho_partial[i*Ny+j] = rho_old[i*Ny+j] - dt/dl * (F[2*Ny+j] - F[i*Ny+j]);
		
	}
	*/
	updatePeriodicBC(rho_partial,Nx,Ny);
	//Calculate y-fluxes everywhere
	for (int i=2; i<Nx-2; i++)
	{
		for (int j=2; j<Ny-2; j++)
		{	/*
		//Calculate interpolated density
			if (Vy[i*Ny+j] >= 0.)
			{
				//Calculate delta_rho's
				delta_rho_plus  = rho_partial[i*Ny+j]   - rho_partial[i*Ny+j-1];
				delta_rho_minus = rho_partial[i*Ny+j-1] - rho_partial[i*Ny+j-2];
				if (delta_rho_plus * delta_rho_minus > 0)
					d_rho = 2. * delta_rho_plus * delta_rho_minus / (delta_rho_plus + delta_rho_minus);
				else
					d_rho = 0.;
				
				//Calculate d_rho
				rho_interp = rho_partial[i*Ny+j-1] + (1.-Vx[i*Ny+j]*dt/dl)*d_rho/2.;

			}
			else
			{
				//Calculate delta_rho's
				delta_rho_plus  = rho_partial[i*Ny+j+1]   - rho_partial[i*Ny+j];
				delta_rho_minus = rho_partial[i*Ny+j] - rho_partial[i*Ny+j-1];
				if (delta_rho_plus * delta_rho_minus > 0)
					d_rho = 2. * delta_rho_plus * delta_rho_minus / (delta_rho_plus + delta_rho_minus);
				else
					d_rho = 0.;
				
				//Calculate d_rho
				rho_interp = rho_partial[i*Ny+j] - (1.+Vx[i*Ny+j]*dt/dl)*d_rho/2.;
			}*/
			rho_interp=vanLeer(rho_partial,Vy[i*Ny+j],i,j,Nx,Ny,dt,dl,1);
			//Calculate F_1
			F[i*Ny+j] = rho_interp * Vy[i*Ny+j];
		
		} 		
	}//end calculating fluxes
	updatePeriodicBC(F,Nx,Ny);
	//Update density
	for (int i=2; i<Nx-2; i++)
	{
		for (int j=2; j<Ny-3; j++)
		{
			rho_new[i*Ny+j] = rho_partial[i*Ny+j] - dt/dl * (F[i*Ny+j+1] - F[i*Ny+j]);
		}
	}
	/*// Boundary Fluxes
	for (int i=2; i<Nx-2; i++)
	{int j=Ny-3;
			rho_new[i*Ny+j] = rho_partial[i*Ny+j] - dt/dl * (F[i*Ny+2] - F[i*Ny+j]);
		
	}
	*/
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

void advectEnergy(double *rho_old,double *rho_new, double *Vx, double *Vy, int Nx, int Ny, double dt, double dl)
{
	double *rho_partial 	= new double[Nx*Ny];
	double *F          		= new double[(Nx)*(Ny)];
	double rho_interp   	= 0.;
	//Calculate x-fluxes everywhere
	for (int i=2; i<Nx-2; i++)
	{
		for (int j=2; j<Ny-2; j++)
		{	/*
		//Calculate interpolated density
			if (Vx[i*Ny+j] > 0.)
			{
				//Calculate delta_rho's
				delta_rho_plus  = rho_old[i*Ny+j]   - rho_old[(i-1)*Ny+j];
				delta_rho_minus = rho_old[(i-1)*Ny+j] - rho_old[(i-2)*Ny+j];
				if (delta_rho_plus * delta_rho_minus > 0)
					d_rho = 2. * delta_rho_plus * delta_rho_minus / (delta_rho_plus + delta_rho_minus);
				else
					d_rho = 0.;
				
				//Calculate d_rho
				rho_interp = rho_old[(i-1)*Ny+j] + (1.-Vx[i*Ny+j]*dt/dl)*d_rho/2.;

			}
			else
			{
				//Calculate delta_rho's
				delta_rho_plus  = rho_old[(i+1)*Ny+j]   - rho_old[i*Ny+j];
				delta_rho_minus = rho_old[i*Ny+j] - rho_old[(i-1)*Ny+j];
				if (delta_rho_plus * delta_rho_minus > 0)
					d_rho = 2. * delta_rho_plus * delta_rho_minus / (delta_rho_plus + delta_rho_minus);
				else
					d_rho = 0.;
				
				//Calculate d_rho
				rho_interp = rho_old[i*Ny+j] - (1.+Vx[i*Ny+j]*dt/dl)*d_rho/2.;
			}*/
			rho_interp = vanLeer(rho_old,Vx[i*Ny+j],i,j,Nx,Ny,dt,dl,0);
			//Calculate F_1
			F[i*Ny+j]  = rho_interp * Vx[i*Ny+j] ;
		
		} 		
	}//end calculating fluxes
	updatePeriodicBC(F,Nx,Ny);
	//Partial update density back into rho_partial
	for (int i=2; i<Nx-2; i++)
	{
		for (int j=2; j<Ny-2; j++)
		{
			rho_partial[i*Ny+j] = rho_old[i*Ny+j] - dt/dl * (F[(i+1)*Ny+j] - F[i*Ny+j]);
		}
	}
	// Boundary Fluxes
	//for (int j=2; j<Nx-2; j++)
	//{int i=Ny-3;
	//		rho_partial[i*Ny+j] = rho_old[i*Ny+j] - dt/dl * (F[2*Ny+j] - F[i*Ny+j]);
		
	//}
	
	updatePeriodicBC(rho_partial,Nx,Ny);
	//Calculate y-fluxes everywhere
	for (int i=2; i<Nx-2; i++)
	{
		for (int j=2; j<Ny-2; j++)
		{	/*
		//Calculate interpolated density
			if (Vy[i*Ny+j] >= 0.)
			{
				//Calculate delta_rho's
				delta_rho_plus  = rho_partial[i*Ny+j]   - rho_partial[i*Ny+j-1];
				delta_rho_minus = rho_partial[i*Ny+j-1] - rho_partial[i*Ny+j-2];
				if (delta_rho_plus * delta_rho_minus > 0)
					d_rho = 2. * delta_rho_plus * delta_rho_minus / (delta_rho_plus + delta_rho_minus);
				else
					d_rho = 0.;
				
				//Calculate d_rho
				rho_interp = rho_partial[i*Ny+j-1] + (1.-Vx[i*Ny+j]*dt/dl)*d_rho/2.;

			}
			else
			{
				//Calculate delta_rho's
				delta_rho_plus  = rho_partial[i*Ny+j+1]   - rho_partial[i*Ny+j];
				delta_rho_minus = rho_partial[i*Ny+j] - rho_partial[i*Ny+j-1];
				if (delta_rho_plus * delta_rho_minus > 0)
					d_rho = 2. * delta_rho_plus * delta_rho_minus / (delta_rho_plus + delta_rho_minus);
				else
					d_rho = 0.;
				
				//Calculate d_rho
				rho_interp = rho_partial[i*Ny+j] - (1.+Vx[i*Ny+j]*dt/dl)*d_rho/2.;
			}*/
			rho_interp=vanLeer(rho_partial,Vy[i*Ny+j],i,j,Nx,Ny,dt,dl,1);
			//Calculate F_1
			F[i*Ny+j] = rho_interp * Vy[i*Ny+j];
		
		} 		
	}//end calculating fluxes
	updatePeriodicBC(F,Nx,Ny);
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
		for (int j=2; j<Ny-2; j++)
		{
			rho_new[i*Ny+j] = rho_partial[i*Ny+j] - dt/dl * (F[i*Ny+j+1] - F[i*Ny+j]);
		}
	}// Boundary Fluxes
	//for (int i=2; i<Nx-2; i++)
	//{int j=Ny-3;
	//		rho_new[i*Ny+j] = rho_partial[i*Ny+j] - dt/dl * (F[i*Ny+2] - F[i*Ny+j]);
		
	//}
	
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
	updatePeriodicBC(rho_new,Nx,Ny);
}

