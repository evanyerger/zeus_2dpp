#include <iostream>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <fstream>
#include <iomanip>

using namespace std;

void advectEnergy(double *e_old, double *e_new, double *rho_old,double *rho_new, double *Vx, double *Vx_new, double *Vy,double *Vy_new,double *Vz,double *Vz_new, int Nx, int Ny, double dt, double dl);

void updatePeriodicBC(double *arr,int Nx, int Ny);

double vanLeerPBV(double V,double ip1j,double im1j,double im2j,double ij,double ijm1,double ijm2,double ijp1,int Nx, int Ny,double dt,double dl,int direction); // PBV = pass by value

int main(int argc, char *argv[])
{
	int Nx = pow(2,6)+4,    Ny = pow(2,6)+4;
	double *rho_old = new double[(Nx)*(Ny)];
	double *rho_new = new double[(Nx)*(Ny)];
	double *e_old   = new double[(Nx)*(Ny)];
	double *e_new   = new double[(Nx)*(Ny)];
	double *Vx_old  = new double[(Nx)*(Ny)];
	double *Vy_old  = new double[(Nx)*(Ny)];
	double *Vx_new  = new double[(Nx)*(Ny)];
	double *Vy_new  = new double[(Nx)*(Ny)];
	double *Vz_old  = new double[(Nx)*(Ny)];
	double *Vz_new  = new double[(Nx)*(Ny)];
	double *temp;
	int steps       = 0;
	double dl       = 1./Nx;
	double dt       = dl/3.;

	if (argc==2)
		sscanf(argv[1],"%d",&steps);

  // setup
	for (int i=2; i<Nx-2; i++)
    {
		for (int j=2; j<Ny-2; j++)
		{

			if (j<Nx/2.  and i<Nx/2.)
			{	rho_old[i*Ny+j] = 1;
				e_old[i*Ny+j] = i;
				Vz_old[i*Ny+j]  = i;
			}
			else{
				rho_old[i*Ny+j] = 1.;
				Vz_old[i*Ny+j]  = 1;
				}
		Vz_old[i*Ny+j] = 10.*exp(-((i-Nx/2.)*(i-Nx/2.)+(j-Nx/2.)*(j-Nx/2.))/400)+1.;
		e_old[i*Ny+j]   = 10*exp(-((i-Nx/2.)*(i-Nx/2.)+(j-Nx/2.)*(j-Nx/2.))/400);
			Vx_old[i*Ny+j]  = .5;
			Vy_old[i*Ny+j]  = 1.5;

		}

    }

	updatePeriodicBC(rho_old,Nx,Ny);
	updatePeriodicBC(e_old,Nx,Ny);
	updatePeriodicBC(Vx_old,Nx,Ny);
	updatePeriodicBC(Vy_old,Nx,Ny);
	updatePeriodicBC(Vz_old,Nx,Ny);
  // end setup

    for (int i=0; i<steps; i++)
    {

    	advectEnergy(e_old, e_new, rho_old, rho_new, Vx_old, Vx_new, Vy_old, Vy_new, Vz_old, Vz_new, Nx, Ny, dt, dl);
    	temp    = rho_old;
    	rho_old = rho_new;
    	rho_new=temp;
    	temp    = e_old;
    	e_old = e_new;
    	e_new=temp;
    	/*temp    = Vx_old;
    	Vx_old = Vx_new;
    	Vx_new=temp;
    	temp    = Vy_old;
    	Vy_old = Vy_new;
    	Vy_new=temp;*/
    	temp    = Vz_old;
    	Vz_old = Vz_new;
    	Vz_new=temp;

    }
    ofstream myfile3 ("output.txt");
    for (int i=0; i<Nx; i++)
    {
		for (int j=0; j<Ny; j++)
		{	myfile3 <<setprecision(2)<<Vz_old[i*Ny+j]<<" ";
		}
		myfile3 <<"\n";
    }
    myfile3.close();

	ofstream myfile5 ("new.txt");
	myfile5<<"Vz_new:" <<"\n";
    for (int j=Nx-1; j>-1; j--)
    {
		for (int i=0; i<Ny; i++)
		{
			myfile5 <<setprecision(2)<<Vz_old[i*Ny+j]<<" ";
		}
		myfile5 <<"\n";
    }myfile5 <<"rho:"<<"\n";
        for (int j=Nx-1; j>-1; j--)
    {
		for (int i=0; i<Ny; i++)
		{	myfile3 <<setprecision(2)<<rho_old[i*Ny+j]<<" ";
		}
		myfile5 <<"\n";
    }
    myfile5.close();

}

void advectEnergy(double *e_old, double *e_new, double *rho_old,double *rho_new, double *Vx,double *Vx_new, double *Vy,double *Vy_new,double *Vz,double *Vz_new, int Nx, int Ny, double dt, double dl)
{	int view=1;
	double *e_partial   	= new double[Nx*Ny];
	double *rho_partial 	= new double[Nx*Ny];
	double *Sx_partial	 	= new double[Nx*Ny];
	double *Sy_partial	 	= new double[Nx*Ny];
	double *Sz_partial	 	= new double[Nx*Ny];
	double *F_e       		= new double[Nx*Ny];
	double *F_rho      		= new double[Nx*Ny];
	double *F_Vx      		= new double[Nx*Ny];
	double *F_Vy      		= new double[Nx*Ny];
	double *F_Vz      		= new double[Nx*Ny];
	double ed_interp   		= 0.;
	double rho_interp   	= 0.;
	double Vx_interp   	    = 0.;
	double Vy_interp   	    = 0.;
	double Vz_interp   	    = 0.;
	double rhoXm2_interp   	= 0.;
	double rhoYm2_interp   	= 0.;

	//Calculate X-fluxes everywhere
	for (int i=2; i<Nx-2; i++)
	{
		for (int j=2; j<Ny-2; j++)
		{
		  	// Calculate interpolated energy density
			ed_interp = vanLeerPBV(Vx[i*Ny+j],e_old[(i+1)*Ny+j]/rho_old[(i+1)*Ny+j] ,e_old[(i-1)*Ny+j]/rho_old[(i-1)*Ny+j],e_old[(i-2)*Ny+j]/rho_old[(i-2)*Ny+j] ,e_old[(i)*Ny+j]/rho_old[(i)*Ny+j],e_old[(i)*Ny+j-1]/rho_old[(i)*Ny+j-1] ,e_old[(i)*Ny+j-2]/rho_old[(i)*Ny+j-2] ,e_old[(i)*Ny+j+1]/rho_old[(i)*Ny+j+1] ,Nx,Ny,dt,dl,0);
			// Calculate interpolated density
			rho_interp = vanLeerPBV(Vx[i*Ny+j],rho_old[(i+1)*Ny+j] ,rho_old[(i-1)*Ny+j] ,rho_old[(i-2)*Ny+j] ,rho_old[i*Ny+j] ,rho_old[i*Ny+j-1] ,rho_old[i*Ny+j-2] ,rho_old[i*Ny+j+1] ,Nx,Ny,dt,dl,0);

			// Calculate X velocity mass fluxes
			Vx_interp = vanLeerPBV(Vx[i*Ny+j],Vx[(i+1)*Ny+j] ,Vx[(i-1)*Ny+j] ,Vx[(i-2)*Ny+j] ,Vx[i*Ny+j] ,Vx[i*Ny+j-1] ,Vx[i*Ny+j-2] ,Vx[i*Ny+j+1] ,Nx,Ny,dt,dl,0);

			rhoXm2_interp = vanLeerPBV(Vx[(i+1)*Ny+j],rho_old[(i+2)*Ny+j] ,rho_old[(i)*Ny+j] ,rho_old[(i-1)*Ny+j] ,rho_old[(i+1)*Ny+j] ,rho_old[(i+1)*Ny+j-1] ,rho_old[(i+1)*Ny+j-2] ,rho_old[(i+1)*Ny+j+1] ,Nx,Ny,dt,dl,0);

			// Calculate Y velocity mass fluxes
			Vy_interp = vanLeerPBV(Vx[i*Ny+j],Vy[(i+1)*Ny+j] ,Vy[(i-1)*Ny+j] ,Vy[(i-2)*Ny+j] ,Vy[i*Ny+j] ,Vy[i*Ny+j-1] ,Vy[i*Ny+j-2] ,Vy[i*Ny+j+1] ,Nx,Ny,dt,dl,0);
			//Fix j-3 term
			rhoYm2_interp = vanLeerPBV(Vx[i*Ny+j-1],rho_old[(i+1)*Ny+j-1] ,rho_old[(i-1)*Ny+j-1] ,rho_old[(i-2)*Ny+j-1] ,rho_old[i*Ny+j-1] ,rho_old[i*Ny+j-2] ,rho_old[i*Ny+j-3] ,rho_old[i*Ny+j] ,Nx,Ny,dt,dl,0);

			// Calculate Z velocity mass fluxes
			Vz_interp = vanLeerPBV(Vx[i*Ny+j],Vz[(i+1)*Ny+j] ,Vz[(i-1)*Ny+j] ,Vz[(i-2)*Ny+j] ,Vz[i*Ny+j] ,Vz[i*Ny+j-1] ,Vz[i*Ny+j-2] ,Vz[i*Ny+j+1] ,Nx,Ny,dt,dl,0);

			//Calculate Fluxes
			F_e[i*Ny+j]   = ed_interp *rho_interp * Vx[i*Ny+j];
			F_rho[i*Ny+j] = rho_interp            * Vx[i*Ny+j];
			F_Vx[i*Ny+j]  = Vx_interp*.5*(rho_interp*Vx[i*Ny+j] + rhoXm2_interp*Vx[(i+1)*Ny+j]);
			F_Vy[i*Ny+j]  = Vy_interp*.5*(rho_interp*Vx[i*Ny+j] + rhoYm2_interp*Vx[i*Ny+j-1]);
			F_Vz[i*Ny+j]  = Vz_interp*(rho_interp*Vx[i*Ny+j]);

		}
	} // End calculating fluxes

	// Enforce Periodic boundary conditions
	updatePeriodicBC(F_e,Nx,Ny);
	updatePeriodicBC(F_rho,Nx,Ny);
	updatePeriodicBC(F_Vx,Nx,Ny);
	updatePeriodicBC(F_Vy,Nx,Ny);
	updatePeriodicBC(F_Vz,Nx,Ny);

    if (view==1){
	ofstream myfile4 ("Fx.txt");
	myfile4<<"F_Vz:" <<"\n";
    for (int j=Nx-1; j>-1; j--)
    {
		for (int i=0; i<Ny; i++)
		{	myfile4 <<setprecision(2)<<F_Vz[i*Ny+j]<<" ";
		}
		myfile4 <<"\n";
    }
    myfile4<<"old Vz:" <<"\n";

    for (int j=Nx-1; j>-1; j--)
    {
		for (int i=0; i<Ny; i++)
		{	myfile4 <<setprecision(2)<<Vz[i*Ny+j]<<" ";
		}
		myfile4 <<"\n";
    }
    myfile4.close();
    }
	//Partial update
	for (int i=2; i<Nx-2; i++)
	{
		for (int j=2; j<Ny-2; j++)
		{
			e_partial[i*Ny+j]   = e_old[i*Ny+j] - dt/dl * (F_e[(i+1)*Ny+j] - F_e[i*Ny+j]);
			rho_partial[i*Ny+j] = rho_old[i*Ny+j] - dt/dl * (F_rho[(i+1)*Ny+j] - F_rho[i*Ny+j]);
			Sx_partial[i*Ny+j]  = ((rho_old[(i-1)*Ny+j]+rho_old[i*Ny+j])/2.*Vx[i*Ny+j] - dt/dl * (F_Vx[(i)*Ny+j] - F_Vx[(i-1)*Ny+j]));
			Sy_partial[i*Ny+j]  = ((rho_old[(i)*Ny+j-1]+rho_old[i*Ny+j])/2.*Vy[i*Ny+j] - dt/dl * (F_Vy[(i+1)*Ny+j] - F_Vy[(i)*Ny+j]));
			Sz_partial[i*Ny+j]  = (rho_old[(i)*Ny+j]*Vz[i*Ny+j] - dt/dl * (F_Vz[(i+1)*Ny+j] - F_Vz[(i)*Ny+j]));
		}
	}
	updatePeriodicBC(e_partial,Nx,Ny);
	updatePeriodicBC(rho_partial,Nx,Ny);
	updatePeriodicBC(Sx_partial,Nx,Ny);
	updatePeriodicBC(Sy_partial,Nx,Ny);
	updatePeriodicBC(Sz_partial,Nx,Ny);

		if (view==1){
	ofstream myfile3 ("partial.txt");
	myfile3<<"Sz_partial:" <<"\n";
    for (int j=Nx-1; j>-1; j--)
    {
		for (int i=0; i<Ny; i++)
		{
			myfile3 <<setprecision(2)<<Sz_partial[i*Ny+j]<<" ";
		}
		myfile3 <<"\n";
    }myfile3 <<"Velocity_z:"<<"\n";
        for (int j=Nx-1; j>-1; j--)
    {
		for (int i=0; i<Ny; i++)
		{	if (i==0)
			myfile3 <<setprecision(2)<<Sz_partial[i*Ny+j]/(rho_partial[i*Ny+j])<<" ";
			else
			myfile3 <<setprecision(2)<<Sz_partial[i*Ny+j]/(rho_partial[i*Ny+j])<<" ";
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
			//Calculate interpolated energy density

			ed_interp=vanLeerPBV(Vy[i*Ny+j],e_partial[(i+1)*Ny+j]/rho_partial[(i+1)*Ny+j] ,e_partial[(i-1)*Ny+j]/rho_partial[(i-1)*Ny+j] ,e_partial[(i-2)*Ny+j]/rho_partial[(i-2)*Ny+j] ,e_partial[(i)*Ny+j]/rho_partial[(i)*Ny+j],e_partial[(i)*Ny+j-1]/rho_partial[(i)*Ny+j-1] ,e_partial[(i)*Ny+j-2]/rho_partial[(i)*Ny+j-2] ,e_partial[(i)*Ny+j+1]/rho_partial[(i)*Ny+j+1] ,Nx,Ny,dt,dl,1);

			//Calculate interpolated density
			rho_interp=vanLeerPBV(Vy[i*Ny+j],rho_partial[(i+1)*Ny+j] ,rho_partial[(i-1)*Ny+j] ,rho_partial[(i-2)*Ny+j] ,rho_partial[i*Ny+j] ,rho_partial[i*Ny+j-1] ,rho_partial[i*Ny+j-2] ,rho_partial[i*Ny+j+1],Nx,Ny,dt,dl,1);

			// Calculate X velocity mass fluxes
			Vx_interp = vanLeerPBV(Vy[i*Ny+j],Vx[(i+1)*Ny+j] ,Vx[(i-1)*Ny+j] ,Vx[(i-2)*Ny+j] ,Vx[i*Ny+j] ,Vx[i*Ny+j-1] ,Vx[i*Ny+j-2] ,Vx[i*Ny+j+1] ,Nx,Ny,dt,dl,1);
			//if (i==2) need to correct for the i-3 term
			//Vx_interp = vanLeerPBV(Vy[i*Ny+j],Sx_partial[(i+1)*Ny+j]/(rho_partial[(i+1)*Ny+j]+rho_partial[(i)*Ny+j])*2. ,Sx_partial[(i-1)*Ny+j]/(rho_partial[(i-1)*Ny+j]+rho_partial[(i-2)*Ny+j])*2  ,Sx_partial[(i-2)*Ny+j]/(rho_partial[(i-2)*Ny+j]+rho_partial[(i-3)*Ny+j])*2 ,Sx_partial[i*Ny+j]/(rho_partial[i*Ny+j]+rho_partial[(i-1)*Ny+j])*2 ,Sx_partial[i*Ny+j-1]/(rho_partial[i*Ny+j-1]+rho_partial[(i-1)*Ny+j-1])*2 ,Sx_partial[i*Ny+j-2]/(rho_partial[i*Ny+j-2]+rho_partial[(i-1)*Ny+j-2])*2 ,Sx_partial[i*Ny+j+1]/(rho_partial[i*Ny+j+1]+rho_partial[(i-1)*Ny+j+1])*2 ,Nx,Ny,dt,dl,1);

			rhoXm2_interp = vanLeerPBV(Vy[(i-1)*Ny+j],rho_partial[(i)*Ny+j] ,rho_partial[(i-2)*Ny+j] ,rho_partial[(i-3)*Ny+j],rho_partial[(i-1)*Ny+j] ,rho_partial[(i-1)*Ny+j-1] ,rho_partial[(i-1)*Ny+j-2] ,rho_partial[(i-1)*Ny+j+1] ,Nx,Ny,dt,dl,1);

			// Calculate Y velocity mass fluxes
			Vy_interp = vanLeerPBV(Vy[i*Ny+j],Vy[(i+1)*Ny+j] ,Vy[(i-1)*Ny+j] ,Vy[(i-2)*Ny+j] ,Vy[i*Ny+j] ,Vy[i*Ny+j-1] ,Vy[i*Ny+j-2] ,Vy[i*Ny+j+1] ,Nx,Ny,dt,dl,1);
			//if (i==2) need to correct for the i-3 term
			//Vx_interp = vanLeerPBV(Vy[i*Ny+j],Sx_partial[(i+1)*Ny+j]/(rho_partial[(i+1)*Ny+j]+rho_partial[(i)*Ny+j])*2. ,Sx_partial[(i-1)*Ny+j]/(rho_partial[(i-1)*Ny+j]+rho_partial[(i-2)*Ny+j])*2  ,Sx_partial[(i-2)*Ny+j]/(rho_partial[(i-2)*Ny+j]+rho_partial[(i-3)*Ny+j])*2 ,Sx_partial[i*Ny+j]/(rho_partial[i*Ny+j]+rho_partial[(i-1)*Ny+j])*2 ,Sx_partial[i*Ny+j-1]/(rho_partial[i*Ny+j-1]+rho_partial[(i-1)*Ny+j-1])*2 ,Sx_partial[i*Ny+j-2]/(rho_partial[i*Ny+j-2]+rho_partial[(i-1)*Ny+j-2])*2 ,Sx_partial[i*Ny+j+1]/(rho_partial[i*Ny+j+1]+rho_partial[(i-1)*Ny+j+1])*2 ,Nx,Ny,dt,dl,1);

			rhoYm2_interp = vanLeerPBV(Vy[i*Ny+j+1],rho_partial[(i+1)*Ny+j+1] ,rho_partial[(i-1)*Ny+j+1] ,rho_partial[(i-2)*Ny+j+1] ,rho_partial[i*Ny+j+1] ,rho_partial[i*Ny+j] ,rho_partial[i*Ny+j-1] ,rho_partial[i*Ny+j+2],Nx,Ny,dt,dl,1);
			// Calculate Z velocity mass fluxes
			Vz_interp = vanLeerPBV(Vy[i*Ny+j],Sz_partial[(i+1)*Ny+j]/rho_partial[(i+1)*Ny+j] ,Sz_partial[(i-1)*Ny+j]/rho_partial[(i-1)*Ny+j] ,Sz_partial[(i-2)*Ny+j]/rho_partial[(i-2)*Ny+j] ,Sz_partial[i*Ny+j]/rho_partial[(i)*Ny+j] ,Sz_partial[i*Ny+j-1]/rho_partial[(i)*Ny+j-1] ,Sz_partial[i*Ny+j-2]/rho_partial[(i)*Ny+j-2] ,Sz_partial[i*Ny+j+1]/rho_partial[(i)*Ny+j+1] ,Nx,Ny,dt,dl,1);
			//Calculate F_1
			F_rho[i*Ny+j] = rho_interp * Vy[i*Ny+j];
			F_e[i*Ny+j]  = ed_interp*rho_interp * Vy[i*Ny+j];
			F_Vx[i*Ny+j]  = Vx_interp*.5*(rho_interp*Vy[i*Ny+j] + rhoXm2_interp*Vy[(i-1)*Ny+j]);
			F_Vy[i*Ny+j]  = Vy_interp*.5*(rho_interp*Vy[i*Ny+j] + rhoYm2_interp*Vy[(i)*Ny+j+1]);
			F_Vz[i*Ny+j]  = Vz_interp*rho_interp*Vy[i*Ny+j];
		}
	}//end calculating fluxes
	updatePeriodicBC(F_e,Nx,Ny);
	updatePeriodicBC(F_rho,Nx,Ny);
	updatePeriodicBC(F_Vx,Nx,Ny);
	updatePeriodicBC(F_Vy,Nx,Ny);
	updatePeriodicBC(F_Vz,Nx,Ny);
			if (view==1){
	ofstream myfile3 ("Fy.txt");
	myfile3<<"F_Vz:" <<"\n";
    for (int j=Nx-1; j>-1; j--)
    {
		for (int i=0; i<Ny; i++)
		{
			myfile3 <<setprecision(2)<<F_Vz[i*Ny+j]<<" ";
		}
		myfile3 <<"\n";
    }
    myfile3.close();
    }
	//Update values
	for (int i=2; i<Nx-2; i++)
	{
		for (int j=2; j<Ny-2; j++)
		{
			e_new[i*Ny+j]   = e_partial[i*Ny+j]   -       dt/dl * (F_e[i*Ny+j+1] - F_e[i*Ny+j]);
			rho_new[i*Ny+j] = rho_partial[i*Ny+j] -   dt/dl * (F_rho[i*Ny+j+1] - F_rho[i*Ny+j]);
			Vx_new[i*Ny+j]  = Sx_partial[i*Ny+j]  - dt/dl * (F_Vx[(i)*Ny+j+1] - F_Vx[(i)*Ny+j]);
			Vy_new[i*Ny+j]  = Sy_partial[i*Ny+j]  - dt/dl * (F_Vy[(i)*Ny+j] - F_Vy[(i)*Ny+j-1]);
			Vz_new[i*Ny+j]  = Sz_partial[i*Ny+j]  - dt/dl * (F_Vz[(i)*Ny+j+1] - F_Vz[(i)*Ny+j]);
		}
	}
	updatePeriodicBC(rho_new,Nx,Ny);
	updatePeriodicBC(e_new,Nx,Ny);
	updatePeriodicBC(Vx_new,Nx,Ny);
	updatePeriodicBC(Vy_new,Nx,Ny);
	updatePeriodicBC(Vz_new,Nx,Ny);

	for (int j=Nx-1; j>-1; j--)
    {
		for (int i=0; i<Ny; i++)
		{	Vz_new[i*Ny+j]=Vz_new[i*Ny+j]/rho_new[i*Ny+j];
			if (i==0)
			{
			    Vx_new[i*Ny+j]=Vx_new[i*Ny+j]/(rho_new[i*Ny+j]+rho_new[(Nx-3)*Ny+j])*2.;
			}
			else
			{
			    if(j==0)
			    {
			        Vy_new[i*Ny+j]=Vy_new[i*Ny+j]/(rho_new[i*Ny+j]+rho_new[(i)*Ny+Nx-3])*2.;
			    }
			    else
			    {
			        Vx_new[i*Ny+j]=Vx_new[i*Ny+j]/(rho_new[i*Ny+j]+rho_new[(i-1)*Ny+j])*2.;
			        Vy_new[i*Ny+j]=Vy_new[i*Ny+j]/(rho_new[i*Ny+j]+rho_new[(i)*Ny+j-1])*2.;
			    }
			}
		}
    }


	delete [] F_e;
	delete [] F_rho;
	delete [] F_Vx;
	delete [] F_Vy;
	delete [] F_Vz;
	delete [] e_partial;
	delete [] rho_partial;
	delete [] Sx_partial;
	delete [] Sy_partial;
	delete [] Sz_partial;
}


void updatePeriodicBC(double *arr,int Nx, int Ny)
{//Implement periodic boundary conditions in Ghost Zones-ooooo
	//Y direction
	for (int i=2; i<Nx-2; i++)
	{
		arr[i*Ny+0]    = arr[i*Ny+Ny-4];        //bottom row
		arr[i*Ny+1]    = arr[i*Ny+Ny-3];
		arr[i*Ny+Ny-2] = arr[i*Ny+2];	        //top row
		arr[i*Ny+Ny-1] = arr[i*Ny+3];
	}
	//X Direction

	for (int j=2; j<Ny-2; j++)
	{
		arr[0*Ny+j]      = arr[(Ny-4)*Ny+j];    //left column
		arr[1*Ny+j]      = arr[(Ny-3)*Ny+j];
		arr[(Nx-2)*Ny+j] = arr[2*Ny+j];	        //right column
		arr[(Nx-1)*Ny+j] = arr[3*Ny+j];
	}
	//BL corner
	arr[0*Ny+0] = arr[(0+Ny-4)*Ny+0+Ny-4];
	arr[0*Ny+1] = arr[(0+Ny-4)*Ny+1+Ny-4];
	arr[1*Ny+0]	= arr[(1+Ny-4)*Ny+0+Ny-4];
	arr[1*Ny+1] = arr[(1+Ny-4)*Ny+1+Ny-4];
	//BR corner
	arr[(0+Ny-2)*Ny+0] = arr[(0+2)*Ny+0+Ny-4];
	arr[(0+Ny-2)*Ny+1] = arr[(0+2)*Ny+1+Ny-4];
	arr[(1+Ny-2)*Ny+0] = arr[(1+2)*Ny+0+Ny-4];
	arr[(1+Ny-2)*Ny+1] = arr[(1+2)*Ny+1+Ny-4];
	//TL corner
	arr[0*Ny+0+Ny-2] = arr[(0+Ny-4)*Ny+0+2];
	arr[0*Ny+1+Ny-2] = arr[(0+Ny-4)*Ny+1+2];
	arr[1*Ny+0+Ny-2] = arr[(1+Ny-4)*Ny+0+2];
	arr[1*Ny+1+Ny-2] = arr[(1+Ny-4)*Ny+1+2];
	//TR corner
	arr[(0+Ny-2)*Ny+0+Ny-2] = arr[(0+2)*Ny+0+2];
	arr[(0+Ny-2)*Ny+1+Ny-2] = arr[(0+2)*Ny+1+2];
	arr[(1+Ny-2)*Ny+0+Ny-2] = arr[(1+2)*Ny+0+2];
	arr[(1+Ny-2)*Ny+1+Ny-2] = arr[(1+2)*Ny+1+2];
}


double vanLeerPBV(double V,double ip1j,double im1j,double im2j,double ij,double ijm1,double ijm2,double ijp1,int Nx, int Ny,double dt,double dl,int direction)
{//Currently just runs Donor Cell method. Uncomment the pointed out lines for vanLeer
	double rho_interp   	= 0.;
	double d_rho			= 0.;
	double delta_rho_plus  	= 0.;
	double delta_rho_minus 	= 0.;
		//Calculate determine what is upwind
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
					rho_interp = im1j;// + (1.-V*dt/dl)*d_rho/2.;  // <--- Uncomment for vanLeer
				else
					rho_interp = ijm1;//  + (1.-V*dt/dl)*d_rho/2.; // <--- Uncomment for vanLeer

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
				rho_interp = ij;// - (1.+V*dt/dl)*d_rho/2.;  //<--- Uncomment for vanLeer
			}
	return rho_interp;
}

