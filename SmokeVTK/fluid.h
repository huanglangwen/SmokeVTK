/* Author: Johannes Schmid, 2006, johnny@grob.org */
#ifndef __FLUID_H
#define __FLUID_H

#include <stdio.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkProperty.h>
#include <vtkXMLImageDataWriter.h>

//#define N 62				// must be N^2-2
//#define SIZE ((N+2)*(N+2)*(N+2))
//#define _I(x,y,z) (((z)<<10)+((y)<<5)+x)
//#define _I(x,y,z) (((z)<<12)+((y)<<6)+x)


#define SWAPFPTR(x,y) {double *t=x;x=y;y=t;}

constexpr int Nx = 200;
constexpr int Ny = 60;
constexpr int Nz = 98;
constexpr int SIZE = (Nx+2)*(Ny+2)*(Nz+2);
constexpr double delta = 0.01612;
constexpr double uBoundary = 0.5f;
constexpr double vBoundary = 0.0f;

typedef double Field[SIZE];

class Fluid
{
protected:
	double buffers[15][SIZE];

public:
	double *T, *T0;         // temperature
	double *d, *d0;			// density
	double *u, *u0;			// velocity in x direction
	double *v, *v0;			// velocity in y direction
	double *w, *w0;			// velocity in z direction
	double Tenv[Nz+2];
	double *curltemp;
	double *dryfallingW;    // temporary w for dry falling of density w-=0.314

protected:
	// simulation methods
		// beware: i changed stam's implementation from a destiation, source ordering
		// of the parameters to source, destiation
	void add_source(double* src, double *dst, double dt);
	void add_buoyancy(double dt);
	void set_bnd(int b, double* x);
	void diffuse(int b, double* x0, double* x, double diff, double dt);
	void advect(int b, double* x0, double* x, double* uu, double* vv, double* ww, double dt);
	void project(void);
	void vorticity_confinement(double dt);

	void vel_step(double dt);
	//void dens_step(double dt);
	//void temp_step(double dt);
	void dens_temp_step(double dt);

	// utility methods
	void clear_buffer(double* x);
	void clear_sources(void);

public:
	Field sd, su, sv, sw, sT;	// sources for density and velocities
	double diffusion, conduction, viscosity, buoyancy, vc_eps;

	int _I(int x, int y, int z);

	Fluid(void);
	~Fluid(void);

	void step(double dt);

	//void store(vtkSmartPointer<vtkImageData> fp);
	void store(int time);

};

#endif

