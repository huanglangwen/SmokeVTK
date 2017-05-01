/* Author: Johannes Schmid, 2006, johnny@grob.org */
#include "fluid.h"
#include <assert.h>
#include <math.h>

inline int Fluid::_I(int x, int y, int z) { return x + y*(Nx+2) + z*(Nx+2)*(Ny+2); };

Fluid::Fluid()
{
	int i;

	for (i=0; i<15; i++)
		clear_buffer(buffers[i]);

	i=0;
	d=buffers[i++]; d0=buffers[i++];
	T=buffers[i++]; T0=buffers[i++];
	u=buffers[i++]; u0=buffers[i++];
	v=buffers[i++]; v0=buffers[i++];
	w=buffers[i++]; w0=buffers[i++];
	curltemp = buffers[i++];
	dryfallingW = buffers[i++];

	//Initialize Tenv
	for (i = 0;i < Nz;i++)Tenv[i] = 0;
	//Set Original Wind Field
	
	for (i = 0;i < SIZE;i++) {
		u[i] = uBoundary;
		v[i] = vBoundary;
	}
	
	clear_sources();
}

Fluid::~Fluid()
{
}

//#define NO_BOUNDARY
inline void Fluid::set_bnd(int b, double* x)
{
#ifndef NO_BOUNDARY
	int i, j, k;

	for (i = 1;i <= Nx;i++) {
		for (j = 1;j <= Ny;j++) {
			x[_I(i, j, 0)] = (b == 3) ? -x[_I(i, j, 1)] : x[_I(i, j, 1)];
			x[_I(i, j, Nz + 1)] = (b == 3) ? -x[_I(i, j, Nz)] : x[_I(i, j, Nz)];
		}
	}

	if (b > 0)
	{
		for (i = 1;i <= Nx;i++) {
			for (k = 1;k <= Nz;k++) {
				x[_I(i, 0, k)] = (b == 2) ? vBoundary - x[_I(i, 1, k)] : x[_I(i, 1, k)];
				x[_I(i, Ny + 1, k)] = (b == 2) ? vBoundary - x[_I(i, Ny, k)] : x[_I(i, Ny, k)];
			}
		}

		for (j = 1;j <= Ny;j++) {
			for (k = 1;k <= Nz;k++) {
				x[_I(0, j, k)] = (b == 1) ? uBoundary - x[_I(1, j, k)] : x[_I(1, j, k)];
				x[_I(Nx + 1, j, k)] = (b == 1) ? uBoundary - x[_I(Nx, j, k)] : x[_I(Nx, j, k)];
			}
		}
	}
	else 
	{
		for (i = 1;i <= Nx;i++) {
			for (k = 1;k <= Nz;k++) {
				x[_I(i, 0, k)] = x[_I(i, 1, k)];
				x[_I(i, Ny + 1, k)] = x[_I(i, Ny, k)];
			}
		}

		for (j = 1;j <= Ny;j++) {
			for (k = 1;k <= Nz;k++) {
				x[_I(0, j, k)] = x[_I(1, j, k)];
				x[_I(Nx + 1, j, k)] = x[_I(Nx, j, k)];
			}
		}
	}

	x[_I(0, 0, 0)] = (x[_I(1, 0, 0)] + x[_I(0, 1, 0)] + x[_I(0, 0, 1)]) / 3;
	x[_I(0, Ny + 1, 0)] = (x[_I(1, Ny + 1, 0)] + x[_I(0, Ny, 0)] + x[_I(0, Ny + 1, 1)]) / 3;
	x[_I(Nx + 1, 0, 0)] = (x[_I(Nx, 0, 0)] + x[_I(Nx + 1, 1, 0)] + x[_I(Nx + 1, 0, 1)]) / 3;
	x[_I(Nx + 1, Ny + 1, 0)] = (x[_I(Nx, Ny + 1, 0)] + x[_I(Nx + 1, Ny, 0)] + x[_I(Nx + 1, Ny + 1, 1)]) / 3;
	x[_I(0, 0, Nz + 1)] = (x[_I(1, 0, Nz + 1)] + x[_I(0, 1, Nz + 1)] + x[_I(0, 0, Nz)]) / 3;
	x[_I(0, Ny + 1, Nz + 1)] = (x[_I(1, Ny + 1, Nz + 1)] + x[_I(0, Ny, Nz + 1)] + x[_I(0, Ny + 1, Nz)]) / 3;
	x[_I(Nx + 1, 0, Nz + 1)] = (x[_I(Nx, 0, Nz + 1)] + x[_I(Nx + 1, 1, Nz + 1)] + x[_I(Nx + 1, 0, Nz)]) / 3;
	x[_I(Nx + 1, Ny + 1, Nz + 1)] = (x[_I(Nx, Ny + 1, Nz + 1)] + x[_I(Nx + 1, Ny, Nz + 1)] + x[_I(Nx + 1, Ny + 1, Nz)]) / 3;

	
	#endif
}

void Fluid::add_source(double* src, double *dst, double dt)
{
	int i;

	for (i=0; i<SIZE; i++)
		dst[i] += src[i]*dt;
}

void Fluid::add_buoyancy(double dt)
{
	int i,j,k;

	//for (i=0; i<SIZE; i++)
		//v[i] += -d[i]*buoyancy*dt;
		//w[i] += d[i] * buoyancy*dt;

	for (k = 0;k < Nz + 2;k++) {
		for (j = 0;j < Ny + 2;j++) {
			for (i = 0;i < Nx + 2;i++) {
				w[_I(i, j, k)] += (T[_I(i, j, k)]-Tenv[k]) / (Tenv[k] + 298)*9.8*dt;//theta'/theta_bar*g
			}
		}
	}
}

inline void Fluid::diffuse(int b, double* x0, double* x, double diff, double dt)
{
	int i, j, k, l;
	double a=dt*diff/(delta*delta);//a=dt*diff*N*N*N
	for (l=0; l<20; l++) 
	{
		for (k=1; k<=Nz; k++)
		{
			for (j=1; j<=Ny; j++)
			{
				for (i=1; i<=Nx; i++)
				{
					x[_I(i,j,k)] = (x0[_I(i,j,k)] + a*(
						x[_I(i-1,j,k)]+x[_I(i+1,j,k)]+
						x[_I(i,j-1,k)]+x[_I(i,j+1,k)]+
						x[_I(i,j,k-1)]+x[_I(i,j,k+1)]))/(1+6*a);
				}
			}
		}
		set_bnd(b,x);
	}
}

inline void Fluid::advect(int b, double* x0, double* x, double* uu, double* vv, double* ww, double dt)
{
	int i, j, k, i0, j0, k0, i1, j1, k1;
	double sx0, sx1, sy0, sy1, sz0, sz1, v0, v1;
	double xx, yy, zz, dtx0, dty0, dtz0;
	//dt0 = dt*N;
	dtx0 = dt / delta;
	dty0 = dt / delta;
	dtz0 = dt / delta;
	for (k=1; k<=Nz; k++)
	{
		for (j=1; j<=Ny; j++)
		{
			for (i=1; i<=Nx; i++)
			{
				xx = i-dtx0*uu[_I(i,j,k)];
				yy = j-dty0*vv[_I(i,j,k)];
				zz = k-dtz0*ww[_I(i,j,k)];
				if (xx<0.5) xx=0.5f; if (xx>Nx+0.5) xx=Nx+0.5f; i0=(int)xx; i1=i0+1;
				if (yy<0.5) yy=0.5f; if (yy>Ny+0.5) yy=Ny+0.5f; j0=(int)yy; j1=j0+1;
				if (zz<0.5) zz=0.5f; if (zz>Nz+0.5) zz=Nz+0.5f; k0=(int)zz; k1=k0+1;
				sx1 = xx-i0; sx0 = 1-sx1;
				sy1 = yy-j0; sy0 = 1-sy1;
				sz1 = zz-k0; sz0 = 1-sz1;
				v0 = sx0*(sy0*x0[_I(i0,j0,k0)]+sy1*x0[_I(i0,j1,k0)])+sx1*(sy0*x0[_I(i1,j0,k0)]+sy1*x0[_I(i1,j1,k0)]);
				v1 = sx0*(sy0*x0[_I(i0,j0,k1)]+sy1*x0[_I(i0,j1,k1)])+sx1*(sy0*x0[_I(i1,j0,k1)]+sy1*x0[_I(i1,j1,k1)]);
				x[_I(i,j,k)] = sz0*v0 + sz1*v1;
			}
		}
	}
	set_bnd(b,d);
}

void Fluid::project(void)
{
	double* p = u0;	double* div = v0;	// temporary buffers, use old velocity buffers
	int i, j, k, l;
	constexpr double h = delta;
	//h = 1.0f/N;
	for (k=1; k<=Nz; k++) {
		for (j=1; j<=Ny; j++) {
			for (i=1; i<=Nx; i++) {
				div[_I(i,j,k)] = -h*(
					u[_I(i+1,j,k)]-u[_I(i-1,j,k)]+
					v[_I(i,j+1,k)]-v[_I(i,j-1,k)]+
					w[_I(i,j,k+1)]-w[_I(i,j,k-1)])/3;
				p[_I(i,j,k)] = 0;
			}
		}
	}
	set_bnd(0,div); set_bnd(0,p);//这里也有边界！！！
	for (l=0; l<20; l++) 
	{
		for (k=1; k<=Nz; k++) {
			for (j=1; j<=Ny; j++) {
				for (i=1; i<=Nx; i++) {
					p[_I(i,j,k)] = (div[_I(i,j,k)]+
						p[_I(i-1,j,k)]+p[_I(i+1,j,k)]+
						p[_I(i,j-1,k)]+p[_I(i,j+1,k)]+
						p[_I(i,j,k-1)]+p[_I(i,j,k+1)])/6;
				}
			}
		}
		set_bnd(0,p);//这里也有！！！
	}
	for (k=1; k<=Nz; k++) {
		for (j=1; j<=Ny; j++) {
			for (i=1; i<=Nx; i++) {
				u[_I(i,j,k)] -= (p[_I(i+1,j,k)]-p[_I(i-1,j,k)])/3/h;
				v[_I(i,j,k)] -= (p[_I(i,j+1,k)]-p[_I(i,j-1,k)])/3/h;
				w[_I(i,j,k)] -= (p[_I(i,j,k+1)]-p[_I(i,j,k-1)])/3/h;
			}
		}
	}
	set_bnd(1,u); set_bnd(2,v);//这里会使u v变化！！！！
}

void Fluid::vorticity_confinement(double dt)
{
	int i,j,k,ijk;
	double *curlx = u0, *curly = v0, *curlz = w0, *curl = curltemp;		// temp buffers 内存泄漏！！！
	double dt0 = dt * vc_eps * delta;
	double x,y,z;


	for (k=1; k<=Nz; k++) {
		for (j=1; j<=Ny; j++) {
			for (i=1; i<=Nx; i++) {
				ijk = _I(i,j,k);
					// curlx = dw/dy - dv/dz
				x = curlx[ijk] = (w[_I(i,j+1,k)] - w[_I(i,j-1,k)]) * 0.5f -
					(v[_I(i,j,k+1)] - v[_I(i,j,k-1)]) * 0.5f;

					// curly = du/dz - dw/dx
				y = curly[ijk] = (u[_I(i,j,k+1)] - u[_I(i,j,k-1)]) * 0.5f -
					(w[_I(i+1,j,k)] - w[_I(i-1,j,k)]) * 0.5f;

					// curlz = dv/dx - du/dy
				z = curlz[ijk] = (v[_I(i+1,j,k)] - v[_I(i-1,j,k)]) * 0.5f -
					(u[_I(i,j+1,k)] - u[_I(i,j-1,k)]) * 0.5f;

					// curl = |curl|
				curl[ijk] = sqrtf(x*x+y*y+z*z);
			}
		}
	}

	for (k=1; k<=Nz; k++) {
		for (j=1; j<=Ny; j++) {
			for (i=1; i<=Nx; i++) {
				ijk = _I(i,j,k);
				double Nx = (curl[_I(i+1,j,k)] - curl[_I(i-1,j,k)]) * 0.5f;
				double Ny = (curl[_I(i,j+1,k)] - curl[_I(i,j-1,k)]) * 0.5f;
				double Nz = (curl[_I(i,j,k+1)] - curl[_I(i,j,k-1)]) * 0.5f;
				double len1 = 1.0f/(sqrtf(Nx*Nx+Ny*Ny+Nz*Nz)+0.0000001f);
				Nx *= len1;
				Ny *= len1;
				Nz *= len1;
				u[ijk] += (Ny*curlz[ijk] - Nz*curly[ijk]) * dt0;
				v[ijk] += (Nz*curlx[ijk] - Nx*curlz[ijk]) * dt0;
				w[ijk] += (Nx*curly[ijk] - Ny*curlx[ijk]) * dt0;
			}
		}
	}
	for (i = 0;i < SIZE;i++) {
		curlx[i] = 0.0f;
		curly[i] = 0.0f;
		curlz[i] = 0.0f;
		curl[i] = 0.0f;
	}
}

#define DIFFUSE
#define ADVECT

void Fluid::vel_step(double dt)
{
	add_source(su, u, dt);
	add_source(sv, v, dt);
	add_source(sw, w, dt);
	add_buoyancy(dt);
	vorticity_confinement(dt);

#ifdef DIFFUSE
	SWAPFPTR(u0, u); SWAPFPTR(v0, v); SWAPFPTR(w0, w);
	diffuse(1, u0, u, viscosity, dt);
	diffuse(2, v0, v, viscosity, dt);
	diffuse(3, w0, w, viscosity, dt);
	project();
#endif
#ifdef ADVECT
	SWAPFPTR(u0, u); SWAPFPTR(v0, v); SWAPFPTR(w0, w);
	advect(1, u0, u, u0, v0, w0, dt);
	advect(2, v0, v, u0, v0, w0, dt);
	advect(3, w0, w, u0, v0, w0, dt);
	project();
#endif
}

void Fluid::dens_temp_step(double dt)
{
	add_source(sd, d, dt);
	add_source(sT, T, dt);
#ifdef DIFFUSE
	SWAPFPTR(d0, d);
	diffuse(0, d0, d, diffusion, dt);
	SWAPFPTR(T0, T);
	diffuse(0, T0, T, conduction, dt);
#endif
#ifdef ADVECT
	SWAPFPTR(d0, d);
	//dry falling d=0.02mm w-=0.0314
	for (int i = 0;i < SIZE;i++)dryfallingW[i] = w[i] - terminalW;
	set_bnd(3, dryfallingW);//设置边界条件
	advect(0, d0, d, u, v, dryfallingW, dt);
	SWAPFPTR(T0, T);
	advect(0, T0, T, u, v, w, dt);
#endif
}

void Fluid::step(double dt)
{
	vel_step(dt);
	dens_temp_step(dt);
}



void Fluid::clear_buffer(double* x)
{
	for (int i=0; i<SIZE; i++) {
		x[i] = 0.0f;
	}
}

void Fluid::clear_sources(void)
{
	for (int i=0; i<SIZE; i++) {
		sT[i] = 0.0f;
		sd[i] = 0.0f;
		su[i] = 0.0f;
		sv[i] = 0.0f;
	}
}

void Fluid::store(int time)
{
	char buff[100];
	snprintf(buff, sizeof(buff), "densityTemperature%d.vti", time);
	std::string densTempName = buff;
	snprintf(buff, sizeof(buff), "Wind%d.vti", time);
	std::string windName = buff;
	auto densTempData = vtkSmartPointer<vtkImageData>::New();
	auto windData = vtkSmartPointer<vtkImageData>::New();
	densTempData->SetDimensions(Nx + 2, Ny + 2, Nz + 2);
	densTempData->AllocateScalars(VTK_DOUBLE, 2);
	windData->SetDimensions(Nx + 2, Ny + 2, Nz + 2);
	windData->AllocateScalars(VTK_DOUBLE, 3);
	int* dims = densTempData->GetDimensions();
	// Fill every entry of the image data with d
	for (int z = 0; z < dims[2]; z++)
	{
		for (int y = 0; y < dims[1]; y++)
		{
			for (int x = 0; x < dims[0]; x++)
			{
				double* pixel1 = static_cast<double*>(densTempData->GetScalarPointer(x, y, z));
				double* pixel2 = static_cast<double*>(windData->GetScalarPointer(x, y, z));
				pixel1[0] = d[_I(x,y,z)];
				pixel1[1] = T[_I(x, y, z)];
				pixel2[0] = u[_I(x, y, z)];
				pixel2[1] = v[_I(x, y, z)];
				pixel2[2] = w[_I(x, y, z)];
			}
		}
	}
	auto writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
#define OUTPUT_DENSTEMP
#ifdef OUTPUT_DENSTEMP
	writer->SetFileName(densTempName.c_str());
	writer->SetInputData(densTempData);
	writer->Write();
#endif
#ifdef OUTPUT_WIND
	writer->SetFileName(windName.c_str());
	writer->SetInputData(windData);
	writer->Write();
#endif
}
