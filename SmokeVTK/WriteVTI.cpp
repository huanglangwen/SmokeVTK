#include "vtkAutoInit.h" 
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
//VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);

#include <vtkVersion.h>

#define RENDER_ON
//#define STORE_ON
#ifdef RENDER_ON
#include <vtkSmartVolumeMapper.h>
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVolumeProperty.h>
#include <vtkCamera.h>
#endif

#include <cstdio>
//#include <Windows.h>
#include "fluid.h"

constexpr int totaltime = 25;
constexpr double dt = 0.1;
constexpr int xpos = 10;
constexpr int ypos = 28;
constexpr int zpos = 3;
constexpr double Tenv_temp[Nz + 2] = { 5.        ,  4.89795918,  4.79591837,  4.69387755,  4.59183673,
4.48979592,  4.3877551 ,  4.28571429,  4.18367347,  4.08163265,
3.97959184,  3.87755102,  3.7755102 ,  3.67346939,  3.57142857,
3.46938776,  3.36734694,  3.26530612,  3.16326531,  3.06122449,
2.95918367,  2.85714286,  2.75510204,  2.65306122,  2.55102041,
2.44897959,  2.34693878,  2.24489796,  2.14285714,  2.04081633,
1.93877551,  1.83673469,  1.73469388,  1.63265306,  1.53061224,
1.42857143,  1.32653061,  1.2244898 ,  1.12244898,  1.02040816,
0.91836735,  0.81632653,  0.71428571,  0.6122449 ,  0.51020408,
0.40816327,  0.30612245,  0.20408163,  0.10204082,  0.        ,
0.        ,   0.40816327,   0.81632653,   1.2244898 ,
1.63265306,   2.04081633,   2.44897959,   2.85714286,
3.26530612,   3.67346939,   4.08163265,   4.48979592,
4.89795918,   5.30612245,   5.71428571,   6.12244898,
6.53061224,   6.93877551,   7.34693878,   7.75510204,
8.16326531,   8.57142857,   8.97959184,   9.3877551 ,
9.79591837,  10.20408163,  10.6122449 ,  11.02040816,
11.42857143,  11.83673469,  12.24489796,  12.65306122,
13.06122449,  13.46938776,  13.87755102,  14.28571429,
14.69387755,  15.10204082,  15.51020408,  15.91836735,
16.32653061,  16.73469388,  17.14285714,  17.55102041,
17.95918367,  18.36734694,  18.7755102 ,  19.18367347,
19.59183673,  20. };//5(0)->0(50)->20(100)

int main(int argc, char *argv[])
{
	Fluid* fluid = NULL;
	fluid = new Fluid();
	fluid->diffusion = 0.00001f;//0.00001f;
	fluid->viscosity = 0.000001f;
	fluid->buoyancy = 1.0f;
	fluid->vc_eps = 300.0f;
	fluid->conduction = 0.000024f;//0.000024f;//¿ÕÆøÈÈÀ©É¢ÂÊ
	fluid->terminalW = 0.08;//m/s for 0.05mm particles
	
#ifdef RENDER_ON
	auto renWin = vtkSmartPointer<vtkRenderWindow>::New();
	auto ren1 = vtkSmartPointer<vtkRenderer>::New();
	ren1->SetBackground(0.45,0.45,0.45);
	renWin->AddRenderer(ren1);

	renWin->SetSize(301, 300); // intentional odd and NPOT  width/height

	auto iren =	vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renWin);

	renWin->Render(); // make sure we have an OpenGL context.

	auto volumeMapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();
	volumeMapper->SetBlendModeToComposite(); // composite first

	auto volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
	volumeProperty->ShadeOff();
	volumeProperty->SetInterpolationType(VTK_LINEAR_INTERPOLATION);

	auto compositeOpacity =	vtkSmartPointer<vtkPiecewiseFunction>::New();
	compositeOpacity->AddPoint(0.0, 0.0);
	compositeOpacity->AddPoint(2.0, 1.0);
	volumeProperty->SetScalarOpacity(compositeOpacity); // composite first.

	auto color = vtkSmartPointer<vtkColorTransferFunction>::New();
	color->AddRGBPoint(0.0, 0.5, 0.5, 0.5);
	color->AddRGBPoint(2.0, 1.0, 1.0, 1.0);
	volumeProperty->SetColor(color);

	auto volume = vtkSmartPointer<vtkVolume>::New();
	volume->SetMapper(volumeMapper);
	volume->SetProperty(volumeProperty);
	ren1->AddViewProp(volume);
	ren1->ResetCamera();

	//volumeMapper->SetRequestedRenderModeToGPU();

#endif
	

	int count = 0;
	for (double t = 0;t < totaltime;t+=dt) 
	{

		int pos;
		for (int k = 0;k < 3;k++)
		{
			for (int i = 0; i < 8; i++)
			{
				for (int j = 0; j < 8; j++)
				{
					double dd = (rand() % 1000) / 1000.0f;
					double T = (rand() % 1000) / 1000.0f * 10;
					//f = genfunc(i,j,8,8,t,gfparams);
					//fluid->d[_I(xpos, i + ypos, j + 28)] = d;
					//fluid->u[_I(xpos,i+ypos,j+28)] = -(1.0f + 2.0f*f);
					//fluid->u[_I(xpos, i + ypos, j + 28)] = -3.0f;
					pos = fluid->_I(i + xpos, j + ypos, k + zpos);
					fluid->d[pos] = dd + 1;
					fluid->w[pos] = 1.5f;
					fluid->T[pos] = T+1;
				}
			}
		}
	/*	
		for (int i = 0;i < Nx + 2;i++) {
			for (int k = 0;k < Nz+2;k++) {
				double vv = (rand() % 1000) / 1000.0f;
				fluid->v[fluid->_I(i, 0, k)] = vv+0.5;
			}
		}

		for (int j = 0;j < Ny + 2;j++) {
			for (int k = 0;k < Nz+2;k++) {
				double uu = (rand() % 1000) / 1000.0f;
				fluid->u[fluid->_I(0, j, k)] = uu+0.5;
			}
		}
		*/
		fluid->step(dt);
		printf("Time: %f\n",t);
		count++;

#ifdef RENDER_ON
		volumeMapper->SetInputData(fluid->getDensTempData(1));
		renWin->Render();
#endif
#ifdef STORE_ON
		if (count % 10 == 0) {
			fluid->store(count/10);
		}
#endif
	}

	



	//system("Pause");
	return EXIT_SUCCESS;
}
