import numpy as np

class Fluid:
    d=0.1
    lx=1
    ly=1
    lz=1
    Nx=int(lx/d)
    Ny=int(ly/d)
    Nz=int(lz/d)
    density=np.zeros([Nx+2,Ny+2,Nz+2],dtype=np.float32)
    u=np.zeros([Nx+2,Ny+2,Nz+2],dtype=np.float32)
    v=np.zeros([Nx+2,Ny+2,Nz+2],dtype=np.float32)
    w=np.zeros([Nx+2,Ny+2,Nz+2],dtype=np.float32)
    diffusion=0.0
    viscosity=0.0
    buoyancy=0.0
    vc_eps=5.0
    MAX_LINEAR_SOLVER_ITERATION=20

    def __init__(self,diff,visc,buoy,vc_eps):
        self.diffusion=diff
        self.viscosity=visc
        self.buoyancy=buoy
        self.vc_eps=vc_eps

    def set_bnd(b,x):
        
    
    def diffuse(self,b,x0,x,diff,dt):
        a=dt*diff/(self.d**3)
        for l in range(self.MAX_LINEAR_SOLVER_ITERATION):
            for i in range(self.Nx):
                for j in range(self.Ny):
                    for k in range(self.Nz):
                        x[i,j,k]=(x0[i,j,k]+a*(
                            x[i-1,j,k]+x[i+1,j,k]+
                            x[i,j-1,k]+x[i,j+1,k]+
                            x[i,j,k-1]+x[i,j,k+1]
                        ))/(1+6*a)
        set_bnd(b,x)
    