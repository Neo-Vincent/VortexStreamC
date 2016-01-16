
import scipy.sparse.linalg as lin
from scipy.sparse import csr_matrix
import numpy as np
Re=400.
sorwei=1.8;
dt=1.e-3


class ns2d():
    def __init__(self,tnx,tny):
        self.nx=tnx
        self.ny=tny
        self.dx= 1./self.nx
        self.dy= 1./self.ny
        self.dx2inv= 1. /(self.dx*self.dx)
        self.dy2inv= 1. /(self.dy*self.dy)
        self.dsinv = 0.5/(self.dx2inv+self.dy2inv)
        self.dxinv= 1./self.dx
        self.dyinv= 1./self.dy
        self.u  =[[0 for i in range(self.nx+1)] for j in range(self.ny+1) ]
        self.v  = [[0 for i in range(self.nx+1)] for j in range(self.ny+1) ]
        self.vor= [[0 for i in range(self.nx+1)] for j in range(self.ny+1) ]
        self.stf= [[0 for i in range(self.nx+1)] for j in range(self.ny+1) ]
        self.cgiter=0
        for i in range(self.nx+1):
            for j in range(self.ny+1):
                self.u  [i][j]= 0.
                self.v  [i][j]= 0.
                self.vor[i][j]= 0.
                self.stf[i][j]= 0.
        for i in range(1,self.nx):
            self.u  [i][self.ny]= 1.

              
        a=-1./self.dsinv
        data=[]
        omega=[]
        indptr=[0]
        indices=[]
        for i in range(self.nx+1):
            for j in range(self.ny+1):
                if i==0 or j==0 or i==self.nx or j==self.ny:
                    omega.append(self.stf[i][j])
                    data.append(1.)
                    indices.append(i*(self.nx+1)+j)
                    indptr.append(len(indices))
                else:
                    data.append(self.dx2inv)
                    indices.append((i-1)*(self.nx+1)+j)
                    data.append(self.dy2inv)
                    indices.append(i*(self.nx+1)+j-1)
                    data.append(a)
                    indices.append(i*(self.nx+1)+j)
                    data.append(self.dy2inv)
                    indices.append(i*(self.nx+1)+j+1)
                    data.append(self.dx2inv)
                    indices.append((i+1)*(self.nx+1)+j)
                    omega.append(self.vor[i][j])
                    indptr.append(len(indices))
        self.operator=csr_matrix((np.array(data),np.array(indices),np.array(indptr)),shape=((self.nx+1)*(self.ny+1),(self.nx+1)*(self.ny+1)))
        self.omega=omega

    def ns_iter(self):
        for i in range(1,self.nx):
            self.vor[i][0] = self.dyinv*( self.u[i][1]  - self.u[i][0] )
            self.vor[i][self.ny]= self.dyinv*( self.u[i][self.ny] - self.u[i][self.ny-1] )
        for j in range(1,self.ny):
            self.vor[0][j] = -self.dxinv*( self.v[1] [j] - self.v[0]   [j] )
            self.vor[self.nx][j]= -self.dxinv*( self.v[self.nx][j] - self.v[self.nx-1][j] )
         # using GMRES to solve the stream function

        for i in range(1,self.nx):
            for j in range(1,self.ny):
                fxr= ( self.vor[i+1][j  ] - self.vor[i  ][j  ] ) * self.dxinv
                fxl= ( self.vor[i  ][j  ] - self.vor[i-1][j  ] ) * self.dxinv
                fyr= ( self.vor[i  ][j+1] - self.vor[i  ][j  ] ) * self.dyinv
                fyl= ( self.vor[i  ][j  ] - self.vor[i  ][j-1] ) * self.dyinv
                flux= -0.5*(( self.u[i][j] - abs(self.u[i][j] ))*fxr + (self.u[i][j] + abs(self.u[i][j]))*fxl )
                fluy= -0.5*(( self.v[i][j] - abs(self.v[i][j] ))*fyr + (self.v[i][j] + abs(self.v[i][j]))*fyl )
            
                self.vor[i][j]= self.vor[i][j]+ dt*( flux + fluy +\
					    1./Re*( self.dx2inv*( self.vor[i+1][j]- 2.*self.vor[i][j]+ self.vor[i-1][j] ) +\
							    self.dy2inv*( self.vor[i][j+1]- 2.*self.vor[i][j]+ self.vor[i][j-1] )  ) )
        omega=[]
        for i in range(self.nx+1):
            for j in range(self.ny+1):
                if i==0 or j==0 or i==self.nx or j==self.ny:
                    omega.append(self.stf[i][j])
                else:
                    omega.append(self.vor[i][j])
        slt=lin.cg(self.operator ,omega)
        self.cgiter=slt[1]
        solution=slt[0]
        for i in range(self.nx):
            for j in range(self.ny+1):
                self.stf[i][j]=solution[i*(self.nx+1)+j]
        for i in range(1,self.nx):
            for j in range(1,self.ny):
                self.u[i][j]= 0.5* self.dyinv * ( self.stf[i][j+1] - self.stf[i][j-1] )
                self.v[i][j]=-0.5* self.dxinv * ( self.stf[i+1][j] - self.stf[i-1][j] )

    def ns2d_solve(self):
        ttime=0.;
        for n in range(5000):
            self.ns_iter()
            print('n=',n,',iter=',self.cgiter,'\n')
            ttime= ttime+dt
            if( n%500==0 ):
                print('n=',n,'\n')
                self.ns2d_output()


    def ns2d_output(self):
        f=open('ns.dat',mode='w')
        f.write(" variables=\"x\",\"y\",\"u\",\"v\",\"vorticity\",\"stream line\"\n")
        f.write(" zone i={0}   j={1} f=point \n".format(self.nx+1,self.ny+1))
        for i in range(self.nx+1):
            for j in range(self.ny+1):
                f.write("{0},{1},{2},{3},{4},{5} \n ".format(self.dx*i,self.dy*j,self.u[i][j],self.v[i][j],self.vor[i][j],self.stf[i][j]))
a=ns2d(100,100)
a.ns2d_solve()
