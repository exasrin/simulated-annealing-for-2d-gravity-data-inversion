# -----------------------------------------------------------------------
#  Program Algoritma Simulated Annealing
#  Refererensi #  Refererensi Won Y. Yang dkk
#  Dibuat Oleh Asrin
#  Kendari  2022
# -----------------------------------------------------------------------

import time

import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib
import pandas as pd

start_time = time.time()
# Calling the input model parameters
itin = pd.read_excel('Fault_Model.xlsx', sheet_name='Sheet1')
tin = itin.head()
#Calling the Gravity Anomaly
grav=np.genfromtxt('fault model.dat')
f=np.array(grav[:,1])
xxx = np.array(grav[:,0])

# Parameter Model Blok 2-D
nlay=int(itin.loc[0,'nb'])   #Banyak grid
nxg=int(itin.loc[0,'nx'])    #Banyak cell lateral
nzg=int(itin.loc[0,'nz'])    #Banyak cell vetikal
dx =int(itin.loc[0,'dx']) #Dimensi cell lateral (m)
dh =int(itin.loc[0,'dz']) #Dimensi cell vetikal (m)
model=np.array([nxg,nzg,dx,dh])


l=itin.loc[1,:] #lower earch of paramater
l=np.array([l[0:int(nlay)]])
u=itin.loc[2,:] #upper search of paramter
u=np.array([u[0:int(nlay)]])

x0=l+(u-l)/2 # initial guess
kmax=1000 #Jumlah Iterasi
q=1 
TolFun=1e-8

def forward_grav(rho,model):
    
    G=6.6732e-11   # Kontanta Gravirtasi damalm Mks
    dx = model[3] #(m)
    dh = model[2] #(m)
    nx = model[0]
    nz = model[1]
    d = dx #station spacing
    h = dh #depth spacing (thickness)

    V = rho
    V = np.reshape(rho,(nz,nx))
    VV = V.transpose().ravel()
    
    x1 = np.arange(0,(nx-1)*dx+dx,dx)
    x = x1+d/2
    z1 = np.transpose([np.arange(0,((nz+1)*h-10),h)])
    z = z1+h/2
    nb = nx*nz

    xx = np.array(numpy.matlib.repmat(x,nz,1))
    xx = xx.transpose().ravel()
    zz = np.array(numpy.matlib.repmat(z,1,nx))
    zz = zz.transpose().ravel()


    AA = np.zeros((nx,nb))

    for i in range(nx):
        for j in range(nb):
            r1 = ((zz[j]-h/2)**2 + (x[i]-xx[j]+d/2)**2)**0.5
            r2 = ((zz[j]+h/2)**2 + (x[i]-xx[j]+d/2)**2)**0.5
            r3 = ((zz[j]-h/2)**2 + (x[i]-xx[j]-d/2)**2)**0.5
            r4 = ((zz[j]+h/2)**2 + (x[i]-xx[j]-d/2)**2)**0.5
            theta1 = np.arctan((x[i]-xx[j]+d/2)/(zz[j]-h/2))
            theta2 = np.arctan((x[i]-xx[j]+d/2)/(zz[j]+h/2))
            theta3 = np.arctan((x[i]-xx[j]-d/2)/(zz[j]-h/2))
            theta4 = np.arctan((x[i]-xx[j]-d/2)/(zz[j]+h/2))
            AA[i,j] = 2*G*((x[i]-xx[j]+d/2)*np.log((r2*r3)/(r1*r4)) 
                    + d*np.log(r4/r3) - (zz[j]+h/2)*(theta4-theta2) 
                    + (zz[j]-h/2)*(theta3-theta1))

    gm = np.dot(AA,VV)*1e5
    return gm

x=x0 
fx=f
xo=x
fo=fx

xx=[]
ff=[]
dff=[]
TT=[]
for k in range(kmax):
    Ti=1000*((0.99)**k)
    mu=10**(100*((k/kmax)**q))
    # print(mu)
    # fungsi mu_inv
    def mu_inv(y, mu):
        x = np.sign(y) * (1.0 / mu) * ((1.0 + mu)**abs(y) - 1.0)
        return x
    dx = mu_inv(2*np.random.rand(np.size(x))-1,mu)*(u-l)
    x1=x+dx
    x1=(x1<l)*l+(l<= x1)*(x1<=u)*x1+(u< x1)*u
    # fungsi forward
    gm = forward_grav(x1,model)
    fx1 = gm
    Nf=len(fx1)
    df=1/Nf*(np.sum(np.sqrt((fx1-fx)**2)))

    if np.any(df)<0 or np.random.rand()<np.any(np.exp(Ti*df/(abs(fx)+np.spacing(fx))/TolFun)):
        x=x1
        fx=fx
    if np.any(fx)<np.any(fo): 
        xo=x
        fo=fx1
    fx11=np.transpose(fx1)
    f11=np.transpose(f)
    dfo=1/Nf*(sum(np.sqrt((fx11-f11)**2)))
    
    dff.append(df)
    xx.append(x1)
    # ff = fx1
    ff.append(fx1)
    TT.append(Ti)
    # print(e_min)

#  Searching the global solution
ff = np.reshape(ff,(kmax,nxg))
e_min=np.array([min(dff)])
ff_min=np.where(dff==e_min)
ff_sol=ff[ff_min,:].ravel()
xsol=np.reshape(np.transpose([xx]),(nlay,kmax))
xx_sol=xsol[:,ff_min]

VV2=xx_sol.ravel().transpose()
VV=np.reshape(xx_sol,(nzg,nxg))
zSA1 = np.arange(0,(nzg-1)*dh,dh)
zSA = zSA1+10/2

plt.figure('Gambar 1',figsize=(9,7.5))
plt.subplot(2,1,1).set(xlim=(0,90),ylim=(0,0.35))
plt.plot(xxx,f,'ob-',label = 'G-obs')
plt.plot(xxx,ff_sol,'or-',label = 'Inversi SA')
# plt.plot(xxx,gas,'og-',label = 'G-cal')
plt.legend()
plt.title('Penampang Hasil Inversi SA 2D Gravitasi')
plt.xlabel('Spasi [M]')
plt.ylabel('Anomali Medan Gravitasi [mGal]')
plt.subplot(2,1,2)
a = plt.imshow(VV,#extent=(x.min()-5, x.max()+5, z.max()+5, z.min()-5),
           aspect='auto')
plt.ylabel('Kedalaman [M]')
plt.colorbar(a, orientation='horizontal')
plt.clim(0,1200)
plt.grid()
plt.figure('Gambar 2',figsize=(9,7.5))
plt.subplot(2,1,1)
plt.plot(TT,'-b',label='Temperature')
plt.title('Kurva Penurunan Temperatur Dan Misfit Rata-Rata')
plt.xlabel('Iterasi')
plt.ylabel('Temperature')
plt.legend()
plt.subplot(2,1,2)
plt.plot(dff,'-r',label='Misfit')
plt.xlabel('Iterasi')
plt.ylabel('Misfit rata-rata [mGal]')
plt.legend()

# plt.subplot(9,1,1)
# plt.plot(ff[:,0],'-b')
# plt.subplot(9,1,2)
# plt.plot(ff[:,1],'-b')
# plt.subplot(9,1,3)
# plt.plot(ff[:,2],'-b')
# plt.subplot(9,1,4)
# plt.plot(ff[:,3],'-b')
# plt.subplot(9,1,5)
# plt.plot(ff[:,4],'-b')
# plt.subplot(9,1,6)
# plt.plot(ff[:,5],'-b')
# plt.subplot(9,1,7)
# plt.plot(ff[:,6],'-b')
# plt.subplot(9,1,8)
# plt.plot(ff[:,7],'-b')
# plt.subplot(9,1,9)
# plt.plot(ff[:,8],'-b')
# plt.plot(np.arange(0,kmax),ff,'.r')
print(VV)
print(ff_sol)
print(e_min)
print("--- %s seconds ---" % (time.time() - start_time))
plt.show()
# plt.plot(np.array([np.arange(9)]),gm,'or')
# plt.show()
# print(gm)

