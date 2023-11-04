# -----------------------------------------------------------------------
#  Program Forward Model Gravitasi 2D
#  Refererensi Last Kubik
#  Difasilitasi oleh PT. Diamond Startup Indonesia
#  Dibuat Oleh Asrin
#  Kendari  2022
# -----------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib as repmat

G=6.6732e-11    # Konstanta gravitasi
dx = 10         # (m)
dh = 10         # (m)
d = dx          # station spacing
h = dh          # depth spacing (thickness)
rho = 1000      # (kg/m^3)
rho2 = 500
nx = 9
nz = 4

# ==================== Model Kotak ==================== 
V = np.zeros((nz,nx))
V[1,3:6] = rho
V[2,3:6] = rho
VV = V.transpose().ravel()
# ==================== Model Cross Homogen ==================== 
# V[1,4] = rho
# V[2,3:6] = rho
# V[3,4] = rho
# ==================== Model Cross Heterogen ==================== 
# V[1,4] = rho
# V[2,3:6] = rho2
# V[3,4] = rho
# ==================== Model Patahan Homogen ==================== 
# V[1:3,2] = rho
# V[2:4,3] = rho
# V[3:5,4] = rho
# ==================== Model Patahan Heterogen ==================== 
# V[1,2] = rho2
# V[2,2] = rho
# V[2,3] = rho2
# V[3,3] = rho
# V[3,4] = rho2
# V[4,4] = rho

x1 = np.arange(0,(nx-1)*dx+dx,dx)
x = x1+d/2
z1 = np.transpose([np.arange(0,((nz+1)*h-10),h)])
z = z1+h/2
nb = nx*nz
xx = np.array(repmat.repmat(x,nz,1))
xx = xx.transpose().ravel()
zz = np.array(repmat.repmat(z,1,nx))
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
GG = np.dot(AA,VV)*1e5

# Visualisasi model penampang gravitasi 2D
plt.figure(figsize=(9,7.5))
plt.subplot(2,1,1)
plt.plot(x,GG,'or-',label = 'G-obs')
plt.legend()
plt.title('Penampang Sintetik 2D Gravitasi')
plt.xlabel('Spasi [M]')
plt.ylabel('Anomali Medan Gravitasi [mGal]')
plt.subplot(2,1,2)
a = plt.imshow(V,extent=(x.min()-5, x.max()+5, z.max()+5, z.min()-5),
           aspect='auto',cmap='viridis')
plt.ylabel('Kedalaman [M]')
plt.colorbar(a, orientation='horizontal')
plt.clim(1,1200)
plt.grid()
plt.show()