import numpy as np
import matplotlib.pyplot as plt
import os
Ncell=15; Nvar=4
dx = 1./Ncell
x = np.arange(0,1,dx); y = x
xi = np.append(x,[1.]); yi = xi
k=0
path = '/media/acharlet/Data/Arthur/Documents/Cours/4A_CRAL/Code/Perso/Hydro/Resultats2D'
os.chdir(path)
for filename in sorted(os.listdir(path)):
    if os.path.isfile(os.path.join(path,filename)) and 'data' in filename:
        with open(filename) as f:
            rho=np.zeros((Ncell,Ncell))
            P=np.zeros((Ncell,Ncell))
            vx=np.zeros((Ncell,Ncell))
            vy=np.zeros((Ncell,Ncell))
            for line in f:
                l=line.split()
                l=[float(z) for z in l]
                i=int(l[0])-1
                j=int(l[1])-1
                rho[i][j]+=l[2]
                vx[i][j]+=l[3]
                vy[i][j]+=l[4]
                P[i][j]+=l[5]
            plt.figure(1)
            plt.subplot(221)
            plt.title('Density')
            plt.pcolormesh(x,y,rho,cmap='OrRd')
            plt.grid(True, zorder=10)
            plt.colorbar()
            plt.subplot(222)
            plt.pcolormesh(x,y,P,cmap='BuPu')
            plt.title('Pressure')
            plt.grid(True, zorder=10)
            plt.colorbar()
            plt.subplot(223)
            plt.streamplot(x,y,vx,vy,density=[0.5, 1])
            plt.title('Speed')
            plt.savefig('Step{0}.png'.format(k))
            plt.close()
            k+=1
k=0
for filename in sorted(os.listdir(path)):
    if os.path.isfile(os.path.join(path,filename)) and 'flux' in filename:
        with open(filename) as f:
            Frhox=np.zeros((Ncell,Ncell))
            Frhoy=np.zeros((Ncell,Ncell))
            FPx=np.zeros((Ncell,Ncell))
            FPy=np.zeros((Ncell,Ncell))
            for line in f:
                l=line.split()
                l=[float(z) for z in l]
                i=int(l[0])-1
                j=int(l[1])-1
                Frhox[i][j]+=l[2]
                Frhoy[i][j]+=l[3]
                FPx[i][j]+=l[4]
                FPy[i][j]+=l[5]
            plt.figure(1)
            plt.streamplot(x,y,Frhox,Frhoy,density=[0.5, 1],color='r')
            plt.title('Density Flux')
            plt.savefig('rFluxStep{0}.png'.format(k))
            plt.close()
            plt.figure(2)
            plt.streamplot(x,y,FPx,FPy,density=[0.5, 1],color='b')
            plt.title('Pressure Flux')
            plt.savefig('PFluxStep{0}.png'.format(k))
            plt.close()
            k+=1

