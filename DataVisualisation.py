import numpy as np
import matplotlib.pyplot as plt
import os

k=0
path = '/media/acharlet/Data/Arthur/Documents/Cours/4A CRAL/Code/Shared/Hydro/Resultats'
os.chdir(path)
maxrho=minrho=maxP=minP=maxv=minv=0.
for filename in os.listdir(path):
    if filename.endswith(".dat"):
        with open(filename) as f:
            rho=[]
            P=[]
            v=[]
            for line in f:
                l=line.split()
                l=[float(i) for i in l]
                rho.append(l[0])
                P.append(l[1])
                v.append(l[2])
                maxrho=max(maxrho,np.amax(rho))
                minrho=min(minrho,np.amin(rho))
                maxP=max(maxP,np.amax(P))
                minP=min(minP,np.amin(P))
                maxv=max(maxv,np.amax(v))
                minv=min(minv,np.amin(v))
max1=max(maxP,maxrho)
min1=min(minP,minrho)
for filename in os.listdir(path):
    if filename.endswith(".dat"):
        with open(filename) as f:
            rho=[]
            P=[]
            v=[]
            for line in f:
                l=line.split()
                l=[float(i) for i in l]
                rho.append(l[0])
                P.append(l[1])
                v.append(l[2])
            fig,ax1 = plt.subplots()
            x=np.arange(0,1,0.1)
            ax1.plot(x,rho,'r',label='density')
            ax1.plot(x,P,'b',label='pressure')
            ax1.set_xlabel('x')
            ax1.set_ylabel('rho, P')
            ax1.set_ylim(1.15*min1,1.15*max1)
            ax2=ax1.twinx()
            ax2.plot(x,v,'g',label='speed')
            ax2.set_ylabel('v',color='g')
            ax2.tick_params('y',colors='g')
            #ax2.set_ylim(1.2*minv,1.2*maxv)
            plt.savefig('Step{0}.png'.format(k))
            k+=1

