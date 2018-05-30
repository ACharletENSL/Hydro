import numpy as np
import matplotlib.pyplot as plt
import os
Ncell=1000; Nvar=3
dx = 1./Ncell
x = np.arange(0,1,dx)
k=0
path = '/home/acharl01/Work/Shared/Resultats2D'
os.chdir(path)
maxrho=minrho=maxP=minP=maxv=minv=0.
#for filename in os.listdir(path):
#    if os.path.isfile(os.path.join(path,filename)) and 'data' in filename:
#        with open(filename) as f:
#            rho=[]
#            P=[]
#            v=[]
#            for line in f:
#                l=line.split()
#                l=[float(i) for i in l]
#                rho.append(l[0])
#                P.append(l[1])
#                v.append(l[2])
#                maxrho=max(maxrho,np.amax(rho))
#                minrho=min(minrho,np.amin(rho))
#                maxP=max(maxP,np.amax(P))
#                minP=min(minP,np.amin(P))
#                maxv=max(maxv,np.amax(v))
#                minv=min(minv,np.amin(v))
#max1=max(maxP,maxrho)
#min1=min(minP,minrho)

for filename in sorted(os.listdir(path)):
    if os.path.isfile(os.path.join(path,filename)) and 'data' in filename:
        with open(filename) as f:
            rho=[]
            P=[]
            v=[]
            for line in f:
                l=line.split()
                l=[float(i) for i in l]
                rho.append(l[0])
                P.append(l[2])
                v.append(l[1])
            fig,ax1 = plt.subplots()
            ax1.plot(x,rho,'r:',label='Density')
            ax1.plot(x,P,'b:',label='Pressure')
            ax1.set_xlabel('x')
            ax1.set_ylabel('rho, P',rotation=0)
            ax1.set_xlim(0,1)
            ax1.set_ylim(0,5)
            ax2=ax1.twinx()
            ax2.plot(x,v,'g:',label='Speed')
            ax2.set_ylabel('v',color='g',rotation=0)
            ax2.tick_params('y',colors='g')
            ax2.set_ylim(-0.2,0.5)
            lines, labels = ax1.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax2.legend(lines + lines2, labels + labels2)
            plt.savefig('Step{0}.png'.format(k))
            plt.close()
            k+=1

