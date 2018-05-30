import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
writer = FFMpegWriter(fps=15, metadata=metadata)

fig = plt.figure()
l, = plt.colormesh(x,y,[],cmap='OrRd')
plt.title('$\rho$')
plt.grid(True, zorder=10)

plt.xlim(0, 49)
plt.ylim(0, 49)

path = '/media/acharlet/Data/Arthur/Documents/Cours/4A_CRAL/Code/Perso/Hydro/Resultats2D'
os.chdir(path)
with writer.saving(fig, "writer_test.mp4", 100):
    for filename in sorted(os.listdir(path)):
        if os.path.isfile(os.path.join(path,filename)) and 'data' in filename:
            with open(filename) as f:
                rho=np.zeros((Ncell,Ncell))
                for line in f:
                    l=line.split()
                    l=[float(z) for z in l]
                    i=int(l[0])-1
                    j=int(l[1])-1
                    rho[i][j]+=l[2]
                    l.set_data(rho)
                    writer.grab_frame()
