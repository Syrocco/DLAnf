import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm 

cmap = matplotlib.cm.get_cmap('cool')



tab=np.loadtxt("tab.txt")



nbParticles=len(tab)


fig, ax = plt.subplots()
fig.set_size_inches(12, 12)
ax.patch.set_facecolor('black')


for i in range(nbParticles):
	circle=plt.Circle((tab[i,0],tab[i,1]),1,color=cmap(i/nbParticles))
	ax.add_artist(circle)


xmax=max(np.abs(np.max(tab[:,0])),np.abs(np.min(tab[:,0])))
ymax=max(np.abs(np.max(tab[:,1])),np.abs(np.min(tab[:,1])))
plt.xlim(-xmax-10,xmax+10)
plt.ylim(-ymax-10,ymax+10)
plt.show()

plt.savefig('test3.png', edgecolor="black",dpi=500)

