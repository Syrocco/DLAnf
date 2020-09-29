import numpy as np
import matplotlib.pyplot as plt

tab=np.loadtxt("tab.txt").astype(int)
taille=tab[-1]
array=np.zeros((taille,taille))


for i in range(taille):
    for j in range(taille):
        array[i,j]=tab[(i*taille)+j]
        

    

cmap=plt.cm.cool
cmap.set_bad(color='black',alpha=1)
array = np.ma.masked_where(array>=0, array)
plt.imsave('test3.png', array, cmap=cmap,dpi=500)
plt.imshow(array,interpolation="none", cmap=cmap)

