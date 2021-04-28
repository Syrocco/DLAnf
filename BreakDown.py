import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

import matplotlib.cm 


#Calcul le (laplacien + tab)/4
def update(tab):
    return 0.25*(np.roll(tab, 1, axis = 0)+np.roll(tab, -1, axis = 0)+np.roll(tab, 1, axis = 1)+np.roll(tab, -1, axis = 1))


#Définit les boundary conditions
def boundary(tab, particles):
    mid = int(len(tab)/2)
    theta = np.linspace(0,2*np.pi,10000)
    x = ((mid-3)*np.cos(theta) + mid).astype(int)
    y = ((mid-3)*np.sin(theta) + mid).astype(int)
    for i in range(len(x)):
            tab[x[i],y[i]] = 1 #Boundary conditions: proba sur un cercle = 1
    
    tab[:,0] = 1
    tab[:,-1] = 1
    tab[-1,:] = 1
    tab[0,:] = 1 #Met les bords à 1 aussi (ce n'est pas vraiment utile mais...)
    
    for i in range(len(particles)):
        tab[particles[i,0], particles[i,1]] = 0 #Met la proba à 0 sur les particules deja présentes

"""
#Boundary pour la géneratio d'un éclair, c'est assez rigolo.    
def boundary(tab, particles):
    
    tab[:,0] = 0
    tab[:,-1] = 1
    tab[-1,:] = 0
    tab[0,:] = 0
    
    for i in range(len(particles)):
        tab[particles[i,0], particles[i,1]] = 0
"""    



def probability(tab, particles):
    possibleLocation = []
    sides = np.array([[-1,0], [1,0], [0,1], [0,-1], [1,1], [-1,1], [-1,-1], [1,-1]])
    for i in range(len(particles)): #Cherche un emplacement possible pour la nouvelle particule
        for j in range(8):
            if not any(np.equal(particles,particles[i]+sides[j]).all(1)):
                possibleLocation.append(particles[i] + sides[j])
    possibleLocation = np.unique(np.array(possibleLocation), axis = 0)      
    nbOfCandidates = len(possibleLocation)
    proba = np.zeros(nbOfCandidates)
    for i in range(nbOfCandidates):
        proba[i] = tab[possibleLocation[i,0], possibleLocation[i,1]]
    proba = proba/np.sum(proba)
    return proba, possibleLocation #Donne les emplacements possibles pour la nouvelle particule et les probas associées à ces emplacements



def newParticle(tab, particles):
    proba, possibleLocation = probability(tab, particles)
    for i in range(1, len(proba)):
        proba[i] += proba[i-1]
    r = np.random.random()
    particles = np.vstack((particles, possibleLocation[np.argmax(proba > r)])) #Ajoute une particule avec une proba donnée par probability()
    return particles 

if __name__ == "__main__":
    
    N = 150
    size = 401
    particles = np.array([[int(size/2), int(size/2)], [int(size/2), int(size/2)+1]])
    #particles = np.array([[int(size/2), 1], [int(size/2), 2]]) #Conditions initiales pour l'éclair
    tab = np.ones((size, size))
    M = 100
    
    boundary(tab, particles)
    
    for j in range(M):
        print(j)
        for i in range(N): #Fait sur-relaxer le potentiel pour trouver Laplacien(tab) = 0
            temp = update(tab)
            tab = (1-0.5)*tab+0.5*temp
            boundary(tab, particles) 
        particles = newParticle(tab, particles) #Rajoute une particule à chaque itération
        
    
    x=np.linspace(0, 1, size)
    y=np.linspace(0, 1, size)
    X,Y=np.meshgrid(x,y)    
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X,Y,tab)
    
    plt.figure()
    
    final = np.zeros((size,size)) #Matrice représentant l'espace contenant les positions des particules
    
    
    for i in range(M):
        final[particles[i,0], particles[i,1]] = i+1
    
    
    cmap = matplotlib.cm.get_cmap('cool')    
    cmap.set_bad(color='white',alpha=1)
    final = np.ma.masked_where(final == 0, final)
    plt.imshow(final.T, cmap = cmap)
