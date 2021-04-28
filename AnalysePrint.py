import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm 

import scipy.optimize as scp
import scipy.integrate as scpi

from BreakDown import boundary, update, probability

def f(x,a,b):
    return a*x+b



#Trouve la dimension fractale à partir des coordonnées des points de l'agrégat
def fractalDim(coord):
    rmax = max(np.sqrt(coord[:,0]**2 + coord[:,1]**2))
    reven = rmax*0.35
    r = np.exp(np.linspace(np.log(10), np.log(reven), 40))
    m = np.arange(len(r))
    for i in range(len(r)):
        m[i] = len(np.where(coord[:,0]**2 + coord[:,1]**2 < r[i]**2)[0])
    return r, m




def plotFractalDim(coord):
    plt.figure()
    r, m = fractalDim(coord)
    x,y=scp.curve_fit(f,np.log(r),np.log(m), (1,1))
    X,Y=x
    string = "D = " + str(round(X,4)) + " ± " + str(round(np.sqrt(y[0][0]),4))
    plt.loglog(r, np.exp(f(np.log(r),X,Y)), linestyle = "--", color = "grey", label = string)
    plt.loglog(r, m, "*")
    plt.legend()
    print(string)
    plt.xlim(10)
    plt.xlabel(r"$l$")
    plt.ylabel(r"$M$")
    plt.grid(linestyle = "--")
    
    

#Convertit les coordonnées cartésiennes en coord polaires
def getAngleAndRadius(coord):
    radius = np.sqrt(coord[:,0]**2+coord[:,1]**2)
    angle = np.arctan2(coord[:,1], coord[:,0])

    radius = np.reshape(radius, (len(radius),1))
    angle = np.reshape(angle, (len(angle),1))

    return np.concatenate((angle, radius), axis = 1)
    
    

#Convertit l'information la donnée par la matrice sur laquelle a poussé l'agrégat en coordonnées de points
def TabToCoord(tab):
    coord = np.zeros((int(np.sum(np.where(tab != 0, 1, tab))),2))
    a = len(tab)
    center = int((a-1)/2)
    
    for i in range(a):
        for j in range(a):
            if tab[i][j] > 0:
                coord[int(tab[i][j]-1)] = [i - center, j - center] #Trie par ordre d'arrivé


    #Décale les coordonnées si necessaire afin que l'arbre soit centré.
    if coord[0,0] != 0:
        coord[:,0] += -coord[0,0]
    if coord[0,1] != 0:
        coord[:,1] += -coord[0,1]   
    
    return coord 

#Affiche l'arbre pour du offLattice
def plotCircles(tab, cmap, nbParticles):
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)
    ax.patch.set_facecolor('black')
    
    for i in range(nbParticles):
    	circle=plt.Circle((tab[i,0],tab[i,1]),1,color=cmap(i/nbParticles))
    	ax.add_artist(circle)
    

    
    
    xmax=max(np.abs(np.max(tab[:,0])),np.abs(np.min(tab[:,0])))
    ymax=max(np.abs(np.max(tab[:,1])),np.abs(np.min(tab[:,1])))
    plt.xlim(-xmax-10,xmax+10)
    plt.ylim(-ymax-10,ymax+10)
    plt.xticks([],[])
    plt.yticks([],[])
    fig.subplots_adjust(bottom = 0)
    fig.subplots_adjust(top = 1)
    fig.subplots_adjust(right = 1)
    fig.subplots_adjust(left = 0)

    
    
def plotAngleAndHist(polarCoord, bins = 200, multipleSize = 3, hist = False):
    
    nbParticles = len(polarCoord)
    

    plt.figure()
    
    Z = np.linspace(0, nbParticles-1, multipleSize+1, dtype = int) 
    for i in range(1,multipleSize+1):
        """
        angleTemp, radiusTemp = sortListWithAnother(polarCoord[:i,0], polarCoord[:i,1])
        radiusMeanTemp = dataSmoothing(angleTemp, radiusTemp)  
        plt.plot(angleTemp, normData(angleTemp, radiusMeanTemp), label = "avec "+str(i)+" particules")"""
        
        if hist:
            plt.hist(polarCoord[Z[i-1]:Z[i],0], bins = bins, density = True, alpha = 0.5, label = "avec "+str(i)+" particules")
        else:
            numberOfData, x = np.histogram(polarCoord[Z[i-1]:Z[i],0], bins = bins)
            x = (x[1:]+x[:-1])/2
            if i == 1:
                string = "1er" + " 1/"+str(multipleSize)
            elif i == 2:
                string = "2nd" + " 1/"+str(multipleSize)
            else:
                string = str(i)+"ème" + " 1/"+str(multipleSize)
            plt.plot(x, normData(x, numberOfData), label = string, linestyle = "-", alpha = 0.8)  #Affiche la proportion de particules arrivant à un angle donné pour une certain portion de la simulation
        plt.xlabel("Angle")
        plt.ylabel("Proportion de particules")
        plt.ylim(0,1)

    plt.legend()
            
            
    #plt.plot(angleMid, normData(angleMid, radiusMid))

        

    plt.xlim(-np.pi-np.pi/bins, np.pi+np.pi/bins)
    

    plt.figure()
    
    for i in np.linspace(nbParticles/multipleSize, nbParticles,multipleSize, dtype = int):
        angleSorted, radiusSorted = sortListWithAnother(polarCoord[:i,0], polarCoord[:i,1])
        angleMid, radiusMid = maxRadForAngle(angleSorted, radiusSorted)
        #plt.polar(angleSorted, radiusSorted,'.')
        plt.polar(np.append(angleMid, angleMid[0]), np.append(radiusMid,radiusMid[0]),"--") #Affiche les bords de l'agrégat pour un nombre de particules croissant
    



#Trie la première liste dans l'ordre de la seconde liste
def sortListWithAnother(listToBeSorted, otherList):
    order = listToBeSorted.argsort()
    listToBeSorted = listToBeSorted[order]
    otherList = otherList[order]
    return listToBeSorted, otherList
    

#Trouve grossièrement les bords de l'agrégat
def maxRadForAngle(angleSorted, radiusSorted):
    division = 100
    step = int(len(angleSorted)/division)
    angleMid = np.zeros(division)
    radiusMax = np.zeros(division)
    
    
    
    for i in range(division):  #Comme je ne suis pas sur que je ne vais pas depasser les bornes de la liste, je prefère utiliser des % et donc ne pas user de i*step:(i+1)*step ou max par exemple
        temp = -1
        for j in range(step):
            angleMid[i] += angleSorted[(i*step+j)%len(angleSorted)]
            if radiusSorted[(i*step+j)%len(angleSorted)]>temp:
                temp = radiusSorted[(i*step+j)%len(angleSorted)]
        radiusMax[i] = temp 
        
    angleMid = angleMid/step

    
    return angleMid, radiusMax


def normData(x,y):
    return y/scpi.trapz(y, x)


#Supprime plus ou moins le bruit sur des données oscillantes
def dataSmoothing(angleSorted, radiusSorted):
    
    step = 0.05
    radiusMean = np.array(radiusSorted)
    steps = np.arange(1-step,0,-step)

    
    
    for i in range(len(steps)): #moyenne
        radiusMean += steps[i]*(np.roll(radiusSorted,i+1)+np.roll(radiusSorted,-i-1))

    
    
    #return radiusSorted
    return radiusMean


#Créer une matrice contenant les probabilités qu'à une particule de se coller pour chaque coordonnée.
def probaFromCoord(probaCoord, size, allEqual = False):
    proba = np.zeros((size, size))
    for i in range(len(probaCoord)):
        proba[probaCoord[i,0], probaCoord[i,1]] += 1
        
    if allEqual:  
        return np.where(proba != 0, 1, proba)
    else:
        return proba/np.sum(proba)
  


#Reshape le premier tableau en fonction de conditions sur le deuxième.
def reshape(tabToReshape, tabToInitShape):
    
    A = np.where(tabToInitShape > 0)
    center = len(tabToInitShape)/2
    
    
    xmin = min(A[0])
    xmax = max(A[0])
    ymin = min(A[1])
    ymax = max(A[1])
    
    xlen = xmax - xmin
    ylen = ymax - ymin
    
    
    finalLen = max([xlen,ylen])*1.5
    
    
    tabReshaped = tabToReshape[int(center-finalLen/2):int(center+finalLen/2),int(center-finalLen/2):int(center+finalLen/2)]
    
    return tabReshaped



#Etablit les probabilités de se coller pour un arbre donné
def probaFromBreakDownModel(coord, size, iterations = 1000):
    
    coordSized = coord + np.array([int((size-1)/2), int((size-1)/2)])
    potential = np.ones((size, size))
    boundary(potential, coordSized)
    
    print("Début du calcul théorique des probabilités, cela va prendre un certain temps (maudit soit mon calcul catastrophique de l'équation de Laplace)")
    for j in range(iterations):
        if j%int(iterations/10)==0:
            print(j,"/", iterations)
        temp = update(potential)
        potential = (1-0.5)*potential+0.5*temp
        boundary(potential, coordSized)
    
    
    
    x=np.linspace(0, 1, size)
    y=np.linspace(0, 1, size)
    X,Y=np.meshgrid(x,y)    
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X,Y,potential)
    
    proba, location = probability(potential, coordSized) #Trouve la probabilité à partir de la forme de l'agrégat et du potentiel trouvé par l'équation de Laplace
      
    
    tabOfProba = np.zeros((size, size))
    for i in range(len(proba)):
        tabOfProba[location[i,0], location[i,1]] = proba[i]
    
    
    
    return tabOfProba/np.sum(tabOfProba)

   
   





tab = np.loadtxt("tab.txt")
a, b = np.shape(tab)
cmap = matplotlib.cm.get_cmap('cool')


if a==b: #OffLattice
    
    
    coord = TabToCoord(tab).astype(int) #On trouve les coordonnées de chaque particule plutot (c'est plus facile de gerer les coordonnées que la matrice NxN entière )
    tabReshaped = reshape(tab, tab)
    
    
    polarCoord = getAngleAndRadius(coord)
    
    probaCoord = np.loadtxt("proba.txt").astype(int)
    if len(probaCoord) != 0: #Regarde si il y a des données sur la distribution de probabilité expérimentale de l'agrégat
        IsThereAProbaDistrib = True
        proba = probaFromCoord(probaCoord, len(tab), False)
        probaReshaped = reshape(proba, tab)
        
        probaTheo = probaFromBreakDownModel(coord, len(tabReshaped))       
        
    else:
        
        IsThereAProbaDistrib = False
        print("L'erreur sur le fichier vide est normale")

    plotAngleAndHist(polarCoord)    
    plotFractalDim(coord)
    
    
    cmap.set_bad(color='black', alpha = 1)  
    tabReshaped = np.ma.masked_where(tabReshaped == 0, tabReshaped)
    plt.imsave('OnLattice.png', tabReshaped, cmap = cmap, dpi = 300)
    

    if IsThereAProbaDistrib:
        alphaProba = np.dstack([np.zeros((len(probaReshaped),len(probaReshaped),3), dtype=np.uint8), (probaReshaped/np.max(probaReshaped)).T])
        alphaProbaTheo = np.dstack([np.zeros((len(probaTheo),len(probaTheo),3), dtype=np.uint8), (probaTheo/np.max(probaTheo)).T])
        alphaProbaSqrt = np.dstack([np.zeros((len(probaReshaped),len(probaReshaped),3), dtype=np.uint8), np.sqrt(probaReshaped/np.max(probaReshaped)).T])
        alphaProbaTheoSqrt = np.dstack([np.zeros((len(probaTheo),len(probaTheo),3), dtype=np.uint8), np.sqrt(probaTheo/np.max(probaTheo)).T])
        
        plt.figure(figsize=(10,10))
        plt.subplot(2,2,1)
        plt.imshow(tabReshaped.T, interpolation = "none", cmap = cmap, alpha = 0.2, origin ="lower")
        plt.imshow(alphaProba, interpolation = "none", origin = "lower")
        plt.ylabel("Probabilité expérimentale")
        plt.xlabel("Probabilité et agrégat")
        plt.subplot(2,2,2)
        plt.imshow(alphaProba, interpolation="none", origin ="lower")
        plt.xlabel("Probabilité")
        plt.subplot(2,2,3)
        plt.imshow(tabReshaped.T, interpolation = "none", cmap = cmap, alpha = 0.2, origin ="lower")
        plt.imshow(alphaProbaTheo, interpolation = "none", origin = "lower")
        plt.ylabel("Probabilité Théorique")
        plt.xlabel("Probabilité et agrégat")
        plt.subplot(2,2,4)
        plt.imshow(alphaProbaTheo, interpolation="none", origin ="lower")
        plt.xlabel("Probabilité")
  
    plt.figure()
    plt.imshow(tabReshaped.T, interpolation="none", cmap = cmap, origin ="lower")
        
else: #OnLattice
    
    if int(tab[0,0]) == 0:
        nbParticles = len(tab)
        coord = tab
        
    else:
        nbParticles = int(tab[0,0])
        tab[0,:] = 0
        coord = tab[:nbParticles,:]
    
    polarCoord = getAngleAndRadius(coord)
    
    plotFractalDim(coord)
    plotAngleAndHist(polarCoord)
    
    print("Affichage de la structure, cela peut prendre un certain temps")
    plotCircles(coord, cmap, nbParticles)
    
    plt.savefig('OffLattice.png', dpi=300)
    
    
    #plt.savefig('OffLattice.pdf')
    
    plt.show()
    
