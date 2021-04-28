#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mersenne.h"

#define SQRT_2_PI 2.5066282746310002

void randomShift(int shift[]);
void boundaryInitialisation(int arraySize, int** array);
void particleInitialisation(int arraySize, int position[], float maxRad);
void randomXandYFromRadius(float rad, int middle, int position[]);
float evolution(int arraySize, int **array, float maxRad, int step, int nbIterations, FILE* proba);
int TestNearby(int **array, int x, int y);
float lenght(int posx, int posy);
int iteration(int arraySize, int **array, int nbIterations, int probaNumber);
void tabToTXT(int arraySize, int **array);
void progressBar(int step, int size);



int main(int argc, char *argv[]){
	srand(time(NULL));
	

	int nbParticles = 2000; 
	int N = 1000; //Taille de la grille
	int nbIterationsForProba = 20000; //Nombre de particules servant à établir la densité de probabilité
	//scanf("%d",&nbParticles);
	

    int **tab = (int **)malloc(N * sizeof(int *)); 
    for (int i=0; i<N; i++) 
         tab[i] = (int *)malloc(N * sizeof(int)); 
	 

	
	if (iteration(N,tab,nbParticles, nbIterationsForProba)){
		printf("Tout s'est bien passe\n");
	}
	else
	{
		printf("L'abre a touche les bords\n");
	}


	tabToTXT(N,tab);
	
	return 0;
	
}



int iteration(int arraySize, int **array, int nbIterations, int nbIterationsForProba){
	float radiusParameter = 0;
	float radiusBuffer;
	FILE *proba;
	proba = fopen("proba.txt","w");
	
	//clock_t begin = clock();
	
	boundaryInitialisation(arraySize,array);
	int ctr = 0;
	for (int i=2;i<nbIterations+2+nbIterationsForProba;i++){ //commence à 2 pour avoir une chronologie correcte
		if (i%(int) (nbIterations/15)==0){
			if (ctr != 14){
				ctr++;
			}
			//clock_t end = clock();
			progressBar(ctr, 15);
			//printf("%lf, ", (double)(end - begin) / CLOCKS_PER_SEC);
			if (i>nbIterations+2)
				printf("\t\t\t\t Phase de calcul de la probabilite en cours...");
		}
		
		radiusBuffer=evolution(arraySize,array,radiusParameter+5,i, nbIterations, proba); //stock la distance au centre de la particule qui vient d'etre ajoutée à l'arbre
		if (radiusBuffer>radiusParameter){
			radiusParameter=radiusBuffer; //Enregistre la distance au centre de la particule la + éloignée
			if (radiusParameter>arraySize/2-15){ //Teste si la particule peut potentiellement être ajoutée en dehors du tableau
				fclose(proba);
				return 0;
			}
		
		}
	}
	fclose(proba);
	return 1;
}


//Fait evoluer la particule jusqu'à ce qu'elle en touche une autre
float evolution(int arraySize, int **array, float maxRad, int step, int nbIterations, FILE* proba){
	int startPosition[2];
	int move[2] = {0,0};
	float radiusPos;
	particleInitialisation(arraySize, startPosition, maxRad);  //initialise l'emplacement de départ d'une particule
	int posx = startPosition[0]; 
	int posy = startPosition[1];
	while (1){
		
		randomShift(move);
		posx += move[0];
		posy += move[1];
		
		
		radiusPos=lenght(posx - arraySize/2, posy - arraySize/2);
		
		if (radiusPos > (maxRad + 5)*2 || posx>=arraySize-1 || posx<=1 || posy>=arraySize-1 || posy<=1){ //Vérifie que la particule n'est pas en dehors de la matrice et qu'elle n'est pas trop loin de la structure

			particleInitialisation(arraySize, startPosition, maxRad);
			posx = startPosition[0];
			posy = startPosition[1];
			
		}
		
		if (radiusPos < maxRad + 4){ //Ne teste pas la position si on est en dehors de l'arbre
			if (TestNearby(array, posx, posy)){ //Teste les collisions
				if (step < nbIterations + 2)
					array[posx][posy] = step;  //Ajoute la particule à l'arbre
				else{
					fprintf(proba, "%d\t%d\n", posx, posy); //Ajoute une particule à la probabilité
				}
				return radiusPos;
			}
		}
	}
}


int TestNearby(int **array, int x, int y){ //teste si une particule est aux alentours de celle qui se balade: 8 cases
	if (array[x+1][y+1]>0 || array[x+1][y]>0 || array[x+1][y-1]>0 || array[x-1][y]>0 || array[x-1][y+1]>0 || array[x-1][y-1]>0 || array[x][y-1]>0 || array[x][y+1]>0)
		return 1;
	else
		return 0;
}

/*
int TestNearby(int **array, int x, int y){ //teste si une particule est aux alentours de celle qui se balade: 4 cases
	if (array[x-1][y]>0 || array[x+1][y]>0 || array[x][y-1]>0 || array[x][y+1]>0)
		return 1;
	else
		return 0;
}
*/




// selectionne un déplacement sur la gauche/bas/haut ou droite
void randomShift(int shift[]){
	int randomNum = (int)(4*genrand64_real1());
	switch (randomNum){
		case 0:
		{
			shift[0] = 0;
			shift[1] = 1;
			break;
		}

		case 1:
		{
			shift[0] = 0;
			shift[1] = -1;
			break;
		}

		case 2:
		{
			shift[0] = 1;
			shift[1] = 0;
			break;
		}

		case 3:
		{
			shift[0] = -1;
			shift[1] = 0;
			break;
		}
	}
}



//Initialise les particules fixées au début de la simulation
void boundaryInitialisation(int arraySize, int **array){
	for (int i = 0; i<arraySize; i++){
		for (int j = 0; j<arraySize; j++){
			array[i][j] = 0;
		}
	}
	array[(arraySize-1)/2][(arraySize-1)/2] = 1; //plante une graine au milieu
}



//Iinitialise la position d'une (et une seule) particule que l'on va faire diffuser
void particleInitialisation(int arraySize, int position[], float maxRad){
	int middle = (arraySize - 1)/2;
	randomXandYFromRadius(maxRad + 5, middle, position);   //initialise une particule aléatoirement, proche des autres.
}




float lenght(int posx, int posy){
	return sqrt(posx*posx+posy*posy);
	
}

//Convertit (r,theta) -> (x,y)
void randomXandYFromRadius(float rad, int middle, int position[]){ //prend une position random sur un cercle de rayon rad et de centre middle
	float theta = ((float)rand()/(float)(RAND_MAX))*2*M_PI;
	int x = floor(rad*cos(theta) + middle);
	int y = floor(rad*sin(theta) + middle);
	position[0] = x;
	position[1] = y;
	
}


void tabToTXT(int arraySize, int **array){
	FILE *fichier;
	fichier = fopen("tab.txt","w");
 
	
	for(int i = 0; i < arraySize; i++){
		for (int j = 0; j < arraySize; j++){
			fprintf(fichier, "%d ", array[i][j]);
		}
		fprintf(fichier, "\n");
    }
	fclose(fichier);
}



void progressBar(int step, int size){
	printf("\033[A\r["); // \033[A monte à la ligne du dessus et \r l'efface, permet de mettre un \n à l'avant dernier print et donc d'avoir un retour à la ligne à la fin du programme... Sous Windows ça ne marche pas.
	
	for (int i = 0; i < size; i++){
		if (i <= step){
			printf("#");
		}
		else{ 
			printf(" ");		
		}
	}
	printf("] - %.2f %%   \n", (float) step/(size-1)*100);
	fflush(stdout); //sinon printf utilise un buffer :/
		
}





