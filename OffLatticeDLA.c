#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#define method 1
#include "mersenne.h"


typedef struct cp gridInfo;
struct cp{
	int number;
	int *tab;	
};

int testNearby(float **array, float position[], gridInfo gridBox);
void savePointToGrid(float x, float y, int index, int i, int j, gridInfo **grid, int boxSize, int boxNumber);
void coordToTXT(int arraySize, float **array, int info);
float** boundaryInitialisation(int arraySize);
gridInfo** gridInitialisation(int boxNumber, int gridSize, int boxSize);
void decompactify(float x1, float y1, float position[]);
void randomPointOnCircle(float rad, float position[]);
float evolution(float **array, float maxRad, int step, gridInfo** grid, int boxSize, int boxNumber, int gridSize);
int iteration(float **array, int nbIterations, gridInfo** grid, int boxSize, int boxNumber, int gridSize);
void coordToIndex(int *i, int *j, float x, float y, int gridSize, int boxSize);
void addPoint(gridInfo *structure, int index);


int main(int argc, char *argv[]){
	
	int nbParticles = 10000; 
	
	
	int boxNumber = 400; //boxNumber^2 = nombre de partitions de l'espace
	int boxSize = 6; //longueur d'une partition
	int gridSize = boxNumber*boxSize;
	
	
	gridInfo **grid = gridInitialisation(boxNumber, gridSize, boxSize); //partition
    float **tab = boundaryInitialisation(nbParticles); //Coordonnées des particules
	
	int infoEnd = iteration(tab, nbParticles, grid, boxSize, boxNumber, gridSize);
	
	if (infoEnd){
		printf("Pas totalement fini, mais pas d'erreur\n");
	}
	else{
		printf("Fini\n");
	}
	
	coordToTXT(nbParticles, tab, infoEnd);
	

	return 0;
	
}



int iteration(float **array, int nbIterations, gridInfo** grid, int boxSize, int boxNumber, int gridSize){
	float radiusParameter = 0;
	float radiusBuffer;
	//clock_t begin = clock();
	for (int i=1;i<nbIterations;i++){ 
		
		if (i%(int) ((nbIterations)/15)==0){
			//clock_t end = clock();
			printf("%d/%d\n",i,nbIterations);
			//printf("%lf, ", (double)(end - begin) / CLOCKS_PER_SEC);
		}
		
		radiusBuffer=evolution(array, radiusParameter+5, i, grid, boxSize, boxNumber, gridSize); //stock la distance au centre de la particule qui vient d'etre ajoutée à l'arbre
		if (radiusBuffer>radiusParameter){
			radiusParameter=radiusBuffer; //Enregistre la distance au centre de la particule la + éloignée
			if (radiusParameter + 30 > gridSize/2){
				printf("%f", radiusParameter);
				return i;
			}
		}
	}

	return 0;
}



//Fait diffuser la particule
float evolution(float **array, float maxRad, int step, gridInfo** grid, int boxSize, int boxNumber, int gridSize){
	float position[2];
	float move[2];
	float rad = maxRad;
	int i,j;
	randomPointOnCircle(maxRad, position);
	while (1){
		
		if (rad>maxRad){
			randomPointOnCircle(rad-maxRad,move); //la sainte optimisation :'u. Permet de faire de grands sauts si la particule est loin de l'agrégat
		}
		else{
			randomPointOnCircle(1,move);
		}
		
		position[0] += move[0];
		position[1] += move[1];
		
		rad=sqrt(position[0]*position[0]+position[1]*position[1]);
		
		

		if (rad > (maxRad + 5)*2){ //Tue la particule si elle est trop loin
			randomPointOnCircle(maxRad, position);
			rad = maxRad;
		}
		
		
		
		if (rad < (maxRad + 5)){ //Checke les collisions uniquement si nécessaire

			
			coordToIndex(&i, &j,  position[0],  position[1],  gridSize,  boxSize);

			
			if (testNearby(array, position, grid[i][j])){ //Teste si la particule entre en collision avec l'agrégat
				array[step][0] = position[0];
				array[step][1] = position[1];
				coordToIndex(&i, &j,  position[0],  position[1],  gridSize,  boxSize); //pour mettre le decompactifier
				savePointToGrid(position[0], position[1], step, i, j, grid, boxSize, boxNumber);

				return rad; //ce n'est pas le bon radius comme j'ai legerement deplacé la particule, mais il ne changera pas de bcp, alors cela reste correct
			}
		}
	}		
}



void decompactify(float x1, float y1, float position[]){ //Déplace la particule pour qu'elle ne soit pas dans celle qu'elle vient de coller.
	
	float x2=position[0];
	float y2=position[1];
	float distance=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	float angle=fabs(asin((x1-x2)/distance));
	
	if ((x1<=x2) && (y1>=y2)){
		position[0]+=(2-distance)*sin(angle);
		position[1]-=(2-distance)*cos(angle);
	}
	
	if ((x1>=x2) && (y1>=y2)){
		position[0]-=(2-distance)*sin(angle);
		position[1]-=(2-distance)*cos(angle);
	}
	
	if ((x1<=x2) && (y1<=y2)){
		position[0]+=(2-distance)*sin(angle);
		position[1]+=(2-distance)*cos(angle);
	}
	
	if ((x1>=x2) && (y1<=y2)){
		position[0]-=(2-distance)*sin(angle);
		position[1]+=(2-distance)*cos(angle);
	}

}



//Initialisation de la partition
gridInfo** gridInitialisation(int boxNumber, int gridSize, int boxSize){
	
	gridInfo** array = malloc(boxNumber*sizeof(gridInfo *));
	for (int i = 0; i < boxNumber; i++){
		array[i] = malloc(boxNumber*sizeof(gridInfo));
	}	
	
	for (int i = 0; i < boxNumber; i++){
		for (int j = 0; j < boxNumber; j++){
			array[i][j].number = 0; //Aucune particule au début
			array[i][j].tab = NULL; //pointeur de longueur nulle (puisqu'il n'y a rien à stocker au début)
		}
	}
	
	int i,j;
	float x = 0.001, y = 0.001;
	
	coordToIndex(&i, &j,  x,  y,  gridSize,  boxSize);  //la fonction SavePointToGrid necessiterait plus de condition pour attribuer correctement les pts x = y = 0, comme ils n'arrivent jamais, je préfère faire une exception ici que modifier savePointToGrid
	savePointToGrid(y, x, 1, i, j, array, boxSize, boxNumber);
	
	return array;
}


//Teste si il y a une collision (à partir des particules se trouvant dans la partition)
int testNearby(float **array, float position[], gridInfo gridBox){
	
	int a = gridBox.number;
	for (int i = 0; i < a; i++){
		if ((position[0] - array[(gridBox.tab)[i]][0])*(position[0] - array[(gridBox.tab)[i]][0]) + (position[1]-array[(gridBox.tab)[i]][1])*(position[1] - array[(gridBox.tab)[i]][1]) <= 4){
			decompactify(array[(gridBox.tab)[i]][0], array[(gridBox.tab)[i]][1], position);
			return 1;
		}
	}
	return 0;
	
}



//Enregistre une particule dans la partition
void savePointToGrid(float x, float y, int index, int i, int j, gridInfo **grid, int boxSize, int boxNumber){

	addPoint(&grid[i][j], index);
	int boundary[2] = {0,0};
	
	if ( ((fabs(fmodf(x, boxSize)) < 1) && (x>=0)) || ((fabs(fmodf(x, boxSize)) > boxSize-1) && (i > 0) && (x<=0)) ){

		addPoint(&grid[i-1][j], index);
		boundary[0] = -1;
		
	}
	if ( ((fabs(fmodf(x, boxSize)) > boxSize-1) && (i < boxNumber-1) && (x>=0)) || ((fabs(fmodf(x, boxSize)) < 1) && (x<=0)) ){
		addPoint(&grid[i+1][j], index);
		boundary[0] = 1;
		
	}
	if ( ((fabs(fmodf(y, boxSize)) < 1) && (y>=0)) || ((fabs(fmodf(y, boxSize)) > boxSize-1) && (j<boxNumber-1) && (y<=0)) ){
		addPoint(&grid[i][j+1], index);
		boundary[1] = 1;
		
	}
	if ( ((fabs(fmodf(y, boxSize)) > boxSize-1) && (j > 0) && (y>=0)) || ((fabs(fmodf(y, boxSize)) < 1) && (y<=0)) ){
		addPoint(&grid[i][j-1], index);
		boundary[1] = -1;
	
	}

	// Le test est trop large et inclue des particules sur les diagonales qu'il ne faudrait pas inclure. Mais c'est pas grave.
	if (boundary[0] != 0 && boundary[1] != 0 ){ 
		addPoint(&grid[i+boundary[0]][j+boundary[1]], index);
	}


}


void addPoint(gridInfo *structure, int index){
	structure->number += 1; //Change le nombre de points
	structure->tab = realloc(structure->tab, structure->number*sizeof(int));  //augmente dynamiquement la mémoire allouée à tab qui contient les index des points (fort heuresement fonctionne comme malloc si (*structure).tab=NULL).
	(structure->tab)[structure->number - 1] = index; //stocke dans le tableau l'index de la particule ajoutée

}



//Donne les coordonnées de la partition correspondant à un (x, y) donné
void coordToIndex(int *i, int *j, float x, float y, int gridSize, int boxSize){
	*i = (gridSize/2 + x)/boxSize;
	*j = (-y + gridSize/2)/boxSize;
}



//Initialise le tableau de coordonnées
float** boundaryInitialisation(int arraySize){
	
	float **array = (float **)malloc(arraySize * sizeof(float *)); 
    for (int i=0; i<arraySize; i++) 
         array[i] = (float *)calloc(2, sizeof(float));
	
	return array;
}


//Sélectionne un point aléatoire sur un cercle de rayon rad
void randomPointOnCircle(float rad, float position[]){
	float theta = (float)(2*M_PI*genrand64_real1());
	position[0]=rad*cos(theta);
	position[1]=rad*sin(theta);
}


//Sauvegarde le tableau de coordonnées
void coordToTXT(int arraySize, float **array, int info){
	FILE *fichier;
	fichier=fopen("tab.txt","w");
 	
	fprintf(fichier,"%d %d\n", info, info);
	
	for(int i=1; i<arraySize; i++){
		for (int j=0; j<2; j++){
			fprintf(fichier,"%f ",array[i][j]);
		}
		fprintf(fichier,"\n");
    }
	fclose(fichier);
}


