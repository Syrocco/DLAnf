#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#define M_PI 3.14159265358979323846



void coordToTXT(int arraySize, float **array);
void boundaryInitialisation(int arraySize, float **array);
void randomPointOnCircle(float rad, float position[]);
int testNearby(float **array, float position[], int step);
float evolution(float **array, float maxRad, int step);
void decompactify(float x1, float y1, float position[]);
void iteration(float **array, int nbIterations);


int main(int argc, char *argv[]){
	srand(time(NULL));
	
	int nbParticles=1000; 
	
	
    float **tab = (float **)malloc(nbParticles * sizeof(float *)); 
    for (int i=0; i<nbParticles; i++) 
         tab[i] = (float *)malloc(2 * sizeof(float));
	 
	 
	 
	iteration(tab, nbParticles);
	coordToTXT(nbParticles, tab);
	
}


void coordToTXT(int arraySize, float **array){
	FILE *fichier;
	fichier=fopen("tab.txt","w");
 
	
	for(int i=0; i<arraySize; i++){
		for (int j=0; j<2; j++){
			fprintf(fichier,"%f ",array[i][j]);
		}
		fprintf(fichier,"\n");
    }
	fclose(fichier);
}


void boundaryInitialisation(int arraySize, float **array){
	
	for (int i = 0; i<arraySize; i++){
		for (int j = 0; j<2; j++){
			array[i][j]=0;
		}
	}
}


void randomPointOnCircle(float rad, float position[]){ //prend une position random sur un cercle de rayon rad et de centre middle
	float theta=((float)rand()/(float)(RAND_MAX))*2*M_PI;
	position[0]=rad*cos(theta);
	position[1]=rad*sin(theta);
	
}


float evolution(float **array, float maxRad, int step){
	float position[2];
	float move[2];
	float rad;
	randomPointOnCircle(maxRad, position);
	while (1){
	
		randomPointOnCircle(1,move);
		position[0]+=move[0];
		position[1]+=move[1];
		
		rad=sqrt(position[0]*position[0]+position[1]*position[1]);
		
		
		
		if (rad>(maxRad+5)*2){
			randomPointOnCircle(maxRad+5, position);
		}
		
		
		
		if (rad<maxRad+4)
			if (testNearby(array, position, step)) {
				array[step][0]=position[0];
				array[step][1]=position[1];
				return rad; //ce n'est pas le bon radius comme j'ai legerement deplacé la particule, mais il ne changera pas de bcp, alors cela reste correct
			}
	}
}

int testNearby(float **array, float position[], int step){
	for (int i=0; i<step; i++){
		if ((position[0]-array[i][0])*(position[0]-array[i][0])+(position[1]-array[i][1])*(position[1]-array[i][1])<=4){
			decompactify(array[i][0], array[i][1], position);
			return 1;
		}
	}
	return 0;

}


void iteration(float **array, int nbIterations){
	float radiusParameter=0;
	float radiusBuffer;
	boundaryInitialisation(nbIterations, array);
	for (int i=1;i<nbIterations;i++){ //commence à 2 pour avoir une chronologie correcte
		if (i%(int) (nbIterations/10)==0){
			printf("%d/%d\n",i,nbIterations);
		}
		
		radiusBuffer=evolution(array,radiusParameter+5,i); //stock la distance au centre de la particule qui vient d'etre ajoutée à l'arbre
		if (radiusBuffer>radiusParameter){
			radiusParameter=radiusBuffer; //Enregistre la distance au centre de la particule la + éloignée
		}
	}
}

void decompactify(float x1, float y1, float position[]){ //La mort de l'inventivité :) :) :) Au moins, ca marche...
	
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





