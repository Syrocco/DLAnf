#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#define M_PI 3.14159265358979323846



void randomShift(int shift[]);
void boundaryInitialisation(int arraySize, int **array);
void particleInitialisation(int arraySize, int position[], float maxRad);
void randomXandYFromRadius(float rad, int middle, int position[]);
void tabPrintTest(int arraySize, int **array);
float evolution(int arraySize, int **array, float maxRad, int step);
int TestNearby(int arraySize, int **array, int x, int y);
float lenght(int posx, int posy);
int iteration(int arraySize, int **array, int nbIterations);
void tabToTXT(int arraySize, int **array);


//gcc -Wl,--stack,4194304 DLA.c -o test.x
//gcc -Ofast -Wl,--stack,30194304 DLA.c -o test.x

int main(int argc, char *argv[]){
	srand(time(NULL));
	
	//int N=2000; //La taille du tableau N * N
	int nbParticles=40000; 
	int N=1500;
	//scanf("%d",&nbParticles);
	

    int **tab = (int **)malloc(N * sizeof(int *)); 
    for (int i=0; i<N; i++) 
         tab[i] = (int *)malloc(N * sizeof(int)); 
	 
	
	//clock_t start, end;
	//double cpu_time_used;
	
	
	//start = clock();
	
	if (iteration(N,tab,nbParticles)){
		printf("Tout s'est bien passe\n");
	}
	else
	{
		printf("L'abre a touche les bords\n");
	}

	//end = clock();
	//cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	//printf("%lf",cpu_time_used);

	tabToTXT(N,tab);
	
	return 0;
	
}





void tabToTXT(int arraySize, int **array){
	FILE *fichier;
	fichier=fopen("tab.txt","w");
 
	
	for(int i=0;i<arraySize;i++){
		for (int j=0;j<arraySize;j++){
			fprintf(fichier,"%d ",array[i][j]);
		}
    }
	fprintf(fichier,"%d ",arraySize);

	fclose(fichier);
}


//parceque j'en ai marre de faire des boucles tout le temps pour debugguer
void tabPrintTest(int arraySize, int **array){
	for (int i = 0; i < arraySize; i++){
    	for (int j = 0; j < arraySize; j++){
    	    printf("%d ", array[i][j]); 
		}
		printf("\n");
	}

}


// selectionne un déplacement sur la gauche/bas/haut ou droite
void randomShift(int shift[]){
	int randomNum;
	randomNum=rand()%4;
	switch (randomNum){
		case 0:
		{
			shift[0]=0;
			shift[1]=1;
			break;
		}

		case 1:
		{
			shift[0]=0;
			shift[1]=-1;
			break;
		}

		case 2:
		{
			shift[0]=1;
			shift[1]=0;
			break;
		}

		case 3:
		{
			shift[0]=-1;
			shift[1]=0;
			break;
		}
	}
}



//Initialise les particules fixées au début de la simulation
void boundaryInitialisation(int arraySize,int **array){
	
	
	
	for (int i = 0; i<arraySize; i++){
		for (int j = 0; j<arraySize; j++){
			array[i][j]=0;
		}
	}
	array[(arraySize-1)/2][(arraySize-1)/2]=1; //plante une graine au milieu
}



//Iinitialise la position d'une (et une seule) particule que l'on va faire diffuser
void particleInitialisation(int arraySize, int position[], float maxRad){
	int middle=(arraySize-1)/2;
	randomXandYFromRadius(maxRad+5,middle,position);   //initialise une particule aléatoirement, proche des autres dans un 
}



void randomXandYFromRadius(float rad, int middle, int position[]){ //prend une position random sur un cercle de rayon rad et de centre middle
	float theta=((float)rand()/(float)(RAND_MAX))*2*M_PI;
	int x=floor(rad*cos(theta)+middle);
	int y=floor(rad*sin(theta)+middle);
	position[0]=x;
	position[1]=y;
	
}


//Fait evoluer la particule jusqu'à ce qu'elle en touche une autre
float evolution(int arraySize, int **array, float maxRad, int step){
	int startPosition[2];
	int move[2];
	float radiusPos;
	particleInitialisation(arraySize, startPosition, maxRad);  //initialise l'emplacement de départ d'une particule
	int posx=startPosition[0]; //stocke ses positions de départ pour reinitialiser la particule si elle va trop loin
	int posy=startPosition[1];
	//int k=0;
	while (1){
		/*
		k++;
		
		
		if (k>10000000 && k<40000000){
		//printf("start pos: %d et %d avec pos: %d et %d et rayon  %f et %f \n",startPosition[0],startPosition[1], posx, posy,radiusPos,maxRad);
		array[posx][posy]=-k;
		}
		
		
		
		if (k>=40000000){
		printf("ok");
		FILE * fp;
		fp = fopen ("tab.txt","w");
	 
		for(int i=0; i < arraySize;i++){
			for (int j=0; j<arraySize;j++){
				fprintf (fp, "%d ",array[i][j]);
		
			}
		
		}
		fprintf (fp, "%d ",arraySize);
		fclose (fp);
		exit(0);
		
		}*/
		randomShift(move);
		posx+=move[0];
		posy+=move[1];
		
		
		radiusPos=lenght(posx-arraySize/2,posy-arraySize/2);
		
		if (radiusPos>(3+maxRad)*2 || posx>=arraySize-1 || posx<=1 || posy>=arraySize-1 || posy<=1){ //marche que pour boundaryType 1 à generaliser après... Quel enfer ! Teste si la particule se trouve encore dans l'array et pas trop loin des autres graines (si elle est deux fois plus loin que la particule la plus eloignée du centre, je reinitialise)

			particleInitialisation(arraySize, startPosition, maxRad);
			posx=startPosition[0];
			posy=startPosition[1];
		}
		
		if (radiusPos<maxRad+4){ //Ne teste pas la position si en dehors de l'arbre*/
			if (TestNearby(arraySize,array,posx,posy)){
				array[posx][posy]=step;  //metre step si on veut avoir la chronologie (première particule = 2, deuxieme particule = 3,...), fixe la position de la particule si elle est proche d'une autre deja fixée.
				return radiusPos;
			}
		}

	}
}


int iteration(int arraySize, int **array, int nbIterations){
	float radiusParameter=0;
	float radiusBuffer;
	
	boundaryInitialisation(arraySize,array);
	for (int i=2;i<nbIterations+2;i++){ //commence à 2 pour avoir une chronologie correcte
		if (i%(int) (nbIterations/100)==0){
			printf("%d/%d\n",i,nbIterations);
		}
		
		radiusBuffer=evolution(arraySize,array,radiusParameter+5,i); //stock la distance au centre de la particule qui vient d'etre ajoutée à l'arbre
		if (radiusBuffer>radiusParameter){
			radiusParameter=radiusBuffer; //Enregistre la distance au centre de la particule la + éloignée
			if (radiusParameter>arraySize/2-15){ //Test si la particule peut potentiellement être ajoutée en dehors du tableau
				
				return 0;
			}
		
		}
	}
	return 1;
}



int TestNearby(int arraySize, int **array, int x, int y){ //test si une particule est aux alentours de celle qui se balade
	if (array[x+1][y+1]>0 || array[x+1][y]>0 || array[x+1][y-1]>0 || array[x-1][y]>0 || array[x-1][y+1]>0 || array[x-1][y-1]>0 || array[x][y-1]>0 || array[x][y+1]>0){
		return 1;
	}
	else
	{
		return 0;
	}
}

float lenght(int posx, int posy){
	return sqrt(posx*posx+posy*posy);
	
}



