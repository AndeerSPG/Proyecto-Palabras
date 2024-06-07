/***************************************************
 AC - OpenMP -- SERIE
 fun_s.c
 rutinas que se utilizan en el modulo grupopal_s.c
****************************************************/
#include <math.h>
#include <stdlib.h>
#include "defineg.h"
#include <float.h>
#include <stdio.h>
 // definiciones

/*******************************************************************
 1 - Funcion para calcular la distancia euclidea entre dos vectores
 Entrada: 2 elementos con NDIM caracteristicas (por referencia)
 Salida:  distancia (double)
********************************************************************/
double gendist (float *vec1, float *vec2){
	// PARA COMPLETAR
	// calcular la distancia euclidea entre dos vectores
	double dist = 0.0;
	for(int i = 0; i< NDIM; i++){
		dist += pow(vec1[i] - vec2[i],2);
	}
	dist = sqrt(dist);
	return dist;
}

/***********************************************************************************
 2 - Funcion para calcular el grupo (cluster) mas cercano (centroide mas cercano)
 Entrada: nvec  numero de vectores, int
          mvec  vectores, una matriz de tamanno MAXV x NDIM, por referencia
          cent  centroides, una matriz de tamanno ngrupos x NDIM, por referencia
 Salida:  popul grupo mas cercano a cada elemento, vector de tamanno MAXV, por ref.
************************************************************************************/
void grupo_cercano (int nvec, float mvec[][NDIM], float cent[][NDIM],int *popul){
	// PARA COMPLETAR
	// popul: grupo mas cercano a cada elemento
	double distmenor,distancia;
	int imenor;
	for(int j=0;j<nvec; j++){
		distmenor = FLT_MAX;
		imenor = 0;
		for(int i = 0; i<ngrupos; i++){
			distancia = gendist(mvec[j],cent[i]);
			if(distancia < distmenor){
				distmenor = distancia;
				imenor = i;
			}
		}
		popul[j]=imenor;
	}

}

/***************************************************************************************
 3 - Funcion para calcular la calidad de la particion de clusteres.
     Ratio entre a y b.
     El termino a corresponde a la distancia intra-cluster.
     El termino b corresponde a la distancia inter-cluster.
 Entrada: mvec    vectores, una matriz de tamanno MAXV x NDIM, por referencia
          listag  vector de ngrupos structs (informacion de grupos generados), por ref.
          cent    centroides, una matriz de tamanno ngrupos x NDIM, por referencia
 Salida:  valor del CVI (double): calidad/bondad de la particion de clusters
****************************************************************************************/
double silhouette_simple(float mvec[][NDIM], struct lista_grupos *listag, float cent[][NDIM],
		float a[]){
   float b[ngrupos];

    // PARA COMPLETAR

    // aproximar a[i] de cada cluster: calcular la densidad de los grupos;
    //		media de las distancia entre todos los elementos del grupo;
    //   	si el numero de elementos del grupo es 0 o 1, densidad = 0
    // ...
    for(int i =0; i< ngrupos ;i++){
	 float sum = 0.0;
	 for (int k=0; k< listag[i].nvecg;k++){
		 for(int j = 0; j < listag[i].nvecg; j++){
			float dist = (float) gendist(mvec[listag[i].vecg[k]],mvec[listag[i].vecg[j]] );
			sum = sum + dist;

		 }
	}
	sum = sum /(float)pow(listag[i].nvecg,2);
	if(listag[i].nvecg < 2){
		sum = 0.0;
		a[i]= 0.0;
	}

	else{
		a[i] = sum;
	}
    }
    // aproximar b[i] de cada cluster
    //...teniendo en cuenta que la distancia con el mismo es 0
    for(int j = 0; j< ngrupos; j++){
    	 float sumc = 0.0;
   	 for( int i = 0; i< ngrupos ; i++){
		float distc =(float)gendist(cent[j],cent[i]);
		sumc = sumc + distc;

	 }
	 sumc = sumc / (float) ngrupos;
    	 b[j]= sumc;
    }

    // calcular el ratio s[i] de cada cluster
    // ...
    float s[ngrupos];
    for (int i=0; i<ngrupos; i++){
	  s[i] = (b[i]-a[i])/  fmax(b[i],a[i]);
    }

    // promedio y devolver
    // ...
    double media = 0.0;
    for(int i=0; i<ngrupos; i++){
    	media =(double) media + s[i];
    }
    media = media / ngrupos;
    return media;
}

//funcion auxiliar que ordena un array de menor a mayor
void ordenar(float *array,int tamano){
	for(int i = 0; i< tamano;i++){
		for(int j=0;j<tamano-1; j++){
			if(array[j+1] < array[j]){
				float aux = 0.0;
				aux = array[j];
				array[j]= array[j+1];
				array[j+1] = aux;
			}
		}
	}
}
/********************************************************************************************
 4 - Funcion para relizar el analisis de campos UNESCO
 Entrada:  listag   vector de ngrupos structs (informacion de grupos generados), por ref.
           mcam     campos, una matriz de tamaño MAXV x NCAM, por referencia
 Salida:   info_cam vector de NCAM structs (informacion del analisis realizado), por ref.
*****************************************************************************************/
void analisis_campos(struct lista_grupos *listag, float mcam[][NCAM],
		struct analisis *info_cam){
	// PARA COMPLETAR
	// Realizar el analisis de campos UNESCO en los grupos:
	//    mediana maxima y el grupo en el que se da este maximo (para cada campo)
	//    mediana minima y su grupo en el que se da este minimo (para cada campo)


	//accede a cada campo 0-23
	for(int i=0;i<NCAM;i++){
	   float mediana[ngrupos];
	   //accede a cada agrupacion de palabras
	   for(int j=0;j<ngrupos;j++){
		float s[listag[j].nvecg];
		//accede a los indices de las palabras del grupo j
		for(int m=0;m< listag[j].nvecg;m++){
			s[m] = mcam[listag[j].vecg[m]][i];
		}
		//ordena de menor a mayor las probabilidades del campo i
		ordenar(s,listag[j].nvecg);

		//la mediana de las probabilidades del grupo j se almacena en mediana[j].
		if(listag[j].nvecg > 1){
		mediana[j] =s[listag[j].nvecg/2];
		}
		if(listag[j].nvecg == 1){
		mediana[j] = s[0];
		}
		if(listag[j].nvecg == 0){
		mediana[j] = 10.0;
		}
	   }

	   //recorre el array donde estan las medianas del campo i para coger el maximo y el minimo.
	   float max = FLT_MIN;
	   float min = FLT_MAX;
	   int grupomin=0;
	   int grupomax=0;
	   for(int s=0; s<ngrupos;s++){
		if(mediana[s] !=10.0){
			if(mediana[s] > max){
				max = mediana[s];
				grupomax = s;
			}
			if(mediana[s] < min){
				min= mediana[s];
				grupomin=s;
			}
		}
	   }
	   info_cam[i].mmax= max;
	   info_cam[i].mmin=min;
	   info_cam[i].gmax = grupomax;
	   info_cam[i].gmin = grupomin;
	}


}



/*************************************
   OTRAS FUNCIONES DE LA APLICACION
**************************************/
void inicializar_centroides(float cent[][NDIM]){
	int i, j;
	float rand_val;
	srand (147);
	for (i=0; i<ngrupos; i++)
		for (j=0; j<NDIM/2; j++){
			rand_val = ((rand() % 10000) / 10000.0)*2-1;
			cent[i][j] = rand_val;
			cent[i][j+(NDIM/2)] = cent[i][j];
		}
}

int nuevos_centroides(float mvec[][NDIM], float cent[][NDIM], int popul[], int nvec){
	int i, j, fin;
	double discent;
	double additions[ngrupos][NDIM+1];
	float newcent[ngrupos][NDIM];

	for (i=0; i<ngrupos; i++)
		for (j=0; j<NDIM+1; j++)
			additions[i][j] = 0.0;

	// acumular los valores de cada caracteristica; numero de elementos al final
	for (i=0; i<nvec; i++){
		for (j=0; j<NDIM; j++) additions[popul[i]][j] += mvec[i][j];
		additions[popul[i]][NDIM]++;
	}

	// calcular los nuevos centroides y decidir si el proceso ha finalizado o no (en funcion de DELTA)
	fin = 1;
	for (i=0; i<ngrupos; i++){
		if (additions[i][NDIM] > 0) { // ese grupo (cluster) no esta vacio
			// media de cada caracteristica
			for (j=0; j<NDIM; j++)
				newcent[i][j] = (float)(additions[i][j] / additions[i][NDIM]);

			// decidir si el proceso ha finalizado
			discent = gendist (&newcent[i][0], &cent[i][0]);
			if (discent > DELTA1){
				fin = 0;  // en alguna centroide hay cambios; continuar
			}

			// copiar los nuevos centroides
			for (j=0; j<NDIM; j++)
				cent[i][j] = newcent[i][j];
		}
	}
	return fin;
}
