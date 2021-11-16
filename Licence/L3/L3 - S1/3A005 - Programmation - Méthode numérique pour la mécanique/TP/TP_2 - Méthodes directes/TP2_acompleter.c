#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 100





void factoriser_LU( float A[MAX][MAX] , float L[MAX][MAX] , float U[MAX][MAX], int n) 
{
	int i,j,k;
	// factorisation à faire
	
}


void resol_trig_inf( float A[MAX][MAX] , float x[MAX], float b[MAX],int n) {
	int i,k;
	
	//la descente à faire
}


void resol_trig_sup( float A[MAX][MAX] , float x[MAX], float b[MAX],int n) {
	int i,k;
	// la remontée à faire
}	

void afficher_vect( float x[MAX],int n) {
	int i;
	for (i=0;i<n;i++) printf("%1f\n",x[i]);
}

void afficher_mat(float A[MAX][MAX] , int n) {
	int i,j;
	
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) {
			printf("%.1f\t",A[i][j]);
		}
		printf("\n");
	}
}

void remplir_vect( float x[MAX],int n) {
	int i;
	for (i=0;i<n;i++) {
		printf("valeur [%d] : ",i+1);
	    if (scanf("%f",&x[i]) !=1) printf("error");
		
	}
		
}

void remplir_mat(float A[MAX][MAX] , int n) {
	int i,j;
	
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) {
		printf("valeur [%d][%d] : ",i+1,j+1);
	    if (scanf("%f",&A[i][j]) !=1) printf("error");
		}
		printf("\n");
	}
}



int main() {
	int n;
	float x[MAX],b[MAX],y[MAX];

	float A[MAX][MAX],L[MAX][MAX],U[MAX][MAX];
	
	printf("Veuillez entrer la valeur de n : ");
	if (scanf("%d",&n) !=1) printf("error");

	printf("Veuillez remplir la matrice A : \n");
	remplir_mat(A , n);
	afficher_mat(A,n);
		
	printf("Veuillez remplir le vecteur b : \n");
	remplir_vect(b,n);
	afficher_vect(b,n);	
	
	/*
	factoriser_LU(A,L,U,n);
	printf("Matrice L : \n");
	afficher_mat(L , n);
	printf("Matrice U : \n");
	afficher_mat(U , n);
	
	printf("Solution de Ly=b: \n");
	resol_trig_inf( L,y,b ,n);
	afficher_vect(y,n);	
	
	
	printf("Solution de Ux=y: \n");
	resol_trig_sup( U,x,y ,n);
	afficher_vect(x,n);
	
	*/
	



	return 0;
}
