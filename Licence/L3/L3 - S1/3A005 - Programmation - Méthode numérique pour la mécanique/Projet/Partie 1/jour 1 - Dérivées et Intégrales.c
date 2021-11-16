// Construction d'outil aidant à la résolution numérique de A et B

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 100 					// Nombre de valeur de temps
#define MAX 255					// Taille des tableaux
#define pi 3.1415

/*
	Jeudi 21 Novembre : difficulté utilisation différences finies; Problématique : comment faire pour utiliser f(i+1) et f(i-1) aux limites.
							--> Solution : définir sois même les valeurs de Theta[-1] et Theta

	Faire de l'allocation dynamique

*/

int main()
{
	float T[MAX] = {0}, Theta[MAX] = {0};			// Initialisation tableaux
	float Pas_de_Temps, T0 = 0, p = 0;
	int k, n;


	Pas_de_Temps = (2*pi)/100;

	T[0] = 0;

// Voir comment utiliser une boucle do
/*
	do
	{
		T[k+1] += Pas_de_Temps;

		k++;

		printf("Pour k = %d, T = %f", k, T[k]);

	}while(k >= 100);
*/

// 1.1

		//a

		printf("\nCreation vecteur Temps : \n\n");

		for(k = 1; k < N + 1; k++)
		{
			T[k] = T0 + Pas_de_Temps;		// Voir comment faire pour T[k] +=
			T0 = T[k];
		}

	//A generaliser
/*
		for(k = 0; k < N + 1; k++)
		{

			printf("T[%d] = %f\n", k, T[k]);
		}
*/

		printf("Vecteur Temps cree\n");

		printf("\n\n");



		//b

		printf("Creation vecteur Theta : \n\n");

		n = 10;		// Nombre étudiant = 3700499

		for(k = 0; k < N + 1; k++)
		{
			Theta[k] = sin(n*T[k]);
		}

	// A generaliser
/*
		for(k = 0; k < N + 1; k++)
		{
			printf("Theta[%f] = %f\n", 2*T[k], Theta[k]); // On affiche de telle sorte qu'on met le théta en evidence, permettant de reconnaitre directement si une erreur dans le sinus est survenu
		}
*/

		printf("Vecteur Theta cree\n");


		printf("\n\n --------------------------------------------------------------------------------------------- \n\n");



// 1.2

		float Derivee_1ere[MAX] = {0}, Derivee_2nde[MAX] = {0}, Verification[MAX]; // Verification permet de s'assurer que les dérivées ce sont bien effectuée
		float h = Pas_de_Temps; // Il faut h<1 ??

		printf("Calcul des derivees premiere et seconde de f(t) = sin(nt) : \n\n");


	/* A FAIRE :
		- justifier pourquoi ces approximations sont bonnes ;
		- justifier que notre h est adapté au calcul ;

		--> l'approximation sur théta (erreur du à approximation sur pi) se répéercute sur la dérivée et la dérivée seconde
				--> lien avec choix de h pour minimiser cela
	*/



	// Différence progressive - dérivée 1ere : f(xi) = [f(xi+h) - f(xi)]/h  --> précision de h

	// Calcul de cosinus pour vérifier nos résultats : dsin(nt) = ncos(nt).

/*
		for(k = 0; k < N + 1; k++)
		{
			Verification[k] = n*cos(n*T[k]);		// On doit avoir Verification[k] = Derivee_1ere[k]
		}
*/

	// Utilisation de gnuplot pour tracer les dérivées 1ere et 2nde
		printf("Ecriture des derivees dans 2 fichiers textes : \n\n");

		FILE *fichier = NULL;
		FILE *fichier_2 = NULL;

    	fichier = fopen("Derivee_1ere.txt", "w");
    	fichier_2 = fopen("Derivee_2nde.txt", "w");


	// On définit nous même Theta(N +2) :  utilisation pour Derivée_1ere et Dérivée_2nde
		Theta[N + 1] = sin(n*(T[N] + Pas_de_Temps));		// On ne peut pas faire sin(nT[N+1]) car la valeur la plus haute de k = N car k < N +1

		for(k = 0; k < N + 1; k++)
		{
			Derivee_1ere[k] = (Theta[k+1] - Theta[k])/h;

			fprintf(fichier, "%f %f\n", T[k], Derivee_1ere[k]);
		}
/*
		for(k = 0; k < N + 1; k++)
		{
			printf("Verification = %f  et  Derivee_1ere[%f]  = %f\n", Verification[k], Theta[k], Derivee_1ere[k]);
		}
*/

		printf("\nDerivee_1ere calculee\n\n");


		// Différence centrée - dérivée 2nde : f''(xi) = [f(xi-1) - 2*f(xi) + f(xi+1)]/2h  --> précision de h²

	// On calcule nous même Theta[-1] afin de l'utiliser pour la dérivée seconde
		Theta[-1] = -sin(Pas_de_Temps);

		for(k = 0; k < N + 1; k++)
		{
			Derivee_2nde[k] = (Theta[k-1] - 2*Theta[k] + Theta[k+1])/(Pas_de_Temps*Pas_de_Temps);

			fprintf(fichier_2, "%f %f\n", T[k], Derivee_2nde[k]);

			// printf("Theta[%f] = %f\n", k + h; Theta[k+h])
		}


		fclose(fichier);
		fclose(fichier_2);

/*
		for(k = 0; k < N + 1; k++)
		{
			printf("Derivee_2nde[%f] = %f\n", Theta[k], Derivee_2nde[k]);
		}
*/

		printf("Derivee_2nde calculee\n\n");

		printf("\nFichiers textes crees\n");


		printf("\n\n --------------------------------------------------------------------------------------------- \n\n");

//1.3

	float Integrale, I0 = 0, I_analytique = 0;
	float Erreur, s1 = 0;

	printf("Calcul de l'integrale - valeur analytique : \n\n");

	//Utilisation de la méthode des trapèzes pour le calcul de l'intégrale entre [0;pi] de Theta

/*
	for(k = 0; k < N + 1; k++)
	{
		Integrale[k] = I0 + (3.1415/2) * (Theta[k+1] - Theta[k]);
		I0 = Integrale[k];
	}


// Voir cette équation et son utilité

	for(k = 0; k < N + 1; k++)
	{
		Integrale += (pi/2) * (Theta[k+1] - Theta[k]);	// On garde la même précision sur pi que pour le calcul de Theta
	}
*/

	for(k = 0; k < N; k++)
    {
        s1 += Theta[k];     // f(a) = f(b) = 0
    }

    Integrale = (pi/n)*s1;

	printf("Integrale = %f\n\n", Integrale);


	// Calcul de la précision - comparasion solution analytique/numérique
	Erreur = fabs(Integrale - I_analytique); 	//Resultat analytique : Integrale = 0 donc pas besoin de la mettre

	printf("L'erreur sur le calcul est : Erreur = %f\n\n", Erreur);	// Voir pourquoi l'erreur est la suivante

	printf("\n\n --------------------------------------------------------------------------------------------- \n\n");

//1.4

	float Matrice_somme_x[MAX][MAX] = {0}, Somme_x[MAX][MAX] = {0};
	int m = 6, i, j;	// numéro étudiant = 3700091

/*	Application de la méthode des moindres carrés : Gm(xi) = Somme(ai*xi), entre [0;m]
	Dans notre cas Gm(x) = a0 + a1x + a2x² + a3x³ + a4⁴ + a5x⁵  On cherche ainsi, a0, ... , a5
	Erreur Quadratique = Somme(Gm(x) - Theta)²

	Résolution d'un système 5x5
*/

	//Création matrice pour système à résoudre

    for (k = 0; k < N + 1; k++)
	{
        for (i = 0; i < m + 1; i++)
        {
            for (j = 0; j < m + 1; j++)
            {
                Somme_x[i][j] += pow(Theta[k], i+j);
                printf("Matrice_somme_x[%d][%d] = %f\n\n", i, j, Matrice_somme_x[i][j]);
            }
        }
	}

// Affichage de la Matrice
    for (i = 0; i < m + 1; i++)
    {
        for (j = 0; j < m + 1; j++)
        {
            Matrice_somme_x[i][j] += pow(Theta[k], i+j);
            printf("Matrice_somme_x[%d][%d] = %f\n\n", i, j, Matrice_somme_x[i][j]);
        }
    }


	return 0;
}
