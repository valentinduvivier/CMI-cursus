// Construction d'outil aidant à la résolution numérique de A et B

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 100 					// Nombre de valeur de temps
#define MAX 255					// Taille des tableaux
#define pi 3.1415
#define g 9.81

/*
	Jeudi 21 Novembre : difficulté utilisation différences finies; Problématique : comment faire pour utiliser f(i+1) et f(i-1) aux limites.
							--> Solution : définir sois même les valeurs de Theta[-1] et Theta


    Dimanche 24 Novembre : difficulté sur la méthode LU. De plus, ne ps avoir d'idée visuelle des courbes à trouver ralentit l'étude.
                            --> Solution : redéfinition des bornes pour les boucles afin de rendre le nombre d'itération clair. Résultat : test de LLt un autre jour

    Jeudi 28 Novembre : difficulté sur l'application de la formule de l'intégrale. On revient sur LLt


	- Faire de l'allocation dynamique
	- Voir le lien entre erreur et méthodes LU - Cholesky
	- Voir interprétation valeur de n : plus précis pour n qui diminue
*/

void Afficher_Vecteur(float x[MAX], int n);
void Afficher_Matrice(float A[MAX][MAX], int n);

void Decomposition_LU(float A[MAX][MAX], float L[MAX][MAX], float U[MAX][MAX], int n);
void Decomposition_LLt(float A[MAX][MAX], float L[MAX][MAX], int n);

void Transpose_Matrice(float At[MAX][MAX], float A[MAX][MAX], int n);

void Resolution_triangulaire_inf(float A[MAX][MAX], float y[MAX], float b[MAX], int n);
void Resolution_triangulaire_sup(float A[MAX][MAX], float x[MAX], float y[MAX], int n);

// int f(float h, int k, float Theta_ini);


int main()
{
        float T[MAX] = {0}, Theta[MAX] = {0};		// Initialisation tableaux de temps et de Théta = sin(T*n)
        float Pas_de_Temps, p = 0;

        int k, n;


printf("---------------------------------- Partie 1 : --------------------------------------- \n\n");

    // 1.1
    printf("1.1\n\n");

		//a
		printf("a\n\n");

        Pas_de_Temps = (2*pi)/100;

		printf(" Creation vecteur Temps : \n\n");

        T[0] = 0;       // Initialisation première valeur de Temps pour qu'elle soit bien à 0

		for(k = 0; k <= N; k++)
		{
			T[k+1] = T[k] + Pas_de_Temps;
		}

//      Afficher_Vecteur(T, N+1);

		printf(" Vecteur Temps cree\n");

		printf("\n\n");


		//b
		printf("b\n\n");

		printf(" Creation vecteur Theta : \n\n");

		n = 10;		// Nombre étudiant = 3700499

		for(k = 0; k <= N; k++)
		{
			Theta[k] = sin(n*T[k]);
		}

//      Afficher_Vecteur(Theta, N+1);

        // A generaliser
/*
        for(k = 0; k < N + 1; k++)
		{
			printf("Theta[10*%f] = %f\n", T[k], Theta[k]); // On affiche de telle sorte qu'on met le théta en evidence, permettant de reconnaitre directement si une erreur dans le sinus est survenu
		}
*/
		printf(" Vecteur Theta cree\n");


		printf("\n\n --------------------------------------------------------------------------------------------- \n\n");



    // 1.2
    printf("1.2 - Calcul des derivees\n\n");

		float Derivee_1ere[MAX] = {0}, Derivee_2nde[MAX] = {0}, Verification[MAX]; // Verification permet de s'assurer que les dérivées ce sont bien effectuée
		float h = pi/N; // Il faut h<1 ??

		printf("Calcul des derivees premiere et seconde de f(t) = sin(nt) : \n\n");


	/* A FAIRE :
		- justifier pourquoi ces approximations sont bonnes ;
		- justifier que notre h est adapté au calcul ;

		--> l'approximation sur théta (erreur du à approximation sur pi) se répéercute sur la dérivée et la dérivée seconde
				--> lien avec choix de h pour minimiser cela
	*/

	// Différence progressive - dérivée 1ere : f(xi) = [f(xi+h) - f(xi)]/h  --> précision de h


	// Utilisation de gnuplot pour tracer les dérivées 1ere et 2nde
		printf("Ecriture des derivees dans 2 fichiers textes : \n\n");

		FILE *fichier = NULL;
		FILE *fichier_2 = NULL;

    	fichier = fopen("Derivee_1ere.txt", "w");
    	fichier_2 = fopen("Derivee_2nde.txt", "w");


	// On définit nous même Theta(N+1) (utilisation pour Derivée_1ere et Dérivée_2nde)
		Theta[N + 1] = sin(n*(T[N] + Pas_de_Temps));		// On ne peut pas faire sin(nT[N+1]) car la valeur la plus haute de k = N car k < N +1

		for(k = 0; k <= N; k++)
		{
			Derivee_1ere[k] = (Theta[k+1] - Theta[k])/h;

			fprintf(fichier, "%f %f\n", T[k], Derivee_1ere[k]);
		}

	// Calcul de cosinus pour vérifier nos résultats : dsin(nt) = ncos(nt)  --> le but est d'avoir un idée visuelle de l'écart entre valeur numérique et valeur analytique

		for(k = 0; k <= N; k++)
		{
			Verification[k] = n*cos(n*T[k]);		// On doit avoir Verification[k] = Derivee_1ere[k]
		}

/*
		for(k = 0; k <= N; k++)
		{
			printf("Verification = %f  et  Derivee_1ere[%f]  = %f\n", Verification[k], Theta[k], Derivee_1ere[k]);
        }

*/

		printf("\nDerivee_1ere calculee\n\n");


		// Différence centrée - dérivée 2nde : f''(xi) = [f(xi-1) - 2*f(xi) + f(xi+1)]/2h  --> précision de h²

        // On calcule nous même Theta[-1] afin de l'utiliser pour la dérivée seconde
		Theta[-1] = -sin(Pas_de_Temps);

		for(k = 0; k <= N; k++)
		{
			Derivee_2nde[k] = (Theta[k-1] - 2*Theta[k] + Theta[k+1])/(pow(h, 2));

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
    printf("1.3 - Calcul integrale\n\n");

        float Integrale = 0, I_analytique = 0, F[MAX] = {0};
        float Erreur_1_3 = 0;

        printf("Calcul de l'integrale - valeur analytique : \n\n");

        //Utilisation de la méthode des trapèzes pour le calcul de l'intégrale entre [0;pi] de Theta

        for(k = 0; k <= N-1; k++)
        {
            F[k] = sin(k*h);
        }

        F[N] = 0;   // On donne la dernière valeur pour F car la boucle précédente ne calcule que jusqu'à F[N-1]
        for(k = 0; k <= N-1; k++)
        {
            Integrale += (h/2)*(F[k] + F[k+1]);
        }

        printf("Integrale = %f\n\n", Integrale);


        // Calcul de la précision - comparasion solution analytique/numérique
        Erreur_1_3 = fabs(Integrale - I_analytique); 	//Resultat analytique : Integrale = 0 donc pas besoin de la mettre

        printf("\nL'erreur sur le calcul est : Erreur Q.3 = %f\n\n", Erreur_1_3);	// Voir pourquoi l'erreur est la suivante

        printf("\n\n --------------------------------------------------------------------------------------------- \n\n");

    //1.4
    printf("1.4 - Calcul du polynome associe Gm(x)\n\n");

        float Matrice_x[MAX][MAX] = {0}, Vect_x[MAX] = {0}, Lt[MAX][MAX] = {0};
        float L[MAX][MAX] = {0}, U[MAX][MAX] = {0}, y[MAX] = {0}, x[MAX] = {0};     // x étant le vecteur que l'on cherche à résoudre dans cette question

        float S1, S2, Gmx[MAX] = {0};
        float Erreur_1_4 = 0;

        int m = 6;      // numéro étudiant = 3700091
        int i, j;

/*
            Application de la méthode des moindres carrés : Gm(xi) = Somme(ai*xi), entre [0;m]
            Dans notre cas Gm(x) = a0 + a1x + a2x² + a3x³ + a4x^4 + a5x^5 + a6x^6  On cherche ainsi, a0, ... , a6
            Erreur Quadratique = Somme(Gm(x) - Theta)²

            on fait des sommes de 0 à N+1  -->  101 points pour approximer la fonction et déterminer le polynôme d'ordre 6

            Résolution d'un système 7x7
*/

// Systeme a resoudre
    //Matrice
        //Création matrice pour système à résoudre
        for (k = 0; k <= N; k++)
        {
            for (i = 0; i <= m; i++)
            {
                for (j = 0; j <= m; j++)
                {
                    Matrice_x[i][j] += pow(T[k], i+j);
                }
            }
        }


//      Afficher_Matrice(Matrice_x, m);

    //Vecteur
        for (k = 0; k <= N; k++)
        {
            for (i = 0; i <= m; i++)
            {
                Vect_x[i] += Theta[k]*pow(T[k], i);
            }
        }


//        Afficher_Vecteur(Vect_x, m);

        printf("Resolution systeme - calcul a0, ..., a6 :\n\n");

//Cholesky

        printf("Decomposition Cholesky\n\n");

//            Decomposition_LLt(Matrice_x, L, m);

//          Afficher_Matrice(L, m);
            Transpose_Matrice(Lt, L, m);
//          Afficher_Matrice(Lt, m);

        printf("Resolution du systeme : \n");

//           Resolution_triangulaire_inf(L, y, Vect_x, m);
//           Resolution_triangulaire_sup(Lt, x, y, m);
        printf("Voici les coefficients pour a1, ..., a6 obtenus :\n");
//            Afficher_Vecteur(x, m);

/*
//LU
        printf("\n    Decomposition LU : \n\n");
            Decomposition_LU(Matrice_x, L, U, m);

//          Afficher_Matrice(L, m);
//          Afficher_Matrice(U, m);

        printf("    Decomposition effectue\n\n");

        printf("\nResolution du systeme : \n\n");
            Resolution_triangulaire_inf(L, y, Vect_x, m);
            Resolution_triangulaire_sup(U, x, y, m);

//          Afficher_Vecteur(y, m);

        printf("Voici les coefficients obtenus pour a1, ..., a6 :\n\n");
        Afficher_Vecteur(x, m);
*/

        // Affichage de la fonction ainsi trouver Gm(x)
        for(k = 0; k <= N; k++)
        {
            for(i = 0; i <= m; i++)
            {
                Gmx[k] += x[i]*pow(T[k], i);
            }
        }

        // A tracer
/*
         for(k = 0; k <= N; k++)
        {
            printf("Gmx[%d] = %f\n", k, Gmx[k]);
        }
*/
        // Calcul d'erreur
        for(k = 0; k <= N; k++)
        {
            Erreur_1_4 += pow(Theta[k] - Gmx[k], 2);
        }

        printf("\nErreur Q.4 = %0.5f\n", Erreur_1_4);


    printf("\n---------------------------------- Partie 2 : Application au pendule --------------------------------------- \n\n");

// 2.1
    printf("2.1 - Evolution du pendule en fonction du temps\n\n");

/*
    - Voir pour une animation du mouvement
    - Faire des typedef pour les listes de donnée sur pendule
*/

        float Theta_lineaire[MAX] = {0};
        float w0, T0, Theta_ini, Vitesse_ini;

        int Longueur;


        // Dimension Pensule
        Longueur = 5;

        // Conditions initiales
        Theta_ini = pi/6;
        Vitesse_ini = 0;

        w0 = sqrt(g/Longueur);
        T0 = (2*pi)/w0;


    /* --------------------------------------------------------------------------------------------- */

        printf("Trace du Theta en temps reel :\n");

        FILE *fichier_3 = NULL;
        fichier_3 = fopen("Theta_lineaire.txt", "w");

        for(k = 0; k <= N; k++)
        {
            Theta_lineaire[k] = Theta_ini*cos(w0*T[k]);
            fprintf(fichier_3, "%f %f\n", T[k], Theta_lineaire[k]);
        }

        fclose(fichier_3);

        printf("\n\n --------------------------------------------------------------------------------------------- \n\n");

// 2.2
    printf("2.2 - Etude de la periode non-lineaire\n\n");

/*
        Méthode quadrature S(f) = ((b-a)/6)*(f(a) + 4*f((a+b)/2) + f(b));
*/

/*
        floatPeriode_theorique;
        float f_int = 0, f[MAX] = {0}, f_int_2;

        float h_2, a = 0, b = Theta_ini;
        float cst = 4*sqrt(Longueur/(2*g));

        float Erreur_2_2 = 0;

        // Evolution de la période en fonction de Theta_ini

        h_2 = (b-a)/N;

        for (k = 0; k <= N-2; k++)
        {
            f[k] = 1/sqrt(cos(a + k*h_2) - cos(Theta_ini));
        }

        for (k = 0; k <= N-2; k++)
        {
            f_int += (h_2/2)*(f[k] + f[k+1]);
        }

        printf("\nL'integrale de la fonction f est : %f\n\n", cst*f_int);

        Periode_theorique = cst*f_int; // Integrale de référence pour theta_ini = pi/6;

        E = 0.02*Periode_theorique; // On veut avoir une bien meilleur précision sur la valeur de la période T

        // Adaptation du pas d'intégration
        do
        {
            for(k = 0; k <= N; k++)
            {
                f_int_2 += (h_2/2)*(f[k] + f[k+1]);
            }

            Erreur_2_2 = fabs(Periode_theorique - T);

        }while(Erreur_2_2 > E);


        on fait une boucle pour i qui augmente avec h0 donné à la fin de la boucle on augmente ce h;
*/


        return 0;
}




















/*
int f(float h, int k, float Theta_ini)
{

    return (1/pow(cos(k*h) - cos(Theta_ini), 0.5));
}
*/






void Afficher_Vecteur(float x[MAX], int n)
{
    int i;

    for(i = 0; i <= n; i++)
    {
        printf("X[%d] = %f\n", i, x[i]);
    }

    printf("\n");
}

void Afficher_Matrice(float A[MAX][MAX], int n)
{
    int i, j;

    for(i = 0; i <= n; i++)
    {
        for(j = 0; j <= n; j++)
        {
            printf("A[%d][%d] = %f\n", i, j, A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void Transpose_Matrice(float At[MAX][MAX], float A[MAX][MAX], int n)
{
    int i, j;

    for (i = 0; i <= n; i++)
    {
        for (j = 0; j <= n; j++)
        {
            At[i][j] = A[j][i];
        }
    }
}


void Decomposition_LU(float A[MAX][MAX], float L[MAX][MAX], float U[MAX][MAX], int n)
{
    int i, j, k;
    float S1, S2;

    for(j = 0; j <= n; j++)
    {
        for(i = 0; i <= j; i++)
        {
            S1 = 0;
            for(k = 0; k <= i-1; k++)
            {
                S1 += L[i][k]*U[k][j];
            }

            U[i][j] = A[i][j] - S1;
        }

        L[j][j] = 1;

        for(i = j+1; i <= n; i++)
        {
            S2 = 0;
            for(k = 0; k <= j-1; k++)
            {
                S2 += L[i][k]*U[k][j];
            }
            L[i][j] = (A[i][j] - S2)/U[j][j];
        }
    }
}

void Decomposition_LLt(float A[MAX][MAX], float L[MAX][MAX], int n)
{
    float S1, S2;
    int i, j, k;

    for(j = 0; j <= n; j++)
    {
        S1 = 0;
        for(k = 0; k <= j-1; k++)
        {
            S1 += pow(L[j][k],2);
        }

        L[j][j] = sqrt(A[j][j] - S1);

        for(i = j; i <= n; i++)
        {
            S2 = 0;
            for(k = 0; k <= j-1; k++)
            {
                S2 += L[i][k]*L[j][k];
            }

            L[i][j] = (A[i][j] - S2)/L[j][j];
        }
	}
}


void Resolution_triangulaire_inf(float A[MAX][MAX], float y[MAX], float b[MAX], int n)
{
    float S1;

    int i, k;

    for(i = 0; i <= n; i++)
    {
        S1 = 0;
        for(k = 0; k <= i-1; k++)
        {
            S1 += A[i][k]*y[k];
        }

        y[i] = (b[i] - S1)/A[i][i];
    }
}

void Resolution_triangulaire_sup(float A[MAX][MAX], float x[MAX], float y[MAX], int n)
{
    float S1;

    int i, k;

    for(i = n; i >= 0; i--)
    {
        S1 = 0;
        for(k = i+1; k <= n; k++)
        {
            S1 += A[i][k]*x[k];
        }

        x[i] = (y[i] - S1)/A[i][i];
    }
}
