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

    Jeudi 28 Novembre : difficulté sur l'application de la formule de l'intégrale. On revient sur LLt  --> on a réalisé la décomposition LLt.

    Dimanche 01 Decembre : difficulté sur l'approximation de l'intégrale de sinus Q 1.3 pour calcul de l'erreur.
                            Il y a une grosse incertitude avec la somme de Taylor
                                --> utilisation pas nécassaire/abandonnée
    Jeudi 05 Décembre : On affiche l'évolution de la période en fonction de théta ini (ou du pas ?)
                        On travail sur eles erreurs et les intégrales :

    Vendredi 06 Décembre : Travail sur la décomposition LU : faire en sorte que ela matrice A soit écraser et remplacé par les valeurs de L et U
                            On affine la précision des résultats de la partie 1 on reviens sur les résultats mtn qu'on a plus de recul sur le sujet.
                            Lien avec Erreurs et Intégrations qui sont les méthodes les plus incertaines dans notre ode à ce jour

    Dimache 08 Décembre : On termine avec les résultats déjà fait en cherchant les dernieres améliorations pour accélerer les calculs et les diminuer par la meme occasion.
                            --> On termine avec les erreurs, les incertitudes et les intégrales. On travail sur l'affichage des erreurs par une évolution graphique.


                          Affichage du polynôme des moindres carrés ssocié à sinus
                          Correction Bornes intégrale sinus Q 1.3
                          Corrections pour des éléments divisants possiblement par zéro.
                            --> On travail sur les valeurs limites :
                                                   - pour le nombre d'itération dans l'optimisation du pas et de Theta_ini
                                                   - Fonction à valeurs interdites comme 1/sqrt(cos...)

	- Faire de l'allocation dynamique
	- Voir le lien entre erreur et méthodes LU - Cholesky
	- Voir interprétation valeur de n : plus précis pour n qui diminue
    - Voir pour faire du contionnement sur Q 1.4, et voir pour les autre Questions
*/

void Afficher_Vecteur(float x[MAX], int n);
void Afficher_Matrice(float A[MAX][MAX], int n);

void Decomposition_LU(float A[MAX][MAX], float L[MAX][MAX], float U[MAX][MAX], int n);
void Decomposition_LU_Ecrase(float A[MAX][MAX], int n);

void Decomposition_LLt(float A[MAX][MAX], float L[MAX][MAX], int n);

// void Transpose_Matrice(float At[MAX][MAX], float A[MAX][MAX], int n);

void Resolution_triangulaire_inf(float A[MAX][MAX], float y[MAX], float b[MAX], int n);
void Resolution_triangulaire_inf_Ecrase(float A[MAX][MAX], float y[MAX], float b[MAX], int n);

void Resolution_triangulaire_sup(float A[MAX][MAX], float x[MAX], float y[MAX], int n);
void Resolution_triangulaire_sup_opti(float A[MAX][MAX], float x[MAX], float y[MAX], int n);

int factorielle(int k);

// int f(float h, int k, float Theta_ini);


int main()
{
        float T[MAX] = {0}, Theta[MAX] = {0};		// Initialisation tableaux de temps et de Théta = sin(T*n)
        float Pas_de_Temps;

        int k, n, i;


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

		float Derivee_1ere[MAX] = {0}, Derivee_2nde[MAX] = {0};
		float Verification_derivee_1ere[MAX], Verification_derivee_2nde[MAX], ERRR[MAX] = {0}; // Verification permet de s'assurer que les dérivées ce sont bien effectuée
		float h = pi/N; // Il faut h<<1 ??

		printf("Calcul des derivees premiere et seconde de f(t) = sin(nt) : \n\n");


	/* A FAIRE :
		- justifier pourquoi ces approximations sont bonnes ;
		- justifier que notre h est adapté au calcul ;

		--> l'approximation sur théta (erreur du à approximation sur pi) se répércute sur la dérivée et la dérivée seconde
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
		Theta[N+1] = sin(n*(T[N] + Pas_de_Temps));		// On ne peut pas faire sin(nT[N+1]) car la valeur la plus haute de k = N car k <= N

		for(k = 0; k <= N; k++)
		{
			Derivee_1ere[k] = (Theta[k+1] - Theta[k])/h;

			fprintf(fichier, "%f %f\n", T[k], Derivee_1ere[k]);
		}

        // Calcul de cosinus pour vérifier nos résultats : dsin(nt) = ncos(nt)  --> le but est d'avoir un idée visuelle de l'écart entre valeur numérique et valeur analytique
		for(k = 0; k <= N; k++)
		{
			Verification_derivee_1ere[k] = n*cos(n*T[k]);		// On doit avoir Verification_derivee_1ere[k] = Derivee_1ere[k]
		}

//      Affichage des valeurs pour comparaison
/*
		for(k = 0; k <= N; k++)
		{
			printf("Verification_derivee_1ere = %f  et  Derivee_1ere[%f]  = %f\n", Verification_derivee_1ere[k], Theta[k], Derivee_1ere[k]);
        }
*/
        // Calcul de l'erreur entre les 2 fonctions. On a ici une erreur qui indique que nous notre résultat sera à 1 près.
        for(k = 0; k <= N; k++)
        {
            if(Verification_derivee_1ere[k] != 0)
            {
                ERRR[k] += fabs(Verification_derivee_1ere[k] - Derivee_1ere[k])/Verification_derivee_1ere[k];
            }
        }


		printf("\nDerivee_1ere calculee\n\n");

        printf("Erreur_deriveee_1ere = %f\n\n", ERRR[N]);

		// Différence centrée - dérivée 2nde : f''(xi) = [f(xi-1) - 2*f(xi) + f(xi+1)]/2h  --> précision jusqu'à h^2

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

        // Calcul de cosinus pour vérifier nos résultats : dsin(nt) = ncos(nt)  --> le but est d'avoir un idée visuelle de l'écart entre valeur numérique et valeur analytique
		for(k = 0; k <= N; k++)
		{
			Verification_derivee_2nde[k] = -(n^2)*sin(n*T[k]);		// On doit avoir Verification_derivee_2nde[k] = Derivee_2nde[k]
		}

		for(k = 0; k <= N; k++)
        {
            ERRR[k] = 0;
        }

        // Calcul de l'erreur entre les 2 fonctions. On a ici une erreur qui indique que nous notre résultat sera à 1 près.
        for(k = 0; k <= N; k++)
        {
            if(Verification_derivee_2nde[k] != 0)
            {
                ERRR[k] += fabs(Verification_derivee_2nde[k] - Derivee_2nde[k])/Verification_derivee_2nde[k];
            }
        }

//      Affichage des valeurs pour comparaison
/*
		for(k = 0; k <= N; k++)
		{
			printf("Verification_derivee_2nde = %f  et  Derivee_2nde[%f]  = %f\n", Verification_derivee_2nde[k], Theta[k], Derivee_2nde[k]);
        }
*/

		printf("\nDerivee_2nde calculee\n\n");

        printf("Erreur_deriveee_2nde = %f\n\n", ERRR[N]);

		printf("\nFichiers textes crees\n");


		printf("\n\n --------------------------------------------------------------------------------------------- \n\n");

//1.3
    printf("1.3 - Calcul de l'integrale\n\n");

        float Integrale_1_3 = 0, F[MAX] = {0};
        float Erreur_1_3 = 0;

        printf("Calcul de l'integrale - valeur analytique : \n\n");

        //Utilisation de la méthode des trapèzes pour le calcul de l'intégrale entre [0;pi] de Theta
/*
    Cours page 2 sur précision :

    La rapidité d’exécution nécessaire pour atteindre ce résultat. De manière générale, toutes les méthodes
    peuvent atteindre de très grandes précisions. Cependant, le temps de calcul augmente avec la précision.
    Ce temps n’augmente pas de la même manière pour toutes les méthodes si bien que certaines s’avèrent plus
    efficaces que d’autres. En particulier, le temps de calcul des méthodes de quadrature est proportionnel
    au nombre de points où la fonction f(x) est évaluée.
*/

        for(k = 0; k <= (N-1)/2; k++) // On s'arrête à PI soit N/2
        {
            F[k] = sin(k*h);
        }

        // On est sensé touvé 1 pour h << 1. Car Integrale = -[(cos(hPI)/h - 1]0,PI
        //F[(N+1)/2] = sin((h*(N+1))/2);
        for(k = 0; k <= (N-2)/2; k++) // (N-2)/2 = n-1, car  pour avoir pi on prend n = 50 --> n-1 = 49 = (N-2)/2
        {
            Integrale_1_3 += (h/2)*(F[k] + F[k+1]);
        }

        printf("Integrale = %f\n\n", Integrale_1_3);


		printf("\n\n --------------------------------------------------------------------------------------------- \n\n");


        // Calcul de la précision - comparasion solution analytique/numérique
         Erreur_1_3 = fabs(Integrale_1_3 - 1); 	//Resultat analytique : Integrale = 0 donc pas besoin de la mettre

        // R(f) = erreur de la méthode des trapèzes = -h^3/12 * f''(a) avec a la borne inférieur et b la borne sup avec b = a+h
/*
        printf("h au cube = %f\n\n", pow(h, 3));
        Erreur_1_3 = (-(pow(h, 3))/12)*(-(n^2)*sin(n*T[0]));
*/
        printf("\nL'erreur sur le calcul est : Erreur Q.3 = %f\n\n", Erreur_1_3);	// Voir pourquoi l'erreur est la suivante

        printf("\n\n --------------------------------------------------------------------------------------------- \n\n");

    //1.4
    printf("1.4 - Calcul du polynome associe Gm(x)\n\n");

/*
        Cette partie, por l'optimisation noua mène à 2 méthodes
            - Cholesky : on a qu'une matrice L déduite de U. On a pas besoin de créer Lt et on a moins de calcul que pour LU
                    --> On a moins d'opération, donc on utilise sans doute moins nos variables qui contiennent des erreurs,
                        d'où le fait qu'il y ait moins d'erreur avec cette méthode

            - LU : On utilise dans notre cas une méthode créant temporairement L et U, pour ensuite les utiliser pour recalculer
                    les composantes de A en écrasant les données précédentes. On a un Matrice de moins que pour Cholesky ??
                        --> Sinon voir comment ne plus avoir L et U dès qu'elle sont inutiles --> Allocation dynamique qui pourrait être très utile


            - IMPORTANT : On peut aussi faire des résolutions par méthode du gradient : fonctionne bien pour A sym definie positive avec cond(A) qui tend vers 1
                            --> Il faut donc évalué si le cond(A) permet d'effectuer cette méthode ou si il y a un moyen d'optimiser cond(A) pour pouvoir faire le gradient
*/

        float Matrice_x[MAX][MAX] = {0}, Vect_x[MAX] = {0};
        float L[MAX][MAX] = {0};

        // A decommenter si on fait LU classique
        // float L[MAX][MAX] = {0}, U[MAX][MAX] = {0};
        float y[MAX] = {0}, x[MAX] = {0};     // x étant le vecteur que l'on cherche à résoudre dans cette question

        float S1, S2, Gmx[MAX] = {0};
        float Erreur_1_4 = 0;

        int m = 6;      // numéro étudiant = 3700091
        int j;      // i et k déjà défini

/*
            Application de la méthode des moindres carrés : Gm(xi) = Somme(ai*xi), entre [0;m]
            Dans notre cas Gm(x) = a0 + a1x + a2x² + a3x³ + a4x^4 + a5x^5 + a6x^6  On cherche ainsi, a0, ... , a6
            Erreur Quadratique = Somme(Gm(x) - Theta)²

            on fait des sommes de 0 à N+1  -->  101 points pour approximer la fonction et déterminer le polynôme d'ordre 6

            Résolution d'un système 7x7

            Cette méthode n'a pas un bon conditionnement : elle n'est pas très fiable au delà de m = 1 ou 2 et ne l'est plus du tout à m = 10
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
/*
        printf("\n    Decomposition Cholesky\n\n");

            Decomposition_LLt(Matrice_x, L, m);

//          Afficher_Matrice(L, m);
            Transpose_Matrice(Lt, L, m);
//          Afficher_Matrice(Lt, m);

        printf("    Decomposition effectue\n\n");

        printf("\nResolution du systeme : \n\n");

        Resolution_triangulaire_inf(L, y, Vect_x, m);
        Resolution_triangulaire_sup(Lt, x, y, m);

        printf("Voici les coefficients pour a1, ..., a6 obtenus :\n");
        Afficher_Vecteur(x, m);
*/

        printf("\n    Decomposition Cholesky\n\n");

        Decomposition_LLt(Matrice_x, L, m);

//          Afficher_Matrice(L, m);
//          Transpose_Matrice(Lt, L, m);
//          Afficher_Matrice(Lt, m);

        printf("    Decomposition effectue\n\n");

        printf("\nResolution du systeme : \n\n");

        Resolution_triangulaire_inf(L, y, Vect_x, m);
        Resolution_triangulaire_sup_opti(L, x, y, m);


        printf("Voici les coefficients pour a1, ..., a6 obtenus :\n");
        Afficher_Vecteur(x, m);

//LU

//        printf("\n    Decomposition LU Ecrase : \n\n");
/*
    2 méthodes possibles :
        - Soit on cré 2 Matrices L et U pour tout le reste d la résolution
        - Soit on cré L et U, puis on écrase les Valeurs de A avec celles de L et U
            Sachant que sur la diagonale on garde A[i][i] = U[i][i] et qu'on pose
            dans la résolution triangulaire inférieure A[i][i] = 1

        Voir si on gagne bien en mémoire système ou si les matrices U et L sont quand même stockées


        La seconde méthode permet de créer les matrices L et U temporairement : on a une seule matrice au lieu de 3
*/

            // A decomenter si on fait la méthode LU  classique
//          Decomposition_LU(Matrice_x, L, U, m);
//            Decomposition_LU_Ecrase(Matrice_x, m);

//          Afficher_Matrice(L, m);
//          Afficher_Matrice(U, m);

//        printf("    Decomposition effectue\n\n");

//        printf("\nResolution du systeme : \n\n");

/*
        Resolution_triangulaire_inf(L, y, Vect_x, m);
        Resolution_triangulaire_sup(U, x, y, m);
*/
/*
        // A decommenter si on fait méthode LU classique
        // Resolution_triangulaire_inf(Matrice_x, y, Vect_x, m);
        Resolution_triangulaire_inf_Ecrase(Matrice_x, y, Vect_x, m);
        Resolution_triangulaire_sup(Matrice_x, x, y, m);

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

/*
        FILE *fichier_3 = NULL;

    	fichier_3 = fopen("Gmx.txt", "w");

        // A tracer
         for(k = 0; k <= N; k++)
        {
            printf("Gmx[%d] = %f\n", k, Gmx[k]);
			fprintf(fichier_3, "%f %f\n", T[k], Gmx[k]);
        }

        fclose(fichier_3);
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
*/

/*
    Utilisation d'un typedef car :
        - on va à de nombreuses reprises redéfinir les valeurs de nos variables : on rappel dans quel contexte elles sont utilisées

            --> cela nous permet de déclarer et redéclarer les variables inhérentes au pendule en mettant en avant le fait que ces
                variables appartienne à ce système précis
*/

        typedef struct Pendule
        {
            float Theta_ini;
            float Theta_ameliore;
            float Theta_lineaire[MAX];

            float Vitesse_ini;

            float w0;
            float T0;
            float Longueur;

        } Pendule;

        // Declaration de la structure
            Pendule pendule = {0, 0, 0, 0, 0, 0};

        // Initialisation Valeurs
            pendule.Theta_ini = pi/6;
            pendule.Longueur = 5;
            pendule.w0 = sqrt(g/pendule.Longueur);
            pendule.T0 = (2*pi)/pendule.w0;

// Test affichage
/*
        printf("Longeur du pendule : %f\n\n", pendule.Longueur);
        printf("Longeur du pendule : %f\n\n", pendule.T0);
        printf("Longeur du pendule : %f\n\n", pendule.w0);
*/


/*
    L'angle initiale restreint énormément la valeur de h possible.
    on a cos(p) - sqrt(3)/2;  p=[0;pi/6]

    On restreint pour se rpprocher d'étude en ca réel : on considère
    qu'on peut linériser le comportement du pendule pour de petits angles
*/


    /* --------------------------------------------------------------------------------------------- */

        printf("Trace du Theta en temps reel :\n");

        FILE *fichier_4 = NULL;
        fichier_4 = fopen("Theta_lineaire.txt", "w");

        for(k = 0; k <= N; k++)
        {
            pendule.Theta_lineaire[k] = (pendule.Theta_ini)*cos((pendule.w0)*T[k]);
            fprintf(fichier_4, "%f %f\n", T[k], pendule.Theta_lineaire[k]);
        }

        fclose(fichier_4);

        printf("\n\n --------------------------------------------------------------------------------------------- \n\n");

// 2.2
    printf("2.2 - Etude de la periode non-lineaire\n\n");

    printf("a - Calcul periode theorique\n\n");

/*
        On se familirise avec la méthode en utilisant la méthode composites des trapèzes :
                Integrale = Somme(i=0 --> n-1) (h/2)[f(a+ih) + f(a+(i+1)h)]

        A voir : Méthode composite de Simpson :
                Integrale = Somme (h/3)*(f(a+(n-2)h) + f(b) + 4*f(a+(n-1)h));
*/


// L'utilisation du Typedef permet de ne pa être perdu dans l'utilisation et la fonctionnalité de chaque variable
        typedef struct Periode
        {
            float f[MAX];
            float f_int;
            float f_int_h_opti;

            float theorique;
            float numerique;

            float h_2;
            float h_3;

            float Erreur_2_2;
            float E;

        } Periode;

        // Declaration de la structure
            Periode periode = {0, 0, 0, 0, 0, 0, 0, 0, 0};

        // Variables qui définissent la période
            periode.theorique = 0;
            periode.numerique = 0;

            periode.f[MAX] = 0;
            periode.f_int = 0;          // Intégrale de base
            periode.f_int_h_opti = 0;   // Intégrale avec h_optimisé

            periode.Erreur_2_2 = 0;     // Erreur sur calculs
            periode.E = 0.15;            // Erreur que l'on fixe pour Theta_opti et h_opti


        float a = 0, b = pendule.Theta_ini;
        float cste = 4*sqrt((pendule.Longueur)/(2*g));

            periode.h_2 = (b-a)/N;      // Pas de base
            periode.h_3 = 0;            // Pas_opti


        for (k = 0; k <= N-1; k++)
        {
            if(cos(a+k*(periode.h_2)) > cos(b)) // On évite de diviser par zéro ou d'avoir une racine négative
            {
                periode.f[k] = 1/sqrt(cos(a+k*(periode.h_2)) - cos(b));
            }

            else
            {
                // printf("\nValeur non definie\n\n");
            }
            // printf("periode.f[%d] = %f et cos(%f) = %f\n", k, periode.f[k], a+k*periode.h_2, cos(a+k*periode.h_2));
        }

        for (k = 0; k <= N-2; k++) // N-1 car trapèze compsite
        {
            periode.f_int += ((periode.h_2)/2)*(periode.f[k] + periode.f[k+1]);
            // printf("\nf_int = %f", f_int);
        }

        periode.theorique = cste*(periode.f_int); // Integrale de référence
        printf("\nPeriode_theorique par methode trapezes composites : %f s\n\n", periode.theorique);


        printf("\n\n --------------------------------------------------------------------------------------------- \n\n");


    printf("b - Adaptation du pas d'integration :\n\n");

        periode.f[k] = 0;
        i = 0;

        do
        {
            i++;

/*
    IMPORTANT : Compendre ce qu'on peut mettre comme pas periode.h_3 += pow(10,-5);
*/

            for (k = 0; k <= N; k++)
            {
                if(cos(a+k*(periode.h_3)) > cos(b)) // On évite de diviser par zéro ou d'avoir une racine négative
                {
                    periode.f[k] = 1/pow((cos(a+k*(periode.h_3)) - cos(b)), 0.5);
                }
                // printf("f[%d] = %f et cos(%f) = %f\n", k, periode.f[k], a+k*periode.h_3, cos(a+k*periode.h_3));
            }

            for (k = 0; k <= N-1; k++)
            {
                periode.f_int_h_opti += ((periode.h_3)/2)*(periode.f[k] + periode.f[k+1]);
            }

            // printf("L'integrale de la fonction f est : %f\n", cste*periode.f_int_h_opti);

            periode.numerique = cste*periode.f_int_h_opti;

/*
            Voir si on fait un écart en pourcent :
                (fabs(periode.numerique - periode.theorique))/periode.numerique
            ou un écart relatif :
                (fabs(periode.numerique - periode.theorique))
*/

            if((fabs(periode.numerique - periode.theorique)) < periode.E)
            {
                break;
            }

            // Trop de précision --> calculs très long !!
            periode.h_3 += pow(10, -5);  // On conseil un pow minimum puissance 4

        }while(periode.h_3 < periode.h_2);

        // En fonction de l'écart maximal souhaité, on voit le h correspondant, en limitant les calculs (h_3 += 0.00001
        printf("\nLe nouveau pas h_3 = %e alors que h_2 = %e. On l'a trouve en %d iterations\n", periode.h_3, periode.h_2, i);
        printf("\nPeriode_numerique avec Theta_ameliore = %f s\n", periode.numerique);


        printf("\n\n --------------------------------------------------------------------------------------------- \n\n");

//      Partie surement à supprimer
/*
        periode.f_int_h_opti = 0;

        printf("c - Etude de l'influence de l'angle initiale :\n\n");*
        for(pendule.Theta_ini = 0; pendule.Theta_ini <= 2*pi; pendule.Theta_ini += pi/4)
        {
            for (k = 0; k <= N; k++)
            {
                if(cos(a+k*(periode.h_2)) > cos(pendule.Theta_ini)) // On évite de diviser par zéro ou d'avoir une racine négative
                {
                    periode.f[k] = 1/sqrt(cos(a+k*(periode.h_2)) - cos(pendule.Theta_ini));
                    // printf("f[%d] = %f\n", k, f[k]);
                }
            }

            for (k = 0; k <= N-1; k++)
            {
                periode.f_int_h_opti = (h/2)*(periode.f[k] + periode.f[k+1]);
                // printf("L'integrale de la fonction f est : %f - f[%d] = %f - f[%d] = %f\n", cste*(periode.f_int_h_opti), k, periode.f[k], k+1, periode.f[k+1]);
            }

            printf("L'integrale de la fonction f est : %f\n", cste*periode.f_int_h_opti);
        }
*/

        printf("\n\n --------------------------------------------------------------------------------------------- \n\n");

        // On réinitialise les valeurs
        pendule.Theta_ini = pi*0.5;
        // pendule.Theta_ameliore = pi/2;
        periode.E = 0.5;
        //periode.f_int_h_opti = 0;
        i = 0;

        int PAS = 5;

        printf("\nd - Calcul de l'angle initial Maximal pour une erreur E donnee :\n\n");
//IMPORTANT        // On fixe ici un pourcentage, pas un écart relatif --> (a-b)/a != a-b

        do
        {
            i++;

            for (k = 0; k <= N; k++)
            {
                if(cos(a+k*(periode.h_2)) > cos(pendule.Theta_ini)) // On évite de diviser par zéro ou d'avoir une racine négative
                {
                    periode.f[k] = 1/sqrt(cos(a+k*(periode.h_2)) - cos(pendule.Theta_ini));
                    // printf("f[%d] = %f\n", k, f[k]);
                }
            }

            for (k = 0; k <= N-1; k++)
            {
                periode.f_int_h_opti += (periode.h_2/2)*(periode.f[k] + periode.f[k+1]);
            }

            periode.numerique = cste*periode.f_int_h_opti;

            // printf("La periode est : %f\n\n", periode.numerique);

            pendule.Theta_ini -= pi/PAS;

            // printf("pendule.Theta_ini = %f\n\n", pendule.Theta_ini);

            // printf("E = %f\n", periode.numerique - periode.theorique);

            if(pendule.Theta_ini <= 0)
            {
                PAS ++;
                pendule.Theta_ini = pi/2;
                for (k = 0; k <= N-1; k++)
                {
                    periode.f[k] = 0;
                    periode.f_int_h_opti = 0;
                }
            }

            else if((fabs(periode.numerique - periode.theorique)) <= periode.E)
            {
                printf("Le PAS %d est adapte\n", PAS);

                // En fonction de l'écart maximal souhaité, on voit le h correspondant, en limitant les calculs (h_3 += 0.00001
                printf("\nLe nouveau Theta_ameliore = %e alors que Theta_ini = %f. Cette valeur a ete atteinte en %d iterations\n", pendule.Theta_ini, pi/2, i);
                printf("La periode de la fonction est : %f\n\n", cste*periode.f_int_h_opti);
            }

        }while(PAS < 50);


        printf("\n\n --------------------------------------------------------------------------------------------- \n\n");

        periode.f[MAX] = 0;
        periode.f_int = 0;

/*
        printf("a = %f\n", a);
        printf("k = %d\n", k);
        printf("h_3 = %f\n", periode.h_3);
        printf("periode.f_int = %f\n", periode.f_int);
        printf("Theta_ameliore = %f\n", pendule.Theta_ameliore);

        for(i = 0; i <= N-1; i++)
        {
            printf("\nf = %f", periode.f[i]);
        }
*/

// Voir pourquoi mauvaise période
/*
        // On recalcule l'intégrale avec h_opti et theta_opti
        for (k = 0; k <= N-1; k++)
        {
            periode.f[k] = 1/sqrt(cos(a+k*(periode.h_3)) - cos(pendule.Theta_ameliore));
            // printf("f[%d] = %f et cos(%f) = %f\n", k, f[k], a+k*h_3, cos(a+k*h_3));
        }

        for (k = 0; k <= N-2; k++)
        {
            periode.f_int += ((periode.h_3)/2)*(periode.f[k] + periode.f[k+1]);
            // printf("\nf_int = %f", f_int);
        }

        printf("\nL'integrale de la fonction f est : %f\n", cste*periode.f_int);
*/


// 2.3
    printf("2.3 - Etude de la periode non-lineaire\n\n");

    printf("a - Calcul periode theorique\n\n");

/*
    Cours sur Matrices creuses :

    Il semblerait que la résolution d'équa diff soit faisable par le biais des matrices bandes/bandes creuses
    cf page 13 cours
    https://fr.wikipedia.org/wiki/Matrice_creuse

    Optimisation des dispositions dans la matrice
    https://fr.wikipedia.org/wiki/Structure_de_donn%C3%A9es

    Si on a ce type de matrice creuse : il existe des méthodes liées pour la résoluion de système (un peu comme LU, Cholesky, etc).
    cf page 23 cours
*/



// PARTIE 2

        // Developpement en série de Maclaurin (Taylor) de la fonction sinus pour approximer l'intégrale
        // On se limite à k = 8 --> factorielle = 2004189184. Au delà on dépasse la place possible pour un float
        // Voir si autre chose qu'un float pourrait convenir

/*
        for(k = 1; k <= 8; k++)
        {
            I_analytique = factorielle(2*k);

            // Valeur approché de l'intégrale
            // Grosse erreur - voir graphique correspondant
            Sommes_Taylor += pow((-1), k+1)*pow(n*pi, 2*k)/I_analytique;

            // Sinon on dit juste que Integrale_analytique = 0

            // printf("Puissance = %f - Factorielle = %f - Sommes Taylor = %f\n", pow(n*T[k], 2*k), I_analytique, Sommes_Taylor);
        }
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


void Decomposition_LU_Ecrase(float A[MAX][MAX], int n)
    {
        float L[MAX][MAX] = {0}, U[MAX][MAX] = {0};

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
                A[i][j] = U[i][j];  // On ecrase
            }

            for(i = j+1; i <= n; i++)
            {
                S2 = 0;
                for(k = 0; k <= j-1; k++)
                {
                    S2 += L[i][k]*U[k][j];
                }
                L[i][j] = (A[i][j] - S2)/U[j][j];
                A[i][j] = L[i][j];  // On ecrase
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


void Resolution_triangulaire_inf_Ecrase(float A[MAX][MAX], float y[MAX], float b[MAX], int n)
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

        y[i] = (b[i] - S1)/1;   // car L[j][j] = 1 et ca nous arrange pour l'écrasement des valeurs de A avec celles de L et U
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

void Resolution_triangulaire_sup_opti(float A[MAX][MAX], float x[MAX], float y[MAX], int n)
{
    float S1;

    int i, k;

    for(i = n; i >= 0; i--)
    {
        S1 = 0;
        for(k = i+1; k <= n; k++)
        {
            S1 += A[k][i]*x[k];
        }

        x[i] = (y[i] - S1)/A[i][i];
    }
}

int factorielle(int k)
{
    int i, factorielle = 1;

    for(i = 1; i <= k; i++)
    {
        factorielle = i*factorielle;
    }

    return factorielle;
}

