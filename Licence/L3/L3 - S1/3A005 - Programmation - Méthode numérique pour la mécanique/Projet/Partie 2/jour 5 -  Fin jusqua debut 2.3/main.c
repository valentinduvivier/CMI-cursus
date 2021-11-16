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

    Dimanche 08 Décembre : On termine avec les résultats déjà fait en cherchant les dernieres améliorations pour accélerer les calculs et les diminuer par la meme occasion.
                            --> On termine avec les erreurs, les incertitudes et les intégrales. On travail sur l'affichage des erreurs par une évolution graphique.


                          Affichage du polynôme des moindres carrés ssocié à sinus
                          Correction Bornes intégrale sinus Q 1.3
                          Corrections pour des éléments divisants possiblement par zéro.
                            --> On travail sur les valeurs limites :
                                                   - pour le nombre d'itération dans l'optimisation du pas et de Theta_ini
                                                   - Fonction à valeurs interdites comme 1/sqrt(cos...)

    Jeudi 12 Décembre : Difficulté sur la compréhension d'exercice pour Theta_numerique et polynôme des moindres carrés Q2.3.1 et 2.3.2
                        DDifficile de définir le caractère précis des résultats.

    Samedi 14 et Dimanche 15 Décembre :
                        --> Correction valeurs de l'erreur Q1.4.
                        --> Correction Méthodes de Taylor
                        --> Compréhension globale du probleme plus precise jusqu'a la Q2.3.2 incluse.

    Lundi 16 Décembre : Travail sur la position avec de grands angles. Pas de difficultés à l'horizon

    Revoir le PAS et cos sinus pour TAYLOR


	- Voir le lien entre erreur et méthodes LU - Cholesky
	- Voir interprétation valeur de n : plus précis pour n qui diminue
*/

// Listing des fonctions qui vont être utiliser lors de ce code

void Afficher_Vecteur(float x[MAX], int n); //Code affichage Vecteurs
void Afficher_Matrice(float A[MAX][MAX], int n); // Code affichage Matrices

void Reinitialisation_Variable(float var1[MAX], int n); // Fonction Reinitialisant une variable de type tableau à zéro

void Decomposition_LU_Ecrase(float A[MAX][MAX], unsigned long int n); // Fonction décomposition LU où l'on écrase les aleurs de L et U au sein de la matrice A
void Decomposition_LLt(float A[MAX][MAX], float L[MAX][MAX], unsigned long int n); // Fonction décomposition Cholesky

void Resolution_triangulaire_inf(float A[MAX][MAX], float y[MAX], float b[MAX], unsigned long int n); // Fonction pour résolution d'un système triangulaire inférieur pour méthode Cholesky
void Resolution_triangulaire_inf_Ecrase(float A[MAX][MAX], float y[MAX], float b[MAX], unsigned long int n); // Fonction pour résolution d'un système triangulaire sinférieur pour méthode LU

void Resolution_triangulaire_sup(float A[MAX][MAX], float x[MAX], float y[MAX], unsigned long int n); // Fonction pour résolution d'un système triangulaire supérieur pour méthode LU
void Resolution_triangulaire_sup_opti(float A[MAX][MAX], float x[MAX], float y[MAX], unsigned long int n); // Fonction pour résolution d'un système triangulaire supérieur Cholesky

int factorielle(int k); // Fonction retournant le vactoriel de k

float F_Taylor(int k, float Var1, float *Var2, float Var3, float Var4, int p); // Fonction pour Q 2.3.2 qui renvoie la dérivée p_ième de la fonction cos(w0*t)

int main()
{
        float T[MAX] = {0}, SIN[MAX] = {0};		// Initialisation tableaux de temps et de SIN = sin(T*n)
        float Pas_de_Temps;

        int k, n, i;   // Variables pour les boucles notamment


printf("---------------------------------- Partie 1 : --------------------------------------- \n\n");

    // 1.1
    printf("1.1\n\n");

		//a
		printf("a : Creation vecteur Temps : \n\n");

        Pas_de_Temps = (2*pi)/100;

        T[0] = 0;       // Initialisation première valeur de Temps pour qu'elle soit bien à 0
		for(k = 0; k <= N; k++)
		{
			T[k+1] = T[k] + Pas_de_Temps;
		}

		printf(" Vecteur Temps cree\n");
		printf("\n\n");


		//b
		printf("b : Creation vecteur SIN : \n\n");

		n = 10;		// Nombre étudiant = 3700499

		for(k = 0; k <= N; k++)
		{
			SIN[k] = sin(n*T[k]);
		}

		printf(" Vecteur SIN cree\n");


		printf("\n\n --------------------------------------------------------------------------------------------- \n\n");


    // 1.2
    printf("1.2 - Calcul des derivees\n\n");

		float Derivee_1ere[MAX] = {0}, Derivee_2nde[MAX] = {0};
		float Verification_derivee_1ere[MAX], Verification_derivee_2nde[MAX], ERRR[MAX] = {0}; // Verification permet de s'assurer que les dérivées ce sont bien effectuée
		float h = (2*pi)/N;  // le pas est égale à l'écart entre 2 valeur de temps

		printf("Calcul des derivees premiere et seconde de SIN(t) = sin(nt) : \n\n");


	// Différence progressive --> précision de h
        // Utilisation de gnuplot pour tracer les dérivées 1ere et 2nde
		printf("Ecriture des derivees dans 2 fichiers textes : \n\n");

		FILE *fichier_1 = NULL;
		FILE *fichier_2 = NULL;

		FILE *fichier_3 = NULL;
		FILE *fichier_4 = NULL;

    	fichier_1 = fopen("Derivee_1ere.txt", "w");
    	fichier_2 = fopen("Derivee_2nde.txt", "w");

    	fichier_3 = fopen("Verification_derivee_1ere.txt", "w");
    	fichier_4 = fopen("Verification_derivee_2nde.txt", "w");


        // On définit nous même SIN(N+1) (utilisation pour Derivée_1ere et Dérivée_2nde)
		SIN[N+1] = sin(n*(T[N] + Pas_de_Temps));		// On ne peut pas faire sin(nT[N+1]) car la valeur la plus haute de k = N car k <= N

		for(k = 0; k <= N; k++)
		{
			Derivee_1ere[k] = (SIN[k+1] - SIN[k])/h;

/*
        Calcul de cosinus pour vérifier nos résultats : dsin(nt) = ncos(nt)  -->  le but est d'avoir une idée visuelle
        de l'écart entre valeur numérique et valeur analytique. On doit avoir Verification_derivee_1ere[k] = Derivee_1ere[k]
*/

			Verification_derivee_1ere[k] = n*cos(n*T[k]);

			fprintf(fichier_1, "%f %f\n", T[k], Derivee_1ere[k]);
            fprintf(fichier_3, "%f %f\n", T[k], Verification_derivee_1ere[k]);

            // Calcul de l'erreur entre les 2 fonctions.
            ERRR[-1] = 0;
            ERRR[k] = ERRR[k-1] + fabs(Verification_derivee_1ere[k] - Derivee_1ere[k]);
		}

		printf("\nDerivee_1ere calculee\n\n");
        printf("Erreur_deriveee_1ere = %f\n", ERRR[N]/N);


        // Reinitialisation de l'erreur pour réutilisation future
        Reinitialisation_Variable(ERRR, N);


    // Différence centrée --> précision de h^2
        // On calcule nous même SIN[-1] afin de l'utiliser pour la dérivée seconde
		SIN[-1] = -sin(Pas_de_Temps);

		for(k = 0; k <= N; k++)
		{
			Derivee_2nde[k] = (SIN[k-1] - 2*SIN[k] + SIN[k+1])/(pow(h, 2));
			Verification_derivee_2nde[k] = -pow(n,2)*sin(n*T[k]);		// On doit avoir Verification_derivee_2nde[k] = Derivee_2nde[k]

			fprintf(fichier_2, "%f %f\n", T[k], Derivee_2nde[k]);
			fprintf(fichier_4, "%f %f\n", T[k], Verification_derivee_2nde[k]);

			ERRR[-1] = 0;
            ERRR[k] = ERRR[k-1] + fabs(Verification_derivee_2nde[k] - Derivee_2nde[k]);
        }

		printf("\nDerivee_2nde calculee\n\n");
        printf("Erreur_deriveee_2nde = %f\n\n", ERRR[N]/N); // On fait la moyenne des erreurs

		printf("\nFichiers textes crees\n");


        fclose(fichier_1);
		fclose(fichier_2);

		fclose(fichier_3);
		fclose(fichier_4);


		printf("\n\n --------------------------------------------------------------------------------------------- \n\n");

//1.3
    printf("1.3 - Calcul de l'integrale\n\n");

        // Variable pour définir l'outil de calcul d'intégration
        float Integrale_1_3 = 0, F[MAX] = {0};
        float Erreur_1_3 = 0;

        h = pi/N; // on a h = (b-a)/N avec b = pi et a  = 0

        printf("Calcul de l'integrale - valeur analytique : \n\n");

        //Utilisation de la méthode des trapèzes pour le calcul de l'intégrale entre [0;pi] de SIN

        for(k = 0; k <= (N-1)/2; k++) // On s'arrête à PI soit (N-1)/2 pour la méthode des trapèzes
        {
            F[k] = sin(k*h);
        }

        for(k = 0; k <= (N-2)/2; k++) // (N-2)/2 = n-1, car  pour la méthode des trapèzes on va de 0 à k_max = n-1 : Soit k_max = pi, on prend n = 50 --> n-1 = 49 = (N-2)/2
        {
            Integrale_1_3 += (h/2)*(F[k+1] + F[k]); // Calcul de l'intégrale Q 1.3
        }

        printf("Integrale = %f\n", Integrale_1_3);


        // Calcul de la précision - comparasion solution analytique/numérique
         Erreur_1_3 = fabs(Integrale_1_3 - ((1 - cos(n*pi))/n)); 	//Resultat analytique : Integrale = 0 car n = 10.


        printf("\nL'erreur sur le calcul est : Erreur Q.3 = %f\n\n", Erreur_1_3);

        printf("\n\n --------------------------------------------------------------------------------------------- \n\n");

    //1.4
    printf("1.4 - Calcul du polynome associe Gm(x)\n\n");

/*
        Cette partie, pour l'optimisation nous mène à 2 méthodes
            - Cholesky : on a qu'une matrice L déduite de U. On a pas besoin de créer Lt et on a moins de calcul que pour LU
                    --> On a moins d'opération, donc on utilise sans doute moins nos variables qui contiennent des erreurs,
                        d'où le fait qu'il y ait moins d'erreur avec cette méthode

            - LU : On utilise dans notre cas une méthode créant temporairement L et U, pour ensuite les utiliser pour recalculer
                    les composantes de A en écrasant les données précédentes. On a un Matrice de moins que la méthode LU classique
*/

        // Définition des variables pour la méthode des moiindres carrés
        float Matrice_x[MAX][MAX] = {0}, Vect_x[MAX] = {0};
        float L[MAX][MAX] = {0};

        float y[MAX] = {0}, x[MAX] = {0};     // x étant le vecteur que l'on cherche à résoudre dans cette question

        float Gmx[MAX] = {0};  // Polynôme Gmx(x) = a0 + a1x + a2x² + ...
        float Erreur_1_4 = 0, S1 = 0;

        int m = 6;      // numéro étudiant = 3700091  --> On aura un polynôme de degré 6
        int j;      // i et k déjà défini

/*
            Application de la méthode des moindres carrés : Gm(xi) = Somme(ai*xi), entre [0;m]
            Dans notre cas Gm(x) = a0 + a1x + a2x² + a3x³ + a4x^4 + a5x^5 + a6x^6  On cherche ainsi, a0, ... , a6
            Erreur Quadratique = Somme(Gm(x) - SIN)²

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


// Afficher_Matrice(Matrice_x, m);
    //Vecteur
        for (k = 0; k <= N; k++)
        {
            for (i = 0; i <= m; i++)
            {
                Vect_x[i] += SIN[k]*pow(T[k], i);
            }
        }


//Cholesky

        printf("\n    Decomposition Cholesky\n\n");

        Decomposition_LLt(Matrice_x, L, m);
        printf("    Decomposition effectue\n\n");

        printf("\nResolution du systeme : \n\n");
        Resolution_triangulaire_inf(L, y, Vect_x, m);
        Resolution_triangulaire_sup_opti(L, x, y, m);

        printf("Voici les coefficients pour a1, ..., a6 obtenus :\n");
        Afficher_Vecteur(x, m);


        // A Decommanter si vous voulez tester méthode LU amélioré avec écrasement de valeur
//LU
//        printf("\n    Decomposition LU Ecrase : \n\n");

/*
        printf("Resolution systeme - calcul a0, ..., a6 :\n\n");

        Decomposition_LU_Ecrase(Matrice_x, m);
        printf("    Decomposition effectue\n\n");

// A décommanter si on souhaite afficher les variables
//          Afficher_Matrice(L, m);
//          Afficher_Matrice(U, m);

        Resolution_triangulaire_inf_Ecrase(Matrice_x, y, Vect_x, m);
        Resolution_triangulaire_sup(Matrice_x, x, y, m);

//          Afficher_Vecteur(y, m);

        printf("Voici les coefficients obtenus pour a0, ..., a6 :\n\n");
        Afficher_Vecteur(x, m);
*/

        FILE *fichier_5 = NULL;
    	fichier_5 = fopen("Gmx.txt", "w"); // Enregistrement des données du polynôme des moindres carrés dans un fichier texte

        for(k = 0; k <= N; k++)
        {
            // Affichage de la fonction ainsi trouver Gm(x)
            for(i = 0; i <= m; i++)
            {
                Gmx[k] += x[i]*pow(T[k],i);
            }

			fprintf(fichier_5, "%f %f\n", T[k], Gmx[k]);

            // Calcul d'erreur quadratique
            S1 += pow((SIN[k] - Gmx[k]), 2);
            Erreur_1_4 += pow(S1, 0.5)/(N-1);
        }

        printf("\nErreur Q.4 = %0.5f\n", Erreur_1_4);

        fclose(fichier_5);


    printf("\n---------------------------------- Partie 2 : Application au pendule --------------------------------------- \n\n");


// 2.1
    printf("2.1 - Evolution du pendule en fonction du temps\n\n");

/*
    Utilisation d'un typedef car :
        - on va à de nombreuses reprises redéfinir les valeurs de nos variables : on rappel dans quel contexte elles sont utilisées

            --> cela nous permet de déclarer et redéclarer les variables inhérentes au pendule en mettant en avant le fait que ces
                variables appartienne à ce système précis
*/

        typedef struct Pendule
        {
            float Theta_ini;
            float Theta_lineaire[MAX];

            float Vitesse_ini;

            float w0;
            float T0;
            float Longueur;

        } Pendule;

            // Declaration de la structure
            Pendule pendule = {0, 0, 0, 0, 0, 0};

        // Initialisation Valeurs
            pendule.Longueur = 5;                   // La longeyur du pendule est de 5m
            pendule.Theta_ini = pi/8;               // Le pendule démarre à pi/8 --> petit angle
            pendule.w0 = sqrt(g/pendule.Longueur);  // Vitesse de rotation du pendule
            pendule.T0 = (2*pi)/pendule.w0;         // Période du pendule issu de la théorie linéarisée

/*
    L'angle initiale restreint énormément la valeur de h possible.

    On restreint pour se rpprocher de l'étude en cas réel : on considère
    qu'on peut linériser le comportement du pendule pour de petits angles
*/


    /* --------------------------------------------------------------------------------------------- */

        printf("Trace du Theta - fonction linearisee :\n");

        FILE *fichier_6 = NULL;
        fichier_6 = fopen("Theta_lineaire.txt", "w");

        for(k = 0; k <= N; k++)
        {
            pendule.Theta_lineaire[k] = (pendule.Theta_ini)*cos((pendule.w0)*T[k]);

            fprintf(fichier_6, "%f %f\n", T[k], pendule.Theta_lineaire[k]);
        }

        fclose(fichier_6);


        printf("La periode pendule.T0 = %f\n", pendule.T0);

        printf("\n\n --------------------------------------------------------------------------------------------- \n\n");


// 2.2
    printf("2.2 - Etude de la periode non-lineaire\n\n");

    printf("a - Calcul periode theorique\n\n");

/*
        On se familirise avec la méthode en utilisant la méthode composites des trapèzes :
                Integrale = Somme(i=0 --> n-1) (h/2)[f(a+ih) + f(a+(i+1)h)]

        A voir : Méthode composite de Simpson :
                Integrale = Somme (h/3)*(f(a+(i-2)h) + f(b) + f(a) + 4*f(a+(i-1)h));
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
            }

            for (k = 0; k <= N; k++) // N-1 car trapèze compsite
            {
                periode.f_int += ((periode.h_2)/2)*(periode.f[k] + periode.f[k+1]);
                //periode.f_int += (periode.h_2/3)*(2*periode.f[2*k-2] + periode.f[3] + periode.f[0] + 4*periode.f[2*k-1]); //Méthode Simpson
            }

        periode.theorique = cste*(periode.f_int); // Integrale de référence
        printf("La periode numerique selon la formule B est : %f\n", periode.theorique);

        periode.Erreur_2_2 = fabs(pendule.T0 - periode.theorique);
        printf("L'erreur sur cette période est de %f\n", periode.Erreur_2_2);

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
            periode.h_3 += pow(10, -5);  // On conseil un pow minimum puissance -5

        }while(periode.h_3 < periode.h_2);

        // En fonction de l'écart maximal souhaité, on voit le h correspondant, en limitant les calculs (h_3 += 0.00001
        printf("\nLe nouveau pas h_3 = %e alors que h_2 = %e. On l'a trouve en %d iterations\n", periode.h_3, periode.h_2, i);
        printf("\nPeriode_numerique avec Theta_ameliore = %f s\n", periode.numerique);


        printf("\n\n --------------------------------------------------------------------------------------------- \n\n");


        float Theta_ini_ameliore[MAX]; // Ces variables vont stocker les Pas idéaux et les périodes associées

        pendule.Theta_ini = pi/2;   // On part avec un angle important et on va décroissant pour trouver Thete_MAX

        i = 0; // i va nous permettre de connaitre le nombre d'itération mais il va surtout permettre de recueillir les theta recherchés


/*
        Ce que nous appelons PAS n'a pas de rapport avec les pas h_2 et h_3
        C'est en fait ce qui va nous aider à voir par quoi diminuer Theta_ini pour aboutir aux précisions souhaité.

        L'utilisation de ce PAS n'est qu'un "Bonus" pour se différencier. On pourrait l'enlever et choisir nous même
        un PAS suffisamment grand.
*/
        periode.E = 0.05; //10e-5; // On ne prend les périodes qui ont un écart de 10e-5 avec la théorie

        int PAS = 1; //Notre pas dans la boucle est tel que les pas PAS = 1 et PAS = 2 ne sont pas utile

        printf("\nd - Calcul de l'angle initial Maximal pour une erreur E donnee :\n\n");

// DECOMMENTER APRES FIN DU CODE
/*
        printf("    Quelle Erreur E souhaitee-vous ?\n");
        scanf("%f", &periode.E);

        printf("\n\n");
*/

// L'utilisation de PAS est du au fait que la boucle n'est pas valide pour tous les PAS dans le cas d'une erreur absolue (E1 - E2)
        do
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
                periode.f_int_h_opti += (periode.h_2/2)*(periode.f[k] + periode.f[k+1]);
            }

            periode.numerique = cste*periode.f_int_h_opti;

            pendule.Theta_ini -= pi/PAS;

            if(pendule.Theta_ini <= 0)
            {
                PAS ++;

                // On réinitialise les variables qui sont soit déclarer en dehors de la boucle ou qui ne se réinitialise pas après chaque itération
                pendule.Theta_ini = pi/2;
                for (k = 0; k <= N; k++)
                {
                    periode.f_int_h_opti = 0;
                }
            }

            else if(fabs(periode.theorique - periode.numerique)/periode.theorique < periode.E) // on fait le rapport des deux pour voir le poucentage de ressemblance
            {
// IMPORTANT A DECOMMENTER EN FIN DE CODE
                // printf("Le PAS %d est adapte\n", PAS); // La periode approchee est a 85 pourcent identique a la periode theorique, avec E = 0.85

                Theta_ini_ameliore[i] = pendule.Theta_ini;
                i++;

                // En fonction de l'écart maximal souhaité, on voit le h correspondant, en limitant les calculs (h_3 += 0.00001
                // printf("\nLe theta_adapte correspondant est = %e alors que Theta_ini = %f. Cette valeur a ete atteinte en %d iterations\n", pendule.Theta_ini, pi/2, i);
                // printf("La periode de la fonction est : %f\n\n", cste*periode.f_int_h_opti);
            }

        }while(PAS <= 100);


/*--------------------------------------------------------------------------------------------- */


/*
        On remarque que les theta obtenus, qui sont stocké dans un tableau, sont quasi
        tout le temps dans un ordre naturellement croissant :
            --> On améliore le pas, ce qui améliore de fait l'angle initiale. Plus le pas
        augmente, plus on a un théta grand. Le temps tant vers la valeur optimale attendue.
        Si on augmente le PAS, alors on remarque que ce théta varie de moins en moins.
*/

        // On ordonne les Théta_ini pour avoir le plus grand en dernier
        float c = 0;

        for(k = 0; k <= i-1; k++)
        {
            for(j = k+1; j <= i-1; j++)
            {
                if(Theta_ini_ameliore[k] > Theta_ini_ameliore[j])
                {
                    c = Theta_ini_ameliore[k];
                    Theta_ini_ameliore[k] = Theta_ini_ameliore[j];
                    Theta_ini_ameliore[j] = c;
                }
            }
        }

        printf("\n  Le theta_ini max pour lequel on atteint la precision souhaitee est : %f\n\n", Theta_ini_ameliore[i-1]);


    printf("\n\n --------------------------------------------------------------------------------------------- \n\n");


// 2.3
    printf("\n2.3 - Etude de l'Equa Diff A\n\n");

    printf("a - Calcul variation angulaire - petits angles\n\n");

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

        typedef struct Numerique
        {
            float position_pendule[MAX];
            float Derivee_2nde_Position[MAX];

            float periode;
            float position_theta[MAX];

            float Theta_ini;
            float Theta[MAX];

            float h;

            float t[MAX];

        } Numerique;

        Numerique numerique = {0, 0, 0, 0, 0, 0, 0, 0};

/*
    On a wo = sqrt(g/l) = 1.400714104
        --> Pour avoir une période il faut donc 2*pi, soit une multiplication par 2*pi/w0 = 4.485701465
*/

        numerique.h = 2*pi/100;
        numerique.Theta_ini = pi/100;

        k = 0;

        // Linearisation de la formule à l'aide de la dérivée seconde : différence finie centrée
        // Il semble y avoir un certain nombre d'itératon avant que le pendule n'ait le comportement
        // oscillant attendu. En deça de cette valeur, le pendule augmente alors qu'il devrait diminuer.
            numerique.Theta[0] = numerique.Theta_ini;
            numerique.Theta[-1] = (2*numerique.Theta[0] - pow(pendule.w0*numerique.h, 2)*sin(numerique.Theta[0]))/2; // = numerique.Theta[1] d'après condition initiale sur la vitesse
            //numerique.Theta[1] = numerique.Theta[-1];

            //numerique.t[1] = numerique.h;

// Code qui fonctionne ??? - Le faire s'arréter après une période : FAIT ??????
        do
        {
            numerique.Theta[k+1] = 2*numerique.Theta[k] - numerique.Theta[k-1] - pow(pendule.w0*numerique.h, 2)*sin(numerique.Theta[k]);

            printf("numerique.Theta[%d] = %f\n", k, numerique.Theta[k]);

            numerique.t[k+1] = numerique.t[k] + numerique.h;
            k++;

        }while(numerique.t[k]*pendule.w0 <= 2*pi); // On s'arrête dès que cos(w0*T) = cos(2pi) soit quand wo*T = 2pi


    printf("\n\n --------------------------------------------------------------------------------------------- \n\n");


    printf("b - Comparaison Moindre-Carres/Taylor pour l'etude agulaire\n\n");

// Moindre Carré

    float O[MAX][MAX] = {0}, P[MAX] = {0};
    float XX[MAX] = {0}, YY[MAX] = {0}, LL[MAX][MAX] = {0};

    int p = 0, S2; // S2 est utilisé pour le calcul de l'erreur  et p qui définira la taille du polynôme pour les moindre carré (p=6)

// Systeme a resoudre
    //Matrice
        //Création matrice pour système à résoudre

        for (p = 0; p <= k-1; p++)  // k-1 = N
        {
            for (i = 0; i <= m; i++)
            {
                for(j = 0; j <= m; j++)
                {
                    O[i][j] += pow(numerique.t[p], i+j);
                }
            }
        }

        // Afficher_Matrice(O, m);

    //Vecteur
        for (p = 0; p <= k-1; p++)  // k-1 = N
        {
            for (i = 0; i <= m; i++)
            {
                P[i] += numerique.Theta[p]*pow(numerique.t[p], i);
            }
        }

        // Afficher_Vecteur(P, m);


        printf("\n    Decomposition Cholesky\n\n");
        Decomposition_LLt(O, LL, m);
        printf("    Decomposition effectue\n\n");

        printf("\nResolution du systeme : \n\n");

        Resolution_triangulaire_inf(LL, YY, P, m);
        Resolution_triangulaire_sup_opti(LL, XX, YY, m);

        //Afficher_Matrice(LL, m);

        printf("Voici les coefficients pour a0, ..., a6 obtenus :\n\n");

        Afficher_Vecteur(XX, m);

    printf("\n\n --------------------------------------------------------------------------------------------- \n\n");


// A decommanter si vous voulez tester la méthode LU optimisée
/*
        printf("\n    Decomposition LU Ecrase : \n\n");

        // A decomenter si on fait la méthode LU  classique
        Decomposition_LU_Ecrase(O, m);
        printf("    Decomposition effectue\n\n");

        printf("\nResolution du systeme : \n\n");

        // Resolution_triangulaire_inf(Matrice_x, y, Vect_x, m);
        Resolution_triangulaire_inf_Ecrase(O, YY, P, m);
        Resolution_triangulaire_sup(O, XX, YY, m);

        printf("Voici les coefficients obtenus pour a0, ..., a6 :\n\n");
        Afficher_Vecteur(XX, m);
*/

//printf("\n\n --------------------------------------------------------------------------------------------- \n\n");


// Taylor à l'ordre 6 de Theta_linearise = Theta_ini*cos(w0*t)
/*
    Reste  = (x-a)^n

    f(x) = Somme (0-->n) dérivée_k_ième__f(a)/(factoriel k) + Reste.
    Le reste est négligeable devant (x-a)^n donc cette formule donne un minorant du Reste


    f(x) = Theta_ini*cos(pendule.wo*T[k])

*/

    numerique.Theta_ini = pi/100;
    numerique.h = 2*pi/100;
    k = 0;

FILE *fichier_8 = NULL;
    fichier_8 = fopen("TAYLOR.txt", "w");

    do
    {
        for(p = 0; p <= m; p++)
        {
            numerique.position_theta[k] += F_Taylor(k, pendule.w0, T, numerique.Theta_ini, numerique.h, p);
        }

        ////numerique.position_theta[k] += pow(b, m);
        //printf("numerique.position_theta[%d] = %f\n", k, numerique.position_theta[k]);
        fprintf(fichier_8, "%f %f\n", numerique.t[k], numerique.position_theta[k]);

        numerique.t[k+1] = numerique.t[k] + numerique.h;
        k++;
    }while(numerique.t[k]*pendule.w0 <= 2*pi); // On s'arrete des que l'on a cow(w0*t) = cos(2*pi) soit qd w0*t = 2*pi

fclose(fichier_8);
//printf("\n\n --------------------------------------------------------------------------------------------- \n\n");


    printf("c - Calcul variation angulaire numerique - grands angles\n\n");

    for(k = 0; k <= N; k++)
    {
        numerique.position_theta[k] = 0;
    }
    numerique.Theta_ini = 2*pi/3;
    numerique.h = 2*pi/100;

    printf("d - Enregistrement des valeurs dans un fichier texte\n\n");
    FILE *fichier_7 = NULL;
    fichier_7 = fopen("Angle_fct_temps.txt", "w");


    k = 0;

    do
    {
        // Calcul des positions pour Theta_ini Grand
        for(p = 0; p <= m; p++)
        {
            numerique.position_theta[k] += F_Taylor(k, pendule.w0, T, numerique.Theta_ini, numerique.h, p);
        }

        ////numerique.position_theta[k] += pow(b, m);
        // Affichage de la position
        // printf("numerique.position_theta[%d] = %f\n", k, numerique.position_theta[k]);

        fprintf(fichier_7, "%f %f\n", numerique.t[k], numerique.position_theta[k]);

        numerique.periode += numerique.h;
        k++;
    }while(numerique.periode*pendule.w0 <= 2*pi);

    printf("    Variations angulaire calculees\n\n");
    printf("Periode numerique pour Thete_grand = %f\n\n", numerique.periode);

    fclose(fichier_7);
    printf("    Valeurs enregistrees dans le fichier texte\n\n");

    return 0;
}



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

void Reinitialisation_Variable(float var1[MAX], int n)
{
    int k;

    for(k = 0; k <= n; k++)
    {
        var1[k] = 0;
    }
}

void Decomposition_LU(float A[MAX][MAX], float L[MAX][MAX], float U[MAX][MAX], unsigned long int n)
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


void Decomposition_LU_Ecrase(float A[MAX][MAX], unsigned long int n)
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


void Decomposition_LLt(float A[MAX][MAX], float L[MAX][MAX], unsigned long int n)
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


void Resolution_triangulaire_inf(float A[MAX][MAX], float y[MAX], float b[MAX], unsigned long int n)
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


void Resolution_triangulaire_inf_Ecrase(float A[MAX][MAX], float y[MAX], float b[MAX], unsigned long int n)
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

void Resolution_triangulaire_sup(float A[MAX][MAX], float x[MAX], float y[MAX], unsigned long int n)
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

void Resolution_triangulaire_sup_opti(float A[MAX][MAX], float x[MAX], float y[MAX], unsigned long int n)
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

float F_Taylor(int k, float Var1, float *Var2, float Var3, float Var4, int p)
{

        // Var1 = pendule.wo
        // Var2 = T[k]
        // Var3 = numerique.Theta_ini
        // Var4 = numerique.h

    if(p==0)
    {
        return Var3*cos(Var1*Var2[k]); // car Theta_ini*cos(w0*T[k]) = Theta_ini*cos(0) = Theta_ini
    }

    else if(p%2 == 0)
    {
        if((p/2)%2 == 0)
        {
            return +(pow(Var1*Var4, p)*Var3*cos(Var1*Var2[k]))/factorielle(p);
        }

        else if((p/2)%2 == 1)
        {
            return -(pow(Var1*Var4, p)*Var3*cos(Var1*Var2[k]))/factorielle(p);
        }
    }

    else if(p%2 == 1)
    {
        if((p/2)%2 == 0)
        {
            return -(pow(Var1*Var4, p)*Var3*sin(Var1*Var2[k]))/factorielle(p);
        }

        else if((p/2)%2 == 1)
        {
            return +(pow(Var1*Var4, p)*Var3*sin(Var1*Var2[k]))/factorielle(p);
        }
    }

    else
    {
        printf("Erreur\n\n");
    }

    return 0.0;
}

int factorielle(int k)
{
    int i;
    unsigned long factorielle = 1;

    if(k == 0)
    {
        return 1;
    }

    for(i = 1; i <= k; i++)
    {
        factorielle = i*factorielle;
    }

    return factorielle;
}

