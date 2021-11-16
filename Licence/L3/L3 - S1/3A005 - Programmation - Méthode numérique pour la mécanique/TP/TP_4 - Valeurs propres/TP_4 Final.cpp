#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 100



void afficher_vect(float x[MAX], int n)
{
	int i;

	for (i = 0; i < n; i++)
    {
        printf("%1f\n", x[i]);
    }
}

void afficher_mat(float A[MAX][MAX], int n)
{
	int i, j;

	for (i = 0; i < n; i++)
    {
		for (j = 0; j < n; j++)
        {
			printf("%.1f\t", A[i][j]);
		}

		printf("\n");
	}
}

void remplir_vect(float x[MAX], int n)
{
	int i;

	for (i = 0; i < n; i++)
    {
		printf("valeur [%d] : ", i + 1);
	    if (scanf("%f", &x[i]) !=1) printf("error");
	}
}

void remplir_mat(float A[MAX][MAX], int n)
{
	int i, j;

	for (i = 0; i < n; i++)
    {
		for (j = 0; j < n; j++)
        {
            printf("valeur [%d][%d] : ", i+1, j+1);
            if (scanf("%f",&A[i][j]) !=1) printf("error");
		}

		printf("\n");
	}
}



void resol_trig_inf( float A[MAX][MAX], float x[MAX], float b[MAX], int n)
{
	int i, k;
	float l;

	for (i = 0; i < n; i++)
    {
		l = 0;

		for (k = 0; k <= i-1; k++)
        {
            l += A[i][k]*x[k];
        }

        x[i] = (b[i] - l)/A[i][i];
    }
}


void resol_trig_sup(float A[MAX][MAX], float x[MAX], float b[MAX], int n)
{
	int i, k;
	float l;

	for (i = n; i >= 0; i--)
    {
		l=0;

		for (k = i+1; k < n; k++)
        {
            l += A[i][k]*x[k];
        }

		x[i] = (b[i] -l)/A[i][i];
	}
}

void factoriser_LU(float A[MAX][MAX], float L[MAX][MAX], float U[MAX][MAX], int n)
{
	int i, j, k;
	float l;

	for (j = 0; j <= n-1; j++)
    {
		for (i = 0; i <= j; i++)
        {
			l = 0;

			for (k = 0; k <= i-1; k++)
            {
                l += L[i][k]*U[k][j];
            }

            U[i][j] = A[i][j] - l;
		}

		L[j][j] = 1;

		for (i = j+1; i <= n-1; i++)
        {
			l = 0;

			for (k = 0; k <= j-1; k++)
            {
                l += L[i][k]*U[k][j];
            }

            L[i][j] = (A[i][j] - l)/U[j][j];
		}
	}
}


void resol_LU(float L[MAX][MAX], float U[MAX][MAX], float u[MAX], float v[MAX], int n)
{
	float y[MAX];

    resol_trig_sup(L, y, v, n);
    resol_trig_inf(U, u, y, n);
}


void Matrice_transp(float TA[MAX][MAX], float A[MAX][MAX], int n)
{
	int i, j;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			TA[i][j] = A[j][i];
		}
	}
}

void Reini_vect(float x0[MAX], int n)
{
	int i;

	for(i = 0; i < n; i++)
	{
		x0[i] = 1;
	}
}


float norme2_vect(float x[MAX], int n)
{
	float var = 0;
	int i;

	for (i = 0; i < n; i++)
    {
        var += (x[i]*x[i]);
    }

	return sqrt(var);
}

float puissance_it(float A[MAX][MAX], float x0[MAX], float x[MAX], float eps, int n)
{
	int i, j, k;

	float v[MAX], L = 1, L0 = 1000000, r, var;

	for(k = 0; k < n; k++)
	{
		for (i = 0; i < n; i++)
		{
			v[i] = x0[i]/norme2_vect(x0, n);
		}

		for (i = 0; i < n; i++)
		{
			x0[i] = 0;

			for(j = 0; j < n; j++)
			{
				x0[i] += A[i][j]*v[j];
			}
		}

		L0 = L;
		L = 0;

		for (i = 0; i < n; i++)
		{
			L += v[i]*x0[i];
		}

		for (i = 0; i < n; i++)
		{
			x[i] = v[i];
		}


		if (fabs(L - L0) < eps)
		{
			break;
		}
	}

	return L;
}


float puissance_inverse(float A[MAX][MAX], float x0[MAX], float x[MAX], float eps, int n)
{
	int i, j, k = 0;

	float L[MAX][MAX], U[MAX][MAX];
	float v[MAX], u[MAX], Lam = 1, L0 = 10000000, r, var;

    factoriser_LU(A, L, U, n);

    for(k = 0; k < n; k++)
	{
	    for (i = 0; i < n; i++)
		{
			u[i] = x0[i];
		}

		for (i = 0; i < n; i++)
		{
			v[i] = u[i]/norme2_vect(u, n);
		}

        resol_LU(L, U, u, v, n);

		L0 = Lam;
		Lam = 0;

		for (i = 0; i < n; i++)
		{
			Lam += v[i]*x0[i];
		}

		for (i = 0; i < n; i++)
		{
			x[i] = v[i];
		}

		if (fabs(Lam - L0) < eps)
		{
			break;
		}
	}

	return Lam;
}


float deflation_it(float A[MAX][MAX], float x0_ini[MAX], float x1[MAX], float x_vect[MAX], float L, int n, float eps)
{

	float M[MAX][MAX], A1[MAX][MAX], y[MAX], var, C;
	float TA[MAX][MAX], PT[MAX][MAX];
	float PS = 0;

	int i, j;

	Matrice_transp(TA, A, n);

    printf("\n\n/*-----------------------------------------------------------------------------------------------------*/\n\n");

    // On cherche à calculer y, le vecteur propre correspondant
    // à la valeur propre Lambda 1 dans le cas de la matrice transposée
	puissance_it(TA, x0_ini, y, eps, n);

	for(i = 0; i < n; i++)
	{
		PS += x1[i]*y[i];
	}

	// Constant
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			PT[i][j] = x1[i]*y[j];
			M[i][j] = PT[i][j]/PS;
   		}
	}

    // Calcul de A1
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			A1[i][j] = A[i][j] - L*M[i][j]; // A1 ne peut varier que en fonction de L
			A[i][j] = A1[i][j];
		}
	}

    // On en déduit la valeur de Lambda
    L = puissance_it(A1, x0_ini, x_vect, eps, n); // x2 = valeur propre 2 et M = A1

    return L;
}


int main()
{
	int n;
	float x0[MAX], x1[MAX], x2[MAX], x3[MAX], x4[MAX];
	float A[MAX][MAX], eps=0.000001, L1, L2, L3, L4, L;


	//printf("Veuillez entrer la valeur de n : ");
	//if (scanf("%d",&n) !=1) printf("error");

	n = 4;

	//Initialisation de A ici matrice de Vandermonde
	A[0][0] = 1;
	A[0][1] = 2;
	A[0][2] = 3;
	A[0][3] = 4;

	A[1][0] = 1;
	A[1][1] = 4;
	A[1][2] = 9;
	A[1][3] = 16;

	A[2][0] = 1;
	A[2][1] = 8;
	A[2][2] = 27;
	A[2][3] = 64;

	A[3][0] = 1;
	A[3][1] = 16;
	A[3][2] = 81;
	A[3][3] = 256;

	afficher_mat(A,n);


	// initialisation de x0
	x0[0] = 1;
	x0[1] = 1;
	x0[2] = 1;
	x0[3] = 1;

	afficher_vect(x0,n);

    /* ------------------------------------------------------------------------------------------------------------------------*/
    /* ------------------------------------------------------------------------------------------------------------------------*/
    /* ------------------------------------------------------------------------------------------------------------------------*/


    // EXO 1 - Puissance itérée

    printf("\n\n QUESTION 1 : Methode de la Puissance iteree pour le calcul des valeurs propres ----------------------------*/\n\n");
	L1 = puissance_it(A, x0, x1, eps, n);
	printf("Lambda 1 = %f \n\n", L1);

	printf("Vecteur propre 1 : \n\n");
	afficher_vect(x1, n);


    /* ------------------------------------------------------------------------------------------------------------------------*/
    /* ------------------------------------------------------------------------------------------------------------------------*/
    /* ------------------------------------------------------------------------------------------------------------------------*/


    // EXO 2 - Déflation

    // L2
    printf("\n\n QUESTION 2 : Methode de Deflation pour le calcul des valeurs propres -------------------------------------*/\n\n");
	L2 = deflation_it(A, x0, x1, x2, L1, n, eps);
	printf("Lambda 2 = %f \n\n", L2);

	printf("Vecteur propre 2: \n\n");
	afficher_vect(x2, n);

    // L3
	L3 = deflation_it(A, x0, x2, x3, L2, n, eps);
	printf("Lambda 3 = %f \n\n", L3);

	printf("Vecteur propre 3: \n\n");
	afficher_vect(x3, n);

    // L4
	L4 = deflation_it(A, x0, x3, x4, L3, n, eps);
	printf("Lambda 4 = %f \n\n", L4);

	printf("Vecteur propre 4: \n\n");
	afficher_vect(x4, n);


    /* ------------------------------------------------------------------------------------------------------------------------*/
    /* ------------------------------------------------------------------------------------------------------------------------*/
    /* ------------------------------------------------------------------------------------------------------------------------*/


    // EXO 3 - Puissance itérée Inverse

    printf("\n\n QUESTION 3 : Methode de la Puissance inverse pour le calcul des valeurs propres ------------------------------*/\n\n");

    // L4
    L4 = puissance_inverse(A, x0, x1, eps, n);
    printf("Lambda 4 = %f \n\n", L4);

	printf("Vecteur propre 4 : \n\n");
	afficher_vect(x1, n);

// Marche PAS
/*
	// L3
    L3 = puissance_inverse(A, x1, x2, L4, eps, n);
    printf("Lambda 3 = %f \n\n", L3);

	printf("Vecteur propre 3 : \n\n");
	afficher_vect(x2, n);

	// L2
    L2 = puissance_inverse(A, x2, x3, L3, eps, n);
    printf("Lambda 2 = %f \n\n", L2);

	printf("Vecteur propre 2 : \n\n");
	afficher_vect(x3, n);

	// L1
    L1 = puissance_inverse(A, x3, x4, L2, eps, n);
    printf("Lambda 1 = %f \n\n", L1);

	printf("Vecteur propre 1 : \n\n");
	afficher_vect(x4, n);
*/
	return 0;
}

