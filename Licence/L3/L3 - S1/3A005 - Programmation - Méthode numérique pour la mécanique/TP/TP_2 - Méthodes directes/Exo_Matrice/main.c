#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define n 3


int main()
{
    int A[n+1][n+1], L[n+1][n+1], U[n+1][n+1];
    float x[n+1], b[n+1] = {9,4,7};

    int i = 1, j = 1, k = 0;
    int s1 = 0, s0 = 0;
    int s2 = 0, s3 = 0;


    printf("\n\n-----------------------------------------\n\n");


    //Initialisation Matrice A : Triangulaire INF
    printf("La matrice A est la suivante : \n");

    for(i = 0; i < n; i++)
    {
        for(j = 0; j <= i; j++)
        {
            printf("Quelles valeurs pour A[%d][%d]\n", i, j);
            scanf("%d", &A[i][j]);
        }
    }


    printf("\n\n-----------------------------------------\n\n");


//Initialisation Ljj et U11
    //Matrice Ljj
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            if (i == j)
                L[i][j] = 1;
        }
    }

    //Valeur U11
    U[1][1] = A[1][1];


//Matrice Uij et Lij
    for(j = 1; j < n; j++)
    {
        for(i = 1; i < j; i++)
        {
            for(k = 1; k < i - 1; k++)
            {
                s1 = s0 + L[i][k]*U[k][j];

                s0 = s1;
            }

            U[i][j] = A[i][j] - s1;
        }

        for(i = j + 1; i < n; i++)
        {
            for(k = 1; k < j - 1; k++)
            {
                s3 = s2 + L[i][k]*U[k][j];
                s2 = s3;
            }

            L[i][j] = (A[i][j] - s3)/U[j][j];
        }
    }

/*
    //Vecteur x

    for(i = 1; i >= 0; i--)
    {
        for(j = i + 1; j < n; j++)
        {
            if(A[i][j] == 0)
            {
                printf("\n - Division par 0 impossible -\n");
                exit(1);
            }

            x[2] = b[2]/A[2][2];
            s1 = A[i][j]*x[j] + s0;

            x[i] = (b[i] - s1)/A[i][i];

            s0 = s1;
        }
    }

    //Affichage Vecteur x
    for(i = 0; i < n; i++)
    {
        printf("x[%d] = %f\n", i + 1, x[i]);
    }
*/
    return 0;
}
