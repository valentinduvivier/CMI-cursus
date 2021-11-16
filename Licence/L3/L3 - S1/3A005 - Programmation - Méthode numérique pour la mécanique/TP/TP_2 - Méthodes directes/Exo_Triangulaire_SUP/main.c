#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define n 3

int main()
{
    int A[n+1][n+1];
    float x[n+1], b[n+1] = {9,4,7};

    int i, j;
    int s1 = 0, s0 = 0;

    printf("La matrice A est la suivante : \n");

//Matrice A triangulaire supérieure
    //Création matrice triangulaire supérieure
    for(j = 0; j < n; j++)
    {
        for(i = 0; i <= j; i++)
        {
            printf("Donner la valeur voulue pour A[%d][%d]\n", i, j);
            scanf("%d", &A[i][j]);
        }
    }

    printf("\n\n-----------------------------------------\n\n");

    //Affichage Matrice A
    for(j = 0; j < n; j++)
    {
        for(i = 0; i <= j; i++)
        {
            printf("A[%d][%d] = %d\n", i + 1, j + 1, A[i][j]);
        }
    }

    printf("\n\n-----------------------------------------\n\n");

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

    return 0;
}
