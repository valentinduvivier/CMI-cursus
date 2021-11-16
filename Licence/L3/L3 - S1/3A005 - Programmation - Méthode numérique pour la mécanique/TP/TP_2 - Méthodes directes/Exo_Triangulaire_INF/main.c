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

//Matrice A triangulaire inférieure
    //Création matrice triangulaire infrieure
    for(i = 0; i < n; i++)
    {
        for(j = 0; j <= i; j++)
        {
            printf("Donner la valeur voulue pour A[%d][%d]\n", i, j);
            scanf("%d", &A[i][j]);
        }
    }

    printf("\n\n-----------------------------------------\n\n");

    //Affichage Matrice A
    for(i = 0; i < n; i++)
    {
        for(j = 0; j <= i; j++)
        {
            printf("A[%d][%d] = %d\n", i + 1, j + 1, A[i][j]);
        }
    }

    printf("\n\n-----------------------------------------\n\n");

//Vecteur x


    for(i = 1; i < n; i++)
    {
        for(j = 0; j < i; j++)
        {
            if(A[j][j] == 0)
            {
                printf("\n - Division par 0 impossible -\n");
                exit(1);
            }

            x[0] = b[0]/A[0][0];
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
