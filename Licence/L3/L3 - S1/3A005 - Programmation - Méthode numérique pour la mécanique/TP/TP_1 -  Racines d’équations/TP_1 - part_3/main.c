#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Iter_Max 50
#define E 10e-6

float F(float x);
float F_prime(float x);
float F_sec(float x);

int main()
{
    float x0 = 0, x1 = 0;
    float x2 = 0, x0_2 = 0;

    int n2 = 0;
    int n = 0;

    printf("Qoui\n");
    scanf("%f", &x0);

    printf("----------------------------------\n\n");

    do
    {
        x1 = x0 - F(x0)/F_prime(x0) - F_sec(x0)*(pow(F(x0),2))/(2*(pow(F_prime(x0),3)));
        n++;

        if(fabs(x1 - x0) < E)
            break;

        x0 = x1;

    }while(n < Iter_Max);

    printf("Qoui\n");
    scanf("%f", &x0_2);

    printf("----------------------------------\n\n");

    do
    {
        x2 = x0_2 - F(x0_2)/F_prime(x0_2);
        n2++;

        if(fabs(x2 - x0_2) < E)
            break;

        x0_2 = x2;

    }while(n2 < Iter_Max);

    printf("La racine est %f et elle a ete trouvee en %d iterations\n", x0, n);
    printf("La racine est %f et elle a ete trouvee en %d iterations\n\n", x1, n);

    printf("La racine est %f et elle a ete trouvee en %d iterations\n", x0_2, n2);
    printf("La racine est %f et elle a ete trouvee en %d iterations\n\n", x2, n2);

    return 0;
}


float F(float x)
{
    return exp(x) - 4*x;
}

float F_prime(float x)
{
    return exp(x) - 4;
}

float F_sec(float x)
{
    return exp(x);
}
