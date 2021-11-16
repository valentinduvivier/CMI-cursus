#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float f(float x) { 
	return exp(x)-4*x;
}

float f1(float x) { 
	return exp(x)-4;
}


float Phi(float x) {
	return (x - f(x)/f1(x));
}

int main() {
	float eps = 0.000001;
	float xn;
	printf("Veuillez entrer la valeur de x0 : ");
	if (scanf("%f",&xn) !=1) printf("error");
	int i=0;
	while ( (Phi(xn) - xn) >= eps || (Phi(xn) - xn) <= -eps ) {
		i++;
		xn = Phi(xn);
	}
	printf("Après %d itérations, la racine trouvée est %f",i,xn);
	return 0;
}
