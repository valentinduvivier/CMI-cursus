#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float Phi1(float x) {
	return exp(x)/4;
}

float Phi2(float x) {
	return log(4*x);
}

float Psi1( float x , float t) {
	return (1-t)*x + t*Phi1(x) ;
}


float Psi2( float x , float t) {
	return (1-t)*x + t*Phi2(x) ;
}



int main() {
	
	
	float eps = 0.000001;
	float xn,x0,t;
	
	printf("Veuillez entrer la valeur 1 de x0 : ");
	if (scanf("%f",&x0) !=1) printf("error");
	
	
	
	for (t = 1; t <=3;t+=0.1) {
	int i=0;
		
	xn = x0;
	
	while ((Psi1(xn,t) - xn) >= eps || (Psi1(xn,t) - xn) <= -eps ) {
		i++;
		xn = Psi1(xn,t);
		
	}
	
	
	printf("Après %d itérations, et pout teta = %f, la racine 1 trouvée est %f\n",i,t,xn);
	
	
    }
    

	printf("Veuillez entrer la valeur 2 de x0 : ");
	if (scanf("%f",&x0) !=1) printf("error");
	
	
	
	for (t = 1; t <=3;t+=0.05) {
	int i=0;
		
	xn = x0;
	
	while ((Psi2(xn,t) - xn) >= eps || (Psi2(xn,t) - xn) <= -eps ) {
		i++;
		xn = Psi2(xn,t);
		
	}
	
	
	printf("Après %d itérations, et pout teta = %f, la racine 1 trouvée est %f\n",i,t,xn);
	
	
    }




	return 0;
}

