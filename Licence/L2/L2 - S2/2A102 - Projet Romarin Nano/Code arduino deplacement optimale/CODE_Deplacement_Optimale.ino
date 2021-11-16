
int throttlePin1 = 0;
int throttlePin2 = 1; //Entrée du joystick axe y 
int throttlePin3 = 2; //Entrée du joystick axe x
 
int val = 0;
int val2 = 0;
int valmhrotation1 = 0;
int valmhrotation2 = 0;
int valmhavance1 = 0;
int valmhavance2 = 0;

void setup(){
    
  Serial.begin(9600);   // Commencer la lecture / ouvrir le moniteur série

  initialisation();     // Fonction réalisant l'initialisation des ESC
  delay(2000);

  int throttle1 = 2047;
  int throttle2 = 2047;
  int throttle3 = 2047;
  delay(500); 
  throttle1 = 0;
  throttle2 = 0;
  throttle3 = 0;
  delay(500);

  }
 
void loop(){

  // Commande de remontée du moteur vertical
  int throttle1 = analogRead(throttlePin1);
  val = analogRead(throttle1);
  val2 = map(throttle1,0, 2047, -127, 126); // Le moteur va aller dans la direction avant quand le potentiomètre va de sa valeur 511 jusqu'à 1023
  analogWrite(9,val2);

  int throttle2 = analogRead(throttlePin2);
  valmhavance1 = analogRead(throttle2);
  valmhavance2 = map(throttle2,0, 1023, 0, 255);

  int throttle3 = analogRead(throttlePin3);
  valmhrotation1 = analogRead(throttle3);
  valmhrotation2 = map(throttle3,0, 1023, 0, 255);


  // test moteur 2 axe X
  // analogWrite(6,valmhrotation2);

  // test moteur 3 axe Y
  // analogWrite(10,valmhavance2);

  //code rotation complet
                   
    if (valmhrotation2 >= 140 && valmhavance2 < 150 ) {
      analogWrite(10,valmhrotation2); //moteur gauche activé
      analogWrite(6,0); //moteur droit désactivé
    } 
    else if (valmhrotation2 < 140 && valmhrotation2 > 130 && valmhavance2 < 150) {
      analogWrite(10,0); //moteur gauche desactivé
      analogWrite(6,0); //moteur droit desactivé
    } 
    
    else if (valmhrotation2 <= 130 && valmhavance2 < 150) {
      int valmhh = 275-valmhrotation2;
      analogWrite(6,valmhh); //moteur droit activé
      analogWrite(10,0); //moteur gauche désactivé
    }   
    // Activer le moteur de gauche pour tourner à droite
    else if (valmhavance2 >= 150) {
      analogWrite(10,valmhavance2); //moteur droit activé
      analogWrite(6,valmhavance2); //moteur gauche activé
    }
    
  //Amélioration Code
   // Moteur X = (throttle2, 0, 1023, -128, 127)    // On affectera aux moteurs droit et gauche les valeurs de moteur X selon les besoins
    // Moteur Y = map(throttle3, 0, 1023, 0, 255);    // On affectera aux moteurs verticaux les valeurs de moteur Y selon les besoins
    
    // Moteur Gauche prendra les valeurs de moteur X pour : 50 < X < 130
    // Moteur Droit prendra les valeurs de moteur X pour : 140 < X < 250
    
    
    // if(moteurX > 130 && moteurX < 140 && moteurY < 10)
            //    analogWrite(7,0);       // On éteint le moteur gauche
            //    analogWrite(8,0);       // On éteint le moteur droit
    
    // else if ( moteurY >= abs(2*moteurX))        // Si pas possible trouver une alternative à absolue. 
                                                  // Si la valeur lue par l'axe verticale est supérieure à y = 2*x, alors on avance. 
                                                  // Voir possibilité et terme exacte pour x et y
                                                  // Normalement X augmente linéairement et donc y aussi et on a pas besoin de fixer de condition sur x                       
    //      analogWrite(7,moteurG)    // on allume le moteur verticale.
    //
    // else if (moteurY <= abs(moteurX/2)          // Si la valeur lue par l'axe verticale est inférieure à y = x/2 (voir possibilité et terme exacte pour x et y), alors :
    //      if (moteurX < 30 && moteurX > 130)    // Si on demande de tourner à gauche
    //          moteurX = 175 - moteurX;          // On affecte à mteur X les valeurs pour lesquelles le moteur fonctionne, soit l plage de fonctionnement du oteur droit ([140,230])
    //          analogWrite(7,moteurX);           // Moteur droit activé afin de tourner à gauche
    //
    //      else if (140 < moteurX && moteurX < 230)
    //          analogWrite(8, moteurX);
    //      else serial.println('erreur')
    
    // else if (abs(moteurX/2) < moteurY && moteurY < abs(moteurX*2))  // Si on se trouve sur une des deux diagonales (tourner en diagonale) : 
    //      if (140 < moteurX)      // Besoin que dune condition sur x ???
                                    // On veut ici tourner en diagonal vers la droite
    //         analogWrite(8, moteurX);   // moteur gauche à plein régime
    //         analogWrite(7, moteurX/2); // moteur droite à une puissance moindre (divis par 2 ici)         
    //      else if (moteurX < 130)
    //         analogWrite(8, moteurX/2);   // moteur gauche à plein régime
    //         analogWrite(7, moteurX); // moteur droite à une puissance moindre (divis par 2 ici)         
    //
    // else serial.println('erreur')           

// Affichage Valeurs Moteurs
      Serial.println(throttle1);
      Serial.println(throttle2);
      Serial.println(val2);
      Serial.println(valmhavance2);
      Serial.println(valmhrotation2); 
}

// Initialisation des Moteurs

void initialisation(int throttle1, int throttle2, int throttle3) // On déclare les différents composants
{
  
  delay(2000);

  int throttle1 = 2047;
  int throttle2 = 2047;
  int throttle3 = 2047;
  delay(500); 
  throttle1 = 0;
  throttle2 = 0;
  throttle3 = 0;
  delay(500);

  return throttle1, throttle2, throttle3
  }
