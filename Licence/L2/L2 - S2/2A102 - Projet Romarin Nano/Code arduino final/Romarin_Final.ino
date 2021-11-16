
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
  valmhavance2 = map(throttle2, 0, 1023, 0, 255);

  int throttle3 = analogRead(throttlePin3);
  valmhrotation1 = analogRead(throttle3);
  valmhrotation2 = map(throttle3, 0, 1023, 0, 255);


  // test moteur 2 axe X
  analogWrite(6,valmhrotation2);
  

  // test moteur 3 axe Y
  analogWrite(10,valmhavance2);

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

//Afichage Valeurs Moteurs (Vérification bon fonctionnement)
    Serial.println(throttle1);
    Serial.println(throttle2);
    Serial.println(val2);
    Serial.println(valmhavance2);
    Serial.println(valmhrotation2); 
}
