#include<SD.h>
//Sensor connection declarations
const int tempSensor=A0;
const int presSensor=A1;
const int lumSensor=A2;
const int hallPin = 7;

double time;
File myFile;
void setup() {
  // put your setup code here, to run once:
Serial.begin(9600);
pinMode(2,OUTPUT);
pinMode(hallPin,INPUT);
//Initialisation
if(SD.begin(10))
  {Serial.println("Initialisation de la carte SD");}
 else
  {Serial.println("Echec lors de l'initialisation de la carte SD"); return;}

//Ecriture
myFile = SD.open("data.txt",FILE_WRITE);
myFile.println("L'écriture du fichier commence maintenant :");
myFile.close();
}

float temp_digit;
float temp_tension; 
float temp_celsius; 
float lum_digit; 
float pres_digit;
float pres_bar;
int sensorValue;

void loop() {
  // put your main code here, to run repeatedly:
time=millis()/1000.0;   //temps (s)


analogRead(tempSensor);
delay(100);
temp_digit = analogRead(tempSensor);

temp_tension = temp_digit*5/1023;
temp_celsius = temp_tension*100;

analogRead(lumSensor);
delay(100);
lum_digit = analogRead(lumSensor);

analogRead(presSensor);
delay(100);
pres_digit = analogRead(presSensor);
pres_bar = (pres_digit + 68)/448;

//lecture du capteur a Effet Hall
sensorValue = digitalRead(hallPin);
sensorValue = not(sensorValue);

if (sensorValue == HIGH){
  digitalWrite(2,HIGH);
  }
else {
  digitalWrite(2,LOW);
  }
 
if (sensorValue == HIGH){
//Affichage écran
  Serial.print("Thermocouple(°C) :"); Serial.print(temp_celsius);
  Serial.print("   Lumiere :"); Serial.print(lum_digit);
  Serial.print("   Presssion :"); Serial.println(pres_bar);
//Ecriture dans un fichier txt
  myFile = SD.open("data.txt",FILE_WRITE);
  myFile.print("Thermocouple(°C) :"); myFile.print(temp_celsius);
  myFile.print("   Lumiere :"); myFile.print(lum_digit);
  myFile.print("   Pression :"); myFile.println(pres_bar);
  myFile.close();
}
else{

}
}
}
