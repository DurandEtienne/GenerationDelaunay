#include <iostream>
#include <fstream>
#include <vector>
#include<string>

#include "lecture.h"
#include "fonctions.h"
using namespace std;

int main()
{

  vector<vector<double>> Res;
  vector<double> Points,PointsInitiaux;
  vector<double> Aretes;
  vector<double> Triangles;

  string nomdufichier = "test_non_convexe.mesh";
  cout << "le nom du fichier est  : " << nomdufichier << endl;
  Res = lecture(nomdufichier);
  Points = Res[0];
  PointsInitiaux  = Res[0]; //Contient les points de la géométrie, ne sera pas modifié, permet de savoir si les triangles créés appartiennet à la géométrie
  Aretes = Res[1];
  //Triangles = Res[2];

  vector <double> point(2);
  point[0]=0.3;
  point[1]=1.67;
  //printf("le point appartient au triangle  %d ",getTriangle(point,Triangles,Aretes,Points));


   //cout << "cercle " << inCircumscribedCircle(point,7,Triangles,Aretes,Points) << endl;

   //deleteEdgesOnCavityAndReconnect(point, Triangles, Aretes, Points);
   buildBoiteEnglobante(Aretes, Points);
   InitializeMeshBoite(Triangles, Aretes, Points);

   int NbofPoints = Points.size()/3 - 4 ; //On enlève les points de la boîte englobante
   cout << "le nombre de points est  : " << NbofPoints << endl;
   vector <double> point2(2);
   for (int i = 1; i < NbofPoints; i++)
   {
     point2[0] = Points[i*3];
     point2[1] = Points[i*3+1];
     deleteEdgesOnCavityAndReconnect(point2, Triangles, Aretes, Points);
   }

   deleteBoiteEnglobante(Triangles, Aretes, Points, PointsInitiaux);


   write (Triangles, Aretes, Points);
  // vector<int> cavityIndex;
  // cavityIndex = getTriangleCavity(point,Triangles,Aretes,Points);
  // int size = cavityIndex.size();
  // cout << "les cavités du point sont " <<endl;
  // for (int i=0; i<size; i++)
  // {
  //   cout << cavityIndex[i] << endl;
  // }


  return 0;
}
