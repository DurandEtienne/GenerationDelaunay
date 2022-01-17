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
  vector<double> Points;
  vector<double> Aretes;
  vector<double> Triangles;

  string nomdufichier = "carre_simple2.mesh";
  cout << "le nom du fichier est  : " << nomdufichier << endl;
  Res = lecture(nomdufichier);
  Points = Res[0];
  Aretes = Res[1];
  Triangles = Res[2];

  vector <double> point(2);
  point[0]=1.99;
  point[1]=1.99;
  printf("le point appartient au triangle  %d ",getTriangle(point,Triangles,Aretes,Points));


   cout << "cercle " << inCircumscribedCircle(point,7,Triangles,Aretes,Points) << endl;

   deleteEdgesOnCavityAndReconnect(point, Triangles, Aretes, Points);
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
