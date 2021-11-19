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

  string nomdufichier = "carre_simple.mesh";
  cout << "le nom du fichier est  : " << nomdufichier << endl;
  Res = lecture(nomdufichier);
  Points = Res[0];
  Aretes = Res[1];
  Triangles = Res[2];
  vector <double> point(2);
  point[0]=1./4.;
  point[1]=1./4.;
  printf("le point appartient au triangle  %d",getTriangle(point,Triangles,Aretes,Points));
  return 0;
}
