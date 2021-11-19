#include <iostream>
#include <math.h>
#include <ctime>
#include <fstream>
#include <string>
#include<vector>

#include "lecture.h"

using namespace std;



vector<vector<double>> lecture(string nomdufichier)
{
  std::vector<double> Points;
  std::vector<double> Aretes;
  std::vector<double> Triangles;



  ifstream fichier(nomdufichier, ios::in);  // on ouvre le fichier en lecture
  if(fichier)  // si l'ouverture a réussi
  {
    // instructions
    string ligne;
     while(getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
     {
       //cout << ligne << endl;  // on l'affiche
       if (ligne.find("#") !=std::string::npos)
       {
           // Ignore this line (comment)
        }
       if (ligne.find("MeshVersionFormatted 1") !=std::string::npos)
       {
             // Ignore this line (comment)
       }
       if (ligne.find("Dimension 2") !=std::string::npos)
       {
           // Ignore this line (comment)
       }
       if (ligne.find("Vertices") !=std::string::npos)
       {
             int nb_vertices;  // notre variable où sera stocké le nombre de points du maillage
             double(vertice_1); double(vertice_2); double(vertice_label);

             fichier >> nb_vertices;
             cout << "Le nombre de points est : " << nb_vertices << endl;  // on l'affiche
             for (int n = 1; n <= nb_vertices; n++) // Boucle en temps
             {
               fichier >> vertice_1 >> vertice_2 >> vertice_label;
               Points.push_back(vertice_1); Points.push_back(vertice_2); Points.push_back(vertice_label);
             }

             // Affichage du vecteur Points
             // int(taille) = Points.size();
             // //cout << taille << endl;  // on l'affiche
             // cout << "Les coordonnées des différents points sont : " << endl;
             // for (int n = 0; n < taille/3; n++) // Boucle en temps
             // {
             //   cout << " Points n° "+to_string(n+1)+" : " << Points[3*n] << "  " << Points[3*n+1] << " et porte le label " << Points[3*n+2] << endl;  // on l'affiche
             // }
       }


       if (ligne.find("Edges") !=std::string::npos)
       {
         int nb_edges; // variable où est stocké le nombre d'arêtes du maillage
         double(label_vertice_1); double(label_vertice_2); double(label_arete);

         fichier >> nb_edges;
         cout << "Le nombre d'aretes est : " << nb_edges << endl;  // on l'affiche
         for (int n = 1; n <= nb_edges; n++) // Boucle en temps
         {
           fichier >> label_vertice_1 >> label_vertice_2 >> label_arete;
           Aretes.push_back(label_vertice_1); Aretes.push_back(label_vertice_2); Aretes.push_back(label_arete);
         }

         // Affichage du vecteur Aretes
         // int(taille) = Aretes.size();
         // //cout << taille << endl;  // on l'affiche
         // cout << "Les extrémités des aretes sont les points : " << endl;
         // for (int n = 0; n < taille/3; n++) // Boucle en temps
         // {
         //   cout << " Arete n° "+to_string(n+1)+" est constitué du point " << Aretes[3*n] << " ainsi que du point  " << Aretes[3*n+1] << " et l'arete porte le label " << Aretes[3*n+2] << endl;  // on l'affiche
         // }

       }


       if (ligne.find("Triangles") !=std::string::npos)
       {
         int nb_triangles; // variable où est stocké le nombre d'arêtes du maillage
         double(label_edge_1); double(label_edge_2); double(label_edge_3); double(label_triangle);

         fichier >> nb_triangles;
         cout << "Le nombre de triangles est : " << nb_triangles << endl;  // on l'affiche
         for (int n = 1; n <= nb_triangles; n++) // Boucle en temps
         {
           fichier >> label_edge_1 >> label_edge_2 >> label_edge_3 >> label_triangle;
           Triangles.push_back(label_edge_1); Triangles.push_back(label_edge_2); Triangles.push_back(label_edge_3); Triangles.push_back(label_triangle);
         }

         // Affichage du vecteur Triangles
         // int(taille) = Triangles.size();
         // //cout << taille << endl;  // on l'affiche
         // cout << "Les aretes du triangles sont : " << endl;
         // for (int n = 0; n < taille/4; n++) // Boucle en temps
         // {
         //   cout << " Triangle n° "+to_string(n+1)+" est constitué des aretes " << Triangles[4*n] << " , " << Triangles[4*n+1] << " et " << Triangles[4*n+2] << " et le triangle porte le label " << Triangles[4*n+3] << endl;  // on l'affiche
         // }
       }
    }

    fichier.close();  // on ferme le fichier
    }

  else  // sinon
    {
      cerr << "Impossible d'ouvrir le fichier !" << endl;
    }


    vector<vector<double>> Res;
    Res = {Points, Aretes, Triangles};
    return Res;
}
