#include <iostream>
#include <math.h>
#include <ctime>
#include <fstream>
#include <string>
#include <vector>

#include "fonctions.h"
#include <bits/stdc++.h>

using namespace std;

vector<double> rotate(vector<double> v, double theta)
{
    vector<double> res(2);
    // convert to radian
    theta *= M_PI / 180.;
    // rotate
    res[0] = cos(theta) * v[0] - sin(theta) * v[1];
    res[1] = sin(theta) * v[0] + cos(theta) * v[1];
    return res;
}
double haveOneCommonEdge(vector<double> triangle1, vector<double> triangle2)
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (triangle1[i] == triangle2[j])
            {
                return triangle1[i];
            }
        }
    }
    return -1;
}

double scalar(vector<double> v1, vector<double> v2)
{
    return v1[0] * v2[0] + v1[1] * v2[1];
}

double norm(vector<double> v)
{
    double s = 0;

    for (int i = 0; i < v.size(); i++)
    {
        s += v[i] * v[i];
    }

    return sqrt(s);
}
double distance(vector<double> v, vector<double> u)
{
    return sqrt((v[0] - u[0]) * (v[0] - u[0]) + (v[1] - u[1]) * (v[1] - u[1]));
}
int getIndexOfElementInArray(vector<double> v, double element)
{
    for (int i = 0; i < v.size(); i++)
    {
        if (v[i] == element)
        {
            return i;
        }
    }
    return -1;
}
void fixEdgesindexing(vector<double> &triangles, double ind)
{
    for (int i = 0; i < triangles.size(); i++)
    {

        if (((int)triangles[i]) % 4 != 3 && triangles[i] > ind)
        {
            triangles[i] -= 1;
        }
    }
}
int checkIfAnEdgeAlreadyExist(vector<double> edges, vector<double> edge)
{
    for (int i = 0; i < edges.size(); i += 3)
    {
        if (edges[i] == edge[0] && edges[i + 1] == edge[1])
        {
            return (i / 3) + 1;
        }
    }
    return -1;
}
int getTriangle(vector<double> point, vector<double> triangles, vector<double> edges, vector<double> vertices)
{
    for (int i = 0; i < triangles.size(); i += 4)
    {
        vector<double> v1(2), v2(2), v3(2);
        // first vertice
        v1[0] = vertices[(triangles[i] - 1) * 3];
        v1[1] = vertices[(triangles[i] - 1) * 3 + 1];
        // second vertice
        v2[0] = vertices[(triangles[i + 1] - 1) * 3];
        v2[1] = vertices[(triangles[i + 1] - 1) * 3 + 1];
        // third vertice
        v3[0] = vertices[(triangles[i + 2] - 1) * 3];
        v3[1] = vertices[(triangles[i + 2] - 1) * 3 + 1];
        // if ((v3[0] == v1[0] && v3[1] == v1[1]) || (v3[0] == v2[0] && v3[1] == v2[1]))
        // {
        //     v3[0] = vertices[edges[triangles[i + 1] * 3 + 1] * 3];
        //     v3[1] = vertices[edges[triangles[i + 1] * 3 + 1] * 3 + 1];
        // }
        vector<double> v1M(2), v1v2(2), v1v3(2);
        v1M[0] = point[0] - v1[0];
        v1M[1] = point[1] - v1[1];
        v1v2[0] = v2[0] - v1[0];
        v1v2[1] = v2[1] - v1[1];
        v1v3[0] = v3[0] - v1[0];
        v1v3[1] = v3[1] - v1[1];
        double lambda2 = scalar(v1M, rotate(v1v3, 90.)) / scalar(v1v2, rotate(v1v3, 90.));
        double lambda3 = scalar(v1M, rotate(v1v2, 90.)) / scalar(v1v3, rotate(v1v2, 90.));
        // checking condition
        if (lambda2 >= 0. && lambda3 >= 0. && lambda2 + lambda3 <= 1.)
        {
            if (lambda2 == 0. || lambda3 == 0. || lambda2 + lambda3 == 1.)
            {
                cout << "You choosed a point in an edge !!! " << endl;
            }

            return i / 4 + 1;
        }
    }
}
vector<int> getTriangleNeighbors(int triangleIndex, vector<double> triangles)
{
    vector<int> res;
    for (int i = 0; i < triangles.size(); i += 4)
    {
        if (haveOneCommonEdge(vector<double>(triangles.begin() + triangleIndex, triangles.begin() + triangleIndex + 2), vector<double>(triangles.begin() + i, triangles.begin() + i + 2)) != -1)
            ; //[triangleIndex:(triangleIndex + 3)], triangles [i:i + 3]))
        {
            res.push_back(i / 4);
        }
    }
    return res;
}
bool inCircumscribedCircle(vector<double> point, int triangleIndex, vector<double> triangles, vector<double> edges, vector<double> vertices)
{
    // Retrieving triangle vertices
    int i = (triangleIndex - 1) * 4;
    vector<double> v1(2), v2(2), v3(2);
    // first vertice
    v1[0] = vertices[(triangles[i] - 1) * 3];
    v1[1] = vertices[(triangles[i] - 1) * 3 + 1];
    // second vertice
    v2[0] = vertices[(triangles[i + 1] - 1) * 3];
    v2[1] = vertices[(triangles[i + 1] - 1) * 3 + 1];
    // third vertice
    v3[0] = vertices[(triangles[i + 2] - 1) * 3];
    v3[1] = vertices[(triangles[i + 2] - 1) * 3 + 1];
    // if ((v3[0] == v1[0] && v3[1] == v1[1]) || (v3[0] == v2[0] && v3[1] == v2[1]))
    // {
    //     v3[0] = vertices[edges[triangles[i + 1] * 3 + 1] * 3];
    //     v3[1] = vertices[edges[triangles[i + 1] * 3 + 1] * 3 + 1];
    // }
    // Creating edges vectors
    vector<double> v1M(2), v1v2(2), v1v3(2);
    v1M[0] = point[0] - v1[0];
    v1M[1] = point[1] - v1[1];
    v1v2[0] = v2[0] - v1[0];
    v1v2[1] = v2[1] - v1[1];
    v1v3[0] = v3[0] - v1[0];
    v1v3[1] = v3[1] - v1[1];

    // Finding weight coefficients of our point on our triangle

    double alpha2 = norm(v1v2) * norm(v1v2);
    double alpha3 = norm(v1v3) * norm(v1v3);
    double alpha23 = scalar(v1v2, v1v3);
    double lambda2 = (alpha3 * (alpha2 - alpha23)) / (2 * (alpha2 * alpha3 - alpha23 * alpha23));
    double lambda3 = (alpha2 * (alpha3 - alpha23)) / (2 * (alpha2 * alpha3 - alpha23 * alpha23));
    // Center vector position
    vector<double> v1C(2);
    v1C[0] = lambda2 * v1v2[0] + lambda3 * v1v3[0];
    v1C[1] = lambda2 * v1v2[1] + lambda3 * v1v3[1];
    // Now the position of M relatively to v1 is known let's use it
    // positon relatively to each other
    vector<double> CM(2);
    CM[0] = v1M[0] - v1C[0];
    CM[1] = v1M[1] - v1C[1];
    // verify condition
    return (norm(CM) <= norm(v1C));
}

vector<int> getTriangleCavity(vector<double> point, vector<double> triangles, vector<double> edges, vector<double> vertices)
{
    int nb_triangles = triangles.size() / 4;
    bool cavity;
    vector<int> cavityIndex;
    // loop over triangles to find the ones that the circumscribed circle contain the point
    for (int i = 1; i <= nb_triangles; i++)
    {
        cavity = inCircumscribedCircle(point, i, triangles, edges, vertices);
        if (cavity == true)
        {
            // adding traingle to the list
            cavityIndex.push_back(i);
        }
    }
    return cavityIndex;
}
// void deleteEdgesOnCavityAndReconnect(vector<double> point, vector<double> &triangles, vector<double> &edges, vector<double> &vertices)
// {
//     vector<int> trianglesInCavity = getTriangleCavity(point, triangles, edges, vertices);
//     vector<double> triangle1(3), triangle2(3);
//     vector<bool> bTriangles(triangles.size() / 4, true), bEdges(edges.size() / 3, true), bVertices(vertices.size() / 3, true);
//     vector<double> pointsToConnect, boundryEdgesToConnect, edgesToDelete;
//     int n = trianglesInCavity.size();
//     int nbOfVertices = vertices.size() / 3;
//     int nbOfEdges = edges.size() / 3;
//     int nbOfTriangles = triangles.size() / 4;
//     int ibegin, jbegin;
//     // deleting elements from cavity and getting vertices of cavity boundary
//     for (int i = 0; i < n; i++)
//     {
//         // constructing first triangle and add edges extremes which are actually the cavity boundry vertices
//         ibegin = (int)(4 * (trianglesInCavity[i] - 1));
//         for (int r = 0; r < 3; r++)
//         {
//             triangle1[r] = triangles[ibegin + r];
//             // adding boundry edges
//             if (getIndexOfElementInArray(boundryEdgesToConnect, triangle1[r])==-1)
//             {
//                 cout << triangle1[r] << endl;
//                 boundryEdgesToConnect.push_back(triangle1[r]);
//             }
//             // getting vertices of cavity boundary
//             int ibeg = (int)(3 * (triangle1[r]-1));
//             // add first boundry vertice if it does not exist already
//             if (getIndexOfElementInArray(pointsToConnect, edges[ibeg])==-1)
//             {
//                 pointsToConnect.push_back(edges[ibeg]);
//             }
//             // add second boundry vertice if it does not exist already
//             if (getIndexOfElementInArray(pointsToConnect, edges[ibeg + 1])==-1)
//             {
//                 pointsToConnect.push_back(edges[ibeg + 1]);
//             }
//         }
//         // second loop for upcomming triangles to make combinations and get common edges
//         for (int j = i + 1; j < n; j++)
//         {
//             // constructiion of second triangle to check if they have a common edge
//             jbegin = (int)(4 * (trianglesInCavity[j] - 1));
//             for (int r = 0; r < 3; r++)
//             {
//                 triangle2[r] = triangles[jbegin + r];
//             }
//             // checking if they have a common edge and erase it in case
//             if (haveOneCommonEdge(triangle1, triangle2)!=-1)
//             {
//                 // getting the beginnig position of the common edge
//                 double communBndIndex = haveOneCommonEdge(triangle1, triangle2);
//                 edgesToDelete.push_back(communBndIndex);
//                 // deleting commun edges from the vector of cavity boundry edges
//                 if (getIndexOfElementInArray(boundryEdgesToConnect, communBndIndex)!=-1)
//                 {
//                     boundryEdgesToConnect.erase(boundryEdgesToConnect.begin() + getIndexOfElementInArray(boundryEdgesToConnect, communBndIndex));
//                 }
//             }
//         }
//     }
//     // erasing triangles after using them
//     int h=0;
//     for (int i = 0; i < n; i++)
//     {
//         ibegin = 4 * (trianglesInCavity[i] - 1)-h*4;
//         triangles.erase(triangles.begin() + ibegin, triangles.begin() + ibegin + 4);
//         h++;
//         nbOfTriangles -= 1;
//     }
//     // adding the new point
//     vertices.push_back(point[0]);
//     vertices.push_back(point[1]);
//     vertices.push_back(0);
//     nbOfVertices++;
//     // Making new edges and triangles
//     int edg1, edg2, edg3;
//     for (int q = 0; q < boundryEdgesToConnect.size(); q++)
//     {
//         edg1 = (int)boundryEdgesToConnect[q];
//         if (checkIfAnEdgeAlreadyExist(edges, {edges[(edg1 - 1) * 3], (double) nbOfVertices})==-1)
//         {
//             edges.push_back(edges[(edg1 - 1) * 3]);
//             edges.push_back(nbOfVertices);
//             edges.push_back(1);
//             nbOfEdges++;
//             edg2 = nbOfEdges;
//         }
//         else
//         {
//             edg2=checkIfAnEdgeAlreadyExist(edges, {edges[(edg1 - 1) * 3],(double) nbOfVertices});
//         }
//         if (checkIfAnEdgeAlreadyExist(edges, {edges[(edg1 - 1) * 3+1], (double)nbOfVertices})==-1)
//         {
//             edges.push_back(edges[(edg1 - 1) * 3+1]);
//             edges.push_back(nbOfVertices);
//             edges.push_back(1);
//             nbOfEdges++;
//             edg3 = nbOfEdges;
//         }
//         else
//         {
//             edg3=checkIfAnEdgeAlreadyExist(edges, {edges[(edg1 - 1) * 3], (double) nbOfVertices});
//         }
//         triangles.push_back(edg1);
//         triangles.push_back(edg2);
//         triangles.push_back(edg3);
//         triangles.push_back(2);
//         nbOfTriangles++;
//     }
//     // finally erase edges
//     for (int p = 0; p < edgesToDelete.size(); p++)
//     {
//         ibegin = (int)(edgesToDelete[p] - 1) * 3-h*3;
//         edges.erase(edges.begin() + ibegin, edges.begin() + ibegin + 3);
//         fixEdgesindexing(triangles, edgesToDelete[p]);
//         fixEdgesindexing(edgesToDelete,edgesToDelete[p]);
//     }
// }
vector<double> measureMeshQuality(vector<double> triangles,vector<double> vertices)
{
    vector<double> res;
    for (int i = 0; i < triangles.size(); i+=4)
    {
        // computing the edges length square and adding them
        double s12x=pow(vertices[3*(triangles[i]-1)]-vertices[3*(triangles[i+1]-1)],2);
        double s13x=pow(vertices[3*(triangles[i]-1)]-vertices[3*(triangles[i+2]-1)],2);
        double s23x=pow(vertices[3*(triangles[i+1]-1)]-vertices[3*(triangles[i+2]-1)],2);
        double s12y=pow(vertices[3*(triangles[i]-1)+1]-vertices[3*(triangles[i+1]-1)+1],2);
        double s13y=pow(vertices[3*(triangles[i]-1)+1]-vertices[3*(triangles[i+2]-1)+1],2);
        double s23y=pow(vertices[3*(triangles[i+1]-1)+1]-vertices[3*(triangles[i+2]-1)+1],2);
        double num= (sqrt(3)/12.)*(s12x+s12y+s13x+s13y+s23x+s23y);
        double a=sqrt(s12x+s12y),b=sqrt(s13x+s13y),c=sqrt(s23x+s23y);
        double p=(a+b+c)/2.;
        // triangle aire with Héron Formula
        double denum=sqrt(p*(p-a)*(p-b)*(p-c));
        res.push_back(num/denum);
    }
    return res;
}
bool edgeAlreadyHere(vector<vector<double>> edgs,vector<double>edg)
{
    for (int i = 0; i < edgs.size(); i++)
    {
        if((edgs[i][0]==edg[0] && edgs[i][1]==edg[1]) || (edgs[i][0]==edg[1] && edgs[i][1]==edg[0]))
        {
            return true;
        }
    }
    return false;
}
bool alreadyConnected(vector<double>Points,double id)
{
    for (int i = 0; i < Points.size(); i++)
    {
        if(Points[i]==id)
        {
            return true;
        }
    }
    return false;
}
bool alreadyDeleted(vector<int> Points, int id)
{
    for (int i = 0; i < Points.size(); i++)
    {
        if (Points[i] == id)
        {
            return true;
        }
    }
    return false;
}
void fixindexing(vector<int> &triangles, int ind)
{
    for (int i = 0; i < triangles.size(); i++)
    {
        if (triangles[i] > ind)
        {
            triangles[i] -= 1;
        }
    }
}
void eraseEdge(vector<double> &edgs,vector<double> edg)
{
    for (int i = 0; i < edgs.size(); i+=3)
    {
        if((edgs[i]==edg[0] && edgs[i+1]==edg[1]) || (edgs[i]==edg[1] && edgs[i+1]==edg[0]))
        {
            edgs.erase(edgs.begin()+i,edgs.begin()+i+3);
        }
    }
}
void deleteEdgesOnCavityAndReconnect(vector<double> point, vector<double> &triangles, vector<double> &edges, vector<double> &vertices)
{
    vector<int> trianglesInCavity = getTriangleCavity(point, triangles, edges, vertices);
    int n = trianglesInCavity.size(), ibeg;
    vector<vector<double>> edgesOnBoundry,edgesToNotKeep;
    vector<double> pointsToConnect;
    // Adding the new point
    vertices.push_back(point[0]);
    vertices.push_back(point[1]);
    vertices.push_back(0);
    // getting its index
    int nv = vertices.size()/3;

    for (int i = 0; i < n; i++)
    {
        ibeg = 4 * (trianglesInCavity[i] - 1);
        // connecting all vertices of triangles on cavity
        if (!alreadyConnected(pointsToConnect, triangles[ibeg]))
        {
            edges.push_back(triangles[ibeg]);
            edges.push_back((double)nv);
            edges.push_back(1);
            pointsToConnect.push_back(triangles[ibeg]);
        }
        if (!alreadyConnected(pointsToConnect, triangles[ibeg+1]))
        {
            edges.push_back(triangles[ibeg+1]);
            edges.push_back((double)nv);
            edges.push_back(1);
            pointsToConnect.push_back(triangles[ibeg+1]);
        }
        if (!alreadyConnected(pointsToConnect, triangles[ibeg+2]))
        {
            edges.push_back(triangles[ibeg+2]);
            edges.push_back((double)nv);
            edges.push_back(1);
            pointsToConnect.push_back(triangles[ibeg+2]);
        }
        if (!edgeAlreadyHere(edgesOnBoundry, {triangles[ibeg], triangles[ibeg + 1],1}))
        {
            edgesOnBoundry.push_back({triangles[ibeg], triangles[ibeg + 1],1});
        }
        else
        {
            // edge that a previous triangle already had , commun one !! we are not going to keep it !!
            eraseEdge(edges, {triangles[ibeg], triangles[ibeg + 1], 1});
            edgesToNotKeep.push_back({triangles[ibeg], triangles[ibeg + 1], 1});
        }
        if (!edgeAlreadyHere(edgesOnBoundry, {triangles[ibeg], triangles[ibeg + 2], 1}))
        {
            edgesOnBoundry.push_back({triangles[ibeg], triangles[ibeg + 2], 1});
        }
        else
        {
            // edge that a previous triangle already had , commun one !! we are not going to keep it !!
            eraseEdge(edges, {triangles[ibeg], triangles[ibeg + 2], 1});
            edgesToNotKeep.push_back({triangles[ibeg], triangles[ibeg + 2], 1});
        }
        if (!edgeAlreadyHere(edgesOnBoundry, {triangles[ibeg+2], triangles[ibeg + 1], 1}))
        {
            edgesOnBoundry.push_back({triangles[ibeg+2], triangles[ibeg + 1], 1});
        }
        else
        {
            // edge that a previous triangle already had , commun one !! we are not going to keep it !!
            eraseEdge(edges, {triangles[ibeg+2], triangles[ibeg + 1], 1});
            edgesToNotKeep.push_back({triangles[ibeg+2], triangles[ibeg + 1], 1});
        }
    }
    // deleting all cavity triangles
    int h = 0;
    for (int i = 0; i < n; i++)
    {
        ibeg = 4 * (trianglesInCavity[i] - 1) - h * 4;
        triangles.erase(triangles.begin() + ibeg, triangles.begin() + ibeg + 4);
        h++;
    }
    // creating new traingles
    for (int i = 0; i < edgesOnBoundry.size(); i++)
    {
        // the two vertices of a boundary edge are connected to new point
        // actually edgesOnBoundry still contain commun edge that has been deleted so we need to filter
        if(!edgeAlreadyHere(edgesToNotKeep,edgesOnBoundry[i]))
        {
            triangles.push_back(edgesOnBoundry[i][0]);
            triangles.push_back(edgesOnBoundry[i][1]);
            triangles.push_back((double)nv);
            triangles.push_back(2);
        }
    }
}

void buildBoiteEnglobante(vector<double> &edges, vector<double> &vertices, vector<double> intialVertices)
{
  // Détermination des minimums et maximums en x et y
  int nbOfVertices = intialVertices.size() / 3;
  int nbOfEdges = edges.size() / 3;

  double x, y, x_min(0.), x_max(0.), y_min(0.), y_max(0.), label(0), label2(1);

  //Parcours de la liste des points pour déterminer les min et max
  for (int p = 0; p < nbOfVertices; p++)
  {

    x = intialVertices[3*p];
    y = intialVertices[3*p+1];
    if (x< x_min)
    {
      x_min = x;
    }
    if (x> x_max)
    {
      x_max = x;
    }
    if (y< y_min)
    {
      y_min = y;
    }
    if (y> y_max)
    {
      y_max = y;
    }
  }

  //Rajout des 4 points de la boite à la fin du vecteur contenant les points
  double moyenne_x = abs(x_max-x_min)/2;
  double moyenne_y = abs(y_max-y_min)/2;
  vertices.push_back(x_min-moyenne_x);
  vertices.push_back(y_min-moyenne_y);
  vertices.push_back(label);
  vertices.push_back(x_max+moyenne_x);
  vertices.push_back(y_min-moyenne_y);
  vertices.push_back(label);
  vertices.push_back(x_max+moyenne_x);
  vertices.push_back(y_max+moyenne_y);
  vertices.push_back(label);
  vertices.push_back(x_min-moyenne_x);
  vertices.push_back(y_max+moyenne_y);
  vertices.push_back(label);


  //Rajout des Aretes de la boite à la fin du vecteur contenant les arêtes
  int q = 0;
  edges.push_back(q+1);
  edges.push_back(q+2);
  edges.push_back(label2);
  edges.push_back(q+2);
  edges.push_back(q+3);
  edges.push_back(label2);
  edges.push_back(q+3);
  edges.push_back(q+4);
  edges.push_back(label2);
  edges.push_back(q+4);
  edges.push_back(q+1);
  edges.push_back(label2);
}

void InitializeMeshBoite(vector<double> &triangles, vector<double> &edges,vector<double> &vertices, vector<double> initialVertices)
{
  int nbOfVertices = vertices.size() / 3;
  double label2,label3;
  label2 = 1;
  label3 = 2;
  // Identification du premier point géométrique
  vector<double> premier_point(2);
  premier_point[0] = initialVertices[0];
  premier_point[1] = initialVertices[1];
  vertices.push_back(premier_point[0]);
  vertices.push_back(premier_point[1]);
  vertices.push_back(0);
  //Liaison de ce premier point avec les 4 coins de la boîte
  //   int q = nbOfVertices - 4; //on enlève les 4 coins de la boîte
  edges.push_back(1);
  edges.push_back(5);
  edges.push_back(label2);
  edges.push_back(2);
  edges.push_back(5);
  edges.push_back(label2);
  edges.push_back(3);
  edges.push_back(5);
  edges.push_back(label2);
  edges.push_back(4);
  edges.push_back(5);
  edges.push_back(label2);

  //Création des 4 premiers triangles de la boîte
  triangles.push_back(5);
  triangles.push_back(1);
  triangles.push_back(2);
  triangles.push_back(label3);
  triangles.push_back(5);
  triangles.push_back(2);
  triangles.push_back(3);
  triangles.push_back(label3);
  triangles.push_back(5);
  triangles.push_back(3);
  triangles.push_back(4);
  triangles.push_back(label3);
  triangles.push_back(1);
  triangles.push_back(4);
  triangles.push_back(5);
  triangles.push_back(label3);
}

// void InitializeMeshBoite(vector<double> point, vector<double> &triangles, vector<double> &edges, vector<double> &vertices)
// {
//   deleteEdgesOnCavityAndReconnect(point, triangles, edges, vertices);
// }

void getBordersBack(vector<double> &triangles, vector<double> &edges, vector<double> vertices)
{
    vector<int> trianglesToDelete;
    vector<vector<double>> edgesToDelete;
    int nbOfVerticesInside=vertices.size()/3-4;
    // From each box vertice we are going to delete triangles that have at least one non boundry edge
    // the loop is stoped at triangles that doesn't contain any new non-boundry edge, which means that we are about to go inside 
    // boundary edges are detected as edges with consecutive vertices index 
    // or with a diffrence equal to the number of vertices inside -1 (end of cycle)
    for (int i = 1; i < 5; i++)
    {
        // loop over triangles comming from a boundary box vertices
        for (int j = 0; j < triangles.size(); j += 4)
        {
            if (i == (int)triangles[j] || i == (int)triangles[j + 1] || i == (int)triangles[j + 2])
            {
                // triangles with  boundary box edges
                if (((triangles[j] < 5) + (triangles[j + 1] < 5) + (triangles[j + 2] < 5)) > 1)
                {
                    if (!alreadyDeleted(trianglesToDelete, j / 4 + 1))
                    {
                        trianglesToDelete.push_back(j / 4 + 1);
                    }
                    continue;
                }
                trianglesToDelete.push_back(j / 4 + 1);
                // identifiying  the third edge and checking if it's a boundry edge 
                int n = 0, m = 0;
                for (int k = 0; k < 3; k++)
                {
                    if ((int)triangles[j + k] != i)
                    {
                        n = (int)triangles[j + k];
                    }
                }
                for (int k = 0; k < 3; k++)
                {
                    if ((int)triangles[j + k] != i && (int)triangles[j + k] != n)
                    {
                        m = (int)triangles[j + k];
                    }
                }
                if (abs(n - m) == 1 || abs(n - m) == nbOfVerticesInside-1)
                {
                    // not a boundry edge!! we are about to go inside !!, move to the other triangle
                    continue;
                }
                else
                {
                    edgesToDelete.push_back({(double)n, (double)m, 1.});
                    bool finished = false;
                    while (!finished)
                    {
                        int nbofTrianglesFound = 0;
                        for (int p = 0; p < triangles.size(); p += 4)
                        {
                            // collecting non-boundary edges and triangle that we should delete 
                            if (edgeAlreadyHere(edgesToDelete, {triangles[p], triangles[p + 1], 1}) || edgeAlreadyHere(edgesToDelete, {triangles[p], triangles[p + 2], 1}) || edgeAlreadyHere(edgesToDelete, {triangles[p + 2], triangles[p + 1], 1}))
                            {
                                if (!alreadyDeleted(trianglesToDelete, p / 4 + 1))
                                {
                                    trianglesToDelete.push_back(p / 4 + 1);
                                    nbofTrianglesFound++;
                                    if (abs(triangles[p] - triangles[p + 1]) != 1 && abs(triangles[p] - triangles[p + 1]) != nbOfVerticesInside-1 )
                                    {
                                        edgesToDelete.push_back({triangles[p], triangles[p + 1], 1});
                                    }
                                    if (abs(triangles[p] - triangles[p + 2]) != 1 && abs(triangles[p] - triangles[p + 2]) != nbOfVerticesInside-1)
                                    {
                                        edgesToDelete.push_back({triangles[p], triangles[p + 2], 1});
                                    }
                                    if (abs(triangles[p + 2] - triangles[p + 1]) != 1 && abs(triangles[p + 2] - triangles[p + 1]) != nbOfVerticesInside-1)
                                    {
                                        edgesToDelete.push_back({triangles[p + 2], triangles[p + 1], 1});
                                    }
                                }
                            }
                        }
                        // if no triangle with new non-boundary edge is detected finish the while loop
                        if (nbofTrianglesFound == 0)
                        {
                            finished = true;
                        }
                    }
                }
            }
        }
    }
    // Deleting triangles that we found
    for (int i = 0; i < trianglesToDelete.size(); i++)
    {
        triangles.erase(triangles.begin() + (trianglesToDelete[i] - 1) * 4, triangles.begin() + (trianglesToDelete[i]) * 4);
        fixindexing(trianglesToDelete, trianglesToDelete[i]);
    }
    // Detecting edges comming from boundry box
    for (int i = 0; i < edges.size(); i+=3)
    {
        if(edges[i]<5 || edges[i+1]<5)
        {
            edgesToDelete.push_back({edges[i],edges[i+1],1});
        }
    }
    // Deleting all edges that we found
    for (int i = 0; i < edgesToDelete.size(); i++)
    {
        eraseEdge(edges, edgesToDelete[i]);
    }
}

// giving labels to boundary edges to make them visible in plots
void LabelizeBorderEdges(std::vector<double> &edges)
{
  for (int i=0; i<edges.size()/3; i++)
  {
    if (abs(edges[i*3+1] - edges[i*3]) == 1 || abs(edges[i*3+1] - edges[i*3]) == 8)
    {
      edges[i*3+2] = 5;
    }
  }
}


void AddPointsInMesh(int IndexTriangle, vector<double> &triangles, vector<double> &edges, vector<double> &vertices)
{
  // adding new points on gravity center to boost the mesh
  int IndexPoint1, IndexPoint2, IndexPoint3;
  IndexPoint1 = triangles[(IndexTriangle-1)*4];
  IndexPoint2 = triangles[(IndexTriangle-1)*4 + 1];
  IndexPoint3 = triangles[(IndexTriangle-1)*4 + 2];
  double x_G,y_G,x1,x2,x3,y1,y2,y3;
  x1 = vertices[(IndexPoint1-1)*3];
  y1 = vertices[(IndexPoint1-1)*3+1];
  x2 = vertices[(IndexPoint2-1)*3];
  y2 = vertices[(IndexPoint2-1)*3+1];
  x3 = vertices[(IndexPoint3-1)*3];
  y3 = vertices[(IndexPoint3-1)*3+1];
  x_G = (x1+x2+x3)/3;
  y_G = (y1+y2+y3)/3;

  vector<double> pointToAdd(2);
  pointToAdd[0] = x_G;
  pointToAdd[1] = y_G;
  deleteEdgesOnCavityAndReconnect(pointToAdd, triangles, edges, vertices);
}

void AddPointsInAllTriangles(vector<double> &triangles, vector<double> &edges, vector<double> &vertices)
{
  int NbOfTrianglesDepart = triangles.size()/4;
   for (int i = 0; i< NbOfTrianglesDepart; i++)
  //for (int i = 0; i< 2; i++)
  {
    AddPointsInMesh(i+1, triangles, edges, vertices);
  }
}


vector<double> getMaxEdgeLengthOfTriangles(vector<double> triangles, vector<double> edges, vector<double> vertices)
{
  vector<double> res;
  for (int i = 0; i < triangles.size(); i+=4)
  {
    double max(0.);
    double lengthEdge1, lengthEdge2, lengthEdge3;
    vector<double> point1(2), point2(2), point3(2);
    point1[0] = vertices[(triangles[i]-1)*3];
    point1[1] = vertices[(triangles[i]-1)*3+1];
    point2[0] = vertices[(triangles[i+1]-1)*3];
    point2[1] = vertices[(triangles[i+1]-1)*3+1];
    point3[0] = vertices[(triangles[i+2]-1)*3];
    point3[1] = vertices[(triangles[i+2]-1)*3+1];
    lengthEdge1 = distance(point1, point2);
    lengthEdge2 = distance(point2, point3);
    lengthEdge3 = distance(point3, point1);
    if (lengthEdge1 >= max)
    {
      max = lengthEdge1;
    }
    if (lengthEdge2 >= max)
    {
      max = lengthEdge2;
    }
    if (lengthEdge3 >= max)
    {
      max = lengthEdge3;
    }
    res.push_back(max);
  }
  return res;
}


void AddPointsInBigTriangles(double minSize, vector<double> Sizes, vector<double> &triangles, vector<double> &edges, vector<double> &vertices)
{
  int NbOfTrianglesDepart = triangles.size()/4;
   for (int i = 0; i< NbOfTrianglesDepart; i++)
  {
    if (Sizes[i] > minSize)
    {
      AddPointsInMesh(i+1, triangles, edges, vertices);
    }
  }
}

double MaximumVector(vector<double> vector)
{
  double max(0.);
  for (int i=0; i<vector.size(); i++)
  {
    if (vector[i] >= max)
    {
      max = vector[i];
    }
  }
}
